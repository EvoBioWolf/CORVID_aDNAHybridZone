#this script plots called genotype datasets with SNPrelate
#also assess the quality of the datasets by looking at missingness and reference bias
#plots genotype PVA from EMU
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager") BiocManager::install("SNPRelate")
library(devtools)
library(ggplot2)
library(SNPRelate)
library(openxlsx")
library(dplyr)
library(psych")
library(ggshadow")
library(ggnewscale")

my_colors <- c(
  "E"  = "#7F3B08",
  "GB" = "#B35806",
  "F"  = "#E08214",
  "D"  = "#FDB863",
  "BE" = "#FEE0B6",
  "NL"=  "#FFDD00",
  "HZ"  = "#E78AC3", #uli's sample: austria, eastgermany,northitaly
  "IRE"= "#D6266C",
  "S" =  "#A6CEE3",
  "PL"=  "#2899B2",
  "RUS"= "#084081",
  "IT" = "#66C2A5",
  "B"  = "#117733",
  "ISR"= "#7A2A7E", 
  "IRQ"= "#54278F")

setwd("/PATH/05_aDNA")
setwd("01_angsd_enrichedallfresh_rescaled")
pop<-read.table("pop.txt", header=FALSE, sep="", col.names = c("pop","sample","age","type"))
# PCA with SNPRelate #

#---------------------------------------------------------------------------------------------------------------------------------------------
# Divergent SNP set -----------------------------------------------------------------------------------------------------------------------------------
vcf.fn <- "enrichedallfresh_rescaled_outlier_111_geno_maxmis_q20_plink.vcf.gz"
#vcf1.fn <- "enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20_dp3_plink.vcf.gz"
snpgdsVCF2GDS(vcf1.fn,"vcf1.gds", method="copy.num.of.ref")
snpgdsSummary("vcf1.gds")
vcf1.gdsfile <- snpgdsOpen("vcf1.gds")

# Count REF and ALT alleles per sample to assess ref bias
geno <- snpgdsGetGeno(vcf1.gdsfile)
ref.counts <- 2 * rowSums(geno == 0, na.rm=TRUE) + rowSums(geno == 1, na.rm=TRUE)
alt.counts <- 2 * rowSums(geno == 2, na.rm=TRUE) + rowSums(geno == 1, na.rm=TRUE)
ref.ratio <- ref.counts / (ref.counts + alt.counts)
sample.stats <- data.frame(REF = ref.counts,ALT = alt.counts,REF_ratio = ref.ratio)
refbias <- cbind(sample.stats,pop)%>%
  filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE", pop != "COR") %>%
  filter(pop != "ori") %>%
  mutate(group= ifelse(type != "present", "ancient", "present"))%>%
  mutate(pop = factor(pop, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","IT","RUS","B","ISR", "IRQ"))) 
ref <- ggplot(refbias, aes(x = factor(group), y = REF_ratio)) +guides(color=guide_legend("Locality", override.aes = list(size=5))) +
  geom_boxplot(alpha = 0.6) + geom_point(aes(color = pop), size = 2, alpha = 0.8) +
  scale_color_manual(values = my_colors) + xlab("Sample Type") + ylab("REF/ALT Ratio") + labs(title="1111 outlier SNPs") + theme_minimal()

vcf.pca <- snpgdsPCA(vcf1.gdsfile, autosome.only=FALSE, missing.rate=NaN) ##use autosome.only=false to include all loci
names(vcf.pca)
pc.percent <- vcf.pca$varprop*100
head(round(pc.percent, 2))
print(pc.percent)
tab <- data.frame(sample.id = vcf.pca$sample.id,
                  EV1 = vcf.pca$eigenvect[,1], # the first eigenvector
                  EV2 = vcf.pca$eigenvect[,2], # the second eigenvector
                  EV3 = vcf.pca$eigenvect[,3], 
                  EV4 = vcf.pca$eigenvect[,4], 
                  EV5 = vcf.pca$eigenvect[,5], 
                  EV6 = vcf.pca$eigenvect[,6], 
                  EV7 = vcf.pca$eigenvect[,7], 
                  EV8 = vcf.pca$eigenvect[,8], 
                  EV9 = vcf.pca$eigenvect[,9], 
                  EV10 = vcf.pca$eigenvect[,10], 
                  stringsAsFactors = FALSE)

dat1 <- data.frame(pop, tab) %>%
  filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE", pop != "COR") %>%
  mutate(pop = factor(pop, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","IT","RUS","B","ISR", "IRQ"))) %>%
  filter(pop != "ori")

#number of biallelic unique SNPs: 35462 / 902 (1111 set)
ggplot() +
  geom_point(data=dat1, aes(x=EV1, y=EV2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat1%>%filter(type=="1500-3000"| type=="more_than_10k"), aes(x=EV1, y=EV2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat1%>%filter(type=="5000-6000" | type=="more_than_16k"), aes(x=EV1, y=EV2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat1%>%filter(pop == "NL" | sample == "LSS003_TE" | age == 5330), aes(x=EV1, y=EV2, fill=pop, shape=type, alpha=type, size=type, color=type))  +
  scale_shape_manual(values = c("100-1500" = 21, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_16k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 3, "1500-3000" = 3, "5000-6000"=3,"more_than_10k"=3,"more_than_16k"=3,"present"=2)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_16k"=1,"present"=1)) +
  scale_color_manual(values = c("100-1500" = "black", "1500-3000" = "black", "5000-6000"="black","more_than_10k"="black","more_than_16k"="black","present"="white")) +
  scale_fill_manual(values=my_colors ) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(21, 22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 0.8), size=c(3, 3, 3, 3, 3, 2))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=3)), color="none", size="none",alpha="none") +
  labs(x=paste("PC1 (",round(pc.percent[1],1),"%)",sep=""), y=paste("PC2 (",round(pc.percent[2],1),"%)",sep="")) +
  scale_x_reverse() +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text.x=element_text (size=10),
        axis.text.y=element_text (size=10),
        legend.position = "right")

#---------------------------------------------------------------------------------------------------------------------------------------------
# Neutral SNP set ----------------------------------------------------------------------------------------------------------------------------
vcf2.fn <- "enrichedallfresh_rescaled_neutral_geno_maxmis_q20_plink.vcf.gz"
#vcf2.fn <- "enrichedallfresh_rescaled_neutral_geno_maxmis_q20_dp3_plink.vcf.gz"
snpgdsVCF2GDS(vcf2.fn,"vcf2.gds", method="copy.num.of.ref")
snpgdsSummary("vcf2.gds")
vcf2.gdsfile <- snpgdsOpen("vcf2.gds")

# Count REF and ALT alleles per sample to assess ref bias
pop<-read.table("pop.txt", header=FALSE, sep="", col.names = c("pop","sample","age","type"))
geno2 <- snpgdsGetGeno(vcf2.gdsfile)
ref.counts <- 2 * rowSums(geno2 == 0, na.rm=TRUE) + rowSums(geno2 == 1, na.rm=TRUE)
alt.counts <- 2 * rowSums(geno2 == 2, na.rm=TRUE) + rowSums(geno2 == 1, na.rm=TRUE)
ref.ratio <- ref.counts / (ref.counts + alt.counts)
sample.stats <- data.frame(REF = ref.counts,ALT = alt.counts,REF_ratio = ref.ratio) 
refbias <- cbind(sample.stats,pop)%>%
  filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE", pop != "COR") %>%
  filter(pop != "ori") %>%
  mutate(group= ifelse(type != "present", "ancient", "present"))%>%
  mutate(pop = factor(pop, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","IT","RUS","B","ISR", "IRQ"))) 
ref1 <- ggplot(refbias, aes(x = factor(group), y = REF_ratio)) +
  geom_boxplot(alpha = 0.6) + geom_point(aes(color = pop), size = 2, alpha = 0.8) + guides(color=guide_legend("Locality", override.aes = list(size=5))) +
  scale_color_manual(values = my_colors) + xlab("Sample Type") + ylab("REF/ALT Ratio") + labs(title="All neutral sites") + theme_minimal()

vcf.pca <- snpgdsPCA(vcf2.gdsfile, autosome.only=FALSE, missing.rate=NaN) ##use autosome.only=false to include all loci
names(vcf.pca)
pc.percent <- vcf.pca$varprop*100
head(round(pc.percent, 2))
print(pc.percent)
tab <- data.frame(sample.id = vcf.pca$sample.id,
                  EV1 = vcf.pca$eigenvect[,1], # the first eigenvector
                  EV2 = vcf.pca$eigenvect[,2], # the second eigenvector
                  EV3 = vcf.pca$eigenvect[,3], 
                  EV4 = vcf.pca$eigenvect[,4], 
                  EV5 = vcf.pca$eigenvect[,5], 
                  EV6 = vcf.pca$eigenvect[,6], 
                  EV7 = vcf.pca$eigenvect[,7], 
                  EV8 = vcf.pca$eigenvect[,8], 
                  EV9 = vcf.pca$eigenvect[,9], 
                  EV10 = vcf.pca$eigenvect[,10], 
                  stringsAsFactors = FALSE)

dat2 <- data.frame(pop, tab) %>%
  filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE", pop != "COR") %>%
  mutate(pop = factor(pop, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","IT","RUS","B","ISR", "IRQ"))) %>%
  filter(pop != "ori")

#The number of biallelic unique SNPs: 177872 / 93431 (dp3)
s1 <- ggplot() +
  geom_point(data=dat2, aes(x=EV1, y=EV2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat2%>%filter(type=="1500-3000"| type=="more_than_10k"), aes(x=EV1, y=EV2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat2%>%filter(type=="5000-6000" | type=="more_than_16k"), aes(x=EV1, y=EV2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat2%>%filter(pop == "NL" | sample == "LSS003_TE" | age == 5330), aes(x=EV1, y=EV2, fill=pop, shape=type, alpha=type, size=type, color=type))  +
  scale_shape_manual(values = c("100-1500" = 21, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_16k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 3, "1500-3000" = 3, "5000-6000"=3,"more_than_10k"=3,"more_than_16k"=3,"present"=2)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_16k"=1,"present"=1)) +
  scale_color_manual(values = c("100-1500" = "black", "1500-3000" = "black", "5000-6000"="black","more_than_10k"="black","more_than_16k"="black","present"="white")) +
  scale_fill_manual(values=my_colors ) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(21, 22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 0.8), size=c(3, 3, 3, 3, 3, 2))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=3)), color="none", size="none",alpha="none") +
  labs(x=paste("PC1 (",round(pc.percent[1],1),"%)",sep=""), y=paste("PC2 (",round(pc.percent[2],1),"%)",sep="")) +
  scale_x_reverse() + theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),
        axis.title.x = element_text(size=11),axis.title.y = element_text(size=11), axis.text.x=element_text (size=10),
        axis.text.y=element_text (size=10),legend.position = "right")

s2 <- ggplot() +
  geom_point(data=dat2, aes(x=EV2, y=EV3, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat2%>%filter(type=="1500-3000"| type=="more_than_10k"), aes(x=EV2, y=EV3, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat2%>%filter(type=="5000-6000" | type=="more_than_16k"), aes(x=EV2, y=EV3, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat2%>%filter(pop == "NL" | sample == "LSS003_TE" | age == 5330), aes(x=EV2, y=EV3, fill=pop, shape=type, alpha=type, size=type, color=type))  +
  scale_shape_manual(values = c("100-1500" = 21, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_16k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 3, "1500-3000" = 3, "5000-6000"=3,"more_than_10k"=3,"more_than_16k"=3,"present"=2)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_16k"=1,"present"=1)) +
  scale_color_manual(values = c("100-1500" = "black", "1500-3000" = "black", "5000-6000"="black","more_than_10k"="black","more_than_16k"="black","present"="white")) +
  scale_fill_manual(values=my_colors ) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(21, 22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 0.8), size=c(3, 3, 3, 3, 3, 2))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=3)), color="none", size="none",alpha="none") +
  labs(x=paste("PC2 (",round(pc.percent[2],1),"%)",sep=""), y=paste("PC3 (",round(pc.percent[3],1),"%)",sep="")) +
  scale_x_reverse() +theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),
        axis.title.x = element_text(size=11),axis.title.y = element_text(size=11), axis.text.x=element_text (size=10),
        axis.text.y=element_text (size=10),legend.position = "right")

plotall2 <- s1 + s2 + 
  plot_layout(widths = c(1, 1)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = c("A", "B")) & theme(plot.margin = unit(c(1, 1, 1, 1), "pt")) 
plotall2 #10x6

#---------------------------------------------------------------------------------------------------------------------------------------------
# Outgroup-ascertained neutral SNP set -------------------------------------------------------------------------------------------------------
vcf3.fn <- "enrichedallfresh_rescaled_neutral_seg72k_geno_maxmis_q20_plink.vcf.gz"
#vcf3.fn <- "enrichedallfresh_rescaled_neutral_seg72k_geno_maxmis_q20_dp3_plink.vcf.gz"
snpgdsVCF2GDS(vcf3.fn,"vcf3.gds", method="copy.num.of.ref")
snpgdsSummary("vcf3.gds")
vcf3.gdsfile <- snpgdsOpen("vcf3.gds")

# Count REF and ALT alleles per sample to assess ref bias
pop<-read.table("pop.txt", header=FALSE, sep="", col.names = c("pop","sample","age","type"))
geno5 <- snpgdsGetGeno(vcf3.gdsfile)
ref.counts <- 2 * rowSums(geno5 == 0, na.rm=TRUE) + rowSums(geno5 == 1, na.rm=TRUE)
alt.counts <- 2 * rowSums(geno5 == 2, na.rm=TRUE) + rowSums(geno5 == 1, na.rm=TRUE)
ref.ratio <- ref.counts / (ref.counts + alt.counts)
sample.stats <- data.frame(REF = ref.counts,ALT = alt.counts,REF_ratio = ref.ratio)
refbias <- cbind(sample.stats, pop) %>%
  filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE", pop != "COR") %>%
  mutate(group= ifelse(type != "present", "ancient", "present"))%>%
  mutate(pop = factor(pop, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","IT","RUS","B","ISR", "IRQ"))) 
ref2 <- ggplot(refbias, aes(x = factor(group), y = REF_ratio)) +
  geom_boxplot(alpha = 0.6) + geom_point(aes(color = pop), size = 2, alpha = 0.8) +
  geom_text_repel(data=refbias %>%filter(REF_ratio<0.1),aes(x = factor(group), y = REF_ratio,label=sample),size=2) +
  scale_color_manual(values = my_colors) + guides(color=guide_legend("Locality", override.aes = list(size=5))) +
  xlab("Sample Type") + ylab("REF/ALT Ratio") + labs(title="Outgroup-asc & transversion")+theme_minimal()
refratio <- ref + ref1 + ref2  + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = c("A", "B", "C")) & theme(plot.margin = unit(c(1, 1, 1, 1), "pt")) #8x6

vcf.pca <- snpgdsPCA(vcf3.gdsfile, autosome.only=FALSE, missing.rate=NaN) ##use autosome.only=false to include all loci
names(vcf.pca)
#variance proportion (%)
pc.percent <- vcf.pca$varprop*100
head(round(pc.percent, 2))
print(pc.percent)
tab <- data.frame(sample.id = vcf.pca$sample.id,
                  EV1 = vcf.pca$eigenvect[,1], # the first eigenvector
                  EV2 = vcf.pca$eigenvect[,2], # the second eigenvector
                  EV3 = vcf.pca$eigenvect[,3], 
                  EV4 = vcf.pca$eigenvect[,4], 
                  EV5 = vcf.pca$eigenvect[,5], 
                  EV6 = vcf.pca$eigenvect[,6], 
                  EV7 = vcf.pca$eigenvect[,7], 
                  EV8 = vcf.pca$eigenvect[,8], 
                  EV9 = vcf.pca$eigenvect[,9], 
                  EV10 = vcf.pca$eigenvect[,10], 
                  stringsAsFactors = FALSE)

dat3 <- data.frame(pop, tab) %>%
  filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE", pop != "COR") %>%
  mutate(pop = factor(pop, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","IT","RUS","B","ISR", "IRQ"))) 

#number of biallelic unique SNPs: 65,794 / 31792 (dp3)
ggplot() +
  geom_point(data=dat3, aes(x=EV2, y=EV3, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat3%>%filter(type=="1500-3000"| type=="more_than_10k"), aes(x=EV2, y=EV3, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat3%>%filter(type=="5000-6000" | type=="more_than_16k"), aes(x=EV2, y=EV3, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat3%>%filter(pop == "NL" | sample == "LSS003_TE" | age == 5330), aes(x=EV2, y=EV3, fill=pop, shape=type, alpha=type, size=type, color=type))  +
  scale_shape_manual(values = c("100-1500" = 21, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_16k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 3, "1500-3000" = 3, "5000-6000"=3,"more_than_10k"=3,"more_than_16k"=3,"present"=2)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_16k"=1,"present"=1)) +
  scale_color_manual(values = c("100-1500" = "black", "1500-3000" = "black", "5000-6000"="black","more_than_10k"="black","more_than_16k"="black","present"="white")) +
  scale_fill_manual(values=my_colors ) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(21, 22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 0.8), size=c(3, 3, 3, 3, 3, 2))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=3)), color="none", size="none",alpha="none") +
  labs(x=paste("PC2 (",round(pc.percent[2],1),"%)",sep=""), y=paste("PC3 (",round(pc.percent[3],1),"%)",sep="")) +
  #scale_x_reverse() +
  scale_y_reverse() +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text.x=element_text (size=10),
        axis.text.y=element_text (size=10),
        legend.position = "right") #6x6

#-------------------------------------------------------------------------------------------------------------------------
# Missing data in neutral sites ------------------------------------------------------------------------------------------
setwd("/PATH/05_aDNA/01_angsd_enrichedallfresh_rescaled")
pop <- read.table("pop.txt") %>%
  mutate(FID = paste(V2,"_",V1, sep=""))
missing <- read.table("enrichedallfresh_rescaled_neutral_geno_maxmis_q20_plink_missing.imiss", header=TRUE, na="N") %>%
#missing <- read.table("enrichedallfresh_rescaled_neutral_geno_maxmis_q20_dp3_plink.imiss", header=TRUE, na="N") %>% mutate(FID=paste(FID,"_",IID,sep="")) %>%
  merge(pop, by="FID") %>%
  filter(V1 != "COR") %>%
  mutate(V1 = factor(V1, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","RUS","IT","B","ISR","IRQ") ))

ggplot(data=missing) +
  geom_point(aes(x=V3, y=F_MISS, color=V1), size=2) +
  geom_text_repel(data=missing%>%filter(V3>1000), aes(x=V3,y=F_MISS,label=V2), size=2.5) +
  geom_text_repel(data=missing%>%filter(FID=="NCP001_TE_B" | FID=="BRW001_TE_PL"), aes(x=V3,y=F_MISS,label=V2), size=2.5) +
  scale_color_manual(values=my_colors) +
  labs(y="Proportion of missing neutral sites", x="Age (years ago)") +
  scale_x_log10(breaks = c(1, 100, 500, 1000, 2000, 5000, 10000,15000)) +
  guides(color=guide_legend("Locality", override.aes = list(size=5))) +
  theme(axis.line = element_line(colour = "black"), panel.background = element_blank(),axis.text.y = element_text(size = 10)) 


#---------------------------------------------------------------------------------------------------------------------------------------------
# Outlier SNPs genotypes at dp 3 for heatmap analysis ----------------------------------------------------------------------------------------
setwd("/PATH/05_aDNA/02_results/uli_pca")
vcf6.fn <- "./01_angsd_ancient/ancient_1111_geno_q20_dp3_plink.vcf"
snpgdsVCF2GDS(vcf6.fn,"vcf6.gds", method="copy.num.of.ref")
snpgdsSummary("vcf6.gds")
vcf6.gdsfile <- snpgdsOpen("vcf6.gds")
vcf7.fn <- "./01_angsd_fresh/fresh_1111_geno_q20_dp3_plink.vcf"
snpgdsVCF2GDS(vcf7.fn,"vcf7.gds", method="copy.num.of.ref")
snpgdsSummary("vcf7.gds")
vcf7.gdsfile <- snpgdsOpen("vcf7.gds")

# Count REF and ALT alleles per sample to assess ref bias
geno6 <- snpgdsGetGeno(vcf6.gdsfile)
ref.counts <- 2 * rowSums(geno6 == 0, na.rm=TRUE) + rowSums(geno6 == 1, na.rm=TRUE)
alt.counts <- 2 * rowSums(geno6 == 2, na.rm=TRUE) + rowSums(geno6 == 1, na.rm=TRUE)
ref.ratio_anc <- ref.counts / (ref.counts + alt.counts)
sample.stats_anc <- data.frame(REF = ref.counts,ALT = alt.counts,REF_ratio = ref.ratio_anc)
geno7 <- snpgdsGetGeno(vcf7.gdsfile)
ref.counts <- 2 * rowSums(geno7 == 0, na.rm=TRUE) + rowSums(geno7 == 1, na.rm=TRUE)
alt.counts <- 2 * rowSums(geno7 == 2, na.rm=TRUE) + rowSums(geno7 == 1, na.rm=TRUE)
ref.ratio_fresh <- ref.counts / (ref.counts + alt.counts)
sample.stats_fresh <- data.frame(REF = ref.counts,ALT = alt.counts,REF_ratio = ref.ratio_fresh)

refbias <- cbind(rbind(sample.stats_anc,sample.stats_fresh), pop) %>%
  filter(!V2 %in% c("BDG002_TE","KCZ001_TE","DSZ007_TE","DVT017_TE"), V1 != "COR") %>%
  mutate(group= ifelse(V4 != "present", "ancient", "present"))%>%
  mutate(V1 = factor(V1, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","IT","RUS","B","ISR", "IRQ"))) 
ref4 <- ggplot(refbias, aes(x = factor(group), y = REF_ratio)) +
  geom_boxplot(alpha = 0.6) + geom_point(aes(color = V1), size = 2, alpha = 0.8) +
  geom_text_repel(data=refbias %>%filter(REF_ratio<0.1),aes(x = factor(group), y = REF_ratio,label=V2),size=2) +
  scale_color_manual(values = my_colors) + guides(color=guide_legend("Locality", override.aes = list(size=5))) +
  xlab("Sample Type") + ylab("REF/ALT Ratio") + labs(title="1111 outliers (dp3)")+theme_minimal()

#---------------------------------------------------------------------------------------------------------------------------------------------
# WGS ----------------------------------------------------------------------------------------------------------------------------------------
setwd("/PATH/05_aDNA")
setwd("01_angsd_wgs_rescaled_outgroup")
pop<-read.table("pop.txt", header=FALSE, sep="", col.names = c("pop","sample","age","type"))

vcf9.fn <- "wgs_rescaled_all_hapconsensus_maxmis_q20_dp3.vcf.gz"
snpgdsVCF2GDS(vcf9.fn,"vcf9.gds", method="copy.num.of.ref")
snpgdsSummary("vcf9.gds")
vcf9.gdsfile <- snpgdsOpen("vcf9.gds") #38,753,604 SNPs, 130 ind

# Count REF and ALT alleles per sample to assess ref bias
pop<-read.table("pop.txt", header=FALSE, sep="", col.names = c("pop","sample","age","type"))
geno2 <- snpgdsGetGeno(vcf9.gdsfile)
ref.counts <- 2 * rowSums(geno2 == 0, na.rm=TRUE) + rowSums(geno2 == 1, na.rm=TRUE)
alt.counts <- 2 * rowSums(geno2 == 2, na.rm=TRUE) + rowSums(geno2 == 1, na.rm=TRUE)
ref.ratio <- ref.counts / (ref.counts + alt.counts)
sample.stats <- data.frame(REF = ref.counts,ALT = alt.counts,REF_ratio = ref.ratio) 
refbias <- cbind(sample.stats,pop)%>%
  #filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE", pop != "COR") %>%
  filter(pop != "COR" & pop != "ori") %>%
  mutate(group= ifelse(type != "present", "ancient", "present"))%>%
  mutate(pop = factor(pop, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","IT","RUS","B","ISR", "IRQ"))) 
ref1 <- ggplot(refbias, aes(x = factor(group), y = REF_ratio)) +
  geom_boxplot(alpha = 0.6) + geom_point(aes(color = pop), size = 2, alpha = 0.8) + guides(color=guide_legend("Locality", override.aes = list(size=5))) +
  scale_color_manual(values = my_colors) + xlab("Sample Type") + ylab("REF/ALT Ratio") + labs(title="All neutral sites") + theme_minimal()

vcf.pca <- snpgdsPCA(vcf9.gdsfile, autosome.only=FALSE, missing.rate=NaN) ##use autosome.only=false to include all loci
names(vcf.pca)
pc.percent <- vcf.pca$varprop*100
head(round(pc.percent, 2))
print(pc.percent)
tab <- data.frame(sample.id = vcf.pca$sample.id,
                  EV1 = vcf.pca$eigenvect[,1], # the first eigenvector
                  EV2 = vcf.pca$eigenvect[,2], # the second eigenvector
                  EV3 = vcf.pca$eigenvect[,3], 
                  EV4 = vcf.pca$eigenvect[,4], 
                  EV5 = vcf.pca$eigenvect[,5], 
                  EV6 = vcf.pca$eigenvect[,6], 
                  EV7 = vcf.pca$eigenvect[,7], 
                  EV8 = vcf.pca$eigenvect[,8], 
                  EV9 = vcf.pca$eigenvect[,9], 
                  EV10 = vcf.pca$eigenvect[,10], 
                  stringsAsFactors = FALSE)

dat2 <- data.frame(pop, tab) %>%
  filter(sample != "DVT017_UDG" & sample != "NCP001_WGS" & sample != "HOC007_WGS" & sample != "DVT022_WGS" & sample != "HOC005_WGS"
         & sample != "KZR002_WGS"  & sample != "NCP002_WGS"  & sample != "DVT016_WGS" & sample != "WMP005_WGS"  & sample != "WMP006_WGS"  & sample != "KCZ003_WGS"
         & sample != "DVT014_WGS"  & sample != "TDN002_WGS"  & sample != "KCZ012_WGS" & sample != "VKP001_WGS") %>%
  mutate(pop = factor(pop, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","IT","RUS","B","ISR", "IRQ"))) %>%
  filter(pop != "ori")

#The number of biallelic unique SNPs: 177872 / 93431 (dp3)
s1 <- ggplot() +
  geom_point(data=dat2, aes(x=EV1, y=EV2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat2%>%filter(type=="1500-3000"| type=="more_than_10k"), aes(x=EV1, y=EV2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat2%>%filter(type=="5000-6000" | type=="more_than_16k"), aes(x=EV1, y=EV2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat2%>%filter(pop == "NL" | sample == "LSS003_TE" | age == 5330), aes(x=EV1, y=EV2, fill=pop, shape=type, alpha=type, size=type, color=type))  +
  geom_text(data=dat2, aes(x=EV1, y=EV2, label=sample)) +
  scale_shape_manual(values = c("100-1500" = 21, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_16k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 3, "1500-3000" = 3, "5000-6000"=3,"more_than_10k"=3,"more_than_16k"=3,"present"=2)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_16k"=1,"present"=1)) +
  scale_color_manual(values = c("100-1500" = "black", "1500-3000" = "black", "5000-6000"="black","more_than_10k"="black","more_than_16k"="black","present"="white")) +
  scale_fill_manual(values=my_colors ) +
  #guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(21, 22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 0.8), size=c(3, 3, 3, 3, 3, 2))), 
  #       fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=3)), color="none", size="none",alpha="none") +
  labs(x=paste("PC1 (",round(pc.percent[1],1),"%)",sep=""), y=paste("PC2 (",round(pc.percent[2],1),"%)",sep="")) +
  scale_x_reverse() + theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),
        axis.title.x = element_text(size=11),axis.title.y = element_text(size=11), axis.text.x=element_text (size=10),
        axis.text.y=element_text (size=10),legend.position = "right")

#---------------------------------------------------------------------------------------------------------------------------------------------
# PCA with EMU -------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
setwd("/PATH/05_aDNA")
setwd("01_angsd_enrichedallfresh_rescaled")
pop<-read.table("pop.txt", header=FALSE, sep="", col.names = c("pop","sample","age","type"))
eigvec <- read.table("enrichedallfresh_rescaled_neutral_seg72k_geno_maxmis_q20_plink.emu.eigvecs", header = FALSE) 
eigval <- scan("enrichedallfresh_rescaled_neutral_seg72k_geno_maxmis_q20_plink.emu.eigvals")

colnames(eigvec) <- c("sample", "pop2", paste0("PC", 1:(ncol(eigvec)-2)))
pca <- eigvec %>%
  left_join(pop, by="sample") %>%
  filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE", pop != "COR") %>%
  mutate(pop = factor(pop, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","IT","RUS","B","ISR", "IRQ")))

f1 <- ggplot() +
  geom_point(data=pca, aes(x=PC1, y=PC2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=pca%>%filter(type=="1500-3000"| type=="more_than_10k"), aes(x=PC1, y=PC2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=pca%>%filter(type=="5000-6000" | type=="more_than_16k"), aes(x=PC1, y=PC2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=pca%>%filter(pop == "NL" | sample == "LSS003_TE" | age == 5330), aes(x=PC1, y=PC2, fill=pop, shape=type, alpha=type, size=type, color=type))  +
  scale_shape_manual(values = c("100-1500" = 21, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_16k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 3, "1500-3000" = 3, "5000-6000"=3,"more_than_10k"=3,"more_than_16k"=3,"present"=2)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_16k"=1,"present"=1)) +
  scale_color_manual(values = c("100-1500" = "black", "1500-3000" = "black", "5000-6000"="black","more_than_10k"="black","more_than_16k"="black","present"="white")) +
  scale_fill_manual(values=my_colors ) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(21, 22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 0.8), size=c(3, 3, 3, 3, 3, 2))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=3)), color="none", size="none",alpha="none") +
  labs(x=paste("PC1 (",round(eigval[1]),"%)",sep=""), y=paste("PC2 (",round(eigval[2]),"%)",sep="")) +
  scale_x_reverse() +
  scale_y_reverse() +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text.x=element_text (size=10),
        axis.text.y=element_text (size=10),
        legend.position = "right") #6x6

f2 <- ggplot() +
  geom_point(data=pca, aes(x=PC2, y=PC3, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=pca%>%filter(type=="1500-3000"| type=="more_than_10k"), aes(x=PC2, y=PC3, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=pca%>%filter(type=="5000-6000" | type=="more_than_16k"), aes(x=PC2, y=PC3, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=pca%>%filter(pop == "NL" | sample == "LSS003_TE" | age == 5330), aes(x=PC2, y=PC3, fill=pop, shape=type, alpha=type, size=type, color=type))  +
  scale_shape_manual(values = c("100-1500" = 21, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_16k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 3, "1500-3000" = 3, "5000-6000"=3,"more_than_10k"=3,"more_than_16k"=3,"present"=2)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_16k"=1,"present"=1)) +
  scale_color_manual(values = c("100-1500" = "black", "1500-3000" = "black", "5000-6000"="black","more_than_10k"="black","more_than_16k"="black","present"="white")) +
  scale_fill_manual(values=my_colors ) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(21, 22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 0.8), size=c(3, 3, 3, 3, 3, 2))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=3)), color="none", size="none",alpha="none") +
  labs(x=paste("PC2 (",round(eigval[2]),"%)",sep=""), y=paste("PC3 (",round(eigval[3]),"%)",sep="")) +
  scale_x_reverse() +
  scale_y_reverse() +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text.x=element_text (size=10),
        axis.text.y=element_text (size=10),
        legend.position = "right") #6x6

f1 + f2  + plot_layout(width = c(1,1)) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = c("A", "B", "C", "D")) & theme(plot.margin = unit(c(1, 1, 1, 1), "pt")) 

