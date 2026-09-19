#install.packages("smartsnp")
library("smartsnp")
library("dplyr")
library("ggrepel")
library("ggplot2")
library("ggshadow")
library(cowplot)

#colorBlind_anc  <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#E69F00","#882255", "#44AA99", "#661100")
#safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#E69F00", "#AA4499", "#6699CC", "#999933", "#882255", "#44AA99", "#661100", "#888888", "#000000", "#E69F00", "#56B4E9")
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

setwd("/PATH/05_aDNA/01_angsd_enrichedallfresh_rescaled")
#dat2 <- read.table(gzfile("enrichedallfresh_rescaled_outlier_hapconsensus_maxmis_q20.haplo.gz"), header=TRUE, na="N")
#dat2 <- read.table(gzfile("enrichedallfresh_rescaled_outlier_1111_hapconsensus_maxmis_q20.haplo.gz"), header=TRUE, na="N")
dat2 <- read.table(gzfile("enrichedallfresh_rescaled_outlier_1111_hapconsensus_maxmis_q20_dp3.haplo.gz"), header=TRUE, na="N")

#basic check:0 are different from major allele, 2 are same as major allele, NA as 9
qc2 <- colSums(!is.na(dat2))/(nrow(dat2))
df2 <- dat2$major == dat2[4:ncol(dat2)]
df2 <- df2*1*2 
df2[is.na(df2)] <- 9
dim(df2) 
#write.table(df2, file = "enrichedallfresh_rescaled_outlier_1111_hapconsensus_maxmis_q20_dp3_smartsnp", col.names = FALSE, row.names = FALSE)

#check ref bias in pseudohaploid
pop <- read.table("pop.txt")
ref.counts <- colSums(df2 == 0, na.rm=TRUE)
alt.counts <- colSums(df2 == 2, na.rm=TRUE)
ref.ratio <- ref.counts / (ref.counts + alt.counts)
sample.stats <- data.frame(REF = ref.counts,ALT = alt.counts,REF_ratio = ref.ratio)
refbias <- cbind(sample.stats,pop)%>%
  filter(!V2 %in% c("BDG002_TE","KCZ001_TE","DSZ007_TE","DVT017_TE")) %>%
  filter(V1 != "ori" , V1 != "COR") %>%
  mutate(group= ifelse(V4 != "present", "ancient", "present"))%>%
  mutate(V1 = factor(V1, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","IT","RUS","B","ISR", "IRQ"))) 
pre <- refbias %>% filter(group == "present")
ref <- ggplot(refbias, aes(x = factor(group), y = REF_ratio)) +guides(color=guide_legend("Locality", override.aes = list(size=5))) +
  geom_boxplot(alpha = 0.6) + geom_point(aes(color = V1), size = 2, alpha = 0.8) +
  #geom_boxplot(data=refbias,aes(x = factor(group), y = REF_ratio), alpha = 0.6) + geom_point(data=refbias,aes(x = factor(group), y = REF_ratio, color = V1), size = 2, alpha = 0.8) +
  geom_text_repel(data=refbias %>%filter(REF_ratio<quantile(pre$REF_ratio, 0) |  REF_ratio>quantile(pre$REF_ratio, 1)),aes(x = factor(group), y = REF_ratio,label=V2),size=2) +
  scale_color_manual(values = my_colors) + xlab("Sample Type") + ylab("REF/ALT Ratio") + labs(title="1111 outliers") + theme_minimal() #12x6

#Read dataSNP and do some filtering
dataSNP2 <- read.table("enrichedallfresh_rescaled_outlier_1111_hapconsensus_maxmis_q20_dp3_smartsnp",sep="")
pop <- read.table("pop.txt")
summary2 <- data.frame(row=c(1:164), country=pop$V1, sample=pop$V2, age=pop$V3, period=pop$V4, miss=(1-qc2[4:length(qc2)])*100)
pass2 <- filter(summary2, miss<100) %>% 
  filter(sample != "HOC007_TE") %>% #this sample gives 0 SNPs in dp3 and cannto be projected
  filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE", country != "COR")
selected2 <- dataSNP2[pass2$row]
#write.table(selected2, file = "enrichedallfresh_rescaled_outlier_1111_hapconsensus_maxmis_q20_dp3_miss100_smartsnp", col.names = FALSE, row.names = FALSE)

### pseudohaploid without projection remaining 200 SNPs for miss40 and 100 SNPs for miss100 as all missing are removed

#sample projection allows miss100
countmodern <- pass2 %>% filter(period == "present")
my_groups <- c(rep("A", nrow(pass2)-nrow(countmodern)),rep("M",nrow(countmodern)))
my_ancient <- c(1:(nrow(pass2)-nrow(countmodern)))
pcaR <- smart_pca(snp_data = "enrichedallfresh_rescaled_outlier_1111_hapconsensus_maxmis_q20_dp3_miss100_smartsnp", 
                  sample_group = my_groups, sample_project = my_ancient, pc_axes=100, missing_impute = "remove")
pcaR_eigen <- pcaR$pca.eigenvalues; dim(pcaR_eigen)
pcaR_load <- pcaR$pca.snp_loadings; dim(pcaR_load) 
pcaR_coord <- pcaR$pca.sample_coordinates; dim(pcaR_coord)
data <- pcaR$pca.sample_coordinates
data$sample <- pass2$sample
data$country <- pass2$country
data$age <- pass2$age
data$miss <- pass2$miss
data$period <- pass2$period
df1 <- data %>% mutate(country = factor(country, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","RUS","IT","B","ISR","IRQ")))

# Normalised eigenvalues
raw_eigenvals <- pcaR$pca.eigenvalues["observed eigenvalues", ]
pc1 <- round(raw_eigenvals[1] / sum(raw_eigenvals) * 100,1)
pc2 <- round(raw_eigenvals[2] / sum(raw_eigenvals) * 100,1)

### Plot pseudohaploid with projection: 29096 SNPs remaining / 801 SNPs for 1111set / 287 SNPs  dp3
proj1 <-ggplot() +
  geom_point(data=df1, aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df1%>%filter(period %in% c("100-1500","undated_100","undated_2k","undated_5k")), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df1%>%filter(period=="1500-3000"| period=="more_than_10k"), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df1%>%filter(period=="5000-6000" | period=="more_than_16k"), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df1%>%filter(country == "NL" | sample == "LSS003_TE" | age == 5330), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period))  +
  #geom_text_repel(data=subset(df1, sample == "LSS003_TE"),aes(PC1,PC2,label="F_LSS003"), size=2.5,box.padding = 1,point.padding = 1,color = "black") +
  scale_shape_manual(values = c("100-1500" = 21, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_16k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 3, "1500-3000" = 3, "5000-6000"=3,"more_than_10k"=3,"more_than_16k"=3,"present"=2)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_16k"=1,"present"=1)) +
  scale_color_manual(values = c("100-1500" = "black", "1500-3000" = "black", "5000-6000"="black","more_than_10k"="black","more_than_16k"="black","present"="white")) +
  scale_fill_manual(values=my_colors ) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(21, 22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 0.8), size=c(3, 3, 3, 3, 3, 2))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=3)), color="none", size="none",alpha="none") +
  labs(x=paste("PC1 (",pc1,"%)"), y=paste("PC2 (",pc2,"%)")) +
  scale_x_reverse() +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x=element_text (size=8),
        axis.text.y=element_text (size=8),
        legend.position = "right") #6x6

#---------------------------------------------------------------------------------------------------------------------------------------------
# All Neutral --------------------------------------------------------------------------------------------------------------------------------
setwd("/PATH/05_aDNA/01_angsd_enrichedallfresh_rescaled")
#dat2 <- read.table(gzfile("enrichedallfresh_rescaled_neutral_hapconsensus_maxmis_q20.haplo.gz"), header=TRUE, na="N")
dat2 <- read.table(gzfile("enrichedallfresh_rescaled_neutral_hapconsensus_maxmis_q20_dp3.haplo.gz"), header=TRUE, na="N")

#basic check
qc2 <- colSums(!is.na(dat2))/(nrow(dat2))
df2 <- dat2$major == dat2[4:ncol(dat2)]
df2 <- df2*1*2 
df2[is.na(df2)] <- 9
dim(df2) 
#write.table(df2, file = "enrichedallfresh_rescaled_neutral_hapconsensus_maxmis_q20_dp3_smartsnp", col.names = FALSE, row.names = FALSE)

#check ref bias in pseudohaploid
pop <- read.table("pop.txt")
ref.counts <- colSums(df2 == 0, na.rm=TRUE)
alt.counts <- colSums(df2 == 2, na.rm=TRUE)
ref.ratio <- ref.counts / (ref.counts + alt.counts)
sample.stats <- data.frame(REF = ref.counts,ALT = alt.counts,REF_ratio = ref.ratio)
refbias <- cbind(sample.stats,pop)%>%
  filter(!V2 %in% c("BDG002_TE","KCZ001_TE","DSZ007_TE","DVT017_TE")) %>%
  filter(V1 != "ori" , V1 != "COR") %>%
  mutate(group= ifelse(V4 != "present", "ancient", "present"))%>%
  mutate(V1 = factor(V1, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","IT","RUS","B","ISR", "IRQ"))) 
pre <- refbias %>% filter(group == "present")
ref1 <- ggplot(refbias, aes(x = factor(group), y = REF_ratio)) +guides(color=guide_legend("Locality", override.aes = list(size=5))) +
  geom_boxplot(alpha = 0.6) + geom_point(aes(color = V1), size = 2, alpha = 0.8) +
  #geom_boxplot(data=refbias,aes(x = factor(group), y = REF_ratio), alpha = 0.6) + geom_point(data=refbias,aes(x = factor(group), y = REF_ratio, color = V1), size = 2, alpha = 0.8) +
  geom_text_repel(data=refbias %>%filter(REF_ratio<quantile(pre$REF_ratio, 0) |  REF_ratio>quantile(pre$REF_ratio, 1)),aes(x = factor(group), y = REF_ratio,label=V2),size=2) +
  scale_color_manual(values = my_colors) + xlab("Sample Type") + ylab("REF/ALT Ratio") + theme_minimal() + labs(title="All neutral SNPs") #+ labs(title="Outgroup-acs & transversion SNPs") #12x6

#Read dataSNP and do some filtering
dataSNP2 <- read.table("enrichedallfresh_rescaled_neutral_seg72k_hapconsensus_maxmis_q20_dp3_smartsnp",sep="")
summary2 <- data.frame(row=c(1:164), country=pop$V1, sample=pop$V2, age=pop$V3, period=pop$V4, miss=(1-qc2[4:length(qc2)])*100)
pass2 <- filter(summary2, miss<100) %>% 
  filter(sample != "HOC007_TE") %>% #problematic in dp3
  filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE" & country != "COR")
selected2 <- dataSNP2[pass2$row]
#write.table(selected2, file = "enrichedallfresh_rescaled_neutral_hapconsensus_maxmis_q20_dp3_miss100_smartsnp", col.names = FALSE, row.names = FALSE)

#sample projection
countmodern <- pass2 %>% filter(period == "present")
my_groups <- c(rep("A", nrow(pass2)-nrow(countmodern)),rep("M",nrow(countmodern)))
my_ancient <- c(1:(nrow(pass2)-nrow(countmodern)))
pcaR <- smart_pca(snp_data = "enrichedallfresh_rescaled_neutral_seg72k_hapconsensus_maxmis_q20_dp3_miss100_smartsnp", 
                  sample_group = my_groups, sample_project = my_ancient, pc_axes=100, missing_impute = "remove")
pcaR_eigen <- pcaR$pca.eigenvalues; dim(pcaR_eigen)
pcaR_load <- pcaR$pca.snp_loadings; dim(pcaR_load) 
pcaR_coord <- pcaR$pca.sample_coordinates; dim(pcaR_coord)
data <- pcaR$pca.sample_coordinates
data$sample <- pass2$sample
data$country <- pass2$country
data$age <- pass2$age
data$miss <- pass2$miss
data$period <- pass2$period
df2 <- data %>% mutate(country = factor(country, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","RUS","IT","B","ISR","IRQ")))

# Normalised eigenvalues
raw_eigenvals <- pcaR$pca.eigenvalues["observed eigenvalues", ]
pc1 <- round(raw_eigenvals[1] / sum(raw_eigenvals) * 100,1)
pc2 <- round(raw_eigenvals[2] / sum(raw_eigenvals) * 100,1)

### Plot pseudohaploid with projection: 129657 SNPs (dp0) / 37701 (dp3) 
proj2 <- ggplot() +
  geom_point(data=df2, aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df2%>%filter(period %in% c("100-1500","undated_100","undated_2k","undated_5k")), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df2%>%filter(period=="1500-3000"| period=="more_than_10k"), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df2%>%filter(period=="5000-6000" | period=="more_than_16k"), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df2%>%filter(country == "NL" | sample == "LSS003_TE" | age == 5330), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period))  +
  #geom_text_repel(data=subset(df2, sample == "LSS003_TE"),aes(PC1,PC2,label="F_LSS003"), size=2.5,box.padding = 2,point.padding = 1,color = "black") +
  scale_shape_manual(values = c("100-1500" = 21, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_16k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 3, "1500-3000" = 3, "5000-6000"=3,"more_than_10k"=3,"more_than_16k"=3,"present"=2)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_16k"=1,"present"=1)) +
  scale_color_manual(values = c("100-1500" = "black", "1500-3000" = "black", "5000-6000"="black","more_than_10k"="black","more_than_16k"="black","present"="white")) +
  scale_fill_manual(values=my_colors ) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(21, 22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 0.8), size=c(3, 3, 3, 3, 3, 2))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=3)), color="none", size="none",alpha="none") +
  labs(x=paste("PC1 (",pc1,"%)"), y=paste("PC2 (",pc2,"%)")) +
  #scale_y_break(c(-450, -50)) +
  scale_y_break(c(-250, -50)) +
  scale_x_reverse() +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size=10,margin = margin(t =0)),
        axis.title.y = element_text(size=10,margin = margin(t =0)),
        axis.text.x=element_text (size=8),
        axis.text.y=element_text (size=8),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.position = "right") #6x6

#---------------------------------------------------------------------------------------------------------------------------------------------
# Neutral outgroup-ascertained ------------------------------------------------------------------------------------------------------------------
setwd("/PATH/05_aDNA/01_angsd_enrichedallfresh_rescaled")
#dat2 <- read.table(gzfile("enrichedallfresh_rescaled_neutral_seg72k_hapconsensus_maxmis_q20.haplo.gz"), header=TRUE, na="N")
dat2 <- read.table(gzfile("enrichedallfresh_rescaled_neutral_seg72k_hapconsensus_maxmis_q20_dp3.haplo.gz"), header=TRUE, na="N")

#basic check
qc2 <- colSums(!is.na(dat2))/(nrow(dat2))
df2 <- dat2$major == dat2[4:ncol(dat2)]
df2 <- df2*1*2 
df2[is.na(df2)] <- 9
dim(df2) 
#write.table(df2, file = "enrichedallfresh_rescaled_neutral_seg72k_hapconsensus_maxmis_q20_dp3_smartsnp", col.names = FALSE, row.names = FALSE)

#check ref bias in pseudohaploid
pop <- read.table("pop.txt")
ref.counts <- colSums(df2 == 0, na.rm=TRUE)
alt.counts <- colSums(df2 == 2, na.rm=TRUE)
ref.ratio <- ref.counts / (ref.counts + alt.counts)
sample.stats <- data.frame(REF = ref.counts,ALT = alt.counts,REF_ratio = ref.ratio)
refbias <- cbind(sample.stats,pop)%>%
  filter(!V2 %in% c("BDG002_TE","KCZ001_TE","DSZ007_TE","DVT017_TE")) %>%
  filter(V1 != "ori" , V1 != "COR") %>%
  mutate(group= ifelse(V4 != "present", "ancient", "present"))%>%
  mutate(V1 = factor(V1, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","IT","RUS","B","ISR", "IRQ"))) 
pre <- refbias %>% filter(group == "present")
ref2 <- ggplot(refbias, aes(x = factor(group), y = REF_ratio)) +guides(color=guide_legend("Locality", override.aes = list(size=5))) +
  geom_boxplot(alpha = 0.6) + geom_point(aes(color = V1), size = 2, alpha = 0.8) +
  #geom_boxplot(data=refbias,aes(x = factor(group), y = REF_ratio), alpha = 0.6) + geom_point(data=refbias,aes(x = factor(group), y = REF_ratio, color = V1), size = 2, alpha = 0.8) +
  geom_text_repel(data=refbias %>%filter(REF_ratio<quantile(pre$REF_ratio, 0) |  REF_ratio>quantile(pre$REF_ratio, 1)),aes(x = factor(group), y = REF_ratio,label=V2),size=2) +
  scale_color_manual(values = my_colors) + xlab("Sample Type") + ylab("REF/ALT Ratio") + theme_minimal() + labs(title="Outgroup-asc & transversion") #+ labs(title="Outgroup-acs & transversion SNPs") #12x6
refratio <- ref + ref1 + ref2  + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = c("A", "B", "C")) & theme(plot.margin = unit(c(1, 1, 1, 1), "pt")) #8x6

#Read dataSNP and do some filtering
dataSNP2 <- read.table("enrichedallfresh_rescaled_neutral_seg72k_hapconsensus_maxmis_q20_dp3_smartsnp",sep="")
summary2 <- data.frame(row=c(1:164), country=pop$V1, sample=pop$V2, age=pop$V3, period=pop$V4, miss=(1-qc2[4:length(qc2)])*100)
pass2 <- filter(summary2, miss<100) %>% 
  filter(sample != "HOC007_TE") %>% #problematic in dp3
  filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE" & country != "COR")
selected2 <- dataSNP2[pass2$row]
#write.table(selected2, file = "enrichedallfresh_rescaled_neutral_seg72k_hapconsensus_maxmis_q20_dp3_miss100_smartsnp", col.names = FALSE, row.names = FALSE)

#sample projection
countmodern <- pass2 %>% filter(period == "present")
my_groups <- c(rep("A", nrow(pass2)-nrow(countmodern)),rep("M",nrow(countmodern)))
my_ancient <- c(1:(nrow(pass2)-nrow(countmodern)))
pcaR <- smart_pca(snp_data = "enrichedallfresh_rescaled_neutral_seg72k_hapconsensus_maxmis_q20_dp3_miss100_smartsnp", 
                  sample_group = my_groups, sample_project = my_ancient, pc_axes=100, missing_impute = "remove")
pcaR_eigen <- pcaR$pca.eigenvalues; dim(pcaR_eigen)
pcaR_load <- pcaR$pca.snp_loadings; dim(pcaR_load) 
pcaR_coord <- pcaR$pca.sample_coordinates; dim(pcaR_coord)
data <- pcaR$pca.sample_coordinates
data$sample <- pass2$sample
data$country <- pass2$country
data$age <- pass2$age
data$miss <- pass2$miss
data$period <- pass2$period
df2 <- data %>% mutate(country = factor(country, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","RUS","IT","B","ISR","IRQ")))
# Normalised eigenvalues
raw_eigenvals <- pcaR$pca.eigenvalues["observed eigenvalues", ]
pc1 <- round(raw_eigenvals[1] / sum(raw_eigenvals) * 100,1)
pc2 <- round(raw_eigenvals[2] / sum(raw_eigenvals) * 100,1)

### Plot pseudohaploid with projection for seg72k: 55607 (dp0) / 15487 (dp3)
proj2 <- ggplot() +
  geom_point(data=df2, aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df2%>%filter(period %in% c("100-1500","undated_100","undated_2k","undated_5k")), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df2%>%filter(period=="1500-3000"| period=="more_than_10k"), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df2%>%filter(period=="5000-6000" | period=="more_than_16k"), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df2%>%filter(country == "NL" | sample == "LSS003_TE" | age == 5330), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period))  +
  #geom_text_repel(data=subset(df2, sample == "LSS003_TE"),aes(PC1,PC2,label="F_LSS003"), size=2.5,box.padding = 2,point.padding = 1,color = "black") +
  scale_shape_manual(values = c("100-1500" = 21, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_16k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 3, "1500-3000" = 3, "5000-6000"=3,"more_than_10k"=3,"more_than_16k"=3,"present"=2)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_16k"=1,"present"=1)) +
  scale_color_manual(values = c("100-1500" = "black", "1500-3000" = "black", "5000-6000"="black","more_than_10k"="black","more_than_16k"="black","present"="white")) +
  scale_fill_manual(values=my_colors ) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(21, 22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 0.8), size=c(3, 3, 3, 3, 3, 2))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=3)), color="none", size="none",alpha="none") +
  labs(x=paste("PC1 (",pc1,"%)"), y=paste("PC2 (",pc2,"%)")) +
  scale_y_break(c(-150, -20)) +
  scale_x_reverse() +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size=10,margin = margin(t =0)),
        axis.title.y = element_text(size=10,margin = margin(t =0)),
        axis.text.x=element_text (size=8),
        axis.text.y=element_text (size=8),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        legend.position = "right") #6x6

#---------------------------------------------------------------------------------------------------------------------------------------------
# Neutral with outgroup -----------------------------------------------------------------------------------------------------------------------------------------
setwd("/PATH/05_aDNA/01_angsd_enrichedallfresh_rescaled_outgroupAM")
dat <- read.table(gzfile("enrichedallfresh_rescaled_outgroupAM_neutral_hapconsensus_maxmis_q20.haplo.gz"), header=TRUE, na="N") #4M requires more MEM

#basic check
qc <- colSums(!is.na(dat))/(nrow(dat))
df <- dat$major == dat[4:ncol(dat)]
df <- df*1*2 
df[is.na(df)] <- 9
dim(df) 
write.table(df, file = "enrichedallfresh_rescaled_outgroupAM_neutral_hapconsensus_maxmis_q20_smartsnp", col.names = FALSE, row.names = FALSE)
pop <- read.table("pop.txt") %>% filter(!V2 %in% c("BDG002_TE","KCZ001_TE","DSZ007_TE","DVT017_TE"))

#Read dataSNP and do some filtering
dataSNP <- read.table("enrichedallfresh_rescaled_outgroupAM_neutral_hapconsensus_maxmis_q20_smartsnp",sep="")
summary <- data.frame(row=c(1:ncol(df)), country=pop$V1, sample=pop$V2, age=pop$V3, period=pop$V4, miss=(1-qc[4:length(qc)])*100)


#---------------------------------------------------------------------------------------------------------------------------------------------
# WGS-----------------------------------------------------------------------------------------------------------------------------------------
setwd("/PATH/05_aDNA/01_angsd_wgs_rescaled_outgroup")
dat2 <- read.table(gzfile("wgs_rescaled_outgroup_all_hapconsensus_maxmis_q20.haplo.gz"), header=TRUE, na="N") #1835510 SNPs (dp3) / 21017638 (dp0)
dat2_nooutgroup <- dat2[ , 1:(ncol(dat2) - 14)]
dat2_nooutgroup <- dat2_nooutgroup[apply(dat2_nooutgroup[, 4:ncol(dat2_nooutgroup)], 1, function(x) {
  alleles <- unique(na.omit(x))
  length(alleles) > 1}), ] #keep variant sites only:787328

#basic check
qc2 <- colSums(!is.na(dat2_nooutgroup))/(nrow(dat2_nooutgroup))
df2 <- dat2_nooutgroup$major == dat2_nooutgroup[4:ncol(dat2_nooutgroup)]
df2 <- df2*1*2 
df2[is.na(df2)] <- 9
dim(df2) 
#write.table(df2, file = "wgs_rescaled_outgroup_all_hapconsensus_maxmis_q20_dp3_smartsnp", col.names = FALSE, row.names = FALSE)
write.table(df2, file = "wgs_rescaled_outgroup_all_hapconsensus_maxmis_q20_noOUTGROUP_smartsnp", col.names = FALSE, row.names = FALSE)
pop <- read.table("pop.txt") %>% 
  filter(V1 != "AM" &  V1 != "MON" &  V1 != "JACKDAW")

#Read dataSNP and do some filtering
dataSNP2 <- read.table("wgs_rescaled_outgroup_all_hapconsensus_maxmis_q20_noOUTGROUP_smartsnp",sep="")
summary2 <- data.frame(row=c(1:ncol(df2)), country=pop$V1, sample=pop$V2, age=pop$V3, period=pop$V4, miss=(1-qc2[4:length(qc2)])*100)
pass2 <- filter(summary2, miss<100) %>% 
  filter(country != "AM" & country != "MON" & country != "JACKDAW") #DVT017_UDG looks different even though missingness was lower than others (~50%)
selected2 <- dataSNP2[pass2$row]
write.table(selected2, file = "wgs_rescaled_outgroup_all_hapconsensus_maxmis_q20_miss100_noOUTGROUP_smartsnp", col.names = FALSE, row.names = FALSE)
write.table(summary2, file = "wgs_rescaled_outgroup_all_hapconsensus_maxmis_q20_miss100_noOUTGROUP.imiss", col.names = TRUE, row.names = TRUE)

#sample projection
countmodern <- pass2 %>% filter(period == "present")
my_groups <- c(rep("A", nrow(pass2)-nrow(countmodern)),rep("M",nrow(countmodern)))
my_ancient <- c(1:(nrow(pass2)-nrow(countmodern)))
pcaR <- smart_pca(snp_data = "wgs_rescaled_outgroup_all_hapconsensus_maxmis_q20_miss100_noOUTGROUP_smartsnp", sample_group = my_groups, sample_project = my_ancient, pc_axes=nrow(countmodern))
pcaR_eigen <- pcaR$pca.eigenvalues; dim(pcaR_eigen) #58358 missing values imputed (dp3) / 218730 missing values imputed (dp0)
pcaR_load <- pcaR$pca.snp_loadings; dim(pcaR_load) 
pcaR_coord <- pcaR$pca.sample_coordinates; dim(pcaR_coord)
data <- pcaR$pca.sample_coordinates
data$sample <- pass2$sample
data$country <- pass2$country
data$age <- pass2$age
data$miss <- pass2$miss
data$period <- pass2$period
df2 <- data %>% mutate(country = factor(country, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","RUS","IT","B","ISR","IRQ")))

# Normalised eigenvalues
raw_eigenvals <- pcaR$pca.eigenvalues["observed eigenvalues", ]
pc1 <- round(raw_eigenvals[1] / sum(raw_eigenvals) * 100,1)
pc2 <- round(raw_eigenvals[2] / sum(raw_eigenvals) * 100,1)

### Plot pseudohaploid with projection: 599482 SNPs from 51 modern samples, ancient samples projected on it (original 1835510 SNPs with outgroup but high missingness in JACKDAW)
ggplot() +
  geom_point(data=df2, aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df2%>%filter(period %in% c("100-1500","undated_100","undated_2k","undated_5k")), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df2%>%filter(period=="1500-3000"| period=="more_than_10k"), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_point(data=df2%>%filter(period=="5000-6000" | period=="more_than_16k"), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_text_repel(data=df2,aes(PC1,PC2,label=sample), size=2.5) +
  scale_shape_manual(values = c("100-1500" = 21, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_16k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 3, "1500-3000" = 3, "5000-6000"=3,"more_than_10k"=3,"more_than_16k"=3,"present"=2)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_16k"=1,"present"=1)) +
  scale_color_manual(values = c("100-1500" = "black", "1500-3000" = "black", "5000-6000"="black","more_than_10k"="black","more_than_16k"="black","present"="white")) +
  scale_fill_manual(values=my_colors ) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(21, 22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 0.8), size=c(3, 3, 3, 3, 3, 2))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=3)), color="none", size="none",alpha="none") +
  labs(x=paste("PC1 (",pc1,"%)"), y=paste("PC2 (",pc2,"%)")) +
  scale_x_reverse() + theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),axis.text.x=element_text (size=10),axis.text.y=element_text (size=10),legend.position = "right") #6x6

### Plot missingness
ggplot() +
  geom_point(data=summary2, aes(x=age, y=miss, fill=country, shape=period, alpha=period, size=period, color=period)) +
  geom_text_repel(data=summary2,aes(x=age, y=miss,label=sample), size=2.5) +
  scale_shape_manual(values = c("100-1500" = 21, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_16k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 3, "1500-3000" = 3, "5000-6000"=3,"more_than_10k"=3,"more_than_16k"=3,"present"=2)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_16k"=1,"present"=1)) +
  scale_color_manual(values = c("100-1500" = "black", "1500-3000" = "black", "5000-6000"="black","more_than_10k"="black","more_than_16k"="black","present"="white")) +
  scale_fill_manual(values=my_colors ) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(21, 22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 0.8), size=c(3, 3, 3, 3, 3, 2))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=3)), color="none", size="none",alpha="none") +
  labs(x="Age (years ago)", y="Missing (%)") + theme_bw(base_size=12) + coord_cartesian(ylim = c(0, 100))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),axis.text.x=element_text (size=10),axis.text.y=element_text (size=10),legend.position = "right") #6x6

  
