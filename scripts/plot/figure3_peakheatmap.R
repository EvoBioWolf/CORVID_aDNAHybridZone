library(openxlsx, lib.loc = "/dss/dsshome1/lxc0E/di67kah/R")
library(gtools, lib.loc = "/dss/dsshome1/lxc0E/di67kah/R")
library(base)
library(gdsfmt)
library(SNPRelate, lib.loc = "/dss/dsshome1/lxc0E/di67kah/R")
library("plotly",lib="/dss/dsshome1/lxc0E/di67kah/R")
library(stats)
library(tidyverse)
library(forcats)
library(patchwork,lib="/dss/dsshome1/lxc0E/di67kah/R")
library(vcfR, lib.loc = "/dss/dsshome1/lxc0E/di67kah/R")
library(memuse, lib.loc = "/dss/dsshome1/lxc0E/di67kah/R")

setwd("/PATH/05_aDNA/02_results/uli_pca")

# (0) Check ancient samples
# (1) Combine all samples
# (2) Filtering and ancestral state
# (3) Create PLINK input
PLINK <- FALSE
# (4) Read data into GDS DB, polarize SNPs and perform PCA
# (5) FST landscape
# (6) Heatmap
# (7) PCA
PCA <- TRUE
# (8) Plumage map
# ------------------------------------------------------------------------------------------------------------------------------------

# (0) FIRST CHECK PCA OF ANCIENT SAMPLES WITH 1111 SNPS ------------------------------------------------------------------------------
path <- "/PATH/05_aDNA/02_results/uli_pca/"
setwd(paste(path,"01_angsd_ancient",sep=""))
vcf.fn <- paste(path,"01_angsd_ancient/ancient_1111_geno_q20_dp3_plink.vcf",sep="")
sampleinfo <- read.xlsx(paste(path,"01_angsd_ancient/data_snps_individuals_ancient.xlsx",sep=""),colNames=TRUE)
sampleinfo <- mutate(sampleinfo,sample.id=paste0(sample_code,"_TE_",country_code))
sample_filter <- sampleinfo %>% dplyr::select(sample_code) %>%
  filter(sample_code !="DSZ007_TE_PL", sample_code !="BDG002_TE_GB", sample_code !="KCZ001_TE_PL", sample_code !="DVT017_TE_B")
snpgdsVCF2GDS(vcf.fn,"vcf.gds", method="copy.num.of.ref")
snpgdsSummary("vcf.gds")
vcf.gdsfile <- snpgdsOpen("vcf.gds")
vcf.pca <- snpgdsPCA(vcf.gdsfile, autosome.only=FALSE, missing.rate=NaN,sample.id = sample_filter$sample.id) ##use autosome.only=false to include all loci
unlink("vcf1.gds", force=TRUE)

names(vcf.pca)
pc.percent <- vcf.pca$varprop*100
head(round(pc.percent, 2))
print(pc.percent)
tab <- data.frame(sample.id = vcf.pca$sample.id,EV1 = vcf.pca$eigenvect[,1],EV2 = vcf.pca$eigenvect[,2],EV3 = vcf.pca$eigenvect[,3], 
                  EV4 = vcf.pca$eigenvect[,4],EV5 = vcf.pca$eigenvect[,5],  EV6 = vcf.pca$eigenvect[,6], EV7 = vcf.pca$eigenvect[,7], 
                  EV8 = vcf.pca$eigenvect[,8], EV9 = vcf.pca$eigenvect[,9], EV10 = vcf.pca$eigenvect[,10], stringsAsFactors = FALSE)

snpplot <- ggplot(tab) +
  geom_point(aes(x=EV1, y=EV2), size=2) +
  geom_text(aes(x=EV1, y=EV2), label=tab$sample, check_overlap=TRUE, size = 2, angle = 45, vjust = -1) +
  scale_shape_manual(name="Age", values=c(15, 17, 18, 16, 13, 8, 9, 10, 11, 12, 13, 5, 15, 16, 17)) +
  #scale_color_manual(name="Population", values=safe_colorblind_palette) +
  labs(x=paste("PC1 ",round(pc.percent[1],2),"%",sep=""), y=paste("PC2 ",round(pc.percent[2],2),"%",sep="")) +
  guides(color=guide_legend("Population"),fill=guide_legend("Population")) +
  ggtitle("") +
  theme_bw(base_size=12) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x=element_text (size=12),axis.text.y=element_text (size=12),
        legend.text=element_text(size=12), legend.title=element_text(size=12)) 
# ------------------------------------------------------------------------------------------------------------------------------------

# (1) COMBINE ANCIENT SAMPLES & NON-CONTACT FRESH SAMPLES WITH ULI'S SAMPLES ---------------------------------------------------------
anc_geno <- read_delim(paste(path,"01_angsd_ancient/ancient_1111_geno_q20_dp3.geno.transposed.txt",sep=""), delim="\t",col_names=TRUE) 
sampleinfo <- read.xlsx(paste(path,"01_angsd_ancient/data_snps_individuals_ancient.xlsx",sep=""),colNames=TRUE)
ancdat <- cbind(sampleinfo, anc_geno) %>%
  filter(sample_code !="DSZ007", sample_code !="BDG002", sample_code !="KCZ001", sample_code !="DVT017") #these samples appeared weird
fresh_geno <- read_delim(paste(path,"01_angsd_fresh/fresh_1111_geno_q20_dp3.geno.transposed.txt",sep=""), delim="\t",col_names=TRUE) 
sampleinfo2 <- read.xlsx(paste(path,"01_angsd_fresh/data_snps_individuals_fresh.xlsx",sep=""),colNames=TRUE)
freshdat <- cbind(sampleinfo2, fresh_geno) %>%
  filter(country_code != "IT", country_code != "PL", country_code != "D", !sample_code %in% c("S02","S04","S06","S07","S09")) #these samples also present in Uli's file but note IT11gg.IT03 and IT13gg.IT05 combined by Uli 

path <- "/PATH/05_aDNA/02_results/uli_pca/"
dat <- read.xlsx(xlsxFile=paste(path,"data_snps_individuals_181130.xlsx",sep=""), sheet="data_individuals_genotypes", colNames=TRUE) 
all_columns <- union(names(ancdat), names(dat))
ancdat[setdiff(all_columns, names(ancdat))] <- "N/N" #these missing scaffolds in anc calling
freshdat[setdiff(all_columns, names(freshdat))] <- "N/N" #these missing scaffolds in add fresh calling
setdiff(all_columns, names(dat)) #"scaffold_7:26556260"
ancdat <- ancdat %>% dplyr::select(-"scaffold_7:26556260")
freshdat <- freshdat %>% dplyr::select(-"scaffold_7:26556260")
#check again for missing columns
all_columns <- union(names(freshdat), names(dat))
setdiff(all_columns, names(dat))
newdat <- rbind(ancdat,dat)
newdat <- rbind(newdat,freshdat)
# ------------------------------------------------------------------------------------------------------------------------------------

# (2) SOME FILTERING ------------------------------------------------------------------------------------------------------------------
newdat <- subset(newdat, !(sample_code %in% c("428V","462V","HF13584","HF73810","HF73830","HF73837","HF73809")))
newdat[newdat=="NA/NA?"] <- "0/0"
newdat[newdat=="N/N"] <- "0/0"
newdat$nestID[which(is.na(dat$nestID))] <- 1:length(which(is.na(dat$nestID)))
snp <- read.xlsx(xlsxFile=paste(path,"data_snps_individuals_181130.xlsx",sep=""), sheet="data_snps", colNames=TRUE)
fst <- read.xlsx(xlsxFile=paste(path,"data_snps_individuals_181130.xlsx",sep=""), sheet="data_FST", colNames=TRUE)
ancestral_uli <- read.xlsx(xlsxFile=paste(path,"data_snps_individuals_181130.xlsx",sep=""), sheet="data_SNPsAncestralState", colNames=TRUE) #10 scaffold_78 SNPs has VarU
#ancestral_uli %>% tidyr::separate(Locus_Name, into = c("Scaffold","Position"), sep = ":") %>% filter(Scaffold == "scaffold_78" & AncestralState == "VarU")

# Adding moneduloides to ancestral - I used 5 inds called by GATK3.8 
mon <- read.table(file=paste(path,"./01_angsd_moneduloides/moneduloides_AA_1111.txt",sep=""),header=TRUE) %>%
  mutate(across(3:ncol(.), ~ case_when(. == "AA" ~ "A", . == "CC" ~ "C", . == "GG" ~ "G", . == "TT" ~ "T",
    . %in% c("CT", "TC", "AC", "CA", "CG", "GC", "TA", "AT", "AG", "GA", "GT", "TG") ~ "Y", . == ".." ~ "N", TRUE ~ .))) %>%
  rowwise() %>%
  mutate(summary_M = {
    alleles <- c_across(3:7)
    observed <- alleles[alleles %in% c("A", "C", "G", "T")]
    counts <- table(observed)
    if (length(observed) >= 1) {
      top_allele <- names(counts)[which.max(counts)]
      if (counts[top_allele] > length(observed) / 2) {top_allele} else {"Y"}} else {"N"}}) %>%
  ungroup() %>%
  mutate(Locus_Name = paste(CHROM,POS,sep=":")) %>%
  dplyr::select(Locus_Name, summary_M)

ancestral <- merge(ancestral_uli,mon, by="Locus_Name") %>%
  relocate(summary_M, .before = AncestralState) %>%
  mutate(AncestralState = case_when(summary_U == summary_M & summary_U %in% c("A","C","G","T") ~ summary_U, 
                                    summary_U == summary_R & summary_U %in% c("A","C","G","T") ~ summary_U, 
                                    summary_U == summary_J & summary_U %in% c("A","C","G","T") ~ summary_U,
                                    summary_U %in% c("R","Y","W","K","N","M","S","B") ~ summary_M, #segregating in U, thus AA inferred to be same as M
                                    summary_U %in% c("A","C","G","T") & summary_U != summary_M ~ summary_U, #fixed in U, thus AA inferred to be same as U
                                    Locus_Name == "scaffold_78:1682876" ~ summary_M, #this locus is identified as the most significant SNP by Uli but it is segregating in AC
                                    TRUE ~ AncestralState)) %>% ungroup()


# Sort according to CHR and POS
snp <- subset(snp, !(Locus_Name %in% c("scaffold_129:218821","scaffold_7:26556260","scaffold_144:237657","scaffold_87:2029")))
snp <- snp[with(snp, mixedorder(POS)), ]
snp <- snp[with(snp, mixedorder(CHR)), ]
newdat <- newdat[ ,c(1:15,match(snp$Locus_Name, colnames(newdat)))]
hindex <- read.xlsx(xlsxFile=paste(path,"data_snps_individuals_181130.xlsx",sep=""), sheet="data_HI_I_C", colNames=TRUE) %>%
  mutate(HI = coalesce(HI.C.18, HI.I.18)) %>%
  dplyr::select(sample_code, HI) 
merged_df <- merge(newdat,hindex, by="sample_code", all.x = TRUE)
merged_df <- merged_df[match(newdat$sample_code, merged_df$sample_code), ]
newdat[, 3] <- merged_df[, ncol(merged_df)] 
colnames(newdat)[3] <- colnames(merged_df)[ncol(merged_df)]
# ------------------------------------------------------------------------------------------------------------------------------------

# (3) Create PLINK input and run PLINK -----------------------------------------------------------------------------------------------
if(PLINK==TRUE) {
  # create .ped file for PLINK -----
  out.ped <- data.frame(matrix(rep(NA,nrow(newdat)*(6+(ncol(newdat)-15)*2)), nrow=nrow(newdat)))
  out.ped[ ,1] <- newdat$hybridcat
  out.ped[ ,2] <- newdat$sample_code
  out.ped[ ,5] <- newdat$Totsexing
  out.ped[ ,c(3,4,6)] <- matrix(rep(0,3*nrow(newdat)),ncol=3)
  for(i in 16:(ncol(newdat))) {
    out.ped[,2*i-25] <- sapply(newdat[,i],function(x) strsplit(x,"/",fixed=TRUE)[[1]][1])
    out.ped[,2*i-24] <- sapply(newdat[,i],function(x) strsplit(x,"/",fixed=TRUE)[[1]][2])
  }
  write.table(out.ped,paste(path,"add_fresh_anc_tmp_data_GG_crows.ped",sep=""), row.names=FALSE, col.names=FALSE, sep=" ", quote=FALSE)
   
  # create .map file for PLINK -----
  snps <- snp[,c("Locus_Name","CHR","POS")]
  out.map <- data.frame(matrix(rep(NA,2*(ncol(newdat)-15)), ncol=2))
  out.map[ ,1] <- colnames(newdat)[16:ncol(newdat)]
  out.map[ ,2] <- rep(0,ncol(newdat)-15)
  out.map <- merge(out.map, snps, by.x="X1", by.y="Locus_Name", sort=FALSE)
  out.map <- out.map[,c(3,1,2,4)]
  write.table(out.map,paste(path,"add_fresh_anc_tmp_data_GG_crows.map",sep=""), row.names=FALSE, col.names=FALSE, sep=" ", quote=FALSE)
  
  # run PLINK to convert .ped to binary .bed -----
  system(paste("/PATH/05_aDNA/02_results/uli_pca/plink --file ",path,"add_fresh_anc_tmp_data_GG_crows --make-bed --allow-extra-chr --chr-set 32 --out ",path,"add_fresh_anc_tmp_data_GG_crows",sep=""))
  
}
# ------------------------------------------------------------------------------------------------------------------------------------

# (4) Read data into GDS DB, polarize SNPs and perform PCA ---------------------------------------------------------------------------
# read data into R for SNPRelate -----
bed.fn <- paste(path,"add_fresh_anc_tmp_data_GG_crows.bed",sep="")
fam.fn <- paste(path,"add_fresh_anc_tmp_data_GG_crows.fam",sep="")
bim.fn <- paste(path,"add_fresh_anc_tmp_data_GG_crows.bim",sep="")
# convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, paste(path,"add_fresh_anc_tmp_data_GG_crows2.gds",sep=""), family=TRUE, cvt.chr="char")
snpgdsSummary(paste(path,"add_fresh_anc_tmp_data_GG_crows2.gds",sep=""))	
### --> 3 SNPs ("scaffold_144:237657","scaffold_87:2029","scaffold_129:218821") appear coded as 0/G or 0/T. This is an artifact (SNPs are fixed and some individuals have not been called). However, this has no consequences, in the genotypes file they are coded correctly as fixed or missing
SNP <- snpgdsOpen(paste(path,"add_fresh_anc_tmp_data_GG_crows2.gds",sep=""), readonly=FALSE)

# open database -----
genos <- read.gdsn(index.gdsn(SNP, "genotype"))
rownames(genos) <- read.gdsn(index.gdsn(SNP, "sample.id"))
colnames(genos) <- read.gdsn(index.gdsn(SNP, "snp.id"))
genos_MAF <- genos

# polarize SNPs (change coding from 0,1,2 copies of the minor allele to 0,1,2 copies of the derived allele) -----
ancestral <- ancestral[order(match(ancestral$Locus_Name,colnames(genos))), ]
xx <- c()
AFs <- data.frame(matrix(rep(NA,5*ncol(genos)),ncol=5))
colnames(AFs) <- c("Locus_Name","X0","X1","X2","X3")
for(i in 1:ncol(genos)) {
  out <- data.frame(table(genos[,i]))
  AFs[i,1] <- colnames(genos)[i]
  if(length(which(out$Var1==0))>0) { AFs[i,2] <- out$Freq[which(out$Var1==0)] }
  if(length(which(out$Var1==1))>0) { AFs[i,3] <- out$Freq[which(out$Var1==1)] }
  if(length(which(out$Var1==2))>0) { AFs[i,4] <- out$Freq[which(out$Var1==2)] }
  if(length(which(out$Var1==3))>0) { AFs[i,5] <- out$Freq[which(out$Var1==3)] }
  y <- data.frame(table(newdat[,which(colnames(newdat)==colnames(genos)[i])]))
  y.Hom <- subset(y, !(Var1 %in% c("A/C","A/G","A/T","C/A","C/G","C/T","G/A","G/C","G/T","T/A","T/C","T/G","0/0")))
  y.Het <- subset(y, Var1 %in% c("A/C","A/G","A/T","C/A","C/G","C/T","G/A","G/C","G/T","T/A","T/C","T/G"))
  if(nrow(y.Hom)==1 & nrow(y.Het)>0) { major.A <- substr(y.Hom$Var1[1],1,1); minor.A <- ifelse(substr(y.Het$Var1[1],1,1)==major.A, substr(y.Het$Var1[1],3,3), substr(y.Het$Var1[1],1,1)) }
  if(nrow(y.Hom)==2 & y.Hom$Freq[1]!=y.Hom$Freq[2]) { major.A <- substr(y.Hom$Var1[which(y.Hom$Freq==max(y.Hom$Freq))],1,1); minor.A <- substr(y.Hom$Var1[which(y.Hom$Freq==min(y.Hom$Freq))],1,1) }
  if(nrow(y.Hom)==1 & nrow(y.Het)==0) { major.A <- substr(y.Hom$Var1[1],1,1); minor.A <- "N" }
  if(nrow(y.Hom)==2 & y.Hom$Freq[1]==y.Hom$Freq[2]) { major.A <- substr(newdat[,which(colnames(newdat)==colnames(genos)[i])][genos[,i]==0][1],1,1); minor.A <- substr(newdat[,which(colnames(newdat)==colnames(genos)[i])][genos[,i]==2][1],1,1) }
  ancestral.A <- ancestral$AncestralState[which(ancestral$Locus_Name==colnames(genos)[i])]
  if(major.A==ancestral.A) { next }
  if(minor.A==ancestral.A) { NAs <- genos[,i]==3; genos[,i] <- abs(genos[,i]-2); genos[NAs,i] <- 3 }
  if(!(ancestral.A %in% c(major.A, minor.A))) { xx <- c(xx,i) }
  if(!(ancestral.A %in% c(major.A, minor.A))) { genos[,i] <- 3 }
}

# create new GDS DB with polarized SNP alleles (be aware that they both are called 'SNP') and open DB -----
setwd(paste(path,sep=""))
snpgdsCreateGeno("add_fresh_anc_tmp_data_GG_crows_polarized2.gds", genmat = genos, snpfirstdim=FALSE, sample.id = rownames(genos), snp.id = colnames(genos), snp.chromosome = read.gdsn(index.gdsn(SNP, "snp.chromosome")),
                 snp.position = read.gdsn(index.gdsn(SNP, "snp.position")), snp.allele = read.gdsn(index.gdsn(SNP, "snp.allele")))
SNP <- snpgdsOpen(paste(path,"add_fresh_anc_tmp_data_GG_crows_polarized2.gds",sep=""), readonly=FALSE)
# ------------------------------------------------------------------------------------------------------------------------------------

# (5) chr18 genomic landscape --------------------------------------------------------------------------------------------------------
# Based on my modern demo paper using 50 kb windows (NOT sliding, so peak value is lower that Poelstra's)
library(CMplot, lib.loc = "/dss/dsshome1/lxc0E/di67kah/R") 
path2mod = "/PATH/05_aDNA/02_results/uli_pca"
y <- read.delim(paste(path2mod,"/hz1_cnx3P_cor2.win.chr18", sep=""), header=TRUE, sep="")
sorted.y <- y %>% filter(scaffold %in% c("scaffold_60", "scaffold_78")) %>%
  mutate(scaffold = factor(scaffold, levels = c("scaffold_78", "scaffold_60"))) %>% arrange(scaffold)
z <- seq(25000, length.out=nrow(sorted.y), by=50000) #dummy position to join the two scaffolds
yz <- cbind(z, sorted.y)
hd <- seq(2525000, 4275000 , by=50000)
yelevated <- data_frame() 
for (region in hd) 
{re <- paste("^",region,"$", sep="")
  elevated <- yz[grepl(re, yz$z),]
  yelevated <- rbind(yelevated, elevated)}
yz_modified <- yz %>%
  group_by(scaffold) %>%
  mutate(Fst = case_when(
    scaffold == "scaffold_78" & row_number(desc(z)) <= 2 ~ NA_real_,
    scaffold == "scaffold_60" & z == min(z) ~ NA_real_,
    TRUE ~ Fst
  )) %>% ungroup()

landscape <- ggplot(yz_modified, aes(x = z, y = Fst)) +
  geom_line(color = "black", size = 1) +  scale_x_reverse() +
  labs(x = "Peak region of chromosome 18 (bp)", y = "Fst") + coord_cartesian(xlim = c(4000000,500000))+
  scale_x_continuous(labels = scientific, position = "top") + theme_minimal() 

# (6) Heat map ----------------------------------------------------------------------------------------------------------------------------
# Combine genotype with available colour PC info
colourpc <- read.xlsx(xlsxFile=paste(path,"data_FIS_FST_PCA_Rec_bgc_GWAS_Phenos_180621.xlsx",sep=""), sheet="Individual_based", colNames=TRUE)
gwasdat <- read.xlsx(xlsxFile=paste(path,"data_FIS_FST_PCA_Rec_bgc_GWAS_Phenos_180621.xlsx",sep=""), sheet="SNP_based", colNames=TRUE) %>%
  dplyr::select("Locus_Name","GWAS_All_chi2.1df")
genos2 <- as.data.frame(genos) %>% 
  mutate(zone = newdat$zone) %>%
  mutate(country_code = newdat$country_code) %>%
  mutate(location_code = newdat$location_code) %>%
  mutate(HI = newdat$HI) %>%
  mutate(HI = ifelse(country_code %in% c("E", "F", "D") & zone %in% c("fresh") & is.na(HI), 0, HI)) %>%
  mutate(HI = ifelse(country_code %in% c("IRQ", "ISR", "B", "S", "IT", "Pl", "PL", "RUS") & zone %in% c("fresh") & is.na(HI), 1, HI)) %>%
  mutate(admix = newdat$newhybrids_P0P1Allo.admix99) %>%
  mutate(classification = newdat$classification)

heatmap <- genos2 %>%
  filter(country_code != "IRE" & country_code != "COR" & country_code != "ori" & country_code != "pec") %>% #IRE hybrid index undetermined & COR subject to island effect
  rownames_to_column(var = "sample_code") %>% 
  mutate(type = ifelse(country_code %in% c("E", "F", "D","BE","GB") & zone =="ancient" & sample_code != "LSS003", "ancient_black", zone)) %>%
  mutate(type = ifelse(type %in% c("South", "East","AllopatricEastCC","AllopatricSouthCC","AllopatricSouthHC","AllopatricEastHC"), "fresh", type)) %>%
  pivot_longer(cols = starts_with("scaffold"), names_to = "scaffold",values_to = "genotype") %>%
  separate(scaffold, into = c("scaff", "pos"), sep=":", remove=FALSE) %>%
  left_join(colourpc %>% dplyr::select("sample_code", "Colour_PC1"), by = "sample_code") %>%
  mutate(Colour_PC1 = ifelse(country_code == "E" & zone %in% c("fresh") & HI <= 0.1 & is.na(Colour_PC1), 4.33, Colour_PC1)) %>% #we are sure the Spanish and Constance samples are pure black
  mutate(Colour_PC1 = ifelse(country_code == "D" & location_code == "Ko" & is.na(Colour_PC1), 4.33, Colour_PC1)) %>% #we are sure these samples are pure black
  mutate(Colour_PC1 = ifelse(country_code %in% c("IRQ", "ISR","S", "Pl") & is.na(Colour_PC1), -3.74, Colour_PC1)) #we are sure these samples are pure grey
heatmap$genotype <- as.factor(heatmap$genotype)

# Calculate FST between pure black vs pure grey
#who are pure black? : country_code %in% c("E", "F", "D") & zone %in% c("ancient", "fresh","AllopatricEastCC","AllopatricSouthCC") & sample_code != "LSS003" #but some AllopatricEastCC (D_Do_C08-C10/C04-C05/C14/C17, D_Kl_C09/C01/C14/C19), FPd01 & AllopatricSouthCC (DKoC62/ DKoC27) are not pure black in HI
#no - only those with PC1> 4.3 (those black ones coded by uli and spain+constance crows) = 69 ind
#who are pure grey? : country_code %in% c("IRQ", "ISR", "B", "S", "IT", "Pl", "PL", "RUS") & zone %in% c("ancient", "fresh","AllopatricEastHC","AllopatricSouthHC")
#no - only those with PC1<- -3.7 (those grey ones by uli and IRQ+ISR+S samples) = 49 ind
#in genotype: count how many 0s, 1s and 2s (3s are missing). count(0)=pp, count(1)=pq, count(2)=qq
pcount <- heatmap %>% 
  #mutate(pure = ifelse(country_code %in% c("E", "F", "D") & zone %in% c("ancient", "fresh","AllopatricEastCC","AllopatricSouthCC") & sample_code != "LSS003", "black", NA)) %>%
  mutate(pure = ifelse(Colour_PC1 > 4.3, "black", NA)) %>%
  #mutate(pure = ifelse(country_code %in% c("IRQ", "ISR", "B", "S", "IT", "Pl", "PL", "RUS") & zone %in% c("ancient", "fresh","AllopatricEastHC","AllopatricSouthHC"), "grey", pure)) %>%
  mutate(pure = ifelse(Colour_PC1 < -3.7, "grey", pure)) %>%
  filter(pure %in% c("black","grey")) %>%
  group_by(scaffold, pure) %>%
  summarise(pp = sum(genotype == 0, na.rm = TRUE), pq = sum(genotype == 1, na.rm = TRUE), qq = sum(genotype == 2 , na.rm = TRUE), ind = sum(genotype %in% c(0, 1, 2), na.rm = TRUE)) %>%
  mutate(p = (pp*2 + pq)/(2*ind), q = 1-p, pq2 = 2*p*q) %>%
  ungroup() %>%
  group_by(scaffold) %>%
  mutate(t_p = (sum(pp)*2+sum(pq))/(sum(ind)*2), t_q =(sum(qq)*2+sum(pq))/(sum(ind)*2), ht = 2*t_p*t_q, hs = sum(pq2)/2, fst = ifelse(ht > 0, pmax((ht - hs) / ht, 0), NA)) %>%
  mutate(q_black = ifelse(pure == "black", q, NA)) %>%
  mutate(q_grey = ifelse(pure == "grey", q, NA)) %>%
  mutate(q_grey = ifelse(is.na(q_grey[1]) & !is.na(q_grey[2]), q_grey[2], q_grey[1])) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(scaffold, ht, hs, fst, q_black, q_grey)
heatmap <- merge(heatmap, pcount, by="scaffold")

# Merge with GWAS results
heatmap <- merge(heatmap, gwasdat, by.x="scaffold", by.y="Locus_Name")

# Define the custom order for country_code and colour
custom_order <- c("E", "F", "GB", "BE", "D","IT","NL","Pl","PL","S","RUS","B","ISR","IRQ")
custom_order2 <- c("ancient_black", "fresh","ancient") 
heatmap$country_code <- ifelse(heatmap$sample_code == "LSS003", "NL", as.character(heatmap$country_code))
heatmap$country_code <- ifelse(heatmap$country_code %in% c("A","IT") & heatmap$zone != "fresh", "D", as.character(heatmap$country_code)) #so that we can sort them by HI not affected by countrycode
heatmap <- heatmap %>% arrange(factor(type, levels = custom_order2),factor(country_code, levels = custom_order), HI, desc(Colour_PC1), as.numeric(classification))
heatmap$sample_code <- forcats::fct_reorder(heatmap$sample_code, heatmap$Colour_PC1,.fun = max, .na_rm = FALSE)
heatmap$sample_code <- factor(heatmap$sample_code, levels = unique(heatmap$sample_code))
heatmap$sample_code <- forcats::fct_rev(heatmap$sample_code)
original_order <- unique(heatmap$sample_code)

genedat <- read.table(file="NCBI_Corvus_cornix_cornix_Annotation_Release_101_peak.CSV", sep=",",header=TRUE) %>%
  mutate(scaffold = ifelse(Accession == "NW_010959954.1","scaffold_78","scaffold_60"))

scaffold_levels <- c("scaffold_60", "scaffold_78")
heatmap_pos <- heatmap %>% 
  filter((scaff == "scaffold_78" & pos >= 1006603 & pos <= 2510417) | (scaff == "scaffold_60" & as.numeric(pos) <= 1000000)) %>%
  group_by(scaff) %>%
  filter(any(genotype != 3)) %>%
  arrange(scaff, desc(as.numeric(pos))) %>%
  ungroup() %>%
  mutate(scaff = factor(scaff, levels = scaffold_levels),
         pos = factor(pos, levels = unique(pos)))

# Include the annotation from snpEff
vcf <- read.vcfR("add_fresh_anc_tmp_data_GG_crows_peaks_modified_7860_annotated.vcf")
tidy_vcf <- vcfR2tidy(vcf)
anndat <- tidy_vcf$fix %>%
  mutate(variant_type = sapply(strsplit(ANN, ","), function(x) strsplit(x[1], "\\|")[[1]][2]), variant_type = gsub("_variant", "", variant_type)) %>%
  dplyr::select(POS, ID, variant_type) %>%
  rename(pos="POS")
anndat$pos <- as.character(anndat$pos)

# First I plot all SNPs in this peak regions -----------------------------------------
# Step 0: Helper function
get_snps_in_gene <- function(gene_start, gene_stop, snp_positions) {
  snp_positions[snp_positions >= gene_start & snp_positions <= gene_stop]}
# Step 1: Prepare SNP positions (numeric and factor)
heatmap_pos_numeric <- as.numeric(as.character(heatmap_pos$pos))
heatmap_pos$pos <- factor(heatmap_pos$pos, levels = sort(unique(heatmap_pos$pos)))

# Step 2: Filter genes that actually overlap with SNP positions
genes_shown <- genedat %>%
  filter(!grepl("^LOC", Gene_symbol)) %>%
  filter(!grepl("^CUNH", Gene_symbol)) %>%
  filter(sapply(1:nrow(.), function(i) {
    gene_range <- as.numeric(Start[i]):as.numeric(Stop[i])
    any(gene_range %in% heatmap_pos_numeric) }))
# Step 3: Get SNPs within each gene
gene_snps <- genes_shown %>%
  rowwise() %>%
  mutate(SNPs_in_gene = list(get_snps_in_gene(as.numeric(Start), as.numeric(Stop), heatmap_pos_numeric))) %>%
  unnest(cols = c(SNPs_in_gene)) %>%
  rename(pos = SNPs_in_gene) %>%
  mutate(pos = factor(pos, levels = levels(heatmap_pos$pos))) %>%
  distinct(Gene_symbol, pos, .keep_all = TRUE)
# Step 4: Ensure we have all SNPs in the plot, even those with no gene
all_snps <- tibble(pos = factor(heatmap_pos$pos, levels = levels(heatmap_pos$pos))) %>% unique
gene_snps_full <- all_snps %>%
  left_join(gene_snps, by = "pos") %>%
  left_join(anndat, by= "pos") %>%
  unique() 


# Finally plot!!!!
gene_snps_full$scaffold <- factor(gene_snps_full$scaffold, levels = c("scaffold_60","scaffold_78"))
gene_snps_full$pos <- factor(gene_snps_full$pos, levels = unique(gene_snps_full$pos))
gene_snps_full$Gene_symbol <- factor(gene_snps_full$Gene_symbol,levels = unique(gene_snps_full$Gene_symbol))
g1 <- ggplot(gene_snps_full, aes(x = pos, y = "Gene", fill = Gene_symbol)) +
  geom_tile() +  
  scale_fill_viridis_d(na.value = "white", option = "D", name="Gene", guide = guide_legend(nrow = 4)) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(),legend.position = "",
        legend.key.size = unit(0.3, "cm"),legend.text = element_text(size = 8))

gene_snps_full$variant_type <- factor(gene_snps_full$variant_type, 
                                      levels = c("5_prime_UTR", "3_prime_UTR", "upstream_gene", "downstream_gene","intergenic_region","intron","intragenic", "synonymous","missense","splice_region&intron" ))
ann1 <- ggplot(gene_snps_full, aes(x = pos, y = "Type", shape = variant_type)) +
  geom_point() +
  scale_shape_manual(name="Type", values=c(0,7,2,6,3,1,4,5,8,11), guide = guide_legend(nrow = 2)) +
  guides(shape = guide_legend(override.aes = list(size = 2), nrow = 2)) +
  #scale_fill_viridis_d(na.value = "white", option = "D", name="Gene", guide = guide_legend(nrow = 2)) +
  theme_minimal() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(),legend.position = "",
        legend.key.size = unit(0.3, "cm"),legend.text = element_text(size = 8))

anctop <- as.character(original_order[1:19])
hanctop1 <- ggplot(heatmap_pos %>% filter(sample_code %in% anctop)) + 
  geom_tile(aes(x=pos, y=sample_code, fill = genotype)) + 
  scale_fill_manual(values=c("0"="#08306B","1"="#1F78B4","2"="#A6CEE3", "3"="white"),
                    labels = c("0" = "Ancestral", "1" = "Heterozygous", "2" = "Derived", "3" = "Missing"), name = "Genotype") + labs(y="Ancient black") +
  theme(axis.text.y = element_text(size = 6), legend.position = "none",
        axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank()) 
  
ancbot <- as.character(tail(original_order, 28))
hancbot1 <- ggplot(heatmap_pos %>% filter(sample_code %in% ancbot )) + 
  geom_tile(aes(x=pos, y=sample_code, fill = genotype)) + 
  scale_fill_manual(values=c("0"="#08306B","1"="#1F78B4","2"="#A6CEE3", "3"="white"),
                    labels = c("0" = "Ancestral", "1" = "Heterozygous", "2" = "Derived", "3" = "Missing"), name = "Genotype") + labs(y="Ancient grey") +
  theme(axis.text.y = element_text(size = 6), legend.position = "none", axis.title.x=element_blank(), 
        axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1), axis.ticks.x=element_blank()) 

modern <- heatmap_pos %>% filter(!sample_code %in% c(anctop,ancbot))
important_samples <- c("D_Rb_Y23", "D_Ne_Y39", "DKoC53", "DLoC21", "D_Lo_C21")
h1 <- ggplot(modern) + 
  geom_tile(aes(x=pos, y=sample_code, fill = genotype)) + 
  scale_fill_manual(values=c("0"="#08306B","1"="#1F78B4","2"="#A6CEE3", "3"="white"),
                    labels = c("0" = "Ancestral", "1" = "Heterozygous", "2" = "Derived", "3" = "Missing"), name = "Genotype") + labs(y="Modern") +
  #scale_y_discrete(breaks = c(as.character(modern$sample_code[seq(1, length(unique(modern$sample_code)), by = 50)]), as.character(modern$sample_code[1]), as.character(modern$sample_code[length(unique(modern$sample_code))]), important_samples)) +
  scale_y_discrete(breaks = c(as.character(modern$sample_code[1]), as.character(modern$sample_code[length(unique(modern$sample_code))]), important_samples)) +
  theme(axis.text.y = element_text(size = 6), legend.position = "top",
        axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank()) 

ndpanctop1 <- ggplot(data = heatmap %>% filter(scaff=="scaffold_2") %>%
                 filter(pos=="41483402" | pos == "41437043" | pos == "41510734") %>%
                 group_by(scaffold) %>%
                 filter(any(genotype != 3)) %>%
                 filter(sample_code %in% c(anctop)), aes(x=pos, y=sample_code, fill = genotype)) + geom_tile() + 
  scale_fill_manual(values=c("0"="#08306B","1"="#1F78B4","2"="#A6CEE3", "3"="white"),
                    labels = c("0" = "Ancestral", "1" = "Heterozygous", "2" = "Derived", "3" = "Missing"), name = "Genotype") +
  theme(axis.text.y = element_blank(), legend.position = "none", axis.title.y=element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank())+ guides(fill = "none")

ndpancbot1 <- ggplot(data = heatmap %>% filter(scaff=="scaffold_2") %>%
                       filter(pos=="41483402" | pos == "41437043" | pos == "41510734") %>%
                       group_by(scaffold) %>%
                       filter(any(genotype != 3)) %>%
                       filter(sample_code %in% c(ancbot)), aes(x=pos, y=sample_code, fill = genotype)) + 
  geom_tile() + 
  scale_fill_manual(values=c("0"="#08306B","1"="#1F78B4","2"="#A6CEE3", "3"="white"),
                    labels = c("0" = "Ancestral", "1" = "Heterozygous", "2" = "Derived", "3" = "Missing"), name = "Genotype") +
  theme(axis.text.y = element_blank(), legend.position = "none", axis.title.y=element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_text(size = 6,angle = 90, vjust = 0.5, hjust = 1), axis.ticks.x=element_blank())+ guides(fill = "none")

ndp1 <- ggplot(data = heatmap %>% filter(scaff=="scaffold_2") %>%
         filter(pos=="41483402" | pos == "41437043" | pos == "41510734") %>%
         group_by(scaffold) %>%
         filter(any(genotype != 3)) %>%
         filter(!sample_code %in% c(anctop,ancbot)), aes(x=pos, y=sample_code, fill = genotype)) + 
  geom_tile() + 
  scale_fill_manual(values=c("0"="#08306B","1"="#1F78B4","2"="#A6CEE3", "3"="white"),
                    labels = c("0" = "Ancestral", "1" = "Heterozygous", "2" = "Derived", "3" = "Missing"), name = "Genotype") +
  theme(axis.text.y = element_blank(), legend.position = "none", axis.title.y=element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank())+ guides(fill = "none")

canctop1 <- ggplot(heatmap_pos %>% filter(sample_code %in% c(anctop))) + 
  geom_tile(aes(x=1, y=sample_code, fill = Colour_PC1)) + 
  scale_fill_gradient(na.value ="white", low = "lightgrey", high = "black") +
  theme(axis.text.y = element_blank(), axis.title.y=element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.position = "")+ guides(fill = "none")

cancbot1 <- ggplot(heatmap_pos %>% filter(sample_code %in% c(ancbot))) + 
  geom_tile(aes(x=1, y=sample_code, fill = Colour_PC1)) + 
  scale_fill_gradient(na.value ="white", low = "lightgrey", high = "black") +
  theme(axis.text.y = element_blank(), axis.title.y=element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.position = "")+ guides(fill = "none")

c1 <- ggplot(heatmap_pos %>% filter(!sample_code %in% c(anctop,ancbot))) + 
  geom_tile(aes(x=1, y=sample_code, fill = Colour_PC1)) + 
  scale_fill_gradient(na.value ="white", low = "lightgrey", high = "black") +
  theme(axis.text.y = element_blank(), axis.title.y=element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.position = "none") + guides(fill = "none")

fst1 <- ggplot(heatmap_pos, aes(x = pos, y = "FST", fill = fst)) +
  geom_tile() + 
  scale_fill_gradientn( na.value="white", colors = c("#fef0d9", "#fdae61", "#8c2d04"),
    values = c(0, 0.5, 1), 
    limits = c(0, 1)) + guides(fill = "none") + theme_minimal() +
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none", legend.text = element_text(size = 8))

derived_black <- ggplot(heatmap_pos, aes(x = pos, y = "DAF_black", fill = q_black)) +
  geom_tile() + 
  scale_fill_gradientn(na.value="black", colors = c("#fef0d9", "#fdae61","#8c2d04"),
    values = c(0, 0.5, 1), 
    limits = c(0, 1)) + guides(fill = "none") + theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8), axis.ticks.y = element_blank(), 
        axis.title.x = element_blank(),  axis.text.x = element_blank(), legend.position = "none", legend.text = element_text(size = 8))

derived_grey <- ggplot(heatmap_pos, aes(x = pos, y = "DAF_grey", fill = q_grey)) +
  geom_tile() + 
  scale_fill_gradientn(na.value="black", colors = c("#fef0d9", "#fdae61","#8c2d04"),
    values = c(0,  0.5, 1), 
    limits = c(0,1)) + guides(fill = "none") + theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none", legend.text = element_text(size = 8))

gwas <- ggplot(heatmap_pos, aes(x = pos, y = "GWAS", fill = GWAS_All_chi2.1df)) + geom_tile() +
  scale_fill_gradientn(na.value="black", colors = c("#fef0d9", "#fdae61","#8c2d04"),
    #values = c(0, quantile(heatmap_pos$GWAS_All_chi2.1df, 0.75), 1), 
    limits = c(0, 100)) + guides(fill = "none") + theme_minimal() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none", legend.text = element_text(size = 8))

row0 <- hanctop1 + ndpanctop1 + canctop1 + plot_layout(widths = c(3, 0.1,0.1))
row1 <- h1 + ndp1 + c1 + plot_layout(widths = c(3, 0.1,0.1))
row1add <- hancbot1 + ndpancbot1 +  cancbot1 + plot_layout(widths = c(3, 0.1,0.1))
row2 <- fst1 + plot_spacer() + plot_spacer() + plot_layout(widths = c(3, 0.1,0.1))
row3 <- derived_black + plot_spacer() + plot_spacer() + plot_layout(widths = c(3, 0.1,0.1))
row4 <- derived_grey + plot_spacer() + plot_spacer() + plot_layout(widths = c(3, 0.1,0.1))
row5 <- gwas + plot_spacer() + plot_spacer() + plot_layout(widths = c(3, 0.1,0.1))
row6 <- g1 + plot_spacer() + plot_spacer() + plot_layout(widths = c(3, 0.1,0.1))
row7 <- ann1 + plot_spacer() + plot_spacer() + plot_layout(widths = c(3, 0.1,0.1))
fst_landscape <-landscape + plot_spacer() + plot_spacer() + plot_layout(widths = c(3, 0.1,0.1))
  
# Combine with relative heights
legend6 <- get_legend(g1 + theme(legend.position = "bottom"))
legend7 <- get_legend(ann1 + theme(legend.position = "bottom"))
combined_legends <- plot_grid(legend6, legend7, ncol = 1,rel_heights = c(1, 1),align = "v")
 
top_legend_plot <- row0 / row1 / row1add + plot_layout(guides = "collect") &
  theme(legend.position = "top", plot.margin = unit(c(0, 0, 0, 0), "pt"))
bottom_legend_plots <- row2 / row3 / row4 / row5 / row6 / row7 &
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))
final_plot <- fst_landscape / top_legend_plot / bottom_legend_plots / combined_legends + plot_layout(heights = c(0.4, 3, 0.6, 0.5))
final_plot #14x12
#png("heatmap_2peaks_with genes_andmore_reverse.png", width =14 , height =12, units = "in", res = 400)
#final_plot 
#dev.off()

#------------------------------------------------------------------------------------------------
# Second I plot only selected SNPs in this peak regions -----------------------------------------
#What sites? Sites with FST>0.8 --> higher differentiated between the pure black and pure grey 
#What sites? Maybe also select sites based on if it's fixed (hom:0/2 state) within SPA/ IRQ+BUL+ISR - anc+fresh: fst>0.8: 73 SNPs, selected: 30 SNPs
sites1 <- heatmap %>% filter(scaff=="scaffold_78") %>%
  filter(pos>=1006603 & pos<=2510417) %>%
  filter(zone == "ancient" | zone == "fresh") %>%
  filter(sample_code != "LSS003") %>%
  filter(sample_code != "FPd01") %>%
  filter(country_code == "E" | country_code == "F" | country_code == "D"  | country_code == "B" | country_code == "ISR" | country_code == "IRQ") %>%
  #filter(country_code == "E" | country_code == "B" | country_code == "ISR" | country_code == "IRQ") %>%
  group_by(scaffold) %>%
  filter(all(genotype!=1) & any(genotype !=3)) %>%
  filter(sum(genotype == 3) < 20) #sites with high missingness are removed
as.data.frame(sort(unique(sites1$scaffold)))

# fixed snps in bulgaria in the second peak region
sites2 <- heatmap %>% filter(scaff=="scaffold_60") %>%
  filter(as.numeric(pos)<=1000000) %>%
  filter(zone == "ancient" | zone == "fresh") %>%
  filter(country_code == "E" | country_code == "F" | country_code == "D" | country_code == "B" | country_code == "ISR" | country_code == "IRQ") %>%
  group_by(scaffold) %>%
  filter(all(genotype!=1) & any(genotype !=3)) %>%
  filter(sum(genotype == 3) < 20) 
as.data.frame(sort(unique(sites2$scaffold)))
selectedsites <- c(unique(sites1$pos), unique(sites2$pos))

scaffold_levels <- c("scaffold_60","scaffold_78")
heatmap_pos <- heatmap %>% 
  filter(fst>0.8) %>% 
  filter(pos %in% selectedsites) %>%
  filter((scaff == "scaffold_78" & pos >= 1006603 & pos <= 2510417) | (scaff == "scaffold_60" & as.numeric(pos) <= 1000000)) %>%
  group_by(scaff) %>%
  filter(any(genotype != 3)) %>%
  arrange(scaff, desc(as.numeric(pos))) %>%
  ungroup() %>%
  mutate(scaff = factor(scaff, levels = scaffold_levels),pos = factor(pos, levels = unique(pos)))

# Step 0: Helper function
get_snps_in_gene <- function(gene_start, gene_stop, snp_positions) {
  snp_positions[snp_positions >= gene_start & snp_positions <= gene_stop]}
# Step 1: Prepare SNP positions (numeric and factor)
heatmap_pos_numeric <- as.numeric(as.character(heatmap_pos$pos))
heatmap_pos$pos <- factor(heatmap_pos$pos, levels = sort(unique(heatmap_pos$pos)))
# Step 2: Filter genes that actually overlap with SNP positions
genes_shown <- genedat %>%
  filter(!grepl("^LOC", Gene_symbol)) %>%
  filter(!grepl("^CUNH", Gene_symbol)) %>%
  filter(sapply(1:nrow(.), function(i) {
    gene_range <- as.numeric(Start[i]):as.numeric(Stop[i])
    any(gene_range %in% heatmap_pos_numeric) }))
# Step 3: Get SNPs within each gene
gene_snps <- genes_shown %>%
  rowwise() %>%
  mutate(SNPs_in_gene = list(get_snps_in_gene(as.numeric(Start), as.numeric(Stop), heatmap_pos_numeric))) %>%
  unnest(cols = c(SNPs_in_gene)) %>%
  rename(pos = SNPs_in_gene) %>%
  mutate(pos = factor(pos, levels = levels(heatmap_pos$pos))) %>%
  distinct(Gene_symbol, pos, .keep_all = TRUE)
# Step 4: Ensure we have all SNPs in the plot, even those with no gene
all_snps <- tibble(pos = factor(heatmap_pos$pos, levels = levels(heatmap_pos$pos))) %>% unique
gene_snps_full <- all_snps %>%
  left_join(gene_snps, by = "pos") %>%
  left_join(anndat, by= "pos") %>%
  unique() 


# Finally plot!!!!
gene_snps_full$scaffold <- factor(gene_snps_full$scaffold, levels = c("scaffold_60","scaffold_78"))
gene_snps_full$pos <- factor(gene_snps_full$pos, levels = unique(gene_snps_full$pos))
gene_snps_full$Gene_symbol <- factor(gene_snps_full$Gene_symbol,levels = unique(gene_snps_full$Gene_symbol))
g2 <- ggplot(gene_snps_full, aes(x = pos, y = "Gene", fill = Gene_symbol)) +
  geom_tile() +
  scale_fill_viridis_d(na.value = "white", option = "D", name="Gene", guide = guide_legend(nrow = 1)) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(),legend.position = "",
        legend.key.size = unit(0.3, "cm"),legend.text = element_text(size = 8))

gene_snps_full$variant_type <- factor(gene_snps_full$variant_type, levels = c("5_prime_UTR", "3_prime_UTR", "upstream_gene", "downstream_gene",
                                                                              "intergenic_region","intron","intragenic", "synonymous","missense","splice_region&intron" ))
ann2 <- ggplot(gene_snps_full, aes(x = pos, y = "Type", shape = variant_type)) +
  geom_point() +
  scale_shape_manual(name="Type", values=c("5_prime_UTR"=0,"3_prime_UTR"=7,"upstream_gene"=2,"downstream_gene"=6,"intergenic_region"=3,"intron"=1,"intragenic"=4,"synonymous"=5,"missense"=8,"splice_region&intron"=11), guide = guide_legend(nrow = 2)) +
  guides(shape = guide_legend(override.aes = list(size = 2), nrow = 1)) +
  #scale_fill_viridis_d(na.value = "white", option = "D", name="Gene", guide = guide_legend(nrow = 2)) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(),legend.position = "",
        legend.key.size = unit(0.3, "cm"),legend.text = element_text(size = 8))

anctop <- as.character(original_order[1:19])
hanctop2 <- ggplot(heatmap_pos %>% filter(sample_code %in% anctop)) + 
  geom_tile(aes(x=pos, y=sample_code, fill = genotype)) +
  scale_fill_manual(values=c("0"="#08306B","1"="#1F78B4","2"="#A6CEE3", "3"="white"),
                    labels = c("0" = "Ancestral", "1" = "Heterozygous", "2" = "Derived", "3" = "Missing"), name = "Genotype") +
  labs(y="Ancient black") +
  theme(axis.text.y = element_text(size = 5), legend.position = "none",
        axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank()) 

ancbot <- as.character(tail(original_order, 28))
hancbot2 <- ggplot(heatmap_pos %>% filter(sample_code %in% ancbot )) + 
  geom_tile(aes(x=pos, y=sample_code, fill = genotype)) +
  scale_fill_manual(values=c("0"="#08306B","1"="#1F78B4","2"="#A6CEE3", "3"="white"),
                    labels = c("0" = "Ancestral", "1" = "Heterozygous", "2" = "Derived", "3" = "Missing"), name = "Genotype") + labs(y="Ancient grey") +
  theme(axis.text.y = element_text(size = 5), legend.position = "none", axis.title.x=element_blank(), 
        axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1), axis.ticks.x=element_blank()) 

modern <- heatmap_pos %>% filter(!sample_code %in% c(anctop,ancbot))
important_samples <- c("D_Rb_Y23", "D_Ne_Y39", "DKoC53", "D_Lo_C21")
h2 <- ggplot(modern) + 
  geom_tile(aes(x=pos, y=sample_code, fill = genotype)) +
  scale_fill_manual(values=c("0"="#08306B","1"="#1F78B4","2"="#A6CEE3", "3"="white"),
                    labels = c("0" = "Ancestral", "1" = "Heterozygous", "2" = "Derived", "3" = "Missing"), name = "Genotype") + labs(y="Modern") +
  scale_y_discrete(breaks = c(as.character(modern$sample_code[seq(1, length(unique(modern$sample_code)), by = 50)]), as.character(modern$sample_code[1]), as.character(modern$sample_code[length(unique(modern$sample_code))]), important_samples)) +
  theme(axis.text.y = element_text(size =5), legend.position = "top",
        axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank()) 

ndpanctop2 <- ggplot(data = heatmap %>% filter(scaff=="scaffold_2") %>%
                       filter(pos=="41483402" | pos == "41437043" | pos == "41510734") %>%
                       group_by(scaffold) %>%
                       filter(any(genotype != 3)) %>%
                       filter(sample_code %in% c(anctop)), aes(x=pos, y=sample_code, fill = genotype)) + geom_tile() +
  scale_fill_manual(values=c("0"="#08306B","1"="#1F78B4","2"="#A6CEE3", "3"="white"),
                    labels = c("0" = "Ancestral", "1" = "Heterozygous", "2" = "Derived", "3" = "Missing"), name = "Genotype") +
  theme(axis.text.y = element_blank(), legend.position = "none", axis.title.y=element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank())+ guides(fill = "none")

ndpancbot2 <- ggplot(data = heatmap %>% filter(scaff=="scaffold_2") %>%
                       filter(pos=="41483402" | pos == "41437043" | pos == "41510734") %>%
                       group_by(scaffold) %>%
                       filter(any(genotype != 3)) %>%
                       filter(sample_code %in% c(ancbot)), aes(x=pos, y=sample_code, fill = genotype)) + geom_tile() +
  scale_fill_manual(values=c("0"="#08306B","1"="#1F78B4","2"="#A6CEE3", "3"="white"),
                    labels = c("0" = "Ancestral", "1" = "Heterozygous", "2" = "Derived", "3" = "Missing"), name = "Genotype") +
  theme(axis.text.y = element_blank(), legend.position = "none", axis.title.y=element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_text(size = 6,angle = 90, vjust = 0.5, hjust = 1), axis.ticks.x=element_blank())+ guides(fill = "none")

ndp2 <- ggplot(data = heatmap %>% filter(scaff=="scaffold_2") %>%
                 filter(pos=="41483402" | pos == "41437043" | pos == "41510734") %>%
                 group_by(scaffold) %>%
                 filter(any(genotype != 3)) %>%
                 filter(!sample_code %in% c(anctop,ancbot)), aes(x=pos, y=sample_code, fill = genotype)) + geom_tile() +
  scale_fill_manual(values=c("0"="#08306B","1"="#1F78B4","2"="#A6CEE3", "3"="white"),
                    labels = c("0" = "Ancestral", "1" = "Heterozygous", "2" = "Derived", "3" = "Missing"), name = "Genotype") +
  theme(axis.text.y = element_blank(), legend.position = "none", axis.title.y=element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank())+ guides(fill = "none")

canctop2 <- ggplot(heatmap_pos %>% filter(sample_code %in% c(anctop))) + 
  geom_tile(aes(x=1, y=sample_code, fill = Colour_PC1)) +
  scale_fill_gradient(na.value ="white", low = "lightgrey", high = "black") +
  theme(axis.text.y = element_blank(), axis.title.y=element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.position = "")+ guides(fill = "none")

cancbot2 <- ggplot(heatmap_pos %>% filter(sample_code %in% c(ancbot))) + 
  geom_tile(aes(x=1, y=sample_code, fill = Colour_PC1)) +
  scale_fill_gradient(na.value ="white", low = "lightgrey", high = "black") +
  theme(axis.text.y = element_blank(), axis.title.y=element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.position = "")+ guides(fill = "none")

c2 <- ggplot(heatmap_pos %>% filter(!sample_code %in% c(anctop,ancbot))) + 
  geom_tile(aes(x=1, y=sample_code, fill = Colour_PC1)) +
  scale_fill_gradient(na.value ="white", low = "lightgrey", high = "black") +
  theme(axis.text.y = element_blank(), axis.title.y=element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.position = "")+ guides(fill = "none")

fst2 <- ggplot(heatmap_pos, aes(x = pos, y = "FST", fill = fst)) +
  geom_tile() +
  scale_fill_gradientn( na.value="white", colors = c("#fef0d9", "#fdae61","#8c2d04"),
                        values = c(0, 0.5, 1), 
                        limits = c(0, 1)) + theme_minimal() +
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "", legend.text = element_text(size = 8))
derived_black2 <- ggplot(heatmap_pos, aes(x = pos, y = "DAF_black", fill = q_black)) +
  geom_tile() +
  scale_fill_gradientn(na.value="black", colors = c("#fef0d9", "#fdae61","#8c2d04"),
                       values = c(0, 0.5, 1), 
                       limits = c(0, 1)) + theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8), axis.ticks.y = element_blank(), 
        axis.title.x = element_blank(),  axis.text.x = element_blank(), legend.position = "", legend.text = element_text(size = 8))+ guides(fill = "none")
derived_grey2 <- ggplot(heatmap_pos, aes(x = pos, y = "DAF_grey", fill = q_grey)) +
  geom_tile() +
  scale_fill_gradientn(na.value="black", colors = c("#fef0d9", "#fdae61","#8c2d04"),
                       values = c(0, 0.5, 1), 
                       limits = c(0, 1)) + theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "", legend.text = element_text(size = 8))+ guides(fill = "none")

gwas2 <- ggplot(heatmap_pos, aes(x = pos, y = "GWAS", fill = GWAS_All_chi2.1df)) +
  geom_tile() +
  scale_fill_gradientn(na.value="black", colors = c("#fef0d9", "#fdae61","#8c2d04"),
                       #values = c(0, quantile(heatmap_pos$GWAS_All_chi2.1df, 0.75), 1), 
                       limits = c(0,100)) + theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "", legend.text = element_text(size = 8))+ guides(fill = "none")

legend8 <- get_legend(fst2 + theme(legend.position = "right"))
row0 <- hanctop2 + ndpanctop2 + canctop2  + plot_layout(widths = c(3, 0.2,0.1))
row1 <- h2 + ndp2 + c2 + plot_layout(widths = c(3, 0.2,0.1))
row1add <- hancbot2 + ndpancbot2 + cancbot2 + plot_layout(widths = c(3, 0.2,0.1))
row2 <- fst2 + plot_spacer() + plot_layout(widths = c(3, 0.3))
row3 <- derived_black2 + plot_spacer() + plot_layout(widths = c(3, 0.3))
row4 <- derived_grey2 + plot_spacer() + plot_layout(widths = c(3, 0.3))
row5 <- gwas2 + plot_spacer() + plot_layout(widths = c(3, 0.3))
row6 <- g2 + plot_spacer() + plot_layout(widths = c(3, 0.3))
row7 <- ann2 + plot_spacer() + plot_layout(widths = c(3, 0.3))
fst_landscape <-landscape + plot_spacer() + plot_layout(widths = c(3, 0.3))

# Combine with relative heights
legend6 <- get_legend(g2 + theme(legend.position = "bottom"))
legend7 <- get_legend(ann2 + theme(legend.position = "bottom"))
combined_legends <- plot_grid(legend6, legend7, ncol = 1,rel_heights = c(1,1),align = "v")
top_legend_plot <- row0 / row1 / row1add + plot_layout(height = c(0.8, 1, 1.2), guides = "collect") &
  theme(legend.position = "top", plot.margin = unit(c(0, 0, 0, 0), "pt"))
bottom_legend_plots <- row2 / row3 / row4 / row5 / row6 / row7 &
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))
final_plot <- fst_landscape / top_legend_plot / bottom_legend_plots / combined_legends + plot_layout(heights = c(0.4, 3, 0.8, 0.5))
final_plot #8x8

modcheck <- h2 + ndp2 + plot_layout(widths = c(1, 0.2))

#png("heatmap_2peaks_with genes_andmore_fst0.8_reversed.png", width =8 , height = 8, units = "in", res = 600)
final_plot 
dev.off()
#png("heatmap_2peaks_with genes_super_stringent_reversed.png", width =8 , height = 8, units = "in", res = 600)
final_plot 
dev.off()

#just to remind myself: 
#the lighter the shade the smaller the value
#for fst & DAF --> the scale is fixed across all plots and from 0 to 1
#for gwas --> the scale is fixed across all plots and from 0 to 100 (chi-square)

# check scaffold 60 - no pattern
heatmap_pos2 <- heatmap %>% 
  filter(scaff == "scaffold_60") %>%
  filter(pos <= 413527) %>% 
  group_by(scaff) %>%
  filter(any(genotype != 3)) %>%
  mutate(pos_num = as.numeric(as.character(pos)), 
    pos = factor(pos, levels = sort(unique(as.numeric(as.character(pos)))))) %>%
  dplyr::select(-pos_num)

ggplot(heatmap_pos2) + 
  geom_tile(aes(x=pos, y=sample_code, fill = genotype)) +
  scale_fill_manual(values=c("0"="#08306B","1"="#1F78B4","2"="#A6CEE3", "3"="white"),
                    labels = c("0" = "Ancestral", "1" = "Heterozygous", "2" = "Derived", "3" = "Missing"), name = "Genotype") +
  theme(axis.text.y = element_text(size = 6), legend.position = "top",axis.text.x = element_text(size = 7,angle = 90, vjust = 0.5, hjust = 1)) 
# -----------------------------------------------------------------------------------------------------------------------------------------

# (7) PLOT PCA (Panel B) ----------------------------------------------------------------------------------------------------------------------------
if(PCA==TRUE) {
  # subsetting samples
  All <- newdat[ ,1:15]$sample_code
  
  # subsetting SNPs
  type <- "both_peaks"
  # SNP.IDs <- subset(snp, CHR!=18 & CHR!="Z")$Locus_Name #"outsidechr18"
  # SNP.IDs <- subset(snp, CHR==18)$Locus_Name #chr18
  # SNP.IDs <- subset(snp, CHR==18 & POS>=2090834 & POS<=3062778)$Locus_Name #first_peak n=62 smallerregion: 2587046
  # SNP.IDs <- subset(snp, CHR==18 & POS>=3085010 & POS<=3594648)$Locus_Name #sec_peak n=50 smallerregion: 3397329
   SNP.IDs <- subset(snp, CHR==18 & POS>=2587046 & POS<=3397329)$Locus_Name #expanded region from 2587046 to 2090834 and 3397329 to 3594648
  
  # PC analysis
  pca_All <- snpgdsPCA(SNP, num.thread=7, maf=0.05, missing.rate=0.05, snp.id=SNP.IDs, sample.id=All)	
 
  # PCs
  eigenvects <- pca_All$eigenvect
  Out_pcs_All <- cbind(pca_All$sample.id, eigenvects)
  Out_pcs_All <- data.frame(Out_pcs_All)
  colnames(Out_pcs_All) <- c("sampleID",paste("PC",seq(1:(ncol(Out_pcs_All)-1)),sep=""))
  
  write.table(Out_pcs_All, paste(path, "add_fresh_anc_out_tmp_PCs_All_pol_smallerregion_", type, ".txt", sep=""), append=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE, eol="\n")
 
  # Loadings
  corr_All <- snpgdsPCACorr(pca_All, SNP, eig.which=1:4, num.thread=6)
 
  Out_corr_All <- data.frame(corr_All$snpcorr)
  colnames(Out_corr_All) <- corr_All$snp.id
  Out_corr_All <- t(Out_corr_All)
  Out_corr_All <- data.frame(Out_corr_All)
  Out_corr_All$snpID <- rownames(Out_corr_All)
  rownames(Out_corr_All) <- NULL
  Out_corr_All$position <- unlist(lapply(strsplit(Out_corr_All$snpID,":",fixed=TRUE), `[[`, 2))
  Out_corr_All <- merge(Out_corr_All, snp[,c("Locus_Name","CHR","POS")], by.x="snpID", by.y="Locus_Name", sort=FALSE)
  
  write.table(Out_corr_All, paste(path, "add_fresh_anc_out_tmp_PCL_All_pol_smallerregion_", type, ".txt", sep=""), append=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE, eol="\n")
  
  # close DB -----
  closefn.gds(SNP)
}

#ALL (both east and south)
dat1 <- read.table("add_fresh_anc_out_tmp_PCs_All_pol_smallerregion_both_peaks.txt", header = TRUE)
df1 <- dat1 %>%
  mutate(identity = case_when(
    newhybrids_P0P1Allo.admix99 %in% c("P0car", "Bxcar") ~ "Carrion",
    newhybrids_P0P1Allo.admix99 %in% c("P1hood", "Bxhood") ~ "Hooded",
    newhybrids_P0P1Allo.admix99 %in% c("F1", "F2") ~ "Hybrid",TRUE ~ NA_character_)) %>%
  #filter(HI >=0.99 | HI<=0.01 | HI>=0.45 & HI<=0.55 | is.na(HI)) %>% # this produce the cleaned plot with 430 inds instead of all 580 inds
  filter(zone != "ancient") %>%
  mutate(identity = factor(identity, levels=c("Carrion", "Hybrid", "Hooded")))

#the sample we chose to long read sequence
to_label <- c("D_Ne_Y39" , "D_Rb_Y23", "D_Lo_C21", "DKoC53")
df_label <- df1 %>% filter(sample_code %in% to_label)

p1 <- ggplot(df1) +
  geom_point(aes(x=PC2, y=PC1, fill=identity, colour=identity), size=1, shape=21) +
  geom_point(data = df_label, aes(x = PC2, y = PC1, fill=identity), size=3, shape=c(24,22,24,22), colour="black") +
  #geom_text(data = df_label,aes(x = PC2, y = PC1, label = gsub("_", "", sampleID)), size = 3,vjust = 0.5, hjust = 1,  nudge_x = 0.015) +
  #geom_text(aes(x=PC1, y=PC2), label=df$sampleID, check_overlap=TRUE, size = 2, vjust=-1) +
  scale_fill_manual(name="", values=c("#E69F00","#009E73", "#56B4E9")) +
  scale_color_manual(name="", values=c("#E69F00","#009E73", "#56B4E9")) +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 4))  ) +
  scale_y_reverse() +
  theme_bw(base_size=10) +
  theme(legend.position="top") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x=element_text (size=10),
        axis.text.y=element_text (size=10, angle=90),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10)) #3x7
# -----------------------------------------------------------------------------------------------------------------------------------------

# (8)  A conclusion map for color distribution across ancient samples (Panel D) -----------------------------------------------------------
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(dplyr)
library(patchwork)

anc_color <- read_delim("/PATH/05_aDNA/pop_map_more.txt", delim="\t",col_names=TRUE) 
coldat <- anc_color %>%
  filter(type == "ancient" & color != "NA")

europe <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")
world <- ne_countries(scale = "medium", returnclass = "sf")
world_outline <- st_union(world)  
carrion <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(subregion == "Western Europe")
carrion2 <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(sovereignt == "Spain") 
carrion3 <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(sovereignt == "Portugal") 
carrion4 <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(sovereignt == "United Kingdom") 
carrion5 <- ne_states(country = "Czech Republic", returnclass = "sf") %>%
  filter(name %in% c("Karlovarský", "Plzeňský", "Ústecký", "Jihočeský", "Středočeský", "Prague"))
carrion6 <- ne_states(country = "Italy", returnclass = "sf") %>%
  filter(name %in% c("Alessandria", "Aoste", "Asti", "Belluno", "Bergamo", "Biella", "Bozen", "Brescia", "Como", "Cuneo", "Genoa", "Gorizia", "La Spezia", "Lecco",  
                     "Mantova", "Milano", "Monza e Brianza", "Novara", "Padova", "Pavia", "Piacenza", "Savona", "Sondrio", "Trento", "Treviso", "Trieste",  
                     "Turin", "Udine", "Varese", "Venezia", "Verbano-Cusio-Ossola", "Vercelli", "Verona", "Vicenza"))
hooded1 <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(sovereignt == c("Italy","Iran"))
hooded2 <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(subregion == "Northern Europe")
hooded3 <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(subregion %in% c("Eastern Europe","Western Asia","Southern Europe","Central Asia"))
hooded8 <- ne_states(country = "United Kingdom", returnclass = "sf") %>%
  filter(name %in% c("Antrim", "Ards", "Armagh", "Ballymena", "Ballymoney", "Banbridge", 
                     "Belfast", "Carrickfergus", "Castlereagh", "Coleraine", "Craigavon", "Derry", "Down", "Dungannon", "Fermanagh", "Larne", "Limavady", 
                     "Lisburn", "Magherafelt", "Mid Ulster", "Moyle", "Newry and Mourne", "Newtownabbey", "North Down", "Omagh", "Strabane"))
hooded9 <- ne_states(country = "Germany", returnclass = "sf") %>%
  filter(name %in% c("Berlin","Brandenburg","Mecklenburg-Vorpommern","Saxony","Saxony-Anhalt","Thuringia"))

ggplot() +
  geom_sf(data=world ,fill="white",color = NA) +
  geom_sf(data=hooded1, fill="#F5F5F5",color = NA) +
  geom_sf(data=hooded2, fill="#F5F5F5",color = NA) +
  geom_sf(data=hooded3, fill="#F5F5F5",color = NA) +
  geom_sf(data=carrion, fill="darkgrey",color = NA) +
  geom_sf(data=carrion2, fill="darkgrey",color = NA) +
  geom_sf(data=carrion3, fill="darkgrey",color = NA) +
  geom_sf(data=carrion4, fill="darkgrey",color = NA) +
  geom_sf(data=hooded8, fill="#F5F5F5",color = NA) +
  geom_sf(data=carrion5, fill="darkgrey",color = NA) +
  geom_sf(data=carrion6, fill="darkgrey",color = NA) +
  geom_sf(data=hooded9, fill="#F5F5F5",color = NA) +
  geom_sf(data = world_outline, fill = NA, color = "black", linewidth = 0.2) +
  geom_point(data=coldat %>% filter(period %in% c("100-1500")), aes(x = as.numeric(longitude), y = as.numeric(latitude),fill = color, shape = as.character(period)), color = "black",size = 4, position = position_jitter(width = 0.5, height = 0.5))+
  geom_point(data=coldat %>% filter(period=="1500-3000"| period=="more_than_10k"), aes(x = as.numeric(longitude), y = as.numeric(latitude),fill = color, shape = as.character(period)), color = "black",size = 4, position = position_jitter(width = 0.5, height = 0.5))+
  geom_point(data=coldat %>% filter(period=="5000-6000" | period=="more_than_16k"), aes(x = as.numeric(longitude), y = as.numeric(latitude),fill = color, shape = as.character(period)), color = "black",size = 4, position = position_jitter(width = 0.5, height = 0.5))+
  scale_fill_manual(values=c("grey"="grey","black"="black")) +
  scale_shape_manual(values = c("100-1500" = 21, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_16k"=25)) +
  annotation_scale(location="br", width_hint = 0.2) +
  #annotate("rect", xmin = -15, xmax = 65, ymin = 30, ymax = 65, fill = NA, color = "black", linewidth = 0.7) +
  coord_sf(xlim = c(-20, 120), ylim = c(36, 61), expand = FALSE) +
  theme_minimal() + labs(shape = "Time (years ago)") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),legend.position = "none") +
  guides(fill="none") #6x8
# -----------------------------------------------------------------------------------------------------------------------------------------

# Sniffles for structural analysis on long reads (supplementary) --------------------------------------------------------------------------
# Here I plot the long read result from Sniffles2
sv <- c("BND", "DEL", "DUP", "INS", "INV")
DKoC21_refDNeY39 <- c(9,249,2,245,0)
DKoC21_refD_RbY23 <- c(9,261,2,235,1)
DKoC53_refDNeY39 <- c(10,150,3,122,0)
DKoC53_refDRbY23 <- c(12,151,4,119,0)
DNeY39_refDKoC21 <- c(18,286,1,203,0)
DNeY39_refDKoC53 <- c(7,233,1,236,0)
DRbY23_refDKoC21 <- c(9,265,1,205,0)
DRbY23_refDKoC53 <- c(5,220,1,229,0)

df <- data.frame(sv,
                 DKoC21_refDNeY39 = c(9,249,2,245,0),
                 DKoC21_refDRbY23 = c(9,261,2,235,1),
                 DKoC53_refDNeY39 = c(10,150,3,122,0),
                 DKoC53_refDRbY23 = c(12,151,4,119,0),
                 DNeY39_refDKoC21 = c(18,286,1,203,0),
                 DNeY39_refDKoC53 = c(7,233,1,236,0),
                 DRbY23_refDKoC21 = c(9,265,1,205,0),
                 DRbY23_refDKoC53 = c(5,220,1,229,0))

df_long <- df %>% 
  pivot_longer(cols = -sv, names_to = "sample_ref",values_to = "count") %>%
  separate(sample_ref, into = c("sample", "ref"), sep = "_ref")
df_long$ref <- factor(df_long$ref,levels = c("DNeY39", "DRbY23", "DKoC21", "DKoC53"))

p2 <- ggplot(df_long, aes(x = sv, y = count, fill = sample)) +
  geom_col(position = "dodge") +
  facet_wrap(~ ref) +
  scale_fill_manual(name="",values = c("DNeY39" = "#1f78b4","DRbY23" = "#a6cee3", "DKoC21"= "#ff7f00", "DKoC53" = "#fdbf6f")) +
  theme_bw() + theme(legend.position = "top") +
  labs(x = "Structural variation",y = "Count",fill = "Sample") +
  theme(legend.text=element_text(size=8.5),legend.key.size = unit(0.4, "cm"))
# -----------------------------------------------------------------------------------------------------------------------------------------

