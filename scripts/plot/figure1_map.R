library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggspatial)
library(dplyr)
library(ggshadow)
library(ggbreak)
library(cowplot)
library(patchwork)
setwd("/PATH")

# Plot Map -------------------------------------------------------------------------------------------------------------------------
colors <- c(
  "SPA"  = "#7F3B08",
  "GB" = "#B35806",
  "EURw1_FR"  = "#E08214",
  "EURw2_DE"  = "#FEE0B6",
  "EURw3_BE" = "#FEE0B6",
  "EURw4_NL"=  "#FEE0B6",
  "IRE"= "#D6266C",
  "HZ"  = "#E78AC3",
  "EURn1_SE" =  "#A6CEE3",
  "EURe1_PL"=  "#2899B2",
  "RUS"= "#084081",
  "COR" = "#66C2A5",
  "EURs1_IT" = "#66C2A5",
  "EURs2_BU"  = "#117733",
  "EURs3_ISR"= "#7A2A7E", 
  "IRQ"= "#54278F")

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
  filter(sovereignt == c("Italy", "Iran"))
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
africa <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(sovereignt == "Morocco" | sovereignt == "Algeria" | sovereignt == "Libya" | sovereignt == "Egypt" | sovereignt == "Sudan" | sovereignt == "Tunisia")
world <- ne_countries(scale = "medium", returnclass = "sf")
others <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(sovereignt == "Afghanistan" | sovereignt == "Pakistan" | sovereignt == "India"| sovereignt == "China")
class(world)
loc <- read.table("pop_map_more.txt",header=TRUE) %>%
  mutate(country = factor(country, levels=c("SPA","GB","EURw1_FR","EURw2_DE","EURw3_BE","EURw4_NL","IRE","HZ","EURn1_SE","EURe1_PL","RUS","COR","EURs1_IT","EURs2_BU","EURs3_ISR","IRQ"))) %>% 
  distinct(longitude, latitude, age, .keep_all = TRUE) %>%
  mutate(color_group = ifelse(period == "present", as.character(country), "black")) 

f2 <- ggplot() + 
  geom_sf(data=world ,fill="white",color = "black", linewidth = 1) +
  geom_sf(data=africa,fill="white",color = NA) +
  geom_sf(data=others,fill="white",color = NA) +
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
  coord_sf(xlim = c(-15, 90), ylim = c(30, 72), expand = FALSE) +
  geom_point(data=loc, aes(x=longitude, y=latitude,  fill = country, color = color_group, shape=period, alpha=period, size=period,stroke = 1)) +
  geom_point(data = filter(loc, period == "present"),aes(x = longitude, y = latitude,color = color_group,shape = period,size = period),fill = "white",stroke = 1) +
  geom_point(data=loc%>%filter(period %in% c("100-1500")), aes(x=longitude, y=latitude, fill = country, color = color_group, shape=period, alpha=period, size=period)) +
  geom_point(data=loc%>%filter(period=="1500-3000"| period=="more_than_10k"), aes(x=longitude, y=latitude, fill = country, color = color_group, shape=period, alpha=period, size=period)) +
  geom_point(data=loc%>%filter(period=="5000-6000" | period=="more_than_20k"), aes(x=longitude, y=latitude, fill = country, color = color_group, shape=period, alpha=period, size=period)) +
  scale_shape_manual(values = c("100-1500" = 22, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_20k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 2, "1500-3000" = 2, "5000-6000"=2,"more_than_10k"=2,"more_than_20k"=2,"present"=1.5)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_20k"=1,"present"=1)) +
  #scale_color_manual(values = c("100-1500" = "black", "1500-3000" = "black", "5000-6000"="black","more_than_10k"="black","more_than_20k"="black","present"="white")) +
  #scale_size(breaks=c(1,5,25), range=c(3,6)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = c(colors, black = "black")) +
  labs(y = "",x="") + 
  annotation_scale(location="br", width_hint = 0.2) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(22,22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 1), size=c(2, 2, 2, 2, 2, 1.5))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=2)), color="none", size="none",alpha="none") +
  theme_bw(base_size=8) +
  theme( axis.text.x=element_text(size = 8),
         axis.title.x = element_text(size = 8),
         axis.ticks.x= element_line(color = "black"),
    legend.position = "right",
    legend.text=element_text(size=8),
    legend.title=element_text(size=9),
    plot.margin = margin(1, 1, 1, 1, "pt"))
#--------------------------------------------------------------------------------------------------------------------------------------------------------

# Age Plot-----------------------------------------------------------------------------------------------------------------------------------------------
f1 <- ggplot() +
  geom_point(data=loc, aes(x=longitude, y=age,  fill = country, color = color_group, shape=period, alpha=period, size=period)) +
  geom_point(data = filter(loc, period == "present"),aes(x = longitude, y = age,color = color_group,shape = period,size = period),fill = "white",stroke = 1) +
  geom_point(data=loc%>%filter(period %in% c("100-1500","undated_100","undated_2k","undated_5k")), aes(x=longitude, y=age, fill = country, color = color_group, shape=period, alpha=period, size=period)) +
  geom_point(data=loc%>%filter(period %in% c("1500-3000","5000-6000","more_than_10k","more_than_20k")), aes(x=longitude, y=age, fill = country, color = color_group, shape=period, alpha=period, size=period)) +
  scale_shape_manual(values = c("100-1500" = 22, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_20k"=25,"present"=21,"undated_100"=21,"undated_2k"=22,"undated_5k"=23)) +
  scale_size_manual(values = c("100-1500" = 2, "1500-3000" = 2, "5000-6000"=2,"more_than_10k"=2,"more_than_20k"=2,"present"=1.5,"undated_100"=2,"undated_2k"=2,"undated_5k"=2)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_20k"=1,"present"=1,"undated_100"=1,"undated_2k"=1,"undated_5k"=1)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 200),breaks = c(0, 1000, 2000, 6000, 16000, 24000)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = c(colors, black = "black")) +
  labs(y = "Age (BP)",x="Longitude") + 
    theme_minimal() +
  theme(panel.grid.major.x = element_blank(),  
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_text(size = 8),
        axis.text.y=element_text(size = 8),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.ticks.x= element_line(color = "black"),
        legend.position = "none", plot.margin = margin(1, 1, 1, 1, "pt"))
#--------------------------------------------------------------------------------------------------------------------------------------------------------

# PCA----------------------------------------------------------------------------------------------------------------------------------------------------
#install.packages("smartsnp")
library("smartsnp")
library("dplyr")
library("ggrepel",lib.loc ="/dss/dsshome1/lxc0E/di67kah/R")
library("ggplot2")
library(ggshadow, lib.loc="/dss/dsshome1/lxc0E/di67kah/R")
library(cowplot, lib.loc="/dss/dsshome1/lxc0E/di67kah/R")

#safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#E69F00", "#AA4499", "#6699CC", "#999933", "#882255", "#44AA99", "#661100", "#888888", "#000000", "#E69F00", "#56B4E9")
my_colors <- c(
  "E"  = "#7F3B08",
  "GB" = "#B35806",
  "F"  = "#E08214",
  "D"  = "#FEE0B6",
  "BE" = "#FEE0B6",
  "NL"=  "#FEE0B6",
  "HZ"  = "#E78AC3", #uli's sample: austria, eastgermany,northitaly
  "IRE"= "#D6266C",
  "S" =  "#A6CEE3",
  "PL"=  "#2899B2",
  "RUS"= "#084081",
  "COR" = "#66C2A5",
  "IT" = "#66C2A5",
  "B"  = "#117733",
  "ISR"= "#7A2A7E", 
  "IRQ"= "#54278F")

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/01_angsd_enrichedallfresh_rescaled")
dat2 <- read.table(gzfile("enrichedallfresh_rescaled_outlier_1111_hapconsensus_maxmis_q20.haplo.gz"), header=TRUE, na="N")

#basic check:0 are different from major allele, 2 are same as major allele, NA as 9
qc2 <- colSums(!is.na(dat2))/(nrow(dat2))
df2 <- dat2$major == dat2[4:ncol(dat2)]
df2 <- df2*1*2 
df2[is.na(df2)] <- 9
dim(df2) 
pop <- read.table("pop.txt")

#Read dataSNP and do some filtering
dataSNP2 <- read.table("enrichedallfresh_rescaled_outlier_1111_hapconsensus_maxmis_q20_smartsnp",sep="")
pop <- read.table("pop.txt")
summary2 <- data.frame(row=c(1:164), country=pop$V1, sample=pop$V2, age=pop$V3, period=pop$V4, miss=(1-qc2[4:length(qc2)])*100)
pass2 <- filter(summary2, miss<100) %>% 
  filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE")
selected2 <- dataSNP2[pass2$row]
#write.table(selected2, file = "enrichedallfresh_rescaled_outlier_1111_hapconsensus_maxmis_q20_miss100_smartsnp_withCOR", col.names = FALSE, row.names = FALSE)

#sample projection allows miss100
countmodern <- pass2 %>% filter(period == "present")
my_groups <- c(rep("A", nrow(pass2)-nrow(countmodern)),rep("M",nrow(countmodern)))
my_ancient <- c(1:(nrow(pass2)-nrow(countmodern)))
pcaR <- smart_pca(snp_data = "enrichedallfresh_rescaled_outlier_1111_hapconsensus_maxmis_q20_miss100_smartsnp_withCOR", sample_group = my_groups, sample_project = my_ancient, pc_axes=100,  missing_impute = "remove")
pcaR_eigen <- pcaR$pca.eigenvalues; dim(pcaR_eigen)
pcaR_load <- pcaR$pca.snp_loadings; dim(pcaR_load) 
pcaR_coord <- pcaR$pca.sample_coordinates; dim(pcaR_coord)
data <- pcaR$pca.sample_coordinates
data$sample <- pass2$sample
data$country <- pass2$country
data$age <- pass2$age
data$miss <- pass2$miss
data$period <- pass2$period
df1 <- data %>% mutate(country = factor(country, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","RUS","COR","IT","B","ISR","IRQ"))) %>%
  mutate(color_group = ifelse(period == "present", as.character(country), "black")) 

save_pc1 <- df1[, c("sample", "PC1", "age", "period")]
#write.table(save_pc1,file = "enrichedallfresh_rescaled_outlier_1111_hapconsensus_maxmis_q20_miss100_smartsnp_withCOR_pc1.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")

# Normalised eigenvalues
raw_eigenvals <- pcaR$pca.eigenvalues["observed eigenvalues", ]
pc1 <- round(raw_eigenvals[1] / sum(raw_eigenvals) * 100,1)
pc2 <- round(raw_eigenvals[2] / sum(raw_eigenvals) * 100,1)

### Plot pseudohaploid with projection: 801 SNPs for 1111set
proj1 <- ggplot() +
  geom_point(data=df1, aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=color_group)) +
  geom_point(data = filter(df1, period == "present"),aes(x = PC1, y = PC2,color = color_group,shape = period,size = period),fill = "white",stroke = 1) +
  geom_point(data=df1%>%filter(period %in% c("100-1500","undated_100","undated_2k","undated_5k")), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color= color_group)) +
  geom_point(data=df1%>%filter(period=="1500-3000"| period=="more_than_10k"), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color= color_group)) +
  geom_point(data=df1%>%filter(period=="5000-6000" | period=="more_than_20k"), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color= color_group)) +
  geom_point(data=df1%>%filter(country == "NL" | sample == "LSS003_TE" | sample == "HOC005_TE"), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color= color_group))  +
  scale_shape_manual(values = c("100-1500" = 22, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_20k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 2, "1500-3000" = 2, "5000-6000"=2,"more_than_10k"=2,"more_than_20k"=2,"present"=1.5)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_20k"=1,"present"=1)) +
  scale_fill_manual(values= my_colors) + coord_fixed() +
  scale_color_manual(values = c(my_colors, black = "black")) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(22, 22, 23, 24, 25,21),alpha = c(1, 1, 1, 1, 1, 1), size=c(2, 2, 2, 2, 2, 2))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=2)), color="none", size="none",alpha="none") +
  labs(x=paste("PC1 (",pc1,"%)"), y=paste("PC2 (",pc2,"%)")) +
  scale_x_reverse() +
  theme_bw(base_size=8) + theme(aspect.ratio = 1) +
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.title.x = element_text(size=9),axis.title.y = element_text(size=9),
        axis.text.x=element_text (size=8),axis.text.y=element_text (size=8),legend.position = "")

################################ NEUTRAL ######################
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/01_angsd_enrichedallfresh_rescaled")
dat2 <- read.table(gzfile("enrichedallfresh_rescaled_neutral_hapconsensus_maxmis_q20.haplo.gz"), header=TRUE, na="N")

#basic check
qc2 <- colSums(!is.na(dat2))/(nrow(dat2))
df2 <- dat2$major == dat2[4:ncol(dat2)]
df2 <- df2*1*2 
df2[is.na(df2)] <- 9
dim(df2) 
pop <- read.table("pop.txt")

#Read dataSNP and do some filtering
dataSNP2 <- read.table("enrichedallfresh_rescaled_neutral_hapconsensus_maxmis_q20_smartsnp",sep="")
summary2 <- data.frame(row=c(1:164), country=pop$V1, sample=pop$V2, age=pop$V3, period=pop$V4, miss=(1-qc2[4:length(qc2)])*100)
pass2 <- filter(summary2, miss<100) %>% 
  filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE")
selected2 <- dataSNP2[pass2$row]
#write.table(selected2, file = "enrichedallfresh_rescaled_neutral_hapconsensus_maxmis_q20_miss100_smartsnp_withCOR", col.names = FALSE, row.names = FALSE)

#sample projection
countmodern <- pass2 %>% filter(period == "present")
my_groups <- c(rep("A", nrow(pass2)-nrow(countmodern)),rep("M",nrow(countmodern)))
my_ancient <- c(1:(nrow(pass2)-nrow(countmodern)))
pcaR <- smart_pca(snp_data = "enrichedallfresh_rescaled_neutral_hapconsensus_maxmis_q20_miss100_smartsnp_withCOR", sample_group = my_groups, sample_project = my_ancient, pc_axes=100, missing_impute = "remove")
pcaR_eigen <- pcaR$pca.eigenvalues; dim(pcaR_eigen)
pcaR_load <- pcaR$pca.snp_loadings; dim(pcaR_load) 
pcaR_coord <- pcaR$pca.sample_coordinates; dim(pcaR_coord)
data <- pcaR$pca.sample_coordinates
data$sample <- pass2$sample
data$country <- pass2$country
data$age <- pass2$age
data$miss <- pass2$miss
data$period <- pass2$period
df2 <- data %>% mutate(country = factor(country, levels=c("E","GB","F","D","BE","NL","IRE","S","PL","RUS","COR","IT","B","ISR","IRQ"))) %>%
  mutate(color_group = ifelse(period == "present", as.character(country), "black")) 

# Normalised eigenvalues
raw_eigenvals <- pcaR$pca.eigenvalues["observed eigenvalues", ]
pc1 <- round(raw_eigenvals[1] / sum(raw_eigenvals) * 100,1)
pc2 <- round(raw_eigenvals[2] / sum(raw_eigenvals) * 100,1)

### Plot pseudohaploid with projection: 129657 SNPs remaining / 55607(seg72k) - exactly same pca
proj2 <- ggplot() +
  geom_point(data=df2, aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=color_group)) +
  geom_point(data = filter(df2, period == "present"),aes(x = PC1, y = PC2,color = color_group,shape = period,size = period),fill = "white",stroke = 1) +
  geom_point(data=df2%>%filter(period %in% c("100-1500","undated_100","undated_2k","undated_5k")), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=color_group)) +
  geom_point(data=df2%>%filter(period=="1500-3000"| period=="more_than_10k"), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=color_group)) +
  geom_point(data=df2%>%filter(period=="5000-6000" | period=="more_than_20k"), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=color_group)) +
  geom_point(data=df2%>%filter(country == "NL" | sample == "LSS003_TE" | sample == "HOC005_TE" ), aes(x=PC1, y=PC2, fill=country, shape=period, alpha=period, size=period, color=color_group))  +
  scale_shape_manual(values = c("100-1500" = 22, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_20k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 2, "1500-3000" = 2, "5000-6000"=2,"more_than_10k"=2,"more_than_20k"=2,"present"=1.5)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_20k"=1,"present"=1)) +
  scale_fill_manual(values= my_colors) +
  scale_color_manual(values = c(my_colors, black = "black")) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(21, 22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 1), size=c(2, 2, 2, 2, 2, 2))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=3)), color="none", size="none",alpha="none") +
  labs(x=paste("PC1 (",pc1,"%)"), y=paste("PC2 (",pc2,"%)")) + coord_fixed() +
  #scale_y_break(c(-480, -50)) +
  scale_x_reverse() + theme_bw(base_size=8) +  theme(aspect.ratio = 1) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),
        axis.title.x = element_text(size=9,margin = margin(t =0)),axis.title.y = element_text(size=9,margin = margin(t =0)),
        axis.text.x=element_text (size=8),axis.text.y=element_text (size=8),
        axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(),axis.line.x.top = element_blank(),legend.position = "") 
#--------------------------------------------------------------------------------------------------------------------------------------------------------

# Combine Plots----------------------------------------------------------------------------------------------------------------------------------------------------
r1 <- f1 +  plot_layout(widths = c(2))
r2 <- f2 + plot_layout(widths = c(2))
r3 <- proj2 +proj1 + plot_layout(widths = c(1,1)) #proj1=801 outlier snps, proj2=129657 neutral snps
fig1 <- r1 / r2 / r3  + 
  plot_layout(heights = c(0.9, 2.1, 2)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = c("A", "B", "C", "D")) & theme(plot.margin = unit(c(1, 1, 1, 1), "pt")) 
fig1 #7x8portrait
