library(devtools)
library(ggplot2)
library(dplyr)
library(ggshadow)

#safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#6699CC", "#999933", "#882255", "#44AA99", "#661100", "#888888", "#000000", "#E69F00", "#56B4E9")
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
pop<-read.table("pop.txt", header=FALSE, sep="", col.names = c("pop","sample","age","type"))


a1 <- as.matrix(read.table("enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20_dp3.cov"))
a2 <- as.matrix(read.table("enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20.cov"))
b1 <- as.matrix(read.table("enrichedallfresh_rescaled_neutral_geno_maxmis_q20_dp3.cov"))
b2 <- as.matrix(read.table("enrichedallfresh_rescaled_neutral_geno_maxmis_q20s.cov"))
C3 <- as.matrix(read.table("enrichedallfresh_rescaled_neutral_seg72k_geno_maxmis_q20_dp3.cov"))
C4 <- as.matrix(read.table("enrichedallfresh_rescaled_neutral_seg72k_geno_maxmis_q20.cov"))

dat1 <- data.frame(pop, C3) %>%
  mutate(pop = factor(pop, levels=c("E","GB","F","D","BE","IRE","NL","S","PL","ITA","B","ISR", "IRQ","RUS"))) %>%
  filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE", pop != "COR") %>%
  filter(sample != "BDG001_TE" & sample != "BDG006_TE") %>% #the UK samples are weird
  filter(pop != "ori")
e <- eigen(C3)
pc.percent <- e$values
print(pc.percent)

t1 <- ggplot() +
  geom_point(data=dat1, aes(x=V1, y=V2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat1%>%filter(type=="1500-3000"| type=="more_than_10k"), aes(x=V1, y=V2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat1%>%filter(type=="5000-6000" | type=="more_than_16k"), aes(x=V1, y=V2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat1%>%filter(pop == "NL" | sample == "LSS003_TE" | age == 5330), aes(x=V1, y=V2, fill=pop, shape=type, alpha=type, size=type, color=type))  +
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
        legend.position = "right")  #6x6 

t2 <- ggplot() +
  geom_point(data=dat1, aes(x=V2, y=V3, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat1%>%filter(type=="1500-3000"| type=="more_than_10k"), aes(x=V2, y=V3, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat1%>%filter(type=="5000-6000" | type=="more_than_16k"), aes(x=V2, y=V3, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat1%>%filter(pop == "NL" | sample == "LSS003_TE" | age == 5330), aes(x=V2, y=V3, fill=pop, shape=type, alpha=type, size=type, color=type))  +
  scale_shape_manual(values = c("100-1500" = 21, "1500-3000" = 22, "5000-6000"=23,"more_than_10k"=24,"more_than_16k"=25,"present"=21)) +
  scale_size_manual(values = c("100-1500" = 3, "1500-3000" = 3, "5000-6000"=3,"more_than_10k"=3,"more_than_16k"=3,"present"=2)) +
  scale_alpha_manual(values = c("100-1500" = 1, "1500-3000" = 1, "5000-6000"=1,"more_than_10k"=1,"more_than_16k"=1,"present"=1)) +
  scale_color_manual(values = c("100-1500" = "black", "1500-3000" = "black", "5000-6000"="black","more_than_10k"="black","more_than_16k"="black","present"="white")) +
  scale_fill_manual(values=my_colors ) +
  guides(shape=guide_legend("Time (years ago)", override.aes = list(shape = c(21, 22, 23, 24, 25, 21),alpha = c(1, 1, 1, 1, 1, 0.8), size=c(3, 3, 3, 3, 3, 2))), 
         fill=guide_legend("Locality", override.aes = list(shape = 21, color = "black", size=3)), color="none", size="none",alpha="none") +
  labs(x=paste("PC2 (",round(pc.percent[2],1),"%)",sep=""), y=paste("PC3 (",round(pc.percent[3],1),"%)",sep="")) +
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

plotall <- t1 + t2 + 
  plot_layout(widths = c(1, 1)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = c("A", "B")) & theme(plot.margin = unit(c(1, 1, 1, 1), "pt")) 
plotall #10x6

#-------------------------------------------------------------------------------------------------------
# wgs---------------------------------------------------------------------------------------------------
setwd("/PATH/05_aDNA/01_angsd_wgs_rescaled")
pop<-read.table("pop.txt", header=FALSE, sep="", col.names = c("pop","sample","age","type"))

a1 <- as.matrix(read.table("wgs_rescaled_all_geno_maxmis_q20.cov"))

dat1 <- data.frame(pop, C2) %>%
  mutate(pop = factor(pop, levels=c("E","GB","F","D","BE","IRE","NL","S","PL","ITA","B","ISR", "IRQ","RUS"))) %>%
  filter(sample != "BDG002_TE" & sample != "KCZ001_TE" & sample != "DSZ007_TE" & sample != "DVT017_TE", pop != "COR") %>%
  filter(sample != "BDG001_TE" & sample != "BDG006_TE") %>% #the UK samples are weird
  filter(pop != "ori")
e <- eigen(C2)
pc.percent <- e$values
print(pc.percent)

t1 <- ggplot() +
  geom_point(data=dat1, aes(x=V1, y=V2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat1%>%filter(type=="1500-3000"| type=="more_than_10k"), aes(x=V1, y=V2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat1%>%filter(type=="5000-6000" | type=="more_than_16k"), aes(x=V1, y=V2, fill=pop, shape=type, alpha=type, size=type, color=type)) +
  geom_point(data=dat1%>%filter(pop == "NL" | sample == "LSS003_TE" | age == 5330), aes(x=V1, y=V2, fill=pop, shape=type, alpha=type, size=type, color=type))  +
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
        legend.position = "right")  #6x6 


