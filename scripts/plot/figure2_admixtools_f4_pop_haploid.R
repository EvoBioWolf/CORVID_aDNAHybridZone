#this is a stricter version of the original admixtools_f4, which excluded 3 German samples: URB001, WMP003, WMP001 due to high missingness and alternate alleles found
#run on haploid datasets
library(usethis)
library(devtools) 
library(igraph)
devtools::install_github("uqrmaie1/admixtools")
library(admixtools)
library(tidyverse)
library(magrittr)
library(dplyr)
library(parallel)
library(ape)
library(ggtree)

setwd("PATH/02_results/admixtools")

# (0) Prepare admixtools input
# (1) f4-stats with SPA0
# (2) f4-stats with IRQ0 (Fig2 in main manuscript)
# (3) f4-ratio
# (4) qpadm for carrion 
# (5) qpadm for hooded
# (6) qpadm map
### supplementary material ###
### validation with wgs samples ###

# ------------------------------------------------------------------------------------------------------------------------------------------------------
# (0) Prepare admixtools input -----------------------------------------------------------------------------------------------------
neutral = "neutral_all_pop"
seg = "neutral_seg72k_pop"
pop = c("SPA1k","SPA0","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k", 
        "IRQ0","EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURse1k","EURe1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k","JACKDAW","MON")
extract_f2(neutral, "neutral_all_pop", auto_only = FALSE, pops=pop,overwrite=TRUE,maxmiss=1, blgsize = 50000) #allow missing data #160053 SNPs and 27 populations / 57475 for no missing
extract_f2(neutral, "neutral_trans_pop", auto_only = FALSE, pops=pop,overwrite=TRUE,maxmiss=1, transitions=FALSE, blgsize = 50000) #104164 SNPs and 27 populations / 38422 for no missing
extract_f2(seg, "neutral_seg72k_pop", auto_only = FALSE, pops=pop,overwrite=TRUE,maxmiss=1, blgsize = 50000) #65190 SNPs and 27 populations / 21705 for no missing

# this is not in manuscript #
# wgs_nooutlier = "wgs_nooutlier_pop"
# wgs_nooutlier_more = "wgs_nooutlier_more_pop"
# pop = c("SPA1k","SPA0","EURw0","EURc0","EURc2k","EURc5k","EURnc1k", "IRQ0","EURse0","EURe2k","EURse1k","EURse2k","EURse10k","EURse20k","JACKDAW","MON","AM")
# extract_f2(wgs_nooutlier, "wgs_nooutlier_pop", auto_only = FALSE, pops=pop,overwrite=FALSE,maxmiss=1, blgsize = 50000) #allow missing data #1828131 SNPs and 17 populations / nomis: 29889 SNPs and 17 populations
# extract_f2(wgs_nooutlier, "wgs_nooutlier_trans_pop", auto_only = FALSE, pops=pop,overwrite=FALSE,maxmiss=1, transitions=FALSE, blgsize = 50000) #104164 SNPs and 27 populations / nomis: 12855 SNPs and 17 populations
# extract_f2(wgs_nooutlier_more, "wgs_nooutlier_more_pop_nomis", auto_only = FALSE, pops=pop,overwrite=FALSE,maxmiss=0, blgsize = 50000) #allow missing data 21017638 SNPs and 17 populations / nomis: 6086141 SNPs 
# extract_f2(wgs_nooutlier_more, "wgs_nooutlier_more_trans_pop_nomis", auto_only = FALSE, pops=pop,overwrite=FALSE,maxmiss=0, transitions=FALSE, blgsize = 50000) #6768085 SNPs and 17 populations / nomis: 1925711 SNPs remain
#Max no. of block = 1x10e9 / 50000 = 22000

# ------------------------------------------------------------------------------------------------------------------------------------------------------
# (1) f4-stats, also known as qpdstat with f4mode=TRUE -----------------------------------------------------------------------------------------------------
# set4: compare EURw to all hoodies across time using SPA0 as pop1
# rmb f4(Pop1, Pop2, Pop3, Pop4) computes f4(Pop2,Pop4;Pop1,Pop3), so pop3 is the outgroup
popused = c("SPA1k","SPA0","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k",
            "IRQ0","EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURse1k","EURe1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k","JACKDAW","MON")
f2_blocks = f2_from_precomp("neutral_all_pop", afprod = TRUE, pops=popused) #Discarding 640 block(s) 
f2_blocks_trans = f2_from_precomp("neutral_trans_pop", afprod = TRUE, pops=popused) #Discarding 1328 block(s) 
f2_blocks_seg = f2_from_precomp("neutral_seg72k_pop", afprod = TRUE, pops=popused) #Discarding 2308  block(s) 

pop1 = c("SPA0")
pop2 = c("SPA1k","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k")
pop3 = c("MON")
pop4 = c("IRQ0","EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURse1k","EURe1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k","JACKDAW")

resd <- f4(f2_blocks, pop1, pop2, pop3, pop4)%>%
  filter(!pop4 %in% c("EURe100","EURse100"))
# write.table(resd%>%arrange(z), file="./f4est_pseudohaploid/neutral_all_pop1_spa.txt", quote=FALSE, row.names = FALSE, sep="\t")
resd_trans <- f4(f2_blocks_trans, pop1, pop2, pop3, pop4)%>%
  filter(!pop4 %in% c("EURe100","EURse100"))
# write.table(resd_trans%>%arrange(z), file="./f4est_pseudohaploid/neutral_trans_AM_pop1_spa.txt", quote=FALSE, row.names = FALSE, sep="\t")
resd_seg <- f4(f2_blocks_seg, pop1, pop2, pop3, pop4) %>%
  filter(!pop4 %in% c("EURe100","EURse100"))
# write.table(resd_seg%>%arrange(z), file="./f4est_pseudohaploid/neutral_outgroupasc_AM_pop1_spa.txt", quote=FALSE, row.names = FALSE, sep="\t")
#Z<3 for all jackdaw and spa1k with IRQ/EURs/EURne

f4st <- resd %>% group_by(pop2) %>% arrange(est, .by_group=TRUE) %>%
  mutate(pop4 = factor(pop4, levels = c("JACKDAW", "IRQ0","EURs0","RUS0","EURse0","EURsw0","EURn0","EURe0","EURse100","EURe100","EURse1k","EURe1k","NL1k","RUS2k","EURse2k","EURe2k","EURw2k_grey","EURse16k","EURse20k"))) %>%
  mutate(pop2 = factor(pop2, levels = c("EURw0","EURc0","SPA1k","EURnc1k","EURw2k","EURnw2k","EURc2k","EURc6k"))) 
d1 <- ggplot(f4st) +
  geom_errorbar(aes(y = pop4, x = est, xmin = est - se, xmax = est + se), width = 0.1, alpha=0.5) +
  geom_point(aes(y = pop4, x = est, fill=pop2,shape=pop2),size=3) +
  scale_fill_manual(values = c("SPA1k"="#7F3B08","EURnc1k"="#FEE0B6", "EURnw2k"="#B35806", "EURw2k"="#E08214", "EURc2k"="#FDB863","EURc6k"="#FDB863")) +
  scale_shape_manual(values = c("SPA1k"=22,"EURw0"=3,"EURc0"=4,"EURnc1k"=22, "EURnw2k"=22, "EURw2k"=22, "EURc2k"=22,"EURc6k"=23)) +
  labs(x="f4-estimates (all) ", y="Hooded crow") +
  guides(shape=guide_legend("Carrion crow", override.aes = list(shape = c(3,4,22,22,22,22,22,23), fill=c("white","white","#7F3B08","#FEE0B6","#E08214","#B35806","#FDB863","#FDB863"))), 
         fill="none") + theme_minimal() + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8)) #8x8

f4st_trans <- resd_trans %>% group_by(pop2) %>% arrange(est, .by_group=TRUE) %>%
  mutate(pop4 = factor(pop4, levels = c("JACKDAW", "IRQ0","EURs0","RUS0","EURse0","EURsw0","EURn0","EURe0","EURse100","EURe100","EURse1k","EURe1k","NL1k","RUS2k","EURse2k","EURe2k","EURw2k_grey","EURse16k","EURse20k"))) %>%
  mutate(pop2 = factor(pop2, levels = c("EURw0","EURc0","SPA1k","EURnc1k","EURw2k","EURnw2k","EURc2k","EURc6k"))) 
d2 <- ggplot(f4st_trans) +
  geom_errorbar(aes(y = pop4, x = est, xmin = est - se, xmax = est + se), width = 0.1, alpha=0.5) +
  geom_point(aes(y = pop4, x = est, fill=pop2,shape=pop2),size=3) +
  scale_fill_manual(values = c("SPA1k"="#7F3B08","EURnc1k"="#FEE0B6", "EURnw2k"="#B35806", "EURw2k"="#E08214", "EURc2k"="#FDB863","EURc6k"="#FDB863")) +
  scale_shape_manual(values = c("SPA1k"=22,"EURw0"=3,"EURc0"=4,"EURnc1k"=22, "EURnw2k"=22, "EURw2k"=22, "EURc2k"=22,"EURc6k"=23)) +
  labs(x="f4-estimates (transversion) ", y="Hooded crow") +
  guides(shape=guide_legend("Carrion crow", override.aes = list(shape = c(3,4,22,22,22,22,22,23), fill=c("white","white","#7F3B08","#FEE0B6","#E08214","#B35806","#FDB863","#FDB863"))), 
         fill="none") + theme_minimal() + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8)) #8x8

f4st_seg <- resd_seg %>% group_by(pop2) %>% arrange(est, .by_group=TRUE) %>%
  mutate(pop4 = factor(pop4, levels = c("JACKDAW", "IRQ0","EURs0","RUS0","EURse0","EURsw0","EURn0","EURe0","EURse100","EURe100","EURse1k","EURe1k","NL1k","RUS2k","EURse2k","EURe2k","EURw2k_grey","EURse16k","EURse20k"))) %>%
  mutate(pop2 = factor(pop2, levels = c("EURw0","EURc0","SPA1k","EURnc1k","EURw2k","EURnw2k","EURc2k","EURc6k"))) 
d3 <- ggplot(f4st_seg) +
  geom_errorbar(aes(y = pop4, x = est, xmin = est - se, xmax = est + se), width = 0.1, alpha=0.5) +
  geom_point(aes(y = pop4, x = est, fill=pop2,shape=pop2),size=3) +
  scale_fill_manual(values = c("SPA1k"="#7F3B08","EURnc1k"="#FEE0B6", "EURnw2k"="#B35806", "EURw2k"="#E08214", "EURc2k"="#FDB863","EURc6k"="#FDB863")) +
  scale_shape_manual(values = c("SPA1k"=22,"EURw0"=3,"EURc0"=4,"EURnc1k"=22, "EURnw2k"=22, "EURw2k"=22, "EURc2k"=22,"EURc6k"=23)) +
  labs(x="f4-estimates (outgroup-asc)", y="Hooded crow") +
  guides(shape=guide_legend("Carrion crow", override.aes = list(shape = c(3,4,22,22,22,22,22,23), fill=c("white","white","#7F3B08","#FEE0B6","#E08214","#B35806","#FDB863","#FDB863"))), 
         fill="none") + theme_minimal() + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))

checkA <- (d1+d2+d3) + plot_layout(width = c(1, 1,1)) + plot_layout(guides = "collect") + plot_annotation(tag_levels = c("A", "B","C")) & theme(plot.margin = unit(c(1, 1, 1, 1), "pt")) #9x6 

#Heatmap matrix
md2 <- ggplot(f4st_trans, aes(x = pop2, y = pop4, fill = est, size=est)) +
  geom_point(shape = 21) +
  scale_fill_viridis_c(option = "C",name="f4-estimates") +
  labs(x = "Carrion crow", y = "Hooded crow", size="", fill="") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  guides(size="none") +
  theme_minimal() + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))

#plot the simplified graph
tree_text <- "((SPA,Carrion_crow),(Hooded_crow, Outgroup));"
tree <- read.tree(text = tree_text)
graph1 <- ggtree(tree, layout = "rectangular") +
  geom_tiplab(angle = 0, hjust = 0.5, size = 4)  + coord_flip() + scale_x_reverse()

r1 <- graph1
r2 <- d2 + md2 + plot_layout(widths = c(1,1))
fig3a <- r1 / r2  + 
  plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = c("A", "B", "C")) & theme(plot.margin = unit(c(1, 1, 1, 1), "pt")) 
fig3a #6x8

# -----------------------------------------------------------------------------------------------------------------------------------------------------
# (2) f4-stats, also known as qpdstat with f4mode=TRUE -----------------------------------------------------------------------------------------------------
# set5: compare EURw to all hoodies across time using IRQ0 as pop1
pop1 = c("IRQ0")
pop2 = c("EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURe1k","EURse1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k")
pop3 = c("MON")
pop4 = c("JACKDAW","SPA1k","SPA0","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k") #jackdaw cannot be in pop2 as it is not closely-related to pop1

rese <- f4(f2_blocks, pop1, pop2, pop3, pop4)%>%
  filter(!pop2 %in% c("EURe100","EURse100"))
# write.table(rese%>%arrange(z), file="./f4est_pseudohaploid/neutral_all_pop1_irq.txt", quote=FALSE, row.names = FALSE, sep="\t")
rese_trans <- f4(f2_blocks_trans, pop1, pop2, pop3, pop4)%>%
  filter(!pop2 %in% c("EURe100","EURse100"))
# write.table(rese_trans%>%arrange(z), file="./f4est_pseudohaploid/neutral_trans_all_pop1_irq.txt", quote=FALSE, row.names = FALSE, sep="\t")
rese_seg <- f4(f2_blocks_seg, pop1, pop2, pop3, pop4)%>%
  filter(!pop2 %in% c("EURe100","EURse100"))
# write.table(rese_seg%>%arrange(z), file="./f4est_pseudohaploid/neutral_outgroupasc_pop1_irq.txt", quote=FALSE, row.names = FALSE, sep="\t")
#z<3 for all jackdaws and between 20k and SPA

f4st <- rese %>% group_by(pop2) %>% arrange(est, .by_group=TRUE) %>%
  mutate(pop2 = factor(pop2, levels = c("EURs0","RUS0","EURse0","EURsw0","EURn0","EURe0","EURse100","EURe100","EURe1k","EURse1k","NL1k","RUS2k","EURse2k","EURe2k","EURse16k","EURse20k"))) %>%
  mutate(pop4 = factor(pop4, levels = c("JACKDAW","SPA0","EURw0","EURc0","SPA1k","EURnc1k","EURw2k","EURnw2k","EURc2k","EURc6k"))) 
e1 <- ggplot(f4st%>%filter(pop4!="JACKDAW")) +
  geom_errorbar(aes(y = pop2, x = est, xmin = est - se, xmax = est + se), width = 0.1, alpha=0.5) +
  geom_point(aes(y = pop2, x = est, fill=pop4, shape=pop4), size = 3) +
  scale_fill_manual(values = c("SPA1k"="#7F3B08","EURnc1k"="#FEE0B6", "EURw2k"="#E08214", "EURnw2k"="#B35806", "EURc2k"="#FDB863","EURc6k"="#FDB863")) +
  scale_shape_manual(values = c("JACKDAW"=8,"SPA0"=10,"SPA1k"=22,"EURw0"=3,"EURc0"=4,"EURnc1k"=22, "EURnw2k"=22, "EURw2k"=22, "EURc2k"=22,"EURc6k"=23)) +
  labs(x="f4-estimates (all)", y="Hooded crow") +
  guides(shape=guide_legend("Carrion crow", override.aes = list(shape = c(10,3,4,22,22,22,22,22,23), fill=c("white","white","white","#7F3B08","#FEE0B6","#E08214","#B35806","#FDB863","#FDB863"))), 
         fill="none") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))

f4st_trans <- rese_trans %>% group_by(pop2) %>% arrange(est, .by_group=TRUE) %>%
  mutate(pop2 = factor(pop2, levels = c("EURs0","RUS0","EURse0","EURsw0","EURn0","EURe0","EURse100","EURe100","EURe1k","EURse1k","NL1k","RUS2k","EURse2k","EURe2k","EURse16k","EURse20k"))) %>%
  mutate(pop4 = factor(pop4, levels = c("JACKDAW","SPA0","EURw0","EURc0","SPA1k","EURnc1k","EURw2k","EURnw2k","EURc2k","EURc6k"))) 
e2 <- ggplot(f4st_trans %>% filter(pop4!="JACKDAW")) +
  geom_errorbar(aes(y = pop2, x = est, xmin = est - se, xmax = est + se), width = 0.1, alpha=0.5) +
  geom_point(aes(y = pop2, x = est, fill=pop4, shape=pop4), size = 3) +
  scale_fill_manual(values = c("SPA1k"="#7F3B08","EURnc1k"="#FEE0B6", "EURw2k"="#E08214", "EURnw2k"="#B35806", "EURc2k"="#FDB863","EURc6k"="#FDB863")) +
  scale_shape_manual(values = c("JACKDAW"=8,"SPA0"=10,"SPA1k"=22,"EURw0"=3,"EURc0"=4,"EURnc1k"=22, "EURnw2k"=22, "EURw2k"=22, "EURc2k"=22,"EURc6k"=23)) +
  labs(x="f4-estimates (transversion)", y="Hooded crow") +
  guides(shape=guide_legend("Carrion crow", override.aes = list(shape = c(10,3,4,22,22,22,22,22,23), fill=c("white","white","white","#7F3B08","#FEE0B6","#E08214","#B35806","#FDB863","#FDB863"))), 
         fill="none") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))

f4st_seg <- rese_seg %>% group_by(pop2) %>% arrange(est, .by_group=TRUE) %>%
  mutate(pop2 = factor(pop2, levels = c("EURs0","RUS0","EURse0","EURsw0","EURn0","EURe0","EURse100","EURe100","EURe1k","EURse1k","NL1k","RUS2k","EURse2k","EURe2k","EURse16k","EURse20k"))) %>%
  mutate(pop4 = factor(pop4, levels = c("JACKDAW","SPA0","EURw0","EURc0","SPA1k","EURnc1k","EURw2k","EURnw2k","EURc2k","EURc6k"))) 
e3 <- ggplot(f4st_seg%>% filter(pop4!="JACKDAW")) +
  geom_errorbar(aes(y = pop2, x = est, xmin = est - se, xmax = est + se), width = 0.1, alpha=0.5) +
  geom_point(aes(y = pop2, x = est, fill=pop4, shape=pop4), size = 3) +
  scale_fill_manual(values = c("SPA1k"="#7F3B08","EURnc1k"="#FEE0B6", "EURw2k"="#E08214", "EURnw2k"="#B35806", "EURc2k"="#FDB863","EURc6k"="#FDB863")) +
  scale_shape_manual(values = c("JACKDAW"=8,"SPA0"=10,"SPA1k"=22,"EURw0"=3,"EURc0"=4,"EURnc1k"=22, "EURnw2k"=22, "EURw2k"=22, "EURc2k"=22,"EURc6k"=23)) +
  labs(x="f4-estimates (outgroup-asc)", y="Hooded crow") +
  guides(shape=guide_legend("Carrion crow", override.aes = list(shape = c(10,3,4,22,22,22,22,22,23), fill=c("white","white","white","#7F3B08","#FEE0B6","#E08214","#B35806","#FDB863","#FDB863"))), 
         fill="none") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))

checkb <- (e1 + e2 +e3) + plot_layout(width = c(1,1,1)) + plot_layout(guides = "collect") + plot_annotation(tag_levels = c("A", "B","C")) & theme(plot.margin = unit(c(1, 1, 1, 1), "pt")) #9x6

# Plot matrix heatmap
me2 <- ggplot(f4st, aes(x = pop4, y = pop2, fill = est, size=est)) +
  geom_point(shape = 21) +
  scale_fill_viridis_c(option = "C",name="f4 estimates") +
  labs(x = "Carrion crow", y = "Hooded crow", size="", fill="") +
  guides(size="none") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))

#plot the simplified graph
tree_text <- "((IRQ,Hooded_crow),(Carrion_crow, Outgroup));"
tree <- read.tree(text = tree_text)
graph2 <- ggtree(tree, layout = "rectangular") +
  geom_tiplab(angle = 0, hjust = 0.5, size = 2)  + coord_flip() + scale_x_reverse()
r1 <- graph2
r2 <- e1 + me2 + plot_layout(widths = c(1,1))
fig3b <- r1 / r2  + 
  plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = c("A", "B", "C")) & theme(plot.margin = unit(c(1, 1, 1, 1), "pt")) 
fig3b #6x8

# e1 with horizontal
rese <- read.csv(file="./f4est_pseudohaploid/neutral_all_pop1_irq.txt", sep="\t")

f4st <- rese %>% group_by(pop2) %>% arrange(est, .by_group=TRUE) %>%
  mutate(pop2 = factor(pop2, levels = rev(c("EURs0","RUS0","EURse0","EURsw0","EURn0","EURe0","EURse100","EURe100","EURe1k","EURse1k","NL1k","RUS2k","EURse2k","EURe2k","EURse16k","EURse20k")))) %>%
  mutate(pop4 = factor(pop4, levels = c("JACKDAW","SPA0","EURw0","EURc0","SPA1k","EURnc1k","EURw2k","EURnw2k","EURc2k","EURc6k"))) 

e1 <- ggplot(f4st%>%filter(pop4!="JACKDAW")) +
  geom_errorbar(aes(y = pop2, x = est, xmin = est - se, xmax = est + se), width = 0.1, alpha=0.5) +
  geom_point(aes(y = pop2, x = est, fill=pop4, shape=pop4), size = 3) +
  scale_fill_manual(values = c("SPA1k"="#7F3B08","EURnc1k"="#FEE0B6", "EURw2k"="#E08214", "EURnw2k"="#B35806", "EURc2k"="#FEE0B6","EURc6k"="#FEE0B6")) +
  scale_shape_manual(values = c("JACKDAW"=8,"SPA0"=10,"SPA1k"=21,"EURw0"=3,"EURc0"=4,"EURnc1k"=21, "EURnw2k"=22, "EURw2k"=22, "EURc2k"=22,"EURc6k"=23)) +
  labs(x="f4-estimates (all)", y="Hooded crow")  + coord_flip() +
  guides(shape=guide_legend("Carrion\n  crow", override.aes = list(shape = c(10,3,4,21,21,22,22,22,23), fill=c("white","white","white","#7F3B08","#FEE0B6","#E08214","#B35806","#FEE0B6","#FEE0B6"))), 
         fill="none") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), axis.text.y = element_text(size = 9), legend.position = "top")
#4x6 landscape

# --------------------------------------------------------------------------------------------------------------------------------------------------------
# (3) f4-ratio -----------------------------------------------------------------------------------------------------------------------------------------------
#α=f4(PO, P5; PX, P2)/f4(PO, P5; P1, P2), PX is admixed, P1 and P2 are sources, PO and P5 are reference pop with no contribution to PX, rmb alpha refer to ancestry from P1!
popused = c("SPA1k","SPA0","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k",
            "IRQ0","EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURe1k","EURse1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k","MON") #EXCLUDE JACKDAW
f2_blocks = f2_from_precomp("neutral_all_pop", afprod = TRUE, pops=popused) #Discarding 428 (vs. 640 with JACKDAW) block(s) 
f2_blocks_trans = f2_from_precomp("neutral_trans_pop", afprod = TRUE, pops=popused) #Discarding 1040 (vs. 1328 with JACKDAW) block(s) 
f2_blocks_seg = f2_from_precomp("neutral_seg72k_pop", afprod = TRUE, pops=popused) #Discarding 1944 (vs. 2308 with JACKDAW)  block(s) 

#first i generate f4-ratio for carrion
PX = c("SPA1k","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k")
num <- f4(f2_blocks, "MON","IRQ0",PX, "EURse0") 
den <- f4(f2_blocks, "MON","IRQ0","SPA0", "EURse0")
alpha <- num$est / den$est
se_alpha <- abs(alpha) * sqrt( (num$se / num$est)^2 + (den$se / den$est)^2 )
z_alpha <- alpha / se_alpha
f4ratio_carrion <- num %>% mutate(alpha=alpha, se_alpha=se_alpha, z_alpha=z_alpha)

PX = c("SPA1k","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k")
num2 <- f4(f2_blocks_trans, "MON","IRQ0",PX, "EURse0") 
den <- f4(f2_blocks_trans, "MON","IRQ0","SPA0", "EURse0")
alpha <- num$est / den$est
se_alpha <- abs(alpha) * sqrt( (num$se / num$est)^2 + (den$se / den$est)^2 )
z_alpha <- alpha / se_alpha
f4ratio_carrion2 <- num %>% mutate(alpha=alpha, se_alpha=se_alpha, z_alpha=z_alpha)

PX = c("SPA1k","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k")
num2 <- f4(f2_blocks_seg, "MON","IRQ0",PX, "EURse0") 
den <- f4(f2_blocks_seg, "MON","IRQ0","SPA0", "EURse0")
alpha <- num$est / den$est
se_alpha <- abs(alpha) * sqrt( (num$se / num$est)^2 + (den$se / den$est)^2 )
z_alpha <- alpha / se_alpha
f4ratio_carrion3 <- num %>% mutate(alpha=alpha, se_alpha=se_alpha, z_alpha=z_alpha)

df_carrion <- f4ratio_carrion  %>%
  mutate(SPA0 = ifelse(alpha < 0, 0, alpha), EURse0 = ifelse(1-alpha > 1, 1, 1-alpha)) %>%
  pivot_longer(cols = c(SPA0, EURse0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(pop3, levels = c("SPA1k","EURw0","EURw2k","EURnw2k","EURnc1k","EURc0","EURc2k","EURc6k"))) %>% group_by(target) %>%
  mutate(Mixture = factor(Mixture, levels = c("SPA0","EURse0"))) %>% 
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

error <- df_carrion %>% filter(Mixture == "SPA0") %>% mutate(error_min = 1-ymax - se_alpha,error_max = 1-ymax + se_alpha)
plot_f4ratiocarrion <- ggplot(df_carrion, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  #geom_text(data = df_carrion,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.05),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + 
  theme_minimal() + labs(x = "Target group (all SNPs)", y = "f4-ratio") + theme(legend.position = "top")

df_carrion2 <- f4ratio_carrion2  %>%
  mutate(SPA0 = ifelse(alpha < 0, 0, alpha), EURse0 = ifelse(1-alpha > 1, 1, 1-alpha)) %>%
  pivot_longer(cols = c(SPA0, EURse0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(pop3, levels = c("SPA1k","EURw0","EURw2k","EURnw2k","EURnc1k","EURc0","EURc2k","EURc6k"))) %>% group_by(target) %>%
  mutate(Mixture = factor(Mixture, levels = c("SPA0","EURse0"))) %>% 
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

error <- df_carrion2 %>% filter(Mixture == "SPA0") %>% mutate(error_min = 1-ymax - se_alpha,error_max = 1-ymax + se_alpha)
plot_f4ratiocarrion2 <- ggplot(df_carrion2, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  #geom_text(data = df_carrion,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.05),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + 
  theme_minimal() + labs(x = "Target group (transversion SNPs)", y = "f4-ratio") + theme(legend.position = "top")

df_carrion3 <- f4ratio_carrion3  %>%
  mutate(SPA0 = ifelse(alpha < 0, 0, alpha), EURse0 = ifelse(1-alpha > 1, 1, 1-alpha)) %>%
  pivot_longer(cols = c(SPA0, EURse0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(pop3, levels = c("SPA1k","EURw0","EURw2k","EURnw2k","EURnc1k","EURc0","EURc2k","EURc6k"))) %>% group_by(target) %>%
  mutate(Mixture = factor(Mixture, levels = c("SPA0","EURse0"))) %>% 
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

error <- df_carrion3 %>% filter(Mixture == "SPA0") %>% mutate(error_min = 1-ymax - se_alpha,error_max = 1-ymax + se_alpha)
plot_f4ratiocarrion3 <- ggplot(df_carrion3, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  scale_fill_manual(values = c("black","grey")) + 
  theme_minimal() + labs(x = "Target group (outgroup ascertained SNPs)", y = "f4-ratio") + theme(legend.position = "top")

f4ratio_carrion_plot <- plot_f4ratiocarrion / plot_f4ratiocarrion2 / plot_f4ratiocarrion3  + 
  plot_layout(heights = c(1, 1,1)) + plot_annotation(tag_levels = c("A", "B")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))

#then i generate f4-ratio for hooded
PX = c("EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURe1k","EURse1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k")
num <- f4(f2_blocks, "MON","SPA0",PX, "EURc0")
den <- f4(f2_blocks, "MON","SPA0","IRQ0", "EURc0")
alpha <- num$est / den$est
se_alpha <- abs(alpha) * sqrt( (num$se / num$est)^2 + (den$se / den$est)^2 )
z_alpha <- alpha / se_alpha
f4ratio_hooded <- num %>% mutate(alpha=alpha, se_alpha=se_alpha, z_alpha=z_alpha)

num <- f4(f2_blocks_trans, "MON","SPA0",PX, "EURc0")
den <- f4(f2_blocks_trans, "MON","SPA0","IRQ0", "EURc0")
alpha <- num$est / den$est
se_alpha <- abs(alpha) * sqrt( (num$se / num$est)^2 + (den$se / den$est)^2 )
z_alpha <- alpha / se_alpha
f4ratio_hooded2 <- num %>% mutate(alpha=alpha, se_alpha=se_alpha, z_alpha=z_alpha)

num <- f4(f2_blocks_seg, "MON","SPA0",PX, "EURc0")
den <- f4(f2_blocks_seg, "MON","SPA0","IRQ0", "EURc0")
alpha <- num$est / den$est
se_alpha <- abs(alpha) * sqrt( (num$se / num$est)^2 + (den$se / den$est)^2 )
z_alpha <- alpha / se_alpha
f4ratio_hooded3 <- num %>% mutate(alpha=alpha, se_alpha=se_alpha, z_alpha=z_alpha)

df_hooded <- f4ratio_hooded  %>%
  mutate(IRQ0 = ifelse(alpha < 0, 0, alpha), EURc0 = ifelse(1-alpha > 1, 1, 1-alpha)) %>%
  pivot_longer(cols = c(IRQ0, EURc0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(pop3, levels = c("EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURe1k","EURse1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k"))) %>% group_by(target) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

error <- df_hooded %>%filter(Mixture == "EURc0") %>% mutate(error_min = pmax(0, 1 - ymax - se_alpha),error_max = pmin(1, 1 - ymax + se_alpha))
plot_f4ratiohooded <- ggplot(df_hooded, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  scale_fill_manual(values = c("black","grey")) + 
  theme_minimal() + labs(x = "Target group (all SNPs)", y = "f4-ratio") + theme(legend.position = "top",axis.text.x = element_text(size = 7,angle = 45, vjust =1, hjust = 1))

df_hooded2 <- f4ratio_hooded2  %>%
  mutate(IRQ0 = ifelse(alpha < 0, 0, alpha), EURc0 = ifelse(1-alpha > 1, 1, 1-alpha)) %>%
  pivot_longer(cols = c(IRQ0, EURc0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(pop3, levels = c("EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURe1k","EURse1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k"))) %>% group_by(target) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

error <- df_hooded2 %>% filter(Mixture == "EURc0") %>% mutate(error_min = 1-ymax - se_alpha,error_max = 1-ymax + se_alpha)
plot_f4ratiohooded2 <- ggplot(df_hooded2, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  scale_fill_manual(values = c("black","grey")) + 
  theme_minimal() + labs(x = "Target group (transversion SNPs)", y = "f4-ratio") + theme(legend.position = "top",axis.text.x = element_text(size = 7,angle = 45, vjust =1, hjust = 1))

df_hooded3 <- f4ratio_hooded3  %>%
  mutate(IRQ0 = ifelse(alpha < 0, 0, alpha), EURc0 = ifelse(1-alpha > 1, 1, 1-alpha)) %>%
  pivot_longer(cols = c(IRQ0, EURc0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(pop3, levels = c("EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURe1k","EURse1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k"))) %>% group_by(target) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

error <- df_hooded3 %>% filter(Mixture == "EURc0") %>% mutate(error_min = 1-ymax - se_alpha,error_max = 1-ymax + se_alpha)
plot_f4ratiohooded3 <- ggplot(df_hooded3, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  scale_fill_manual(values = c("black","grey")) + 
  theme_minimal() + labs(x = "Target group (outgroup ascertained SNPs)", y = "f4-ratio") + theme(legend.position = "top",axis.text.x = element_text(size = 7,angle = 45, vjust =1, hjust = 1))

f4ratio_hooded_plot <- plot_f4ratiohooded / plot_f4ratiohooded2 /plot_f4ratiohooded3  + 
  plot_layout(heights = c(1, 1,1)) + plot_annotation(tag_levels = c("A", "B")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt")) #8x6 

# -----------------------------------------------------------------------------------------------------------------------------------------------------
# (4) qpadm for carrion: how target relates to left and right ---------------------------------------------------------------------------------------------------------
# Note: if only 2 outgroup are used, then p-value cannot be generated for rank drop (how well the 2-sources model fit), We rely on pnested to show if 2-sources is significantly better than 1 source
# if >3 outgroups, then p-value is calculated. p_value>0.05 means the model fit well (i cannot reject 2-sources model). 
popused = c("SPA1k","SPA0","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k",
            "IRQ0","EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURse1k","EURe1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k","JACKDAW","MON")
f2_blocks = f2_from_precomp("neutral_all_pop", afprod = TRUE, pops=popused) #Discarding 640 block(s) 
f2_blocks_trans = f2_from_precomp("neutral_trans_pop", afprod = TRUE, pops=popused) #Discarding 1328 block(s) 
f2_blocks_seg = f2_from_precomp("neutral_seg72k_pop", afprod = TRUE, pops=popused) #Discarding 2308  block(s) 

popa = c("SPA1k","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k")
target=popa
adm1 <-data.frame(target=popa, SPA0=1, EURse0=0, se=0,z1=0,z2=0,p_rankdrop=0)
adm2 <-data.frame(target=popa, SPA0=1, EURse0=0, se=0,z1=0,z2=0,p_rankdrop=0)
adm3 <-data.frame(target=popa, SPA0=1, EURse0=0, se=0,z1=0,z2=0,p_rankdrop=0)
left =c("SPA0","EURse0")
right=c("MON","IRQ0","JACKDAW") #no AM for 3 outgroup
#qpadm(f2_blocks, left, right, "SPA1k")

for (i in 1:length(target)) { 
  adm1[i,2] <- qpadm(f2_blocks, left, right, target[i])$weights$weight[1]
  adm1[i,3] <- qpadm(f2_blocks, left, right, target[i])$weights$weight[2]
  adm1[i,4] <- qpadm(f2_blocks, left, right, target[i])$weights$se[1]
  adm1[i,5] <- qpadm(f2_blocks, left, right, target[i])$weights$z[1]
  adm1[i,6] <- qpadm(f2_blocks, left, right, target[i])$weights$z[2]
  adm1[i,7] <- qpadm(f2_blocks, left, right, target[i])$rankdrop$p[1]} #rankdrop$p_nested[1]
# write.table(adm1, file="./f4adm_pseudohaploid/neutral_all_target_carrion_3outgrp.txt", quote=FALSE, row.names = FALSE, sep="\t")
for (i in 1:length(target)) { 
  adm2[i,2] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$weight[1]
  adm2[i,3] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$weight[2]
  adm2[i,4] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$se[1]
  adm2[i,5] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$z[1]
  adm2[i,6] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$z[2]
  adm2[i,7] <- qpadm(f2_blocks_trans, left, right, target[i])$rankdrop$p[1]}
# write.table(adm2, file="./f4adm_pseudohaploid/neutral_trans_target_carrion_3outgrp.txt", quote=FALSE, row.names = FALSE, sep="\t")
for (i in 1:length(target)) { 
  adm3[i,2] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$weight[1]
  adm3[i,3] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$weight[2]
  adm3[i,4] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$se[1]
  adm3[i,5] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$z[1]
  adm3[i,6] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$z[2]
  adm3[i,7] <- qpadm(f2_blocks_seg, left, right, target[i])$rankdrop$p[1]}
# write.table(adm3, file="./f4adm_pseudohaploid/neutral_outgroupasc_target_carrion_3outgrp.txt", quote=FALSE, row.names = FALSE, sep="\t")

adm1<- read_table(file="./f4adm_pseudohaploid/neutral_all_target_carrion_3outgrp.txt",col_name=TRUE)
adm2 <- read_table(file="./f4adm_pseudohaploid/neutral_trans_target_carrion_3outgrp.txt",col_name=TRUE)
adm3 <- read_table(file="./f4adm_pseudohaploid/neutral_outgroupasc_target_carrion_3outgrp.txt",col_name=TRUE)
df_long <- adm1 %>%
  mutate(SPA0 = ifelse(SPA0 < 0, 0, SPA0), EURse0 = ifelse(EURse0 > 1, 1, EURse0)) %>%
  pivot_longer(cols = c(SPA0, EURse0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("SPA1k","EURw0","EURw2k","EURnw2k","EURnc1k","EURc0","EURc2k","EURc6k"))) %>% group_by(target) %>%
  mutate(Mixture = factor(Mixture, levels = c("SPA0","EURse0"))) %>% 
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

df_long2 <- adm2 %>%
  mutate(SPA0 = ifelse(SPA0 < 0, 0, SPA0), EURse0 = ifelse(EURse0 > 1, 1, EURse0)) %>%
  pivot_longer(cols = c(SPA0, EURse0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("SPA1k","EURw0","EURw2k","EURnw2k","EURnc1k","EURc0","EURc2k","EURc6k"))) %>% group_by(target) %>%
  mutate(Mixture = factor(Mixture, levels = c("SPA0","EURse0"))) %>% 
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

df_long3 <- adm3 %>%
  mutate(SPA0 = ifelse(SPA0 < 0, 0, SPA0), EURse0 = ifelse(EURse0 > 1, 1, EURse0)) %>%
  pivot_longer(cols = c(SPA0, EURse0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("SPA1k","EURw0","EURw2k","EURnw2k","EURnc1k","EURc0","EURc2k","EURc6k"))) %>% group_by(target) %>%
  mutate(Mixture = factor(Mixture, levels = c("SPA0","EURse0"))) %>% 
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

error <- df_long %>% filter(Mixture == "SPA0") %>% mutate(error_min = 1-ymax - se,error_max = 1-ymax + se)
x1 <- ggplot(df_long, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.05),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + 
  theme_minimal() + labs(x = "Target group (all SNPs)", y = "qpAdm") + theme(legend.position = "top")

error <- df_long2 %>% filter(Mixture == "SPA0") %>% mutate(error_min = 1-ymax - se,error_max = 1-ymax + se)
x1_trans <- ggplot(df_long2, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long2,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.05),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + 
  theme_minimal() + labs(x = "Target group (transversion SNPs)", y = "qpAdm") + theme(legend.position = "top")

error <- df_long3 %>% filter(Mixture == "SPA0") %>% mutate(error_min = 1-ymax - se,error_max = 1-ymax + se)
x1_outasc <- ggplot(df_long3, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long3,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.05),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + 
  theme_minimal() + labs(x = "Target group (outgroup ascertained SNPs)", y = "qpAdm") + theme(legend.position = "top")

#plot_grid(x1, x1_trans, x1_outasc , labels = c("A", "B","C"), ncol = 1)
adm1_plot <- x1 / x1_trans / x1_outasc  + 
  plot_layout(heights = c(1, 1, 1)) + plot_annotation(tag_levels = c("A", "B", "C")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", plot.margin = unit(c(0, 0, 0, 0), "pt")) #8.5x6 protrait

# -----------------------------------------------------------------------------------------------------------------------------------------------------
# (5) qpadm for hooded: how target relates to left and right ---------------------------------------------------------------------------------------------------------
popb = c("EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURe1k","EURse1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k")
target=popb
adm4 <-data.frame(target=popb, IRQ0=1, EURc0=0,se=0,z1=0,z2=0,p_rankdrop=0)
adm5 <-data.frame(target=popb, IRQ0=1, EURc0=0,se=0,z1=0,z2=0,p_rankdrop=0)
adm6 <-data.frame(target=popb, IRQ0=1, EURc0=0,se=0,z1=0,z2=0,p_rankdrop=0)
left =c("IRQ0","EURc0")
right=c("MON","SPA0","JACKDAW")
#qpadm(f2_blocks_seg, left, right, target)

for (i in 1:length(target)) { 
  adm4[i,2] <- qpadm(f2_blocks, left, right, target[i])$weights$weight[1]
  adm4[i,3] <- qpadm(f2_blocks, left, right, target[i])$weights$weight[2] 
  adm4[i,4] <- qpadm(f2_blocks, left, right, target[i])$weights$se[1]
  adm4[i,5] <- qpadm(f2_blocks, left, right, target[i])$weights$z[1]
  adm4[i,6] <- qpadm(f2_blocks, left, right, target[i])$weights$z[2]
  adm4[i,7] <- qpadm(f2_blocks, left, right, target[i])$rankdrop$p[1]
}
# write.table(adm4, file="./f4adm_pseudohaploid/neutral_all_target_hooded_3outgrp.txt", quote=FALSE, row.names = FALSE, sep="\t")
for (i in 1:length(target)) { 
  adm5[i,2] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$weight[1]
  adm5[i,3] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$weight[2]
  adm5[i,4] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$se[1]
  adm5[i,5] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$z[1]
  adm5[i,6] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$z[2]
  adm5[i,7] <- qpadm(f2_blocks_trans, left, right, target[i])$rankdrop$p[1]}
# write.table(adm5, file="./f4adm_pseudohaploid/neutral_trans_target_hooded_3outgrp.txt", quote=FALSE, row.names = FALSE, sep="\t")
for (i in 1:length(target)) { 
  adm6[i,2] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$weight[1]
  adm6[i,3] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$weight[2] 
  adm6[i,4] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$se[1]
  adm6[i,5] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$z[1]
  adm6[i,6] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$z[2]
  adm6[i,7] <- qpadm(f2_blocks_seg, left, right, target[i])$rankdrop$p[1]}
# write.table(adm6, file="./f4adm_pseudohaploid/neutral_outgroupasc_target_hooded_3outgrp.txt", quote=FALSE, row.names = FALSE, sep="\t")

adm4<- read_table(file="./f4adm_pseudohaploid/neutral_all_target_hooded_3outgrp.txt",col_name=TRUE)
adm5 <- read_table(file="./f4adm_pseudohaploid/neutral_trans_target_hooded_3outgrp.txt",col_name=TRUE)
adm6 <- read_table(file="./f4adm_pseudohaploid/neutral_outgroupasc_target_hooded_3outgrp.txt",col_name=TRUE)
df_long4 <- adm4 %>%
  mutate(IRQ0 = ifelse(IRQ0 < 0, 0, IRQ0), EURc0 = ifelse(EURc0 > 1, 1, EURc0)) %>%
  pivot_longer(cols = c(IRQ0, EURc0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("JACKDAW", "IRQ0","EURs0","EURsw0","EURn0","RUS0","RUS2k","NL1k","EURe0","EURe100","EURe1k","EURe2k","EURse0","EURse100","EURse1k","EURse2k","EURse16k","EURse20k"))) %>% group_by(target) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

df_long5 <- adm5 %>%
  mutate(IRQ0 = ifelse(IRQ0 < 0, 0, IRQ0), EURc0 = ifelse(EURc0 > 1, 1, EURc0)) %>%
  pivot_longer(cols = c(IRQ0, EURc0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("JACKDAW", "IRQ0","EURs0","EURsw0","EURn0","RUS0","RUS2k","NL1k","EURe0","EURe100","EURe1k","EURe2k","EURse0","EURse100","EURse1k","EURse2k","EURse16k","EURse20k"))) %>% group_by(target) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

df_long6 <- adm6 %>%
  mutate(IRQ0 = ifelse(IRQ0 < 0, 0, IRQ0), EURc0 = ifelse(EURc0 > 1, 1, EURc0)) %>%
  pivot_longer(cols = c(IRQ0, EURc0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("JACKDAW", "IRQ0","EURs0","EURsw0","EURn0","RUS0","RUS2k","NL1k","EURe0","EURe100","EURe1k","EURe2k","EURse0","EURse100","EURse1k","EURse2k","EURse16k","EURse20k"))) %>% group_by(target) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

error <- df_long4 %>%filter(Mixture == "EURc0") %>% mutate(error_min = pmax(0, 1 - ymax - se),error_max = pmin(1, 1 - ymax + se))
x2 <- ggplot(df_long4, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long4,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.01),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0.00, 1.00))+
  theme_minimal() + labs(x = "Target group (all SNPs)", y = "qpAdm") + theme(legend.position = "top",axis.text.x = element_text(size = 7,angle = 45, vjust = 1, hjust = 1))

error <- df_long5 %>% filter(Mixture == "EURc0") %>% mutate(error_min = 1-ymax - se,error_max = 1-ymax + se)
x2_trans <- ggplot(df_long5, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long5,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.01),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0, 1))+
  theme_minimal() + labs(x = "Target group (transversion SNPs)", y = "qpAdm") + theme(legend.position = "top",axis.text.x = element_text(size = 7,angle = 45, vjust = 1, hjust = 1))

error <- df_long6 %>% filter(Mixture == "EURc0") %>% mutate(error_min = 1-ymax - se,error_max = 1-ymax + se)
x2_outasc <- ggplot(df_long6, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long6,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.01),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0, 1))+
  theme_minimal() + labs(x = "Target group (outgroup ascertained SNPs)", y = "qpAdm") + theme(legend.position = "top",axis.text.x = element_text(size = 7,angle = 45, vjust =1, hjust = 1))

bul <- c("EURse0", "EURse1k","EURse2k","EURse10k","EURse20k")
error <- df_long4 %>% filter(target %in% bul) %>% filter(Mixture == "EURc0") %>% mutate(error_min = pmax(0, 1 - ymax - se),error_max = pmin(1, 1 - ymax + se))

x2_outasc_bul <- ggplot(df_long4 %>% filter(target %in% bul), aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long4%>% filter(target %in% bul),aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.01),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0, 1))+ labs(fill = "Ancestry") +
  theme_minimal() + labs(x = "Target group", y = "qpAdm") + theme(legend.position = "right",axis.text.x = element_text(size = 7,angle = 0, vjust =0, hjust = 0.5))

adm2_plot <- x2 / x2_trans / x2_outasc  + 
  plot_layout(heights = c(1, 1, 1)) + plot_annotation(tag_levels = c("A", "B", "C")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", plot.margin = unit(c(0, 0, 0, 0), "pt")) #8x10.5

#only bulgaria
df_long4bul <- adm4 %>%
  mutate(IRQ0 = ifelse(IRQ0 < 0, 0, IRQ0), EURc0 = ifelse(EURc0 > 1, 1, EURc0)) %>%
  pivot_longer(cols = c(IRQ0, EURc0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("JACKDAW", "IRQ0","EURs0","EURsw0","EURn0","RUS0","RUS2k","NL1k","EURe0","EURe100","EURe1k","EURe2k","EURse20k","EURse16k","EURse2k","EURse1k","EURse100","EURse0"))) %>% group_by(target) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

error <- df_long4bul %>% filter(Mixture == "EURc0") %>% 
  mutate(error_min = pmax(0, 1 - ymax - se),error_max = pmin(1, 1 - ymax + se)) %>% 
  filter(target %in% c("EURse0", "EURse1k", "EURse2k","EURse16k","EURse20k"))

x2bul <- ggplot(df_long4bul%>% filter(target %in% c("EURse0", "EURse1k", "EURse2k","EURse16k","EURse20k")), aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long4bul %>% filter(target %in% c("EURse0", "EURse1k", "EURse2k","EURse16k","EURse20k")),aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.01),vjust = 0, size = 4, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0.00, 1.04)) + 
  theme_minimal() + labs(x = "Target group (all SNPs)", y = "qpAdm") + theme(legend.position = "right",axis.text.x = element_text(size = 10,angle = 0, vjust = 1, hjust = 1))
#8x3 landscape
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# qpadm for test ---------------------------------------------------------------------------------------------------------
popc = c("SPA1k","EURw0","EURw2k","EURc0","EURc2k","EURc5k","EURnw2k","EURnc1k","EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURse1k","EURe2k","EURse2k","EURse10k","EURse20k","NL1k")
#fitting is bad throughout 
target=popc
adm7 <-data.frame(target=popc, SPA0=0, IRQ0=1,se=0,z1=0,z2=0,p_rankdrop=0)
adm8 <-data.frame(target=popc, SPA0=0, IRQ0=1,se=0,z1=0,z2=0,p_rankdrop=0)
adm9 <-data.frame(target=popc, SPA0=0, IRQ0=1,se=0,z1=0,z2=0,p_rankdrop=0)
left =c("SPA0","IRQ0")
right=c("MON","JACKDAW")
#qpadm(f2_blocks_seg, left, right, "SPA1k")

for (i in 1:length(target)) { 
  adm7[i,2] <- qpadm(f2_blocks, left, right, target[i])$weights$weight[1]
  adm7[i,3] <- qpadm(f2_blocks, left, right, target[i])$weights$weight[2] 
  adm7[i,4] <- qpadm(f2_blocks, left, right, target[i])$weights$se[1]
  adm7[i,5] <- qpadm(f2_blocks, left, right, target[i])$weights$z[1]
  adm7[i,6] <- qpadm(f2_blocks, left, right, target[i])$weights$z[2]
  adm7[i,7] <- qpadm(f2_blocks, left, right, target[i])$rankdrop$p_nested[1]} #p_nested for 2 outgroup
#write.table(adm7, file="./f4adm_pseudohaploid/neutral_all_target_SPAIRQ_2outgrp.txt", quote=FALSE, row.names = FALSE, sep="\t")
for (i in 1:length(target)) { 
  adm8[i,2] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$weight[1]
  adm8[i,3] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$weight[2]
  adm8[i,4] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$se[1]
  adm8[i,5] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$z[1]
  adm8[i,6] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$z[2]
  adm8[i,7] <- qpadm(f2_blocks_trans, left, right, target[i])$rankdrop$p_nested[1]}
#write.table(adm8, file="./f4adm_pseudohaploid/neutral_trans_target_SPAIRQ_2outgrp.txt", quote=FALSE, row.names = FALSE, sep="\t")
for (i in 1:length(target)) { 
  adm9[i,2] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$weight[1]
  adm9[i,3] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$weight[2] 
  adm9[i,4] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$se[1]
  adm9[i,5] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$z[1]
  adm9[i,6] <- qpadm(f2_blocks_seg, left, right, target[i])$weights$z[2]
  adm9[i,7] <- qpadm(f2_blocks_seg, left, right, target[i])$rankdrop$p_nested[1]}
#write.table(adm9, file="./f4adm_pseudohaploid/neutral_outgroupasc_target_SPAIRQ_2outgrp.txt", quote=FALSE, row.names = FALSE, sep="\t")
adm7<- read_table(file="./f4adm_pseudohaploid/neutral_all_target_SPAIRQ_3outgrp.txt",col_name=TRUE)
adm8 <- read_table(file="./f4adm_pseudohaploid/neutral_trans_target_SPAIRQ_3outgrp.txt",col_name=TRUE)
adm9 <- read_table(file="./f4adm_pseudohaploid/neutral_outgroupasc_target_SPAIRQ_3outgrp.txt",col_name=TRUE)

df_long7 <- adm7 %>%
  mutate(IRQ0 = ifelse(IRQ0 < 0, 0, IRQ0), SPA0 = ifelse(SPA0 > 1, 1, SPA0)) %>%
  pivot_longer(cols = c(IRQ0, SPA0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("SPA0","SPA1k","EURw0","EURw2k","EURnw2k","EURnc1k","EURc0","EURc2k","EURc5k","JACKDAW","EURs0","EURsw0","EURn0","RUS0","RUS2k","NL1k","EURe0","EURe100","EURe2k","EURse0","EURse100","EURse1k","EURse2k","EURse10k","EURse20k","IRQ0"))) %>% group_by(target) %>%
  mutate(Mixture = factor(Mixture, levels = c("SPA0","IRQ0"))) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

df_long8 <- adm8 %>%
  mutate(IRQ0 = ifelse(IRQ0 < 0, 0, IRQ0), SPA0 = ifelse(SPA0 > 1, 1, SPA0)) %>%
  pivot_longer(cols = c(IRQ0, SPA0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("SPA0","SPA1k","EURw0","EURw2k","EURnw2k","EURnc1k","EURc0","EURc2k","EURc5k","JACKDAW","EURs0","EURsw0","EURn0","RUS0","RUS2k","NL1k","EURe0","EURe100","EURe2k","EURse0","EURse100","EURse1k","EURse2k","EURse10k","EURse20k","IRQ0"))) %>% group_by(target) %>%
  mutate(Mixture = factor(Mixture, levels = c("SPA0","IRQ0"))) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

df_long9 <- adm9 %>%
  mutate(IRQ0 = ifelse(IRQ0 < 0, 0, IRQ0), SPA0 = ifelse(SPA0 > 1, 1, SPA0)) %>%
  pivot_longer(cols = c(IRQ0, SPA0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("SPA0","SPA1k","EURw0","EURw2k","EURnw2k","EURnc1k","EURc0","EURc2k","EURc5k","JACKDAW","EURs0","EURsw0","EURn0","RUS0","RUS2k","NL1k","EURe0","EURe100","EURe2k","EURse0","EURse100","EURse1k","EURse2k","EURse10k","EURse20k","IRQ0"))) %>% group_by(target) %>%
  mutate(Mixture = factor(Mixture, levels = c("SPA0","IRQ0"))) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

error <- df_long7 %>% filter(Mixture == "SPA0") %>% mutate(error_min = 1-ymax - se,error_max = 1-ymax + se)
z2 <- ggplot(df_long7, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") +
  #geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long7,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.0),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0.00, 1.00))+
  theme_minimal() + labs(x = "Target group (all SNPs)", y = "qpAdm") + theme(legend.position = "top",axis.text.x = element_text(size = 7,angle = 45, vjust = 1, hjust = 1))

error <- df_long8 %>% filter(Mixture == "SPA0") %>% mutate(error_min = 1-ymax - se,error_max = 1-ymax + se)
z2_trans <- ggplot(df_long8, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  #geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long8,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.0),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0, 1))+
  theme_minimal() + labs(x = "Target group (transversion SNPs)", y = "qpAdm") + theme(legend.position = "top",axis.text.x = element_text(size = 7,angle = 45, vjust = 1, hjust = 1))

error <- df_long9 %>% filter(Mixture == "SPA0") %>% mutate(error_min = 1-ymax - se,error_max = 1-ymax + se)
z2_outasc <- ggplot(df_long9, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  #geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long9,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.0),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0, 1))+
  theme_minimal() + labs(x = "Target group (outgroup ascertained SNPs)", y = "qpAdm") + theme(legend.position = "top",axis.text.x = element_text(size = 7,angle = 45, vjust =1, hjust = 1))

adm3_plot <- z2 / z2_trans / z2_outasc  + 
  plot_layout(heights = c(1, 1, 1)) + plot_annotation(tag_levels = c("A", "B", "C")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", plot.margin = unit(c(0, 0, 0, 0), "pt"))

# -----------------------------------------------------------------------------------------------------------------------------------------------------
# (6) qpadm map -------------------------------------------------------------------------------------------------------------------------------------------
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
#devtools::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)
library(ggplot2)
library(scatterpie) 
library(dplyr)
library(patchwork)
#library("pastclim", lib.loc="/dss/dsshome1/lxc0E/di67kah/R")
#library(terra)

locality_carrion <- data.frame(target=c("SPA1k","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k"),
                               age = c(1000,0,0,2000,2000,6000,2000,1000),
                               longitude=c(1.672744, 3.105836, 9.198, 3.908588, 9.1571, 9.001, 0.063, 5.783),
                               latitude=c(41.209305, 45.777968, 47.676, 43.566198, 49.2285, 48.8859, 50.8404,50.133))

adm3ed <- read_table(file="./f4adm_pseudohaploid/neutral_all_target_carrion_3outgrp.txt",col_name=TRUE) %>%
  rename(carrion = SPA0, hooded =  EURse0) %>%
  merge(locality_carrion, by="target")

locality_hooded <- data.frame(target=c("EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURe1k","EURse1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k"),
                              age = c(0,0,0,0,0,0,100,100,2000,1000,1000,2000,2000,16000,20000,1000),
                              longitude=c(34.74,49.669,12.46,17.643,21.012,23.326,16.013303,26.503967,49.061816,16.013303,26.819494,19.551111,25.610631,24.895556,24.895556,5.9145),
                              latitude=c(32.067,58.597,41.91,59.858,52.229,42.679,50.687106,43.670692,55.877009,50.687106,43.164038,50.464444,43.219839,43.236389,43.236389,52.558))

adm6ed <- read_table(file="./f4adm_pseudohaploid/neutral_all_target_hooded_3outgrp.txt",col_name=TRUE)  %>%
  rename(carrion = EURc0, hooded =  IRQ0) %>%
  merge(locality_hooded, by="target")

admall <- rbind(adm3ed,adm6ed) %>%
  bind_rows(tibble(target="SPA0", carrion=1, hooded=0,se=0,z1=0,z2=0,p_rankdrop=NA, age=0,latitude=42.588,longitude=-5.477),
            tibble(target="SPA1k", carrion=1, hooded=0,se=0,z1=0,z2=0,p_rankdrop=NA, age=1000,latitude=42.588,longitude=-5.477),
            tibble(target="SPA2k", carrion=1, hooded=0,se=0,z1=0,z2=0,p_rankdrop=NA, age=2000,latitude=42.588,longitude=-5.477),
            tibble(target="SPA6k", carrion=1, hooded=0,se=0,z1=0,z2=0,p_rankdrop=NA, age=6000,latitude=42.588,longitude=-5.477),
            tibble(target="SPA16k", carrion=1, hooded=0,se=0,z1=0,z2=0,p_rankdrop=NA, age=16000,latitude=42.588,longitude=-5.477),
            tibble(target="SPA20k", carrion=1, hooded=0,se=0,z1=0,z2=0,p_rankdrop=NA, age=20000,latitude=42.588,longitude=-5.477),
            tibble(target="IRQ0", carrion=0, hooded=1,se=0,z1=0,z2=0,p_rankdrop=NA, age=0,latitude=33.94,longitude=44.3614),
            tibble(target="IRQ1k", carrion=0, hooded=1,se=0,z1=0,z2=0,p_rankdrop=NA, age=1000,latitude=33.94,longitude=44.3614),
            tibble(target="IRQ2k", carrion=0, hooded=1,se=0,z1=0,z2=0,p_rankdrop=NA, age=2000,latitude=33.94,longitude=44.3614),
            tibble(target="IRQ6k", carrion=0, hooded=1,se=0,z1=0,z2=0,p_rankdrop=NA, age=6000,latitude=33.94,longitude=44.3614),
            tibble(target="IRQ16k", carrion=0, hooded=1,se=0,z1=0,z2=0,p_rankdrop=NA, age=16000,latitude=33.94,longitude=44.3614),
            tibble(target="IRQ20k", carrion=0, hooded=1,se=0,z1=0,z2=0,p_rankdrop=NA, age=20000,latitude=33.94,longitude=44.3614)) %>% 
  mutate(longitude = as.numeric(longitude),latitude = as.numeric(latitude))

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

hooded4 <- ne_states(country = "France", returnclass = "sf") %>% filter(region == "Corse")
hooded8 <- ne_states(country = "United Kingdom", returnclass = "sf") %>%
  filter(name %in% c("Antrim", "Ards", "Armagh", "Ballymena", "Ballymoney", "Banbridge", 
                     "Belfast", "Carrickfergus", "Castlereagh", "Coleraine", "Craigavon", "Derry", "Down", "Dungannon", "Fermanagh", "Larne", "Limavady", 
                     "Lisburn", "Magherafelt", "Mid Ulster", "Moyle", "Newry and Mourne", "Newtownabbey", "North Down", "Omagh", "Strabane"))
hooded9 <- ne_states(country = "Germany", returnclass = "sf") %>%
  filter(name %in% c("Berlin","Brandenburg","Mecklenburg-Vorpommern","Saxony","Saxony-Anhalt","Thuringia"))

plot_admixture_map <- function(df, title) {
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
    geom_sf(data=hooded4, fill="#F5F5F5",color = NA) +
    geom_sf(data=hooded9, fill="#F5F5F5",color = NA) +
    geom_sf(data = world_outline, fill = NA, color = "black", linewidth = 0.2) +
    geom_scatterpie(data = df,aes(x = longitude, y = latitude, group = factor(target)),cols = c("carrion", "hooded"),color = "black", 
                    pie_scale = 1.5) +
    geom_point(aes(x = -5.477, y = 42.588),shape = 22,size = 4,fill = "black", color = "black")+
    geom_point(aes(x = 44.3614, y = 33.94),shape = 22,size = 4,fill = "white", color = "black")+
    annotate("rect", xmin = -15, xmax = 55, ymin = 30, ymax = 65,fill = NA, color = "black", linewidth = 0.7) +
    coord_sf(xlim = c(-15, 55), ylim = c(30, 65), expand = FALSE) +
    scale_fill_manual(values = c("carrion" = "black", "hooded" = "white")) +
    theme_minimal()  +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none") }

adm_plot <- admall %>% filter(age == "0")
map0 <- plot_admixture_map(adm_plot, "Present")
adm_plot <- admall %>% filter(age == "1000")
map1 <- plot_admixture_map(adm_plot, "1 kya")
adm_plot <- admall %>% filter(age == "2000")
map2 <- plot_admixture_map(adm_plot, "2 kya")
adm_plot <- admall %>% filter(age == "6000")
map5 <- plot_admixture_map(adm_plot, "6 kya")
adm_plot <- admall %>% filter(age == "16000")
map10 <- plot_admixture_map(adm_plot, "16 kya")
adm_plot <- admall %>% filter(age == "20000")
map20 <- plot_admixture_map(adm_plot, "20 kya")

admmap <- map20 + map10 + map5 + map2 + map1 + map0 +
  plot_layout(ncol=3) &
  theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "pt")) #3x6, edit size of piecharts later

### supplementary material ###
# ------------------------------------------------------------------------------------------------------------------------------------------------------
# f4-stats, also known as qpdstat with f4mode=TRUE -----------------------------------------------------------------------------------------------------
# set6: compare EURc0 to all hoodies across time using SPA0 as pop1
results_list <- list(
# between each pair of IRQ and another hooded crow pop, another pop always has higher allele sharing with carrion crows (-3<Z<3)
f4(f2_blocks_seg, "SPA0", "EURc0", "IRQ0", "EURs0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "IRQ0", "EURse0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "IRQ0", "EURsw0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "IRQ0", "RUS0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "IRQ0", "EURn0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "IRQ0", "EURe0"), 

# between each pair of EURs and another hooded crow pop, another pop always has slightly significantly higher allele sharing with carrion crows except in RUS, Sweden and Poland
f4(f2_blocks_seg, "SPA0", "EURc0", "EURs0", "EURse0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "EURs0", "EURsw0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "EURs0", "RUS0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "EURs0", "EURn0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "EURs0", "EURe0"), 

# between each pair of hooded crow pop below, no significant difference in allele sharing with carrion crows (-3<Z<3)
f4(f2_blocks_seg, "SPA0", "EURc0", "EURse0", "EURsw0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "EURse0", "RUS0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "EURse0", "EURn0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "EURse0", "EURe0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "EURsw0", "RUS0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "EURsw0", "EURn0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "EURsw0", "EURe0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "RUS0", "EURn0"), 
f4(f2_blocks_seg, "SPA0", "EURc0", "RUS0", "EURe0"),
f4(f2_blocks_seg, "SPA0", "EURc0", "EURn0", "EURe0"), 

# between each pair of Bulgarian pop across time, there is no significant difference in allele sharing with carrion crows (-3<Z<3)
f4(f2_blocks_seg, "SPA0", "EURc0", "EURse20k","EURse0"),
f4(f2_blocks_seg, "SPA0", "EURc0", "EURse16k","EURse0"),
f4(f2_blocks_seg, "SPA0", "EURc0", "EURse2k","EURse0"),
f4(f2_blocks_seg, "SPA0", "EURc0", "EURse1k","EURse0")
)
results_hooded <- do.call(rbind, results_list) %>%
  mutate(pop3 = factor(pop3, levels = c("IRQ0","EURs0","RUS0","EURse0","EURsw0","EURn0","EURe0","EURse1k","EURse2k","EURse16k","EURse20k")),
         pop3_pop4 = paste(pop3, pop4, sep = " & "),
         pop3_pop4 = factor(pop3_pop4, levels = unique(pop3_pop4)),
         sig = z > 3 ) 

#plot results_hooded
hoodies <- ggplot(results_hooded) +
  geom_errorbar(aes(y = pop3_pop4, x = est, xmin = est - se, xmax = est + se), width = 0.1, alpha=0.5) +
  geom_point(aes(y = pop3_pop4, x = est, color = sig),size = 2) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  labs(x="f4-estimates (all)", y="Genetic affinity to carrion crows", color = "Z-score > 3") +
   theme_minimal() + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))

hoodies2 <- ggplot(results_hooded) +
  geom_errorbar(aes(y = pop3_pop4, x = est, xmin = est - se, xmax = est + se), width = 0.1, alpha=0.5) +
  geom_point(aes(y = pop3_pop4, x = est, color = sig),size = 2) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  labs(x="f4-estimates (transversion)", y="Genetic affinity to carrion crows", color = "Z-score > 3") +
  theme_minimal() + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))

hoodies3 <- ggplot(results_hooded) +
  geom_errorbar(aes(y = pop3_pop4, x = est, xmin = est - se, xmax = est + se), width = 0.1, alpha=0.5) +
  geom_point(aes(y = pop3_pop4, x = est, color = sig),size = 2) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  labs(x="f4-estimates (outgroup ascertained)", y="Genetic affinity to carrion crows", color = "Z-score > 3") +
  theme_minimal() + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))

f4stat_hooded <- hoodies + hoodies2 + hoodies3  + plot_layout(guides = "collect") +
  plot_layout(width = c(1, 1,1)) + plot_annotation(tag_levels = c("A", "B")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt")) #10x6 

# ------------------------------------------------------------------------------------------------------------------------------------------------------
# f4-stats, also known as qpdstat with f4mode=TRUE -----------------------------------------------------------------------------------------------------
# set7: compare EURse0 to all carrion across time using IRQ0 as pop1

results_list2 <- list(
  # between each pair of SPA and another carrion crow pop, another pop always has higher allele sharing with hooded crows (-3<Z<3)
  f4(f2_blocks, "IRQ0", "EURse0", "SPA0", "EURw0"),
  f4(f2_blocks, "IRQ0", "EURse0", "SPA0", "EURc0"),
  f4(f2_blocks, "IRQ0", "EURse0", "SPA0", "EURnc1k"),
  f4(f2_blocks, "IRQ0", "EURse0", "SPA0", "EURnw2k"),
  
  # between each pair of carrion crow pop below, no significant difference in allele sharing with hooded crows except france with germany and belgium
  f4(f2_blocks, "IRQ0", "EURse0", "EURw0", "EURc0"),
  f4(f2_blocks, "IRQ0", "EURse0", "EURw0", "EURnc1k"),
  f4(f2_blocks, "IRQ0", "EURse0", "EURw0", "EURnw2k"),
  f4(f2_blocks, "IRQ0", "EURse0", "EURc0", "EURnc1k"),
  f4(f2_blocks, "IRQ0", "EURse0", "EURc0", "EURnw2k"),
  f4(f2_blocks, "IRQ0", "EURse0", "EURnc1k", "EURnw2k"),
  
  # between each pair of west european pop across time, there is no significant difference in allele sharing with carrion crows (-3<Z<3)
  f4(f2_blocks, "IRQ0", "EURse0", "EURw2k","EURw0"),
  f4(f2_blocks, "IRQ0", "EURse0", "EURc2k", "EURc0"),
  f4(f2_blocks, "IRQ0", "EURse0", "EURc6k", "EURc0"),
  f4(f2_blocks, "IRQ0", "EURse0", "EURc6k", "EURc2k"))
results_carrion <- do.call(rbind, results_list2) %>%
  mutate(pop3 = factor(pop3, levels = c("SPA0","EURw0","EURc0","EURnc1k","EURnw2k","EURw2k","EURc2k","EURc6k")),
         pop3_pop4 = paste(pop3, pop4, sep = " & "),
         pop3_pop4 = factor(pop3_pop4, levels = unique(pop3_pop4)),
         sig = abs(z) > 3 ) 

#plot results_carrion
carrion <- ggplot(results_carrion) +
  geom_errorbar(aes(y = pop3_pop4, x = est, xmin = est - se, xmax = est + se), width = 0.1, alpha=0.5) +
  geom_point(aes(y = pop3_pop4, x = est, color = sig),size = 2) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  labs(x="f4-estimates (all)", y="Genetic affinity to hooded crows", color = "Z-score > 3") +
  theme_minimal() + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))

carrion2 <- ggplot(results_carrion) +
  geom_errorbar(aes(y = pop3_pop4, x = est, xmin = est - se, xmax = est + se), width = 0.1, alpha=0.5) +
  geom_point(aes(y = pop3_pop4, x = est, color = sig),size = 2) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  labs(x="f4-estimates (transversion)", y="Genetic affinity to hooded crows", color = "Z-score > 3") +
  theme_minimal() + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))

carrion3 <- ggplot(results_carrion) +
  geom_errorbar(aes(y = pop3_pop4, x = est, xmin = est - se, xmax = est + se), width = 0.1, alpha=0.5) +
  geom_point(aes(y = pop3_pop4, x = est, color = sig),size = 2) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  labs(x="f4-estimates (outgroup ascertained)", y="Genetic affinity to hooded crows", color = "Z-score > 3") +
  theme_minimal() + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))

f4stat_carrion <- carrion + carrion2 + carrion3  + plot_layout(guides = "collect") +
  plot_layout(width = c(1, 1,1)) + plot_annotation(tag_levels = c("A", "B")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt")) #10x6 

# --------------------------------------------------------------------------------------------------------------------------------------------------------
# f3-stats -----------------------------------------------------------------------------------------------------------------------------------------------
pop1 = c("SPA0")
pop2 = c("EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k")
pop3 = c("IRQ0","EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","RUS2k","EURse1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k","JACKDAW")
pop4 = c("IRQ0")
f3_1_seg <- f3(f2_blocks_seg,pop2,pop1,pop3) #none of the est <0
f3_2_seg <- f3(f2_blocks_seg,pop3,pop4,pop2) #none of the est <0
#e.g the below f2 equation give the same output as f3
1/2*(f2(f2_blocks_seg,"EURc0","EURse0")$est + f2(f2_blocks_seg,"EURc0","SPA0")$est - f2(f2_blocks_seg,"EURse0","SPA0")$est)
f3(f2_blocks_seg,"EURc0","EURse0","SPA0")
# why f3 is not negative in any case? Even though f3<0 is expected when the target pop has admixture with the sources, post-admixture drift can cause positive f3

### validation with wgs samples ###
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# f4-ratio with wgs samples ------------------------------------------------------------------------------------------------------------------------------
OUT <- c("AM")
#popused = c("SPA1k","SPA0","EURw0","EURc0","EURc2k","EURnc1k", "IRQ0","EURse0","EURe2k","EURse1k","EURse2k","EURse10k",OUT) #removed EURse20k,"EURc5k"
#f2_blocks = f2_from_precomp("wgs_nooutlier_more_pop", afprod = TRUE, pops=popused) #Discarding 7145 block / Discarding 509 block(s) w/o EURse20k / Discarding 41 block(s) w/o EURc5k
#f2_blocks_trans = f2_from_precomp("wgs_nooutlier_more_trans_pop", afprod = TRUE, pops=popused) #Discarding 11128 block(s) / Discarding 1332 block(s) w/o EURse20k / Discarding 67 block(s) w/o EURc5k
popused = c("SPA1k","SPA0","EURw0","EURc0","EURc2k","EURnc1k", "IRQ0","EURse0","EURe2k","EURse1k","EURse2k","EURse10k","EURse20k","EURc5k",OUT)
f2_blocks = f2_from_precomp("wgs_nooutlier_more_pop", afprod = TRUE, pops=popused) #discarding 37 blocks (8 if removed EURse20k and EURc5k)
f2_blocks_trans = f2_from_precomp("wgs_nooutlier_more_trans_pop", afprod = TRUE, pops=popused) #discarding 80 blocks (14 if removed EURse20k and EURc5k)

#first i generate f4-ratio for carrion
PX = c("SPA1k","EURw0","EURc0","EURc2k","EURnc1k","EURc5k")
num <- f4(f2_blocks, OUT,"IRQ0",PX, "EURse0") 
den <- f4(f2_blocks, OUT,"IRQ0","SPA0", "EURse0")
alpha <- num$est / den$est
se_alpha <- abs(alpha) * sqrt( (num$se / num$est)^2 + (den$se / den$est)^2 )
z_alpha <- alpha / se_alpha
f4ratio_carrion <- num %>% mutate(alpha=alpha, se_alpha=se_alpha, z_alpha=z_alpha)

num2 <- f4(f2_blocks_trans, OUT,"IRQ0",PX, "EURse0") 
den <- f4(f2_blocks_trans, OUT,"IRQ0","SPA0", "EURse0")
alpha <- num$est / den$est
se_alpha <- abs(alpha) * sqrt( (num$se / num$est)^2 + (den$se / den$est)^2 )
z_alpha <- alpha / se_alpha
f4ratio_carrion2 <- num %>% mutate(alpha=alpha, se_alpha=se_alpha, z_alpha=z_alpha)

df_carrion <- f4ratio_carrion  %>%
  mutate(SPA0 = ifelse(alpha < 0, 0, alpha), EURse0 = ifelse(1-alpha > 1, 1, 1-alpha)) %>%
  pivot_longer(cols = c(SPA0, EURse0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(pop3, levels = c("SPA1k","EURw0","EURw2k","EURnw2k","EURnc1k","EURc0","EURc2k","EURc5k"))) %>% group_by(target) %>%
  mutate(Mixture = factor(Mixture, levels = c("SPA0","EURse0"))) %>% 
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)


error <- df_carrion %>% filter(Mixture == "SPA0") %>% mutate(error_min = 1-ymax - se_alpha,error_max = 1-ymax + se_alpha)
plot_f4ratiocarrion <- ggplot(df_carrion, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  #geom_text(data = df_carrion,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.05),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0, 1))+
  theme_minimal() + labs(x = "Target group (all SNPs)", y = "f4-ratio") + theme(legend.position = "top")

df_carrion2 <- f4ratio_carrion2  %>%
  mutate(SPA0 = ifelse(alpha < 0, 0, alpha), EURse0 = ifelse(1-alpha > 1, 1, 1-alpha)) %>%
  pivot_longer(cols = c(SPA0, EURse0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(pop3, levels = c("SPA1k","EURw0","EURw2k","EURnw2k","EURnc1k","EURc0","EURc2k","EURc5k"))) %>% group_by(target) %>%
  mutate(Mixture = factor(Mixture, levels = c("SPA0","EURse0"))) %>% 
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

error <- df_carrion2 %>% filter(Mixture == "SPA0") %>% mutate(error_min = 1-ymax - se_alpha,error_max = 1-ymax + se_alpha)
plot_f4ratiocarrion2 <- ggplot(df_carrion2, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  #geom_text(data = df_carrion,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.05),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0, 1))+
  theme_minimal() + labs(x = "Target group (transversion SNPs)", y = "f4-ratio") + theme(legend.position = "top")

f4ratio_carrion_plot <- plot_f4ratiocarrion / plot_f4ratiocarrion2  + 
  plot_layout(heights = c(1, 1,1)) + plot_annotation(tag_levels = c("A", "B")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))

#then i generate f4-ratio for hooded
PX = c("EURse0","EURe2k","EURse1k","EURse2k","EURse10k","EURse20k")
num <- f4(f2_blocks, OUT,"SPA0",PX, "EURc0")
den <- f4(f2_blocks, OUT,"SPA0","IRQ0", "EURc0")
alpha <- num$est / den$est
se_alpha <- abs(alpha) * sqrt( (num$se / num$est)^2 + (den$se / den$est)^2 )
z_alpha <- alpha / se_alpha
f4ratio_hooded <- num %>% mutate(alpha=alpha, se_alpha=se_alpha, z_alpha=z_alpha)

num <- f4(f2_blocks_trans, OUT,"SPA0",PX, "EURc0")
den <- f4(f2_blocks_trans, OUT,"SPA0","IRQ0", "EURc0")
alpha <- num$est / den$est
se_alpha <- abs(alpha) * sqrt( (num$se / num$est)^2 + (den$se / den$est)^2 )
z_alpha <- alpha / se_alpha
f4ratio_hooded2 <- num %>% mutate(alpha=alpha, se_alpha=se_alpha, z_alpha=z_alpha)

df_hooded <- f4ratio_hooded  %>%
  mutate(IRQ0 = ifelse(alpha < 0, 0, alpha), EURc0 = ifelse(1-alpha > 1, 1, 1-alpha)) %>%
  pivot_longer(cols = c(IRQ0, EURc0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(pop3, levels = c("MON", "IRQ0","EURs0","EURsw0","EURn0","RUS0","RUS2k","NL1k","EURe0","EURe100","EURe2k","EURse0","EURse100","EURse1k","EURse2k","EURse10k","EURse20k"))) %>% group_by(target) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

error <- df_hooded %>%filter(Mixture == "EURc0") %>% mutate(error_min = pmax(0, 1 - ymax - se_alpha),error_max = pmin(1, 1 - ymax + se_alpha))
plot_f4ratiohooded <- ggplot(df_hooded, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0, 1))+
  theme_minimal() + labs(x = "Target group (all SNPs)", y = "f4-ratio") + theme(legend.position = "top",axis.text.x = element_text(size = 7,angle = 45, vjust =1, hjust = 1))

df_hooded2 <- f4ratio_hooded2  %>%
  mutate(IRQ0 = ifelse(alpha < 0, 0, alpha), EURc0 = ifelse(1-alpha > 1, 1, 1-alpha)) %>%
  pivot_longer(cols = c(IRQ0, EURc0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(pop3, levels = c("MON", "IRQ0","EURs0","EURsw0","EURn0","RUS0","RUS2k","NL1k","EURe0","EURe100","EURe2k","EURse0","EURse100","EURse1k","EURse2k","EURse10k","EURse20k"))) %>% group_by(target) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

error <- df_hooded2 %>% filter(Mixture == "EURc0") %>% mutate(error_min = 1-ymax - se_alpha,error_max = 1-ymax + se_alpha)
plot_f4ratiohooded2 <- ggplot(df_hooded2, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0, 1))+
  theme_minimal() + labs(x = "Target group (transversion SNPs)", y = "f4-ratio") + theme(legend.position = "top",axis.text.x = element_text(size = 7,angle = 45, vjust =1, hjust = 1))

f4ratio_hooded_plot <- plot_f4ratiohooded / plot_f4ratiohooded2  + 
  plot_layout(heights = c(1, 1,1)) + plot_annotation(tag_levels = c("A", "B")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))

# -----------------------------------------------------------------------------------------------------------------------------------------------------
# qpadm with wgs samples ------------------------------------------------------------------------------------------------------------------------------
popused = c("SPA1k","SPA0","EURw0","EURc0","EURc2k","EURnc1k", "IRQ0","EURse0","EURe2k","EURse1k","EURse2k","EURse10k","MON","JACKDAW","AM","EURse20k","EURc5k") #removed EURse20k,"EURc5k"
f2_blocks = f2_from_precomp("wgs_nooutlier_more_pop_nomis", afprod = TRUE, pops=popused) #Discarding 62 block(s) | 39 blocks (more SNPs)
f2_blocks_trans = f2_from_precomp("wgs_nooutlier_more_trans_pop_nomis", afprod = TRUE, pops=popused) #Discarding 105 block(s) | 83 blocks (more SNPs)

popb = c("EURse0","EURe2k","EURse1k","EURse2k","EURse10k","EURse20k")
target=popb
adm4 <-data.frame(target=popb, IRQ0=1, EURc0=0,se=0,z1=0,z2=0,p_rankdrop=0)
adm5 <-data.frame(target=popb, IRQ0=1, EURc0=0,se=0,z1=0,z2=0,p_rankdrop=0)
adm6 <-data.frame(target=popb, IRQ0=1, EURc0=0,se=0,z1=0,z2=0,p_rankdrop=0)
left =c("IRQ0","EURc0")
right=c("MON","SPA0","JACKDAW") #w/o AM for 3 outgroup

for (i in 1:length(target)) { 
  adm4[i,2] <- qpadm(f2_blocks, left, right, target[i])$weights$weight[1]
  adm4[i,3] <- qpadm(f2_blocks, left, right, target[i])$weights$weight[2] 
  adm4[i,4] <- qpadm(f2_blocks, left, right, target[i])$weights$se[1]
  adm4[i,5] <- qpadm(f2_blocks, left, right, target[i])$weights$z[1]
  adm4[i,6] <- qpadm(f2_blocks, left, right, target[i])$weights$z[2]
  adm4[i,7] <- qpadm(f2_blocks, left, right, target[i])$rankdrop$p[1]
}
write.table(adm4, file="./f4adm_pseudohaploid/wgs_nooutlier_target_hooded_3outgrp_moreSNPs_nomis.txt", quote=FALSE, row.names = FALSE, sep="\t")
for (i in 1:length(target)) { 
  adm5[i,2] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$weight[1]
  adm5[i,3] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$weight[2]
  adm5[i,4] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$se[1]
  adm5[i,5] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$z[1]
  adm5[i,6] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$z[2]
  adm5[i,7] <- qpadm(f2_blocks_trans, left, right, target[i])$rankdrop$p[1]}
write.table(adm5, file="./f4adm_pseudohaploid/wgs_nooutlier_trans_target_hooded_3outgrp_moreSNPs_nomis.txt", quote=FALSE, row.names = FALSE, sep="\t")

#adm4<- read_table(file="./f4adm_pseudohaploid/wgs_nooutlier_target_hooded_4outgrp.txt",col_name=TRUE)
#adm5 <- read_table(file="./f4adm_pseudohaploid/wgs_nooutlier_trans_target_hooded_4outgrp.txt",col_name=TRUE)

df_long4 <- adm4 %>%
  mutate(IRQ0 = ifelse(IRQ0 < 0, 0, IRQ0), EURc0 = ifelse(EURc0 > 1, 1, EURc0)) %>%
  pivot_longer(cols = c(IRQ0, EURc0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("JACKDAW", "IRQ0","EURs0","EURsw0","EURn0","RUS0","RUS2k","NL1k","EURe0","EURe100","EURe2k","EURse0","EURse100","EURse1k","EURse2k","EURse10k","EURse20k","EURse10k1","EURse10k2","EURse20k2","EURse20k1"))) %>% group_by(target) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

df_long5 <- adm5 %>%
  mutate(IRQ0 = ifelse(IRQ0 < 0, 0, IRQ0), EURc0 = ifelse(EURc0 > 1, 1, EURc0)) %>%
  pivot_longer(cols = c(IRQ0, EURc0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("JACKDAW", "IRQ0","EURs0","EURsw0","EURn0","RUS0","RUS2k","NL1k","EURe0","EURe100","EURe2k","EURse0","EURse100","EURse1k","EURse2k","EURse10k","EURse20k","EURse10k1","EURse10k2","EURse20k2","EURse20k1"))) %>% group_by(target) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)

error <- df_long4 %>%filter(Mixture == "EURc0") %>% mutate(error_min = pmax(0, 1 - ymax - se),error_max = pmin(1, 1 - ymax + se))
x2 <- ggplot(df_long4, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long4,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.01),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0.00, 1.00))+
  theme_minimal() + labs(x = "Target group (all SNPs)", y = "qpAdm") + theme(legend.position = "top",axis.text.x = element_text(size = 7,angle = 45, vjust = 1, hjust = 1))

error <- df_long5 %>% filter(Mixture == "EURc0") %>%  mutate(error_min = pmax(0, 1 - ymax - se),error_max = pmin(1, 1 - ymax + se))
x2_trans <- ggplot(df_long5, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long5,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.01),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0.00, 1.00))+
  theme_minimal() + labs(x = "Target group (transversion SNPs)", y = "qpAdm") + theme(legend.position = "top",axis.text.x = element_text(size = 7,angle = 45, vjust = 1, hjust = 1))

adm2_plot <- x2 / x2_trans  + 
  plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = c("A", "B", "C")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", plot.margin = unit(c(0, 0, 0, 0), "pt")) #8x10.5

popa = c("SPA1k","EURw0","EURc0","EURc2k","EURnc1k","EURc5k")
target=popa
adm1 <-data.frame(target=popa, SPA0=1, EURse0=0, se=0,z1=0,z2=0,p_rankdrop=0)
adm2 <-data.frame(target=popa, SPA0=1, EURse0=0, se=0,z1=0,z2=0,p_rankdrop=0)
adm3 <-data.frame(target=popa, SPA0=1, EURse0=0, se=0,z1=0,z2=0,p_rankdrop=0)
left =c("SPA0","EURse0")
right=c("MON","IRQ0","JACKDAW") #no AM for 3 outgroup

for (i in 1:length(target)) { 
  adm1[i,2] <- qpadm(f2_blocks, left, right, target[i])$weights$weight[1]
  adm1[i,3] <- qpadm(f2_blocks, left, right, target[i])$weights$weight[2]
  adm1[i,4] <- qpadm(f2_blocks, left, right, target[i])$weights$se[1]
  adm1[i,5] <- qpadm(f2_blocks, left, right, target[i])$weights$z[1]
  adm1[i,6] <- qpadm(f2_blocks, left, right, target[i])$weights$z[2]
  adm1[i,7] <- qpadm(f2_blocks, left, right, target[i])$rankdrop$p[1]} 
write.table(adm1, file="./f4adm_pseudohaploid/wgs_nooutlier_target_carrion_3outgrp_moreSNPs_nomis.txt", quote=FALSE, row.names = FALSE, sep="\t")
for (i in 1:length(target)) { 
  adm2[i,2] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$weight[1]
  adm2[i,3] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$weight[2]
  adm2[i,4] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$se[1]
  adm2[i,5] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$z[1]
  adm2[i,6] <- qpadm(f2_blocks_trans, left, right, target[i])$weights$z[2]
  adm2[i,7] <- qpadm(f2_blocks_trans, left, right, target[i])$rankdrop$p[1]}
write.table(adm2, file="./f4adm_pseudohaploid/wgs_nooutlier_trans_target_carrion_3outgrp_moreSNPs_nomis.txt", quote=FALSE, row.names = FALSE, sep="\t")

#adm1<- read_table(file="./f4adm_pseudohaploid/wgs_nooutlier_target_carrion_3outgrp.txt",col_name=TRUE)
#adm2 <- read_table(file="./f4adm_pseudohaploid/wgs_nooutlier_trans_target_carrion_3outgrp.txt",col_name=TRUE)

df_long <- adm1 %>%
  mutate(SPA0 = ifelse(SPA0 < 0, 0, SPA0), EURse0 = ifelse(EURse0 > 1, 1, EURse0)) %>%
  pivot_longer(cols = c(SPA0, EURse0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("SPA1k","EURw0","EURw2k","EURnw2k","EURnc1k","EURc0","EURc2k","EURc5k"))) %>% group_by(target) %>%
  mutate(Mixture = factor(Mixture, levels = c("SPA0","EURse0"))) %>% 
  arrange(target, Mixture) %>%  mutate(Proportion = ifelse(Proportion>1, 1, Proportion)) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion) 

df_long2 <- adm2 %>%
  mutate(SPA0 = ifelse(SPA0 < 0, 0, SPA0), EURse0 = ifelse(EURse0 > 1, 1, EURse0)) %>%
  pivot_longer(cols = c(SPA0, EURse0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("SPA1k","EURw0","EURw2k","EURnw2k","EURnc1k","EURc0","EURc2k","EURc5k"))) %>% group_by(target) %>%
  mutate(Mixture = factor(Mixture, levels = c("SPA0","EURse0"))) %>% 
  arrange(target, Mixture) %>%  mutate(Proportion = ifelse(Proportion>1, 1, Proportion)) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion) 

error <- df_long %>% filter(Mixture == "SPA0") %>%  mutate(error_min = pmax(0, 1 - ymax - se),error_max = pmin(1, 1 - ymax + se))
x1 <- ggplot(df_long, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.01),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0.00, 1.00))+
  theme_minimal() + labs(x = "Target group (all SNPs)", y = "qpAdm") + theme(legend.position = "top")

error <- df_long2 %>% filter(Mixture == "SPA0") %>% mutate(error_min = pmax(0, 1 - ymax - se),error_max = pmin(1, 1 - ymax + se))
x1_trans <- ggplot(df_long2, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long2,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.01),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0.00, 1.00))+
  theme_minimal() + labs(x = "Target group (transversion SNPs)", y = "qpAdm") + theme(legend.position = "top")

#plot_grid(x1, x1_trans, x1_outasc , labels = c("A", "B","C"), ncol = 1)
adm1_plot <- x1 / x1_trans  + 
  plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = c("A", "B", "C")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", plot.margin = unit(c(0, 0, 0, 0), "pt")) #8.5x6 protrait
