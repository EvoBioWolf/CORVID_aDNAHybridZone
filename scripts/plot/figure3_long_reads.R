library("plotly")
library("dplyr")

# Panel A-----------------------------------------------------------------------------------------------------------------------------------
setwd("PATH/05_aDNA/02_results/uli_pca")
dat1 <- read.table("add_fresh_anc_out_tmp_PCs_All_pol_smallerregioncleaned_both_peaks.txt", header = TRUE)

df1 <- dat1 %>%
  mutate(identity = case_when(
    newhybrids_P0P1Allo.admix99 %in% c("P0car", "Bxcar") ~ "carrion",
    newhybrids_P0P1Allo.admix99 %in% c("P1hood", "Bxhood") ~ "hooded",
    newhybrids_P0P1Allo.admix99 %in% c("F1", "F2") ~ "hybrid",TRUE ~ NA_character_)) %>%
  filter(HI >=0.99 | HI<=0.01 | HI>=0.45 & HI<=0.55 | is.na(HI)) %>% 
  filter(zone != "ancient") 

#the sample we chose to long read sequence
to_label <- c("D_Ne_Y39" , "D_Rb_Y23", "DKoC21", "DKoC53")
df_label <- df1 %>% filter(sample_code %in% to_label)

p1 <- ggplot(df1) +
  geom_point(aes(x=PC1, y=PC2, color=identity), size=1) +
  geom_text(data = df_label,aes(x = PC1, y = PC2, label = gsub("_", "", sampleID)),
    size = 3,vjust = 0.5, hjust = 1,  nudge_x = 0.015) +
  #geom_text(aes(x=PC1, y=PC2), label=df$sampleID, check_overlap=TRUE, size = 2, vjust=-1) +
  scale_color_manual(name="", values=c("#E69F00", "#56B4E9","#009E73")) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
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
        axis.text.y=element_text (size=10),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10))
# -----------------------------------------------------------------------------------------------------------------------------------------

# Panel C ---------------------------------------------------------------------------------------------------------------------------------
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

# Panel D & F (carrion crows)--------------------------------------------------------------------------------------------------------------
# How chr18 change through time due to recombination (genotype >=3)
library(vcfR)
library(tidyverse)
library(data.table)
library(dplyr)
library(slider)

setwd("PATH/05_aDNA/01_angsd_enriched_allsites_carrion")
SPA0 <- c("E01_E", "E02_E", "E03_E", "E04_E", "E05_E", "E06_E", "E07_E") #picked 7 out of 15 to match with hooded crow
EURse0 <- c("B01_B", "B02-b_B", "B03_B", "B04_B", "B05_B", "B06_B", "B07-b_B")
EURc0  <- c("D03_D", "D04_D", "D11_D", "D12_D") # reduced to four of the poorest quality
EURc2k <- c("ENG001_TE_D", "RTT001_TE_D", "WMP005_TE_D", "WMP006_TE_D")
EURc5k <- c("HOC005_TE_D")
EURw2k<- c("GAB001_TE_F", "LSS004_TE_F")
groups <- list(SPA0 = SPA0, EURse0 = EURse0, EURc0 = EURc0, EURc2k = EURc2k, EURc5k = EURc5k)
chromosomes <- c(1:15, 17:24, 26:28, "1A", "4A")
bin_size <- 1000
min_snps_per_bin <- 1 
eps   <- 1e-9    #to avoid exact 0/1 that gives INF/-INF in LL and deltaLL
alpha <- 0.00     # pseudocount numerator
beta  <- 0.00     # pseudocount denominator
summary_list <- list()
ll_all_bins_list <- list()
ancient_pops <- c("EURc0","EURc2k","EURc5k")
ref_pops <- c(carrion = "SPA0", hooded = "EURse0")

for(chr in chromosomes){
  cat("Processing chromosome:", chr, "\n")
  vcf_file <- paste0("enriched_allsites_carrion_geno_q20_dp3_chr", as.character(chr), ".vcf.gz")
  vcf <- read.vcfR(vcf_file)
  CHR <- vcf@fix[,"CHROM"]
  POS <- as.numeric(vcf@fix[,"POS"])
  meta <- data.frame(CHR = CHR, POS = POS, stringsAsFactors = FALSE)
  gt <- extract.gt(vcf)
  
  # keep SNPs where at least one EURc2k sample is called 
  eurc2k_called <- rowSums(!is.na(gt[, groups$EURc2k, drop = FALSE])) > 0
  gt <- gt[eurc2k_called, , drop = FALSE]
  meta <- meta[eurc2k_called, , drop = FALSE]
  
  # convert genotypes to ALT counts
  gt_to_altcount <- function(gt_vec){
    gt_vec <- str_replace_all(gt_vec, "\\.", NA_character_)
    alleles <- str_split_fixed(gt_vec, "/", 2)
    apply(alleles, 1, function(x){
      if(any(is.na(x))) return(NA_real_)
      sum(as.numeric(x))
    })
  }
  
  alt_counts <- lapply(groups, function(samples){
    sub_gt <- gt[, samples, drop = FALSE]
    ac <- apply(sub_gt, 2, gt_to_altcount)
    if(is.vector(ac)) ac <- matrix(ac, ncol=1)
    alt_sum <- rowSums(ac, na.rm = TRUE)
    called <- apply(ac, 1, function(x) 2*sum(!is.na(x)))
    data.frame(alt = alt_sum, called = called)
  })
  
  # --- calculate allele frequencies per group ---
  bin_counts <- list()
  for(pop in names(alt_counts)){
    df <- data.frame(bin = floor(meta$POS / bin_size) * bin_size,
                     alt = alt_counts[[pop]]$alt,
                     called = alt_counts[[pop]]$called)
    bin_counts[[pop]] <- df %>%
      group_by(bin) %>%
      summarise(alt = sum(alt, na.rm = TRUE),
                called = sum(called, na.rm = TRUE)) }
  
  bin_af <- data.frame(bin = bin_counts[[1]]$bin)
  bin_called <- data.frame(bin = bin_counts[[1]]$bin)
  for(pop in names(bin_counts)){
    #bin_af[[pop]] <- bin_counts[[pop]]$alt / bin_counts[[pop]]$called #previous 
    #bin_called[[pop]] <- bin_counts[[pop]]$called #previous
    alt    <- bin_counts[[pop]]$alt 
    called <- bin_counts[[pop]]$called 
    af <- (alt + alpha) / (called + alpha + beta) #smoothed AF as fixed AF (0/1) creates -INF/INF in LL computation when complete mismatch or match
    af <- pmin(pmax(af, eps), 1 - eps) #avoid 0/1 to prevent INF/-INF
    bin_af[[pop]] <- af 
    bin_called[[pop]] <- called 
  }
  
  # --- binomial log-likelihood per bin ---
  ll_plot_df <- lapply(ancient_pops, function(anc){
    # skip if no data
    if(all(is.na(bin_af[[anc]])) || all(is.na(bin_called[[anc]]))){
      return(NULL) }
    alt_anc <- bin_af[[anc]] * bin_called[[anc]]
    n_called_anc <- bin_called[[anc]]
    p_carrion <- bin_af[[ref_pops["carrion"]]]
    p_hooded <- bin_af[[ref_pops["hooded"]]]
    
    # skip bins with zero calls
    valid <- n_called_anc > 0 & !is.na(alt_anc) & !is.na(p_carrion) & !is.na(p_hooded)
    if(sum(valid) == 0) return(NULL)
    ll_carrion <- rep(NA, length(alt_anc))
    ll_hooded <- rep(NA, length(alt_anc))
    ll_carrion[valid] <- dbinom(round(alt_anc[valid]), size = n_called_anc[valid], prob = p_carrion[valid], log = TRUE)
    ll_hooded[valid] <- dbinom(round(alt_anc[valid]), size = n_called_anc[valid], prob = p_hooded[valid], log = TRUE)
    
    data.frame(bin = bin_af$bin[valid], 
               delta_ll = ll_carrion[valid] - ll_hooded[valid], 
               ancient = anc,
               chromosome = as.character(chr))  }) %>% bind_rows()
  
  # --- categorize by assignment ---
  if(nrow(ll_plot_df) > 0){
    summary_bins <- ll_plot_df %>%
      mutate(assignment = case_when(delta_ll > 1  ~ "carrion",
                                    delta_ll < -1 ~ "hooded",
                                    TRUE ~ "neutral")) %>%
      filter(assignment != "neutral") %>%
      group_by(ancient) %>%
      summarise(n_bins_total = n(),
                n_carrion_like = sum(delta_ll > 1),
                prop_carrion = n_carrion_like / n_bins_total,
                n_hooded_like = sum(delta_ll < -1),
                prop_hooded = n_hooded_like / n_bins_total,
                chromosome = as.character(chr),
                .groups = "drop")
  } else {
    summary_bins <- data.frame(
      ancient = character(0),
      n_bins_total = integer(0),
      n_carrion_like = integer(0),
      prop_carrion = numeric(0),
      n_hooded_like = integer(0),
      prop_hooded = numeric(0),
      chromosome = character(0)
    )
  }
  
  summary_list[[as.character(chr)]] <- summary_bins
  ll_all_bins_list[[as.character(chr)]] <- ll_plot_df
}

all_summary_carrion <- bind_rows(summary_list)
all_bins_carrion <- bind_rows(ll_all_bins_list) 
all_summary_carrion$chromosome[is.na(all_summary_carrion$chromosome) | all_summary_carrion$chromosome==""] <- "empty_chr"  
all_summary_carrion <- all_summary_carrion %>%
  mutate(chromosome = factor(chromosome, levels = unique(chromosome)),
         ancient = factor(ancient, levels = c("EURc0","EURc2k","EURc5k","EURw2k"))) 
all_summary_carrion$chromosome <- factor(all_summary_carrion$chromosome, levels = c(1,"1A",2:4,"4A",5:28))

#Significance test
sig_df <- all_summary_carrion %>%
  filter(ancient %in% c("EURc0", "EURc2k")) %>%
  dplyr::select(chromosome, ancient, n_carrion_like, n_bins_total) %>%
  tidyr::pivot_wider(names_from = ancient,values_from = c(n_carrion_like, n_bins_total)) %>%
  rowwise() %>%
  mutate(pval = ifelse(is.na(n_carrion_like_EURc0) | is.na(n_carrion_like_EURc2k) |
                         n_bins_total_EURc0 == 0 | n_bins_total_EURc2k == 0,NA_real_,
                       prop.test(x = c(n_carrion_like_EURc0, n_carrion_like_EURc2k), n = c(n_bins_total_EURc0, n_bins_total_EURc2k))$p.value),
         star = case_when(is.na(pval) ~ "", pval < 0.001 ~ "***", pval < 0.01 ~ "**", pval < 0.05 ~ "*", TRUE~ "")) %>%
  ungroup()

star_df <- all_summary_carrion %>%
  group_by(chromosome) %>%
  summarise(y = max(prop_carrion, na.rm = TRUE) + 0.05) %>%
  left_join(sig_df %>% dplyr::select(chromosome, star), by = "chromosome")

# Proportion SPA per chromosome
p_grouped_bar <- ggplot(all_summary_carrion%>% filter(ancient != "EURc5k"), aes(x = chromosome, y = prop_carrion, colour = ancient, group = ancient, fill=ancient), linetype = ancient) + #  
  geom_col(position = position_dodge(width = 0.8)) +
  geom_line(size = 0.8) +
  #geom_text(aes(label = n_bins_total),vjust = -0.5,position = position_dodge(width = 1),angle = 45,hjust = 0, size=2.5, color = "black") +
  geom_text(data = star_df, aes(x = chromosome, y = y, label = star), angle = 90, color = "black", size = 5, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(0.25,0.8), labels = scales::percent_format()) +
  scale_colour_manual(name="", values = c("EURc0" = "#FFD580", "EURc2k" = "#FFA500", "EURc5k" = "#FF8C00" )) + 
  scale_fill_manual(name="", values = c("EURc0" = "#FFD580", "EURc2k" = "#FFA500", "EURc5k" = "#FF8C00" )) + 
  scale_linetype_manual(name="",values = c("EURc0" = "solid", "EURc2k" = "dashed", "EURc5k"="dotted")) +
  labs(x = "Chromosome", y = "Proportion of SPA bins") +
  theme_minimal(base_size = 10) +
  theme(legend.position="top",
        axis.text.x = element_text(angle = 45, hjust = 1))

p_grouped_bar_intro <- ggplot(all_summary_carrion%>% filter(ancient != "EURc5k"), aes(x = chromosome, y = prop_hooded, colour = ancient, group = ancient, fill=ancient), linetype = ancient) + #  
  geom_col(position = position_dodge(width = 0.8)) +
  geom_line(size = 1) +
  geom_text(aes(label = n_bins_total),vjust = -0.5,position = position_dodge(width = 1),angle = 45,hjust = 0, size=2.5, color = "black") +
  geom_text(data = star_df, aes(x = chromosome, y = y, label = star), angle = 90, color = "black", size = 6, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(0.25,0.8), labels = scales::percent_format()) +
  scale_colour_manual(name="", values = c("EURc0" = "#FFD580", "EURc2k" = "#FFA500", "EURc5k" = "#FF8C00" )) + 
  scale_fill_manual(name="", values = c("EURc0" = "#FFD580", "EURc2k" = "#FFA500", "EURc5k" = "#FF8C00" )) + 
  scale_linetype_manual(name="",values = c("EURc0" = "solid", "EURc2k" = "dashed", "EURc5k"="dotted")) +
  labs(x = "Chromosome", y = "Proportion of introgressed bins") +
  theme_minimal(base_size = 10) +
  theme(legend.position="top",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Scatter plot of chr18
bins18 <- all_bins_carrion %>% filter(chromosome == 18) %>%
  filter(delta_ll > 1 | delta_ll < -1) %>%
  mutate(assignment = case_when(delta_ll > 1  ~ "carrion",delta_ll < -1 ~ "hooded")) 
bins18$ancient <- factor(bins18$ancient, levels = c("EURc0", "EURc2k", "EURc5k"))

# Scatter plot: bin 6665000 is an outlier for EURc2k and 5k (not shown in plot)
p <- ggplot(bins18%>% filter(ancient != "EURc5k"), aes(x = bin, y = delta_ll, color = assignment)) +
  geom_point(alpha = 1, size = 1) +
  geom_vline(xintercept = 10000000, linetype = "dotted", color = "red", linewidth=0.6) +
  geom_vline(xintercept = 8500000, linetype = "dotted", color = "red", linewidth=0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("carrion" = "black", "hooded" = "grey"), name = "") +
  labs( x = "Chromosome 18 (bp)",  y = expression(Delta*" log-likelihood")) +
  facet_wrap(~ ancient, ncol = 1, scales = "free_y") +  # separate panels per population
  guides(color = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(limits = c(-500,500)) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 10, face = "bold"),
        panel.spacing = unit(0.2, "lines"),
        plot.title = element_text(hjust = 0),
        legend.position = "top",
        panel.grid.minor = element_blank())

# -----------------------------------------------------------------------------------------------------------------------------------------

# Panel E & G (hooded crows)---------------------------------------------------------------------------------------------------------------
setwd("PATH/05_aDNA/01_angsd_enriched_allsites_hooded")
SPA0red <- c("E01_E", "E02_E", "E03_E", "E04_E", "E05_E")
IRQ0 <- c("IRQ641313_IRQ", "IRQ641314_IRQ", "IRQ641315_IRQ", "IRQ641341_IRQ", "IRQ645931_IRQ")
EURc0  <- c("D03_D", "D04_D", "D11_D", "D12_D") # reduced to four of the poorest quality
EURse0red <- c("B01_B", "B02-b_B")
EURse2k<- c("NCP001_TE_B", "NCP002_TE_B")
EURse10k<- c("DVT014_TE_B", "DVT016_TE_B")
EURse20k<- c("DVT022_TE_B", "KRZ002_TE_B")
groups <- list(EURc0 = EURc0,IRQ0=IRQ0, EURse2k=EURse2k, EURse10k=EURse10k, EURse20k=EURse20k, EURse0red=EURse0red)
chromosomes <- c(1:15, 17:24, 26:28, "1A", "4A")
bin_size <- 1000
min_snps_per_bin <- 1 
eps   <- 1e-9    #to avoid exact 0/1 that gives INF/-INF in LL and deltaLL
alpha <- 0.00     # pseudocount numerator
beta  <- 0.00     # pseudocount denominator
summary_list <- list()
ll_all_bins_list <- list()
ancient_pops <- c("EURse0red","EURse2k","EURse10k","EURse20k")
ref_pops <- c(carrion = "EURc0", hooded = "IRQ0")

for(chr in chromosomes){
  cat("Processing chromosome:", chr, "\n")
  #vcf_file <- paste0("enriched_allsites_carrion_geno_q20_dp3_chr", as.character(chr), ".vcf.gz")
  vcf_file <- paste0("enriched_allsites_hooded_geno_q20_dp3_chr", as.character(chr), ".vcf.gz")
  vcf <- read.vcfR(vcf_file)
  CHR <- vcf@fix[,"CHROM"]
  POS <- as.numeric(vcf@fix[,"POS"])
  meta <- data.frame(CHR = CHR, POS = POS, stringsAsFactors = FALSE)
  gt <- extract.gt(vcf)

  # keep SNPs where at least one EURse2k sample is called 
  eurse2k_called <- rowSums(!is.na(gt[, groups$EURse2k, drop = FALSE])) > 0
  gt <- gt[eurse2k_called, , drop = FALSE]
  meta <- meta[eurse2k_called, , drop = FALSE]
  
  # convert genotypes to ALT counts
  gt_to_altcount <- function(gt_vec){
    gt_vec <- str_replace_all(gt_vec, "\\.", NA_character_)
    alleles <- str_split_fixed(gt_vec, "/", 2)
    apply(alleles, 1, function(x){
      if(any(is.na(x))) return(NA_real_)
      sum(as.numeric(x))
    })
  }
  
  alt_counts <- lapply(groups, function(samples){
    sub_gt <- gt[, samples, drop = FALSE]
    ac <- apply(sub_gt, 2, gt_to_altcount)
    if(is.vector(ac)) ac <- matrix(ac, ncol=1)
    alt_sum <- rowSums(ac, na.rm = TRUE)
    called <- apply(ac, 1, function(x) 2*sum(!is.na(x)))
    data.frame(alt = alt_sum, called = called)
  })
  
  # --- calculate allele frequencies per group ---
  bin_counts <- list()
  for(pop in names(alt_counts)){
    df <- data.frame(bin = floor(meta$POS / bin_size) * bin_size,
                     alt = alt_counts[[pop]]$alt,
                     called = alt_counts[[pop]]$called)
    bin_counts[[pop]] <- df %>%
      group_by(bin) %>%
      summarise(alt = sum(alt, na.rm = TRUE),
                called = sum(called, na.rm = TRUE)) }
  
  bin_af <- data.frame(bin = bin_counts[[1]]$bin)
  bin_called <- data.frame(bin = bin_counts[[1]]$bin)
  for(pop in names(bin_counts)){
    alt    <- bin_counts[[pop]]$alt 
    called <- bin_counts[[pop]]$called 
    af <- (alt + alpha) / (called + alpha + beta) #smoothed AF as fixed AF (0/1) creates -INF/INF in LL computation when complete mismatch or match
    af <- pmin(pmax(af, eps), 1 - eps) #avoid 0/1 to prevent INF/-INF
    bin_af[[pop]] <- af 
    bin_called[[pop]] <- called 
  }
  
  # --- binomial log-likelihood per bin ---
  ll_plot_df <- lapply(ancient_pops, function(anc){
    # skip if no data
    if(all(is.na(bin_af[[anc]])) || all(is.na(bin_called[[anc]]))){
      return(NULL) }
    alt_anc <- bin_af[[anc]] * bin_called[[anc]]
    n_called_anc <- bin_called[[anc]]
    p_carrion <- bin_af[[ref_pops["carrion"]]]
    p_hooded <- bin_af[[ref_pops["hooded"]]]
    
    # skip bins with zero calls
    valid <- n_called_anc > 0 & !is.na(alt_anc) & !is.na(p_carrion) & !is.na(p_hooded)
    if(sum(valid) == 0) return(NULL)
    ll_carrion <- rep(NA, length(alt_anc))
    ll_hooded <- rep(NA, length(alt_anc))
    ll_carrion[valid] <- dbinom(round(alt_anc[valid]), size = n_called_anc[valid], prob = p_carrion[valid], log = TRUE)
    ll_hooded[valid] <- dbinom(round(alt_anc[valid]), size = n_called_anc[valid], prob = p_hooded[valid], log = TRUE)
    
    data.frame(bin = bin_af$bin[valid], 
               delta_ll = ll_carrion[valid] - ll_hooded[valid], 
               ancient = anc,
               chromosome = as.character(chr))  }) %>% bind_rows()
  
  # --- categorize by assignment ---
  if(nrow(ll_plot_df) > 0){
    summary_bins <- ll_plot_df %>%
      mutate(assignment = case_when(delta_ll > 1  ~ "carrion",
                                    delta_ll < -1 ~ "hooded",
                                    TRUE ~ "neutral")) %>%
      filter(assignment != "neutral") %>%
      group_by(ancient) %>%
      summarise(n_bins_total = n(),
                n_carrion_like = sum(delta_ll > 1),
                prop_carrion = n_carrion_like / n_bins_total,
                n_hooded_like = sum(delta_ll < -1),
                prop_hooded = n_hooded_like / n_bins_total,
                chromosome = as.character(chr),
                .groups = "drop")
  } else {
    summary_bins <- data.frame(
      ancient = character(0),
      n_bins_total = integer(0),
      n_carrion_like = integer(0),
      prop_carrion = numeric(0),
      n_hooded_like = integer(0),
      prop_hooded = numeric(0),
      chromosome = character(0)
    )
  }
  
  summary_list[[as.character(chr)]] <- summary_bins
  ll_all_bins_list[[as.character(chr)]] <- ll_plot_df
}

all_summary_hooded <- bind_rows(summary_list)
all_bins_hooded <- bind_rows(ll_all_bins_list) 
all_summary_hooded$chromosome[is.na(all_summary_hooded$chromosome) | all_summary_hooded$chromosome==""] <- "empty_chr"  
all_summary_hooded <- all_summary_hooded %>%
  mutate(chromosome = factor(chromosome, levels = unique(chromosome)),
         ancient = factor(ancient, levels = c("EURse0red","EURse2k","EURse10k","EURse20k"))) 
all_summary_hooded$chromosome <- factor(all_summary_hooded$chromosome, levels = c(1,"1A",2:4,"4A",5:28))

sig_df <- all_summary_hooded %>%
  filter(ancient %in% c("EURse0red", "EURse2k")) %>%
  dplyr::select(chromosome, ancient, n_hooded_like, n_bins_total) %>%
  tidyr::pivot_wider(names_from = ancient,values_from = c(n_hooded_like, n_bins_total)) %>%
  rowwise() %>%
  mutate(pval = ifelse(is.na(n_hooded_like_EURse0red) | is.na(n_hooded_like_EURse2k) |
                         n_bins_total_EURse0red == 0 | n_bins_total_EURse2k == 0,NA_real_,
                       prop.test(x = c(n_hooded_like_EURse0red, n_hooded_like_EURse2k), n = c(n_bins_total_EURse0red, n_bins_total_EURse2k))$p.value),
         star = case_when(is.na(pval) ~ "", pval < 0.001 ~ "***", pval < 0.01 ~ "**", pval < 0.05 ~ "*", TRUE~ "")) %>%
  ungroup()

star_df <- all_summary_hooded %>%
  group_by(chromosome) %>%
  summarise(y = max(prop_hooded, na.rm = TRUE) + 0.05) %>%
  left_join(sig_df %>% dplyr::select(chromosome, star), by = "chromosome")

# Proportion IRQ per chromosome
p_grouped_bar2 <- ggplot(all_summary_hooded%>% filter(!ancient %in% c("EURse10k","EURse20k")), aes(x = chromosome, y = prop_hooded, colour = ancient, group = ancient, fill=ancient), linetype = ancient) + #  
  geom_col(position = position_dodge(width = 0.8)) +
  geom_line(size = 0.8) +
  #geom_text(aes(label = n_bins_total),vjust = -0.5,position = position_dodge(width = 1),angle = 45,hjust = 0, size=2.5, color = "black") +
  geom_text(data = star_df, aes(x = chromosome, y = y, label = star), angle = 90, color = "black", size = 5, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(0.25,0.8), labels = scales::percent_format()) +
  scale_colour_manual(name="", values = c("EURse0red" = "#89CFF0", "EURse2k" = "#1F78B4", "EURse10k" = "#08519C", "EURse20k" = "#08306B")) + 
  scale_fill_manual(name="", values = c("EURse0red" = "#89CFF0", "EURse2k" = "#1F78B4", "EURse10k" = "#08519C", "EURse20k" = "#08306B")) + 
  scale_linetype_manual(name="",values = c("EURse0red" = "solid", "EURse2k" = "dashed", "EURse10k"="dotted", "EURse20k"="dotdash")) +
  labs(x = "Chromosome", y = "Proportion of IRQ bins") +
  theme_minimal(base_size = 10) +
  theme(legend.position="top",
        axis.text.x = element_text(angle = 45, hjust = 1))

p_grouped_bar_intro2 <- ggplot(all_summary_hooded%>% filter(!ancient %in% c("EURse10k","EURse20k")), aes(x = chromosome, y = prop_carrion, colour = ancient, group = ancient, fill=ancient), linetype = ancient) + #  
  geom_col(position = position_dodge(width = 0.8)) +
  geom_line(size = 1) +
  #geom_text(aes(label = n_bins_total),vjust = -0.5,position = position_dodge(width = 1),angle = 45,hjust = 0, size=2.5, color = "black") +
  geom_text(data = star_df, aes(x = chromosome, y = y, label = star), angle = 90, color = "black", size = 6, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(0.25,0.8), labels = scales::percent_format()) +
  scale_colour_manual(name="", values = c("EURse0red" = "#89CFF0", "EURse2k" = "#1F78B4", "EURse10k" = "#08519C", "EURse20k" = "#08306B")) + 
  scale_fill_manual(name="", values = c("EURse0red" = "#89CFF0", "EURse2k" = "#1F78B4", "EURse10k" = "#08519C", "EURse20k" = "#08306B")) + 
  scale_linetype_manual(name="",values = c("EURse0red" = "solid", "EURse2k" = "dashed", "EURse10k"="dotted", "EURse20k"="dotdash")) +
  labs(x = "Chromosome", y = "Proportion of introgressed bins") +
  theme_minimal(base_size = 10) +
  theme(legend.position="top",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Scatter plot of chr18
bins18 <- all_bins_hooded %>% filter(chromosome == 18) %>%
  filter(delta_ll > 1 | delta_ll < -1) %>%
  mutate(assignment = case_when(delta_ll > 1  ~ "carrion",delta_ll < -1 ~ "hooded")) 
bins18$ancient <- factor(bins18$ancient, levels = c("EURse0red","EURse2k","EURse10k","EURse20k"))

q <- ggplot(bins18%>% filter(!ancient %in% c("EURse10k","EURse20k")), aes(x = bin, y = delta_ll, color = assignment)) +
  geom_point(alpha = 1, size = 1) +
  geom_vline(xintercept = 10000000, linetype = "dotted", color = "red", linewidth=0.6) +
  geom_vline(xintercept = 8500000, linetype = "dotted", color = "red", linewidth=0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("carrion" = "black", "hooded" = "grey"), name = "") +
  labs( x = "Chromosome 18 (bp)",  y = expression(Delta*" log-likelihood")) +
  facet_wrap(~ ancient, ncol = 1, scales = "free_y") +  # separate panels per population
  guides(color = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(limits = c(-500,500)) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 10, face = "bold"),
        panel.spacing = unit(0.2, "lines"),
        plot.title = element_text(hjust = 0),
        legend.position = "top",
        panel.grid.minor = element_blank())
# -----------------------------------------------------------------------------------------------------------------------------------------

# Combine plots ---------------------------------------------------------------------------------------------------------------------------
labeled_spacer <- ggplot() + theme_void() + labs(title = "")
t1 <- p1 + labeled_spacer +p2 + plot_layout(widths = c(0.8, 1.1, 1.1))
t2 <- p + q + plot_layout(widths = c(0.8, 0.8 + 0.1))
t3 <- p_grouped_bar + plot_layout(widths = c(3))
t4 <- p_grouped_bar2 + plot_layout(widths = c(3))
t1/t2/t3/t4  + plot_layout(heights = c(0.4, 0.8, 0.4,0.4)) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D")) & theme(plot.margin = unit(c(1, 1, 1, 1), "pt")) 


