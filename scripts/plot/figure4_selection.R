library(vcfR)
library(tidyverse)
library(data.table)
library(dplyr)
library(slider)

# allele frequency based method to decide which reference explains the target population better based on log-likelihood
# Panel A & B (carrion crows)--------------------------------------------------------------------------------------------------------------
# How each chromosome changed between present and 2kya using AF assignment
setwd("/PATH/05_aDNA/01_angsd_enriched_allsites_carrion")
SPA0 <- c("E01_E", "E02_E", "E03_E", "E04_E", "E05_E", "E06_E", "E07_E") #picked 7 out of 15 to match with EURse0
EURse0 <- c("B01_B", "B02-b_B", "B03_B", "B04_B", "B05_B", "B06_B", "B07-b_B")
EURc0  <- c("D03_D", "D04_D", "D11_D", "D12_D") # reduced to four of the poorest quality
EURc2k <- c("ENG001_TE_D", "RTT001_TE_D", "WMP005_TE_D", "WMP006_TE_D")
EURc5k <- c("HOC005_TE_D")
EURw2k<- c("GAB001_TE_F", "LSS004_TE_F")
groups <- list(SPA0 = SPA0, EURse0 = EURse0, EURc0 = EURc0, EURc2k = EURc2k, EURc5k = EURc5k)
chromosomes <- c(1:15, 17:24, 26:28, "1A", "4A")
bin_size <- 1000
eps   <- 1e-9    #to avoid exact 0/1 that gives INF/-INF in LL and deltaLL
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
    af <- (alt) / (called) #smoothed AF as fixed AF (0/1) creates -INF/INF in LL computation when complete mismatch or match
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

# Prepare plot
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

# Scatter plot of chr18
bins18 <- all_bins_carrion %>% filter(chromosome == 18) %>%
  filter(delta_ll > 1 | delta_ll < -1) %>%
  mutate(assignment = case_when(delta_ll > 1  ~ "carrion",delta_ll < -1 ~ "hooded")) 
bins18$ancient <- factor(bins18$ancient, levels = c("EURc2k","EURc0", "EURc5k"))

# Scatter plot: bin 6665000 is an outlier for EURc2k and 5k (not shown in plot)
p <- ggplot(bins18%>% filter(ancient != "EURc5k"), aes(x = bin, y = delta_ll, color = assignment)) +
  geom_point(alpha = 1, size = 1) +
  #geom_vline(xintercept = 10000000, linetype = "dotted", color = "red", linewidth=0.6) +
  #geom_vline(xintercept = 8500000, linetype = "dotted", color = "red", linewidth=0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("carrion" = "black", "hooded" = "grey"), name = "") +
  labs( x = "Chromosome 18 (bp)",  y = expression(Delta*" log-likelihood")) +
  facet_wrap(~ ancient, ncol = 1, scales = "free_y") +  # separate panels per population
  guides(color = guide_legend(override.aes = list(size = 4))) +
  #scale_y_continuous(limits = c(-500,500)) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 10, face = "bold"),
        panel.spacing = unit(0.2, "lines"),
        plot.title = element_text(hjust = 0),
        legend.position = "top",
        panel.grid.minor = element_blank())

# mark chromosomes that have at least one star (pval<0.05)
chrom_highlight <- star_df %>%
  group_by(chromosome) %>%
  summarise(highlight = any(star != ""))  # TRUE if at least one "*"
all_summary_carrion2 <- all_summary_carrion %>%
  left_join(chrom_highlight, by = "chromosome") %>%
  mutate(highlight=  factor(highlight,levels=c(FALSE,TRUE))) %>% arrange(highlight)
all_summary_carrion2$ancient <- factor(all_summary_carrion2$ancient, levels = c("EURc2k","EURc0", "EURc5k"))
present <- all_summary_carrion2 %>% filter(ancient == "EURc0", chromosome != 18)
old <-  all_summary_carrion2 %>% filter(ancient == "EURc2k", chromosome != 18)

chart1 <- ggplot(all_summary_carrion2 %>% filter(ancient != "EURc5k"),aes(x = ancient, y = prop_hooded, fill = ancient)) +
  geom_line(data = all_summary_carrion2 %>% filter(highlight == FALSE, ancient != "EURc5k"),
            aes(group = chromosome),color = "lightgrey",size = 0.8) +
  geom_line( data = all_summary_carrion2 %>% filter(highlight == TRUE, ancient != "EURc5k"),
             aes(group = chromosome),color = "black", size = 0.8) +
  geom_line( data = all_summary_carrion2 %>% filter(highlight == TRUE, ancient != "EURc5k", chromosome == 18),
             aes(group = chromosome),color = "black", size = 0.8) +
  geom_point(data = present,aes(x = "EURc0", y = mean(prop_hooded)), color="#FFA500", size=2) +
  geom_point(data = old,aes(x = "EURc2k", y = mean(prop_hooded)), color="#FFA500",size=2) +
  geom_text(data = all_summary_carrion2 %>% filter(ancient == "EURc0", highlight == TRUE, chromosome != 18),
            aes(y = prop_hooded,label = chromosome), x = "EURc0", hjust = 1.5,size = 3,color = "black") +
  geom_text(data = all_summary_carrion2 %>% filter(ancient == "EURc2k", highlight == TRUE, chromosome == 18),
            aes(y = prop_hooded,label = chromosome), x = "EURc2k", hjust = -0.1,size = 3.5,color = "black") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0.2, 0.7)) +
  labs(x = "", y = "Proportion of introgressed bins") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(hjust = 1),legend.position = "none")
# -----------------------------------------------------------------------------------------------------------------------------------------

# Panel B & C (hooded crows)---------------------------------------------------------------------------------------------------------------
setwd("/PATH/05_aDNA/01_angsd_enriched_allsites_hooded")
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
summary_list <- list()
ll_all_bins_list <- list()
ancient_pops <- c("EURse0red","EURse2k","EURse10k","EURse20k")
ref_pops <- c(carrion = "EURc0", hooded = "IRQ0")

for(chr in chromosomes){
  cat("Processing chromosome:", chr, "\n")
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
    af <- (alt) / (called) #smoothed AF as fixed AF (0/1) creates -INF/INF in LL computation when complete mismatch or match
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

#Prepare plot
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

# Scatter plot of chr18
bins18 <- all_bins_hooded %>% filter(chromosome == 18) %>%
  filter(delta_ll > 1 | delta_ll < -1) %>%
  mutate(assignment = case_when(delta_ll > 1  ~ "carrion",delta_ll < -1 ~ "hooded")) 
bins18$ancient <- factor(bins18$ancient, levels = c("EURse2k","EURse0red","EURse10k","EURse20k"))

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

# mark chromosomes that have at least one star (pval<0.05)
chrom_highlight <- star_df %>%
  group_by(chromosome) %>%
  summarise(highlight = any(star != ""))  # TRUE if at least one "*"
all_summary_hooded2 <- all_summary_hooded %>%
  left_join(chrom_highlight, by = "chromosome") %>%
  mutate(highlight=  factor(highlight,levels=c(FALSE,TRUE))) %>% 
  filter(ancient %in% c("EURse0red","EURse2k"))
all_summary_hooded2$ancient <- factor(all_summary_hooded2$ancient, levels = c("EURse2k","EURse0red","EURse10k","EURse20k"))
present <- all_summary_hooded2 %>% filter(ancient == "EURse0red", chromosome != 18)
old <-  all_summary_hooded2 %>% filter(ancient == "EURse2k", chromosome != 18)

chart2 <- ggplot(all_summary_hooded2, aes(x = ancient, y = prop_carrion, fill = ancient)) +
  geom_line(data = all_summary_hooded2 %>% filter(highlight == FALSE),
            aes(group = chromosome),color = "lightgrey",size = 0.8) +
  geom_line( data = all_summary_hooded2 %>% filter(highlight == TRUE),
             aes(group = chromosome),color = "black", size = 0.8) +
  geom_line( data =all_summary_hooded2 %>% filter(highlight == TRUE, chromosome == 18),
             aes(group = chromosome),color = "black", size = 0.8) +
  geom_point(data = present,aes(x = "EURse0red", y = mean(prop_carrion)),colour="#1F78B4", size = 2) +
  geom_point(data = old,aes(x = "EURse2k", y = mean(prop_carrion)),colour="#1F78B4", size = 2) +
  geom_text(data = all_summary_hooded2 %>% filter(ancient == "EURse2k", highlight == TRUE, chromosome != 18),
            aes(y = prop_carrion,label = chromosome), x = "EURse2k", hjust = -0.1,size = 3,color = "black") +
  geom_text(data = all_summary_hooded2 %>% filter(ancient == "EURse2k", highlight == TRUE, chromosome == 18),
            aes(y = prop_carrion,label = chromosome), x = "EURse2k", hjust = -0.1,size = 3.5,color = "black") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0.2, 0.7)) +
  labs(x = "", y = "Proportion of introgressed bins") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(hjust = 1),legend.position = "none")
# -----------------------------------------------------------------------------------------------------------------------------------------

# Combine plots ---------------------------------------------------------------------------------------------------------------------------
t1 <- p + chart1 + plot_layout(widths = c(1, 0.5))
t2 <- q + chart2 + plot_layout(widths = c(1, 0.5))
t1/t2  + plot_layout(heights = c(1,1)) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D")) & theme(plot.margin = unit(c(1, 1, 1, 1), "pt")) 

# HWE based method to decide which reference explains the target population better based on log-likelihood
# Panel A & B (carrion crows)--------------------------------------------------------------------------------------------------------------
# How each chromosome changed between present and 2kya using HWE assignment
setwd("/PATH/05_aDNA/01_angsd_enriched_allsites_carrion")
SPA0 <- c("E01_E", "E02_E", "E03_E", "E04_E", "E05_E", "E06_E", "E07_E")
EURse0 <- c("B01_B", "B02-b_B", "B03_B", "B04_B", "B05_B", "B06_B", "B07-b_B")
EURc0  <- c("D03_D", "D04_D", "D11_D", "D12_D")
EURc2k <- c("ENG001_TE_D", "RTT001_TE_D", "WMP005_TE_D", "WMP006_TE_D")
EURc5k <- c("HOC005_TE_D")
EURw2k <- c("GAB001_TE_F", "LSS004_TE_F")

groups <- list(SPA0 = SPA0, EURse0 = EURse0, EURc0 = EURc0, EURc2k = EURc2k, EURc5k = EURc5k)
chromosomes <- c(1:15, 17:24, 26:28, "1A", "4A")
bin_size <- 1000
eps <- 1e-9

ancient_pops <- c("EURc0","EURc2k","EURc5k")
ref_pops <- c(carrion = "SPA0", hooded = "EURse0")

summary_list <- list()
ll_all_bins_list <- list()

# ----------------------------
# Helper function: convert genotype strings to ALT counts
# ----------------------------
gt_to_altcount <- function(gt_vec){
  gt_vec <- str_replace_all(gt_vec, "[/|]", "/")   # handle phased genotypes
  gt_vec <- str_replace_all(gt_vec, "\\.", NA_character_)
  alleles <- str_split_fixed(gt_vec, "/", 2)
  apply(alleles, 1, function(x){
    if(any(is.na(x))) return(NA_real_)
    sum(as.numeric(x))
  })
}

# ----------------------------
# HWE log-likelihood function
# ----------------------------
hwe_loglik <- function(gt_matrix, p_ref, q_ref){
  loglik_mat <- matrix(NA_real_, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix))
  for(i in 1:nrow(gt_matrix)){
    for(j in 1:ncol(gt_matrix)){
      g <- gt_matrix[i,j]
      if(is.na(g) || g == "."){
        loglik_mat[i,j] <- NA
      } else {
        alleles <- as.numeric(str_split_fixed(g, "/", 2))
        g_sum <- sum(alleles)
        if(g_sum == 0) loglik_mat[i,j] <- log(p_ref^2 + eps)
        else if(g_sum == 1) loglik_mat[i,j] <- log(2*p_ref*q_ref + eps)
        else if(g_sum == 2) loglik_mat[i,j] <- log(q_ref^2 + eps)
      }
    }
  }
  sum(loglik_mat, na.rm = TRUE)
}

# ----------------------------
# Main chromosome loop
# ----------------------------
for(chr in chromosomes){
  cat("Processing chromosome:", chr, "\n")
  
  # Read VCF
  vcf_file <- paste0("enriched_allsites_carrion_geno_q20_dp3_chr", chr, ".vcf.gz")
  vcf <- read.vcfR(vcf_file)
  CHR <- vcf@fix[,"CHROM"]
  POS <- as.numeric(vcf@fix[,"POS"])
  meta <- data.frame(CHR = CHR, POS = POS, stringsAsFactors = FALSE)
  gt <- extract.gt(vcf)
  
  # Filter SNPs: at least one EURc2k called
  eurc2k_called <- rowSums(!is.na(gt[, groups$EURc2k, drop = FALSE])) > 0
  gt <- gt[eurc2k_called, , drop = FALSE]
  meta <- meta[eurc2k_called, , drop = FALSE]
  
  # Convert genotypes to ALT counts per group
  alt_counts <- lapply(groups, function(samples){
    sub_gt <- gt[, samples, drop = FALSE]
    ac <- apply(sub_gt, 2, gt_to_altcount)
    if(is.vector(ac)) ac <- matrix(ac, ncol=1)
    alt_sum <- rowSums(ac, na.rm = TRUE)
    called <- apply(ac, 1, function(x) 2*sum(!is.na(x)))
    data.frame(alt = alt_sum, called = called)
  })
  
  # Bin-wise allele frequencies
  bin_counts <- list()
  for(pop in names(alt_counts)){
    df <- data.frame(
      bin = floor(meta$POS / bin_size) * bin_size,
      alt = alt_counts[[pop]]$alt,
      called = alt_counts[[pop]]$called
    )
    bin_counts[[pop]] <- df %>%
      group_by(bin) %>%
      summarise(
        alt = sum(alt, na.rm = TRUE),
        called = sum(called, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  bin_af <- data.frame(bin = bin_counts[[1]]$bin)
  for(pop in names(bin_counts)){
    alt <- bin_counts[[pop]]$alt
    called <- bin_counts[[pop]]$called
    af <- (alt) / (called)
    af <- pmin(pmax(af, eps), 1 - eps)
    bin_af[[pop]] <- af
  }
  
  # ----------------------------
  # HWE-based log-likelihood per bin
  # ----------------------------
  ll_plot_df_list <- list()
  for(anc in ancient_pops){
    gt_anc <- gt[, groups[[anc]], drop = FALSE]
    
    # Ensure full chromosome coverage
    all_bins <- seq(0, max(meta$POS, na.rm = TRUE), by = bin_size)
    delta_ll_vec <- rep(NA_real_, length(all_bins))
    
    for(k in seq_along(all_bins)){
      b <- all_bins[k]
      snp_idx <- which(floor(meta$POS / bin_size) * bin_size == b)
      
      if(length(snp_idx) == 0){
        delta_ll_vec[k] <- NA
        next
      }
      
      # Reference allele frequencies
      p_carrion <- mean(bin_af[[ref_pops["carrion"]]][snp_idx], na.rm = TRUE)
      p_hooded <- mean(bin_af[[ref_pops["hooded"]]][snp_idx], na.rm = TRUE)
      q_carrion <- 1 - p_carrion
      q_hooded <- 1 - p_hooded
      
      # Compute HWE log-likelihood
      ll_carrion <- hwe_loglik(gt_anc[snp_idx, , drop = FALSE], p_carrion, q_carrion)
      ll_hooded <- hwe_loglik(gt_anc[snp_idx, , drop = FALSE], p_hooded, q_hooded)
      
      delta_ll_vec[k] <- ll_carrion - ll_hooded
    }
    
    ll_plot_df_list[[anc]] <- data.frame(
      bin = all_bins,
      delta_ll = delta_ll_vec,
      ancient = anc,
      chromosome = as.character(chr),
      stringsAsFactors = FALSE
    )
  }
  
  ll_plot_df <- bind_rows(ll_plot_df_list)
  
  # ----------------------------
  # Bin assignment based on delta_ll
  # ----------------------------
  if(nrow(ll_plot_df) > 0){
    summary_bins <- ll_plot_df %>%
      mutate(assignment = case_when(
        delta_ll > 1  ~ "carrion",
        delta_ll < -1 ~ "hooded",
        TRUE ~ "neutral"
      )) %>%
      filter(assignment != "neutral") %>%
      group_by(ancient) %>%
      summarise(
        n_bins_total = n(),
        n_carrion_like = sum(assignment == "carrion"),
        prop_carrion = n_carrion_like / n_bins_total,
        n_hooded_like = sum(assignment == "hooded"),
        prop_hooded = n_hooded_like / n_bins_total,
        chromosome = as.character(chr),
        .groups = "drop"
      )
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
