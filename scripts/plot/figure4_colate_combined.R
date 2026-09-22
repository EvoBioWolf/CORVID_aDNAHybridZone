library(tidyverse)
library(stringr)
library(ggplot2)
library(patchwork)
# ============================================================
# SETTINGS
# ============================================================

generation_time <- 5.79
REGIONS <- c(peak_chr18 = "peak_chr18", all_chr = "all_chr")
excluded_samples <- c("VKP001", "BRW001","DVT017") #non-udg samples 
excluded_age_groups <- c("more_than_20k", "1000-1500") #no wgs samples of such

# ============================================================
# PATHS
# ============================================================

dat <- "PATH"
TYP <- "wgs"
colate_dir <- file.path(dat, "02_results", "colate")
heatmap_dir <- file.path(colate_dir, "heatmap", TYP, "hooded")
pop_file <- file.path(colate_dir, "pop.txt")
summary_dir <- file.path(colate_dir, "plots", "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PLOT SETTINGS
# ============================================================

age_colours <- c(
  "present" = "#0072B2",
  "200-1000" = "#FDD49E",
  "1000-1500" = "#C2410C",
  "1500-2500" = "#F16913",
  "more_than_10k" = "#8C2D04",
  "more_than_20k" = "#542005"
)

age_shapes <- c(
  "present" = 16,
  "200-1000" = 17,
  "1000-1500" = 15,
  "1500-2500" = 18,
  "more_than_10k" = 8,
  "more_than_20k" = 4
)

hooded <- c(
  "KCZ003"  = "lightblue",
  "KCZ012"  = "blue",
  "NCP001"  = "lightgreen",
  "DVT014"=  "darkgreen")
# ============================================================
# READ SAMPLE METADATA
# ============================================================

pop_data <- read_tsv(pop_file, col_names = c("population", "sample", "age_years", "age_group"),
                     show_col_types = FALSE) %>%
  mutate(population = as.character(population),
         sample = as.character(sample),
         age_years = as.numeric(age_years),
         age_group = as.character(age_group),
         age_generations = age_years / generation_time,
         hooded_age_group = factor(age_group, levels = names(age_colours)))

# ============================================================
# FIND .COAL FILES
# ============================================================

coal_file_table <- imap_dfr(REGIONS,
                            function(chr_label, region_name) {
                              files <- list.files(heatmap_dir, pattern = paste0("^pairwise_E0[1-5]_vs_.*_", chr_label, "\\.coal$"),
                                                  full.names = TRUE)
                              tibble(region = region_name, chr_label = chr_label, source_file = files)
                            }
)

cat("Files found by region:\n")
coal_file_table %>%
  count(region) %>%
  print()

# ============================================================
# PARSE SAMPLE IDS FROM FILENAME
# ============================================================

get_pair_ids <- function(file, chr_label) {
  pair_name <- basename(file) %>%
    str_remove("\\.coal$") %>%
    str_remove("^pairwise_") %>%
    str_remove(paste0("_", chr_label, "$"))
  pair_ids <- str_split(pair_name, "_vs_", simplify = TRUE)
  if (ncol(pair_ids) != 2) {
    stop(
      "Could not parse filename:\n",
      basename(file)
    )
  }
  tibble(
    carrion_ref = pair_ids[, 1],
    hooded_sample = pair_ids[, 2]
  )
}

# ============================================================
# READ ONE .COAL FILE
# ============================================================

read_pairwise_coal <- function(file, region, chr_label) {
  lines <- readLines(file)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 3) {
    warning("Skipping incomplete file: ", basename(file))
    return(tibble())
  }
  pair_ids <- get_pair_ids(file, chr_label)
  generation <- suppressWarnings(as.numeric(strsplit(trimws(lines[2]), "\\s+")[[1]]))
  
  if (any(!is.finite(generation)) || length(generation) < 2) {
    warning("Invalid generation line: ", basename(file))
    return(tibble())
  }
  bootstrap_lines <- lines[3:length(lines)]
  map_dfr(seq_along(bootstrap_lines),
          function(i) {
            values <- suppressWarnings(as.numeric(strsplit(trimws(bootstrap_lines[i]), "\\s+")[[1]]))
            expected_length <- length(generation) + 2
            
            if (length(values) != expected_length) {
              warning("Skipping malformed row in ", basename(file),"; row ", i)
              return(tibble())
            }
            
            tibble(region = region, carrion_ref = pair_ids$carrion_ref, hooded_sample = pair_ids$hooded_sample,
                   bootstrap = values[2],
                   column_index = seq_along(generation),
                   generation = generation,
                   coal_rate = values[-c(1, 2)],
                   source_file = basename(file)
            )
          }
  )
}

# ============================================================
# READ ALL FILES
# ============================================================

pairwise_raw <- pmap_dfr(
  coal_file_table,
  function(region, chr_label, source_file) {
    read_pairwise_coal(
      file = source_file,
      region = region,
      chr_label = chr_label
    )
  }
)

cat("\nRows read by region:\n")

pairwise_raw %>%
  count(region) %>%
  print()

# ============================================================
# CONSTRUCT COLATE EPOCHS
# ============================================================

pairwise_epochs <- pairwise_raw %>%
  group_by(region,
           carrion_ref,
           hooded_sample,
           bootstrap) %>%
  arrange(column_index, .by_group = TRUE) %>%
  mutate(epoch_start = generation,
         epoch_end = lead(generation),
         epoch_duration = epoch_end - epoch_start) %>%
  ungroup() %>%
  filter(!is.na(epoch_end), epoch_duration > 0)

# ============================================================
# DEFINE COMMON TIME BINS
# ============================================================

standard_boundaries <- c(0,
                         273.729,
                         433.832,
                         687.577,
                         1089.74,
                         1727.12,
                         2737.29,
                         4338.32,
                         6875.77,
                         10897.4,
                         17271.2,
                         27372.9,
                         43383.2,
                         68757.7,
                         108974,
                         172712,
                         273729,
                         433832,
                         687577,
                         1.08974e+06,
                         1.72712e+06,
                         1.72712e+07)

standard_bins <- tibble(bin_start = head(standard_boundaries, -1),
                        bin_end = tail(standard_boundaries, -1)) %>%
  mutate(bin_duration = bin_end - bin_start,
         bin_label = paste0(round(bin_start, 1), "–", round(bin_end, 1)),
         bin_start_years = bin_start * generation_time,
         bin_end_years = bin_end * generation_time,
         bin_start_kya = bin_start_years / 1000,
         bin_end_kya = bin_end_years / 1000,
         bin_year_label = case_when(bin_start == 0 & bin_end_kya < 1 ~ paste0("0 –", round(bin_end_years)),
                                    bin_start_kya < 1 & bin_end_kya >= 1 ~ paste0(round(bin_start_years), " – ", round(bin_end_kya, 2)),
                                    bin_end_kya < 1 ~ paste0(round(bin_start_years), "–", round(bin_end_years)), TRUE ~
                                      paste0(round(bin_start_kya, 2)," – ",round(bin_end_kya, 2))))

# Only plot relatively recent time bins.
bins_to_plot <- standard_bins %>%filter(bin_start <= 4338.32)

# ============================================================
# ADD AGE INFORMATION
# ============================================================

pairwise_epochs <- pairwise_epochs %>%
  left_join(pop_data %>%
              transmute(hooded_sample = sample,
                        hooded_age_years = age_years,
                        hooded_age_generations = age_generations,
                        hooded_age_group = hooded_age_group,
                        hooded_population = population),by = "hooded_sample")

# ============================================================
# REMOVE PRE-SAMPLING TIME FOR ANCIENT SAMPLES
# ============================================================

pairwise_epochs_valid <- pairwise_epochs %>%
  mutate(valid_epoch_start = pmax(
    epoch_start,hooded_age_generations),
    valid_epoch_end = epoch_end,
    valid_epoch_duration = valid_epoch_end - valid_epoch_start) %>%
  filter(!is.na(hooded_age_generations),valid_epoch_duration > 0,is.finite(coal_rate),coal_rate >= 0)

# ============================================================
# PROJECT ONTO COMMON BINS
# ============================================================

pair_epoch_bins <- pairwise_epochs_valid %>%
  cross_join(standard_bins) %>%
  mutate(overlap_start = pmax(valid_epoch_start, bin_start),
         overlap_end = pmin(valid_epoch_end, bin_end),
         overlap_duration = overlap_end - overlap_start) %>%
  filter(overlap_duration > 0) %>%
  group_by(region,
           carrion_ref,
           hooded_sample,
           hooded_age_years,
           hooded_age_generations,
           hooded_age_group,
           bootstrap,
           bin_start,
           bin_end,
           bin_label) %>%
  summarise(coal_rate_bin =csum(coal_rate * overlap_duration, na.rm = TRUE) /
              sum(overlap_duration, na.rm = TRUE),
            duration_used_generations = sum(overlap_duration), .groups = "drop")

# ============================================================
# SUMMARISE ACROSS BOOTSTRAPS
# ============================================================

pair_epoch_summary <- pair_epoch_bins %>%
  group_by(region,
           carrion_ref,
           hooded_sample,
           hooded_age_years,
           hooded_age_generations,
           hooded_age_group,
           bin_start,
           bin_end,
           bin_label) %>%
  summarise(mean_coal_rate = mean(coal_rate_bin, na.rm = TRUE),
            median_coal_rate = median(coal_rate_bin, na.rm = TRUE),
            lower_bootstrap = quantile(coal_rate_bin,0.025,na.rm = TRUE),
            upper_bootstrap = quantile(coal_rate_bin,0.975,na.rm = TRUE),
            n_bootstraps = sum(is.finite(coal_rate_bin)),
            duration_used_generations = first(duration_used_generations), .groups = "drop") %>%
  filter(is.finite(mean_coal_rate), mean_coal_rate > 0)

cat("\nSummary rows by region:\n") 
pair_epoch_summary %>% count(region) %>%
  print()

# ============================================================
# SAVE REGION-SPECIFIC SUMMARY TABLES
# ============================================================

write_tsv(pair_epoch_summary, file.path(summary_dir,"pairwise_colate_summary_peak_chr18_and_all_chr.tsv"))
write_tsv(pair_epoch_summary %>% filter(region == "peak_chr18"), file.path(summary_dir,"pairwise_colate_summary_peak_chr18.tsv"))
write_tsv(pair_epoch_summary %>% filter(region == "all_chr"),file.path(summary_dir, "pairwise_colate_summary_all_chr.tsv"))

# ============================================================
# MAIN PLOT
# MEDIAN ACROSS CARRION REFERENCES PER ANCIENT HOODED SAMPLE
# ============================================================

plot_max_kya <- 30
individual_pairwise_native <- pairwise_epochs_valid %>%
  filter(hooded_age_group != "present",
         region %in% c("all_chr", "peak_chr18"),
         !hooded_sample %in% excluded_samples,
         !hooded_age_group %in% excluded_age_groups,
         !is.na(hooded_age_group),
         is.finite(coal_rate),
         coal_rate > 0) %>%
  group_by(region, carrion_ref, hooded_sample,
           hooded_age_years, hooded_age_generations, hooded_age_group,
           valid_epoch_start, valid_epoch_end) %>%
  summarise(mean_coal_rate = mean(coal_rate, na.rm = TRUE),
            n_bootstraps = n_distinct(bootstrap), .groups = "drop") %>%
  mutate(epoch_start_kya = valid_epoch_start * generation_time / 1000,
         epoch_end_kya = valid_epoch_end * generation_time / 1000,
         x_start = pmax(epoch_start_kya, 0),
         x_end = pmin(epoch_end_kya, plot_max_kya),
         region_label = case_when(
           region == "all_chr" ~ "All chromosomes",
           region == "peak_chr18" ~ "Peak chr18")) %>%
  filter(x_end > x_start, is.finite(mean_coal_rate), mean_coal_rate > 0)

# Median and confidence interval
hooded_sample_summary_native <- individual_pairwise_native %>%
  group_by(region,region_label,
           hooded_sample,
           hooded_age_years,
           hooded_age_generations,
           hooded_age_group,
           valid_epoch_start,
           valid_epoch_end,
           epoch_start_kya,
           epoch_end_kya,
           x_start,
           x_end) %>%
  summarise(median_coal_rate = median(mean_coal_rate, na.rm = TRUE),
            lower_coal_rate = quantile(mean_coal_rate, probs = 0.025,na.rm = TRUE),
            upper_coal_rate = quantile(mean_coal_rate, probs = 0.975, na.rm = TRUE),
            n_carrion_refs = n_distinct(carrion_ref),.groups = "drop") %>%
  mutate(x_mid = (x_start + x_end) / 2, hooded_sample = factor(hooded_sample,levels = names(hooded)),
         region_label = factor(region_label,levels = c("All chromosomes", "Peak chr18"))) %>%
  arrange(region, hooded_sample, x_start)

#Step plot
hooded_sample_step_data <- hooded_sample_summary_native %>%
  group_by(region, hooded_sample) %>%
  group_modify(~ {dat <- .x %>% arrange(x_start)
  final_row <- dat %>%
    slice_tail(n = 1) %>%
    mutate(x_start = x_end)
  bind_rows(dat, final_row)}) %>%
  ungroup() %>%
  rename(x_step = x_start) %>%
  arrange(region, hooded_sample, x_step)

allchr_sample_data <- hooded_sample_summary_native %>%
  filter(region == "all_chr")
allchr_step_data <- hooded_sample_step_data %>%
  filter(region == "all_chr")
peak18_step_data <- hooded_sample_step_data %>%
  filter(region == "peak_chr18")

p_hooded_allchr_peak18 <- ggplot() +
  geom_rect(data = allchr_sample_data, aes(xmin = x_start,xmax = x_end,
                                           ymin = lower_coal_rate, ymax = upper_coal_rate, fill = hooded_sample), inherit.aes = FALSE, alpha = 0.1, colour = NA) +
  geom_step(data = allchr_step_data, aes(x = x_step, y = median_coal_rate,
                                         colour = hooded_sample, group = hooded_sample),
            direction = "hv", linewidth = 0.9, linetype = "solid") +
  geom_step(data = peak18_step_data, aes(x = x_step, y = median_coal_rate,
                                         colour = hooded_sample, group = hooded_sample), direction = "hv", linewidth = 0.9, linetype = "dashed") +
  #geom_point(data = allchr_sample_data,aes(x = x_mid,y = median_coal_rate,colour = hooded_sample),size = 1.8) +
  #geom_point(data = hooded_sample_summary_native %>%filter(region == "peak_chr18"),aes(x = x_mid,y = median_coal_rate, colour = hooded_sample),shape = 1,size = 1.9,stroke = 0.7) +
  scale_colour_manual(values = hooded,drop = FALSE,name = "Anc hooded") +
  scale_fill_manual(values = hooded,drop = FALSE,guide = "none") +
  scale_x_continuous(name = "Time before present (kya)", limits = c(0, plot_max_kya),
                     breaks = seq(0, plot_max_kya, by = 5),
                     expand = c(0, 0)) +
  scale_y_log10(name = "Pairwise coalescence rate between\n hooded and carrion crows") +
  labs() +
  theme_classic() +
  theme(legend.position = "right")

print(p_hooded_allchr_peak18)

ggsave(filename = file.path(summary_dir,"coalescence_rates_all_chr_chr18.pdf"),
       plot = p_hooded_allchr_peak18, width = 10,height = 4,units = "in")

########### PCA ##############

# ============================================================
# SETTINGS
# ============================================================

generation_time <- 5.79
PCA_REGION <- "all_chr" #inc selfing
CHR_LABEL <- "all_chr"
pca_boundaries_kya <- c(0,5,10,20,30,150,600)
pca_time_bins <- tibble(bin_start_kya = head(pca_boundaries_kya, -1),
                        bin_end_kya = tail(pca_boundaries_kya, -1)) %>%
  mutate(pca_bin = paste0(bin_start_kya,"-",bin_end_kya," kya"),
         pca_bin = factor(pca_bin,levels = pca_bin),
         bin_start_gen = bin_start_kya * 1000 / generation_time,
         bin_end_gen = bin_end_kya * 1000 / generation_time)
pca_bin_names <- pca_time_bins$pca_bin

excluded_samples <- c("VKP001","BRW001", "SLN001", "DVT017") #non-udg 

dat <- "PATH"
TYP <- "wgs"
colate_dir <- file.path(dat, "02_results", "colate")
heatmap_dir <- file.path(colate_dir, "heatmap", TYP)
pop_file <- file.path(colate_dir, "pop.txt")
summary_dir <- file.path(colate_dir,"plots", "summary")

# ============================================================
# AGE-GROUP PLOT SETTINGS
# ============================================================

colors <- c(
  "E"  = "#7F3B08",
  "D"  = "#FDB863",
  "BE" = "#FEE0B6",
  "PL"=  "#2899B2",
  "B"  = "#117733",
  "IRQ"= "#54278F")

# ============================================================
# READ SAMPLE METADATA
# ============================================================

pop_data <- read_tsv(pop_file, col_names = c("population", "sample", "age_years", "age_group"),
                     show_col_types = FALSE) %>%
  mutate(sample = as.character(sample),
         population = as.character(population),
         age_years = as.numeric(age_years),
         age_group = factor(as.character(age_group),levels = names(age_colours)))

# ============================================================
# FIND ALL PAIRWISE .COAL FILES
# ============================================================

coal_files <- list.files(heatmap_dir, pattern = paste0("^pairwise_.*_vs_.*_", CHR_LABEL, "\\.coal$"), full.names = TRUE)
cat("Number of .coal files found:", length(coal_files), "\n")
print(head(basename(coal_files), 10))

# ============================================================
# PARSE SAMPLE IDs FROM FILENAME
# ============================================================

get_pair_ids <- function(file, chr_label) {
  pair_name <- basename(file) %>%
    str_remove("\\.coal$") %>%
    str_remove("^pairwise_") %>%
    str_remove(paste0("_", chr_label, "$"))
  
  ids <- str_split(pair_name,"_vs_",simplify = TRUE)
  
  if (ncol(ids) != 2) {
    stop("Could not parse: ", basename(file))
  }
  
  tibble(sample1 = ids[, 1], sample2 = ids[, 2])
}

# ============================================================
# READ ONE .COAL FILE
# ============================================================

read_pairwise_coal <- function(file, chr_label) {
  lines <- readLines(file)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 3) {
    warning("Skipping incomplete file: ", basename(file))
    return(tibble())
  }
  
  ids <- get_pair_ids(file, chr_label)
  generation <- suppressWarnings(as.numeric(strsplit(trimws(lines[2]), "\\s+")[[1]]))
  
  if (length(generation) < 2 || any(!is.finite(generation))) {
    warning("Bad generation line: ", basename(file))
    return(tibble())
  }
  
  bootstrap_lines <- lines[3:length(lines)]
  map_dfr(seq_along(bootstrap_lines),
          function(i) {
            v <- suppressWarnings(as.numeric(strsplit(trimws(bootstrap_lines[i]), "\\s+")[[1]]))
            if (length(v) != length(generation) + 2) {
              warning(
                "Skipping malformed row ",
                i,
                " in ",
                basename(file)
              )
              return(tibble())
            }
            
            tibble(sample1 = ids$sample1, sample2 = ids$sample2,
                   bootstrap = v[2], generation = generation,
                   column_index = seq_along(generation),
                   coal_rate = v[-c(1, 2)],
                   source_file = basename(file)
            )
          }
  )
}

# ============================================================
# READ ALL FILES
# ============================================================

pairwise_raw <- map_dfr(coal_files, ~ read_pairwise_coal( file = .x, chr_label = CHR_LABEL))
cat("Rows read:", nrow(pairwise_raw), "\n")

# ============================================================
# MAKE NATIVE COLATE EPOCHS
# ============================================================

pairwise_epochs <- pairwise_raw %>%
  group_by(sample1, sample2, bootstrap) %>%
  arrange(column_index, .by_group = TRUE) %>%
  mutate(epoch_start_gen = generation,
         epoch_end_gen = lead(generation),
         epoch_duration_gen = epoch_end_gen - epoch_start_gen) %>%
  ungroup() %>%
  filter(!is.na(epoch_end_gen), epoch_duration_gen > 0, is.finite(coal_rate), coal_rate > 0)

# ============================================================
# ADD GENERATION BOUNDARIES TO PCA BINS
# ============================================================

pca_time_bins <- pca_time_bins %>%
  mutate(bin_start_gen = bin_start_kya * 1000 / generation_time,
         bin_end_gen = bin_end_kya * 1000 / generation_time,
         pca_bin = as.character(pca_time_bins$pca_bin))

# ============================================================
# PROJECT RAW EPOCHS INTO 0-5 / 5-10 / 10-15 / 15-30 KYA
# Each RAW epoch contributes to the amount that overlaps the broad time interval.
# ============================================================

pairwise_broad_bins_bootstrap <- pairwise_epochs %>%
  cross_join(pca_time_bins) %>%
  mutate(overlap_start_gen = pmax(epoch_start_gen, bin_start_gen),
         overlap_end_gen = pmin(epoch_end_gen, bin_end_gen),
         overlap_duration_gen =overlap_end_gen - overlap_start_gen) %>%
  filter(overlap_duration_gen > 0) %>%
  group_by(sample1, sample2, bootstrap,
           bin_start_kya, bin_end_kya, pca_bin) %>%
  summarise(coal_rate_broad = sum(coal_rate * overlap_duration_gen, na.rm = TRUE) /
              sum(overlap_duration_gen, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(coal_rate_broad), coal_rate_broad > 0)

# ============================================================
# MEAN ACROSS COLATE BOOTSTRAPS
# One mean rate per sample pair per broad time bin.
# ============================================================

pairwise_broad_summary <- pairwise_broad_bins_bootstrap %>%
  group_by(sample1, sample2, bin_start_kya, bin_end_kya, pca_bin) %>%
  summarise(mean_coal_rate = mean(coal_rate_broad, na.rm = TRUE),
            n_bootstraps = n_distinct(bootstrap), .groups = "drop") %>%
  filter(is.finite(mean_coal_rate), mean_coal_rate > 0) %>%
  mutate(pair_a = pmin(sample1, sample2), pair_b = pmax(sample1, sample2)) %>% # Convert A/B and B/A into the same canonical pair.
  group_by(pair_a, pair_b, bin_start_kya, bin_end_kya, pca_bin) %>%
  summarise(mean_coal_rate = mean(mean_coal_rate, na.rm = TRUE), .groups = "drop") %>%
  rename(sample1 = pair_a, sample2 = pair_b)

write_tsv(pairwise_broad_summary, file.path(summary_dir, "exploratory_PCA_pairwise_all_chr_0_fig4e.tsv"))

# ============================================================
# FUNCTION: BUILD PAIRWISE MATRIX 
# Self-comparison files (e.g. E01_vs_E01) are retained on the diagonal
# Missing pairwise comparisons remain NA - checked NO MISSING 
# ============================================================

build_rate_matrix <- function(pair_data, one_bin) {
  bin_data <- pair_data %>%
    filter(pca_bin == one_bin,
           !sample1 %in% excluded_samples,
           !sample2 %in% excluded_samples)
  
  samples <- sort(unique(c(bin_data$sample1, bin_data$sample2)))
  
  if (length(samples) < 2) {
    stop("Too few samples in bin: ", one_bin)
  }
  
  rate_matrix <- matrix(NA_real_, nrow = length(samples),
                        ncol = length(samples), dimnames = list(samples, samples))
  
  for (i in seq_len(nrow(bin_data))) {
    s1 <- bin_data$sample1[i]
    s2 <- bin_data$sample2[i]
    value <- bin_data$mean_coal_rate[i]
    rate_matrix[s1, s2] <- value
    rate_matrix[s2, s1] <- value
  }
  
  missing_index <- which(is.na(rate_matrix), arr.ind = TRUE)
  missing_pairs <- tibble(sample1 = rownames(rate_matrix)[missing_index[, "row"]],
                          sample2 = colnames(rate_matrix)[missing_index[, "col"]]) %>%
    mutate(pca_bin = one_bin, is_self_comparison = sample1 == sample2,
           pair_a = pmin(sample1, sample2), pair_b = pmax(sample1, sample2)) %>%
    distinct(pca_bin, pair_a, pair_b, .keep_all = TRUE) %>%
    arrange(is_self_comparison, pair_a, pair_b)
  
  matrix_qc <- tibble(pca_bin = one_bin,
                      n_samples = nrow(rate_matrix),
                      total_matrix_cells = length(rate_matrix),
                      missing_cells = sum(is.na(rate_matrix)),
                      observed_cells = sum(!is.na(rate_matrix)),
                      missing_fraction = mean(is.na(rate_matrix)),
                      n_missing_unique_pairs = nrow(missing_pairs),
                      n_missing_self_comparisons = sum(missing_pairs$is_self_comparison),
                      n_missing_between_sample_pairs = sum(!missing_pairs$is_self_comparison))
  
  cat("Time bin: ", one_bin, "\n", sep = "")
  cat("Samples: ", nrow(rate_matrix), "\n", sep = "")
  cat("Missing matrix cells: ", sum(is.na(rate_matrix)), "\n", sep = "")
  cat("Missing unique pairs: ", nrow(missing_pairs), "\n",sep = "")
  
  if (any(is.na(rate_matrix))) {
    warning("PCA not run for ", one_bin, " because it contains missing values.")
    return(list(rate_matrix = rate_matrix,
                missing_pairs = missing_pairs,
                matrix_qc = matrix_qc,
                pca_result = NULL,
                scores = tibble(),
                variance_explained = NA_real_))
  }
  
  list(
    rate_matrix = rate_matrix,
    missing_pairs = missing_pairs,
    matrix_qc = matrix_qc
  )
}

# ============================================================
# BUILD MATRICES FOR ALL TIME BINS
# ============================================================

pca_bin_names <- as.character(pca_time_bins$pca_bin)
matrix_results <- map(pca_bin_names,
                      ~ build_rate_matrix(
                        pair_data = pairwise_broad_summary,
                        one_bin = .x))
names(matrix_results) <- pca_bin_names

walk(names(matrix_results),
     function(one_bin) {
       safe_name <- one_bin %>%
         str_replace_all(" ", "_") %>%
         str_replace_all("-", "_")
       matrix_out <- as.data.frame(matrix_results[[one_bin]]$rate_matrix) %>%
         rownames_to_column(var = "sample")
       write_tsv(matrix_out, file.path(summary_dir,paste0("exploratory_rate_matrix_all_chr_fig4e",safe_name,".tsv")))
     }
)

# ============================================================
# PLOT PCA with log10 rates
# ============================================================

rate_0_5 <- matrix_results[["0-5 kya"]]$rate_matrix
rate_5_10 <- matrix_results[["5-10 kya"]]$rate_matrix
rate_10_20 <- matrix_results[["10-20 kya"]]$rate_matrix
rate_20_30 <- matrix_results[["20-30 kya"]]$rate_matrix

plot_rate_pca <- function(rate_matrix, time_label) {
  pca_result <- prcomp(rate_matrix, center = TRUE, scale. = TRUE)
  out <- as.data.frame(pca_result$x)
  out <- as.data.frame(pca_result$x) %>%
    rownames_to_column("sample") %>%
    left_join(pop_data %>% select(sample, population, age_group,age_years), by = "sample") %>%
    mutate(population=factor(population, levels=c("E", "D", "BE","PL","B","IRQ")))
   if (mean(out$PC1[out$population == "E"], na.rm = TRUE) > 0) {
    out$PC1 <- -out$PC1
  }
  if (mean(out$PC2[out$population == "E"], na.rm = TRUE) > 0) {
    out$PC2 <- -out$PC2
  variance_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
  ggplot(out, aes(x = PC1, y = PC2, colour = population, shape = age_group)) +
    geom_point(size = 3) +
    scale_shape_manual(values = c("1500-2500" = 16, "200-1000" = 16, "more_than_10k"=16,"1000-1500"=16,"present"=21)) +
    scale_color_manual(values = colors) +
    geom_text(data = out %>%filter(age_group != "present"),aes(label = sample), vjust = -0.7, size = 2, show.legend = FALSE) +
    coord_fixed(ratio=1, xlim=c(-5,4), ylim=c(-4, 3), expand=F) +
    theme_classic() +
    guides(shape="none") +
    labs(title = time_label,x = paste0("PC1 (",round(variance_explained[1], 1),"%)"),
         y = paste0("PC2 (",round(variance_explained[2], 1),"%)"),colour = "Population") +
    theme(plot.title = element_text(hjust = 0.5, size=12), legend.position = "top")
}

p_pca_0_5 <- plot_rate_pca(rate_matrix = rate_0_5,time_label = "0-5 kya")
p_pca_5_10 <- plot_rate_pca(rate_matrix = rate_5_10,time_label = "5-10 kya")
p_pca_10_20 <- plot_rate_pca(rate_matrix = rate_10_20,time_label = "10-20 kya")

p_pca_all <- (p_pca_0_5 + p_pca_5_10 + p_pca_10_20) + plot_layout(ncol = 3,guides = "collect") & theme(legend.position = "right") &
  guides(colour = guide_legend(ncol = 1,byrow = TRUE))
print(p_pca_all)

ggsave(filename = file.path(summary_dir,"PCA_all_chr_pairwise_rate_matrices_fig4.pdf"), plot = p_pca_all, width = 10, height = 4, units = "in")

