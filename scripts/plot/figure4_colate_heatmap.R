library(tidyverse)

# ============================================================
# SETTINGS
# ============================================================

generation_time <- 5.79
CHR_LABEL <- "all_chr"
WINDOW_START_YEARS <- 0
WINDOW_END_YEARS <- 20000

# ============================================================
# PATHS
# ============================================================

dat <- "/PATH
TYP = "wgs"
colate_dir <- file.path(dat, "02_results", "colate")
heatmap_dir <- file.path(colate_dir, "heatmap", TYP)
pop_file <- file.path(heatmap_dir, "colate_samples.txt")

# ============================================================
# READ SAMPLE METADATA
# ============================================================

pop_data <- read_tsv(pop_file,col_types = cols(
    sample = col_character(),
    age_years = col_double(),
    population = col_character()
  )) %>%
  mutate(age_generations = age_years / generation_time)

print(pop_data, n = Inf)

# ============================================================
# FIND PAIRWISE .COAL FILES
# ============================================================

coal_files <- list.files(heatmap_dir, pattern = paste0("^pairwise_.*_vs_.*_", CHR_LABEL, "\\.coal$"),
  full.names = TRUE)

cat("Number of pairwise .coal files:", length(coal_files), "\n")
head(coal_files)

get_pair_ids <- function(file, chr_label = CHR_LABEL) {
  pair_name <- basename(file) %>%
    str_remove("\\.coal$") %>%
    str_remove("^pairwise_") %>%
    str_remove(paste0("_", chr_label, "$"))
  pair_ids <- str_split(pair_name, "_vs_", simplify = TRUE)
  
  if (ncol(pair_ids) != 2) {
    stop(
      "Could not parse target/reference samples from filename:\n",
      basename(file))
  }
  
  tibble(target = pair_ids[, 1],reference = pair_ids[, 2])
}

# example
get_pair_ids(coal_files[1])
read_pairwise_coal <- function(file) {
  lines <- readLines(file)
  lines <- lines[nzchar(trimws(lines))]
  pair_ids <- get_pair_ids(file)
  
  if (length(lines) < 3) {
    warning(
      "Skipping incomplete file with fewer than 3 non-empty lines: ",
      basename(file)
    )
    return(tibble())
  }
  
  generation <- as.numeric(
    strsplit(trimws(lines[2]), "\\s+")[[1]]
  )
  
  if (any(!is.finite(generation)) || length(generation) < 2) {
    warning(
      "Skipping file with invalid generation line: ",
      basename(file)
    )
    return(tibble())
  }
  
  bootstrap_lines <- lines[3:length(lines)]
  
  map_dfr(
    seq_along(bootstrap_lines),
    function(i) {
      
      x <- bootstrap_lines[i]
      
      v <- suppressWarnings(
        as.numeric(
          strsplit(trimws(x), "\\s+")[[1]]
        )
      )
      
      # Require the two metadata fields plus one rate per generation point.
      expected_length <- length(generation) + 2
      
      if (length(v) != expected_length) {
        warning(
          paste0(
            "Skipping malformed bootstrap row.\n",
            "File: ", basename(file), "\n",
            "Bootstrap-line index: ", i, "\n",
            "Expected fields: ", expected_length, "\n",
            "Observed fields: ", length(v), "\n",
            "Line: ", x
          )
        )
        return(tibble())
      }
      
      bootstrap_id <- v[2]
      coal_values <- v[-c(1, 2)]
      
      if (length(coal_values) != length(generation)) {
        warning(
          "Skipping generation/rate mismatch in: ",
          basename(file)
        )
        return(tibble())
      }
      
      tibble(
        target = pair_ids$target,
        reference = pair_ids$reference,
        bootstrap = bootstrap_id,
        column_index = seq_along(generation),
        generation = generation,
        coal_rate = coal_values,
        source_file = basename(file)
      )
    }
  )
}

pairwise_raw <- map_dfr(coal_files, read_pairwise_coal)

cat("Rows:", nrow(pairwise_raw), "\n")
cat("Pairs:", n_distinct(paste(pairwise_raw$target, pairwise_raw$reference)), "\n")

pairwise_epochs <- pairwise_raw %>%
  group_by(target, reference, bootstrap) %>%
  arrange(column_index, .by_group = TRUE) %>%
  mutate(
    epoch_start = generation,
    epoch_end = lead(generation),
    epoch_duration = epoch_end - epoch_start
  ) %>%
  ungroup() %>%
  filter(!is.na(epoch_end), epoch_duration > 0)

pairwise_epochs <- pairwise_epochs %>%
  left_join(pop_data %>%
      select(target = sample,
        target_age_years = age_years,
        target_age_generations = age_generations,
        target_population = population), by = "target") %>%
  left_join(pop_data %>%
      select(reference = sample,
        reference_age_years = age_years,
        reference_age_generations = age_generations,
        reference_population = population), by = "reference")

if (any(is.na(pairwise_epochs$target_age_generations))) {
  warning("Some target IDs were not found in colate_samples.txt")
}

if (any(is.na(pairwise_epochs$reference_age_generations))) {
  warning("Some reference IDs were not found in colate_samples.txt")
}

window_start_gen <- WINDOW_START_YEARS / generation_time
window_end_gen <- WINDOW_END_YEARS / generation_time

pair_bootstrap_summary <- pairwise_epochs %>%
  mutate(pair_entry_time = pmax(target_age_generations,reference_age_generations),
    valid_start = pmax(epoch_start,pair_entry_time,window_start_gen),
    valid_end = pmin(epoch_end,window_end_gen),
    overlap_generations = valid_end - valid_start) %>%
  filter(overlap_generations > 0,is.finite(coal_rate)) %>%
  group_by(target,
    reference,
    target_population,
    reference_population,
    target_age_years,
    reference_age_years,
    bootstrap) %>%
  summarise(weighted_rate = sum(coal_rate * overlap_generations) / sum(overlap_generations),
    duration_generations = sum(overlap_generations),
    .groups = "drop") %>%
  mutate(duration_years = duration_generations * generation_time,
    duration_kya = duration_years / 1000)

pair_summary <- pair_bootstrap_summary %>%
  group_by(target,
    reference,
    target_population,
    reference_population,
    target_age_years,
    reference_age_years) %>%
  summarise(mean_coal_rate = mean(weighted_rate, na.rm = TRUE),
    median_coal_rate = median(weighted_rate, na.rm = TRUE),
    lower_bootstrap = quantile(weighted_rate,0.025,na.rm = TRUE),
    upper_bootstrap = quantile(
      weighted_rate, 0.975, na.rm = TRUE),
    n_bootstraps = sum(!is.na(weighted_rate)),
    duration_kya_used = first(duration_kya),
    .groups = "drop")


# symmetric matrix
heatmap_data <- bind_rows(pair_summary %>%
    transmute(sample_x = target,
      sample_y = reference,
      coal_rate = mean_coal_rate,
      lower_bootstrap,
      upper_bootstrap,
      duration_kya_used,
      pop_x = target_population,
      pop_y = reference_population,
      age_x = target_age_years,
      age_y = reference_age_years),
  pair_summary %>%
    transmute(
      sample_x = reference,
      sample_y = target,
      coal_rate = mean_coal_rate,
      lower_bootstrap,
      upper_bootstrap,
      duration_kya_used,
      pop_x = reference_population,
      pop_y = target_population,
      age_x = reference_age_years,
      age_y = target_age_years))

diagonal_data <- pop_data %>%
  transmute(sample_x = sample,
    sample_y = sample,
    coal_rate = NA_real_,
    lower_bootstrap = NA_real_,
    upper_bootstrap = NA_real_,
    duration_kya_used = NA_real_,
    pop_x = population,
    pop_y = population,
    age_x = age_years,
    age_y = age_years)

heatmap_data <- bind_rows(heatmap_data, diagonal_data)

# sample order
sample_order <- pop_data %>%
  arrange(population, age_years, sample) %>%
  pull(sample)

heatmap_data <- heatmap_data %>%
  mutate(sample_x = factor(sample_x, levels = sample_order),sample_y = factor(sample_y, levels = rev(sample_order))) 

# plot
ggplot(heatmap_data, aes(x = sample_x, y = sample_y, fill = coal_rate)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_viridis_c(trans = "log10",option = "cividis",direction = 1,na.value = "grey90", name = "Mean\ncoal. rate") +
  coord_fixed() +
  labs(x = NULL, y = NULL, title = paste0("Pairwise Colate coalescence rate: ",
      WINDOW_START_YEARS / 1000,
      "–",WINDOW_END_YEARS / 1000," kya"),
    subtitle = paste0("Duration-weighted mean across native Colate epochs; ",CHR_LABEL)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7), axis.ticks = element_blank(), panel.grid = element_blank())

#saveRDS(pair_summary,file.path(heatmap_dir, "pair_summary_all_chr.rds"))
#saveRDS(pair_summary,file.path(heatmap_dir, "pair_summary_peak_chr18.rds"))

### plot both together ###
pair_allchr <- readRDS(file.path(heatmap_dir, "pair_summary_all_chr.rds")) %>%
  mutate(dataset = "All chromosomes") %>%
  filter(!target %in% c("BRW001","SLN001","VKP001")) %>%
  filter(!reference %in% c("BRW001","SLN001","VKP001"))

pair_chr18 <- readRDS(file.path(heatmap_dir, "pair_summary_peak_chr18.rds")) %>%
  mutate(dataset = "chr18 peak") %>%
  filter(!target %in% c("BRW001","SLN001","VKP001")) %>%
  filter(!reference %in% c("BRW001","SLN001","VKP001"))

make_symmetric_heatmap <- function(pair_summary, dataset_name) {
  bind_rows(pair_summary %>%
      transmute(sample_x = target,
        sample_y = reference,
        coal_rate = mean_coal_rate,
        lower_bootstrap,
        upper_bootstrap,
        duration_kya_used,
        dataset = dataset_name),
    pair_summary %>%
      transmute(sample_x = reference,
        sample_y = target,
        coal_rate = mean_coal_rate,
        lower_bootstrap,
        upper_bootstrap,
        duration_kya_used,
        dataset = dataset_name))
}

heatmap_allchr <- make_symmetric_heatmap(pair_allchr,"All chromosomes") 
heatmap_chr18 <- make_symmetric_heatmap(pair_chr18,"chr18 peak") 
pop_data <- pop_data %>%
  filter(!sample %in% c("BRW001","SLN001","VKP001"))

sample_order <- pop_data %>%
  mutate(
    sort_group = case_when(
      age_years == 0 & str_detect(sample, "^E")   ~ 1,
      age_years == 0 & str_detect(sample, "^D")   ~ 2,
      age_years == 0 & str_detect(sample, "^IRQ") ~ 3,
      age_years == 0 & str_detect(sample, "^B")   ~ 4,
      age_years == 0                               ~ 5,
      TRUE                                         ~ 10
    )
  ) %>%
  arrange(population,sort_group, age_years, sample) %>%
  pull(sample)

sample_index <- tibble(
  sample = sample_order,
  sample_number = seq_along(sample_order)
)

heatmap_allchr <- heatmap_allchr %>%
  left_join(sample_index %>% rename(sample_x = sample, x_index = sample_number),by = "sample_x") %>%
  left_join(sample_index %>% rename(sample_y = sample, y_index = sample_number),by = "sample_y") %>%
  filter(x_index < y_index)

heatmap_chr18 <- heatmap_chr18 %>%
  left_join(sample_index %>% rename(sample_x = sample, x_index = sample_number),by = "sample_x") %>%
  left_join(sample_index %>% rename(sample_y = sample, y_index = sample_number),by = "sample_y") %>%
  filter(x_index > y_index)

diagonal_data <- tibble(
  sample_x = sample_order,
  sample_y = sample_order,
  coal_rate = NA_real_,
  lower_bootstrap = NA_real_,
  upper_bootstrap = NA_real_,
  duration_kya_used = NA_real_,
  dataset = "Diagonal"
)

heatmap_combined <- bind_rows(heatmap_allchr,heatmap_chr18,diagonal_data) %>%
  mutate(sample_x = factor(sample_x, levels = sample_order),
    sample_y = factor(sample_y, levels = rev(sample_order)))

ggplot(heatmap_combined, aes(x = sample_x, y = sample_y, fill = coal_rate)) +
  geom_tile(colour = "white",linewidth = 0.15) +
  scale_fill_viridis_c(trans = "log10",option = "cividis",na.value = "grey85",name = "Mean coal. rate to 20kya") +
  coord_fixed() +
  labs(x = NULL,y = NULL) +
    #title = paste0("Pairwise Colate coalescence rates: ",WINDOW_START_YEARS / 1000,"–",WINDOW_END_YEARS / 1000," kya"),
    #subtitle = paste0("Lower triangle: all chromosomes; upper triangle: chr18 peak")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5,size = 7),
    axis.text.y = element_text(size = 7),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top"
  ) +
  guides(
    fill = guide_colourbar(
      direction = "horizontal",
      barwidth = grid::unit(12, "cm"),
      barheight = grid::unit(0.35, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )
  )
  