# =============================================================================
# MAPK3 (ERK1) nuclear intensity analysis
# -----------------------------------------------------------------------------
# hfLMNA WT vs R377L nuclei, across three mechanical conditions
# (Rigid-unstrained / Soft-unstrained / Soft-strained), quantified from
# per-cell segmentation output (e.g. QuPath annotation export).
#
# Outputs (written to ./figures):
#   1. main_violin.pdf         - per-cell intensity, WT vs R377L across substrates
#   2. intensity_bins_count.pdf  \_ stacked bar charts of intensity-bin
#   3. intensity_bins_pct.pdf    /  distribution per condition
#   4. above10_barchart.pdf    - modeled proportion of cells with intensity > 10
#   5. supplement_area_vs_intensity.pdf       (6-panel, Genotype x Substrate)
#   6. supplement_perimeter_vs_intensity.pdf  (6-panel, Genotype x Substrate)
#
# Stats: cells are nested within images, nested within wells. All comparisons
# therefore use mixed-effects models with a (1 | ReplicateID/ImageID) random
# term rather than treating individual cells as independent replicates, and
# p-values are Benjamini-Hochberg (BH) adjusted across the reported contrasts.
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)
library(scales)

# ── 0. Config ─────────────────────────────────────────────────────────────
data_file  <- "Fig.6 MAPK3 analysis.xlsx"   # <- update path as needed
figure_dir <- "/figures"
dir.create(figure_dir, showWarnings = FALSE)

condition_levels <- c(
  "hfLMNA_WT Rigid-unstrained", "hfLMNA_WT Soft-unstrained", "hfLMNA_WT Soft-strained",
  "hfLMNA_R377L Rigid-unstrained", "hfLMNA_R377L Soft-unstrained", "hfLMNA_R377L Soft-strained"
)

condition_colors <- c(
  "hfLMNA_WT Rigid-unstrained"    = "chartreuse4",
  "hfLMNA_WT Soft-unstrained"     = "chartreuse3",
  "hfLMNA_WT Soft-strained"       = "chartreuse2",
  "hfLMNA_R377L Rigid-unstrained" = "orangered3",
  "hfLMNA_R377L Soft-unstrained"  = "orangered2",
  "hfLMNA_R377L Soft-strained"    = "orangered1"
)

# ── 1. Load data ──────────────────────────────────────────────────────────
# Reads one sheet per condition and standardizes columns by name (robust to
# column reordering), keeping only true per-cell annotation rows.
read_intensity <- function(tab, cond_label) {
  df <- read_excel(data_file, sheet = tab)
  names(df) <- trimws(names(df))
  
  int_col   <- grep("Intensity subtracted", names(df), value = TRUE)[1]
  area_col  <- grep("^Area", names(df), value = TRUE)[1]
  peri_col  <- grep("^Perimeter", names(df), value = TRUE)[1]
  roi_col   <- grep("^ROI$", names(df), value = TRUE)[1]        # Polygon / Rectangle
  otype_col <- grep("Object type", names(df), value = TRUE)[1]
  
  df %>%
    fill(Well, Picture, .direction = "down") %>%   # guard against merged/blank cells
    transmute(
      Well      = as.character(Well),
      Picture   = as.character(Picture),
      Image     = as.character(Image),              # unique per image already
      ROI       = .data[[roi_col]],
      ObjType   = .data[[otype_col]],
      Area      = suppressWarnings(as.numeric(.data[[area_col]])),
      Perimeter = suppressWarnings(as.numeric(.data[[peri_col]])),
      Intensity = suppressWarnings(as.numeric(.data[[int_col]]))
    ) %>%
    filter(
      ObjType == "Annotation",        # drop summary rows (e.g. "average negative control")
      ROI != "Rectangle",             # drop whole-image background ROI
      !is.na(Well), !is.na(Image), !is.na(Intensity),
      Area <= 1000                    # remove debris / oversized objects
    ) %>%
    mutate(
      Condition   = cond_label,
      ReplicateID = paste(cond_label, Well, sep = "_"),   # well = biological replicate
      ImageID     = Image
    )
}

dat <- bind_rows(
  read_intensity("LMNA WT_rigid_MAPK3",          "hfLMNA_WT Rigid-unstrained"),
  read_intensity("LMNA WT_soft_MAPK3",            "hfLMNA_WT Soft-unstrained"),
  read_intensity("LMNA WT_soft+stretch_MAPK3",    "hfLMNA_WT Soft-strained"),
  read_intensity("LMNA R377L_rigid_MAPK3",        "hfLMNA_R377L Rigid-unstrained"),
  read_intensity("LMNA R377L_soft_MAPK3",         "hfLMNA_R377L Soft-unstrained"),
  read_intensity("LMNA R377L_soft+stretch_MAPK3", "hfLMNA_R377L Soft-strained")
) %>%
  mutate(Condition = factor(Condition, levels = condition_levels))

# Sanity check before running any stats
dat %>%
  group_by(Condition) %>%
  summarise(n_wells = n_distinct(Well), n_images = n_distinct(ImageID), n_cells = n(),
            .groups = "drop") %>%
  print()

# ── Helper: BH-adjusted pairwise contrasts -> plotting-ready data frame ────
# Matches an emmeans::pairs() result to a requested list of comparisons and
# returns group1/group2/p.value/stars, in the order requested (works for both
# "A - B" contrasts from lmer models and "A / B" contrasts from glmer models
# fit on the response scale).
get_significance_df <- function(pairwise_obj, comparisons) {
  res <- as.data.frame(pairwise_obj) %>%
    mutate(
      contrast_clean = gsub("[()]", "", contrast),
      stars = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  
  match_one <- function(pair) {
    candidates <- c(
      paste(pair[1], "-", pair[2]), paste(pair[2], "-", pair[1]),
      paste(pair[1], "/", pair[2]), paste(pair[2], "/", pair[1])
    )
    hit <- res %>% filter(contrast_clean %in% candidates)
    if (nrow(hit) == 0) {
      return(data.frame(group1 = pair[1], group2 = pair[2], p.value = NA_real_, stars = "ns"))
    }
    data.frame(group1 = pair[1], group2 = pair[2], p.value = hit$p.value[1], stars = hit$stars[1])
  }
  
  bind_rows(lapply(comparisons, match_one))
}

# =============================================================================
# 1. Main violin plot — per-cell MAPK3 intensity
# =============================================================================
model <- lmer(Intensity ~ Condition + (1 | ReplicateID/ImageID), data = dat)
cat("\n=== Mixed-effects model: Intensity ~ Condition ===\n")
print(summary(model))

emm <- emmeans(model, ~ Condition)
pairwise_results <- pairs(emm, adjust = "BH")
cat("\n=== Pairwise contrasts (BH-adjusted) ===\n")
print(pairwise_results)

comparisons_main <- list(
  c("hfLMNA_WT Rigid-unstrained", "hfLMNA_WT Soft-strained"),
  c("hfLMNA_WT Rigid-unstrained", "hfLMNA_WT Soft-unstrained"),
  c("hfLMNA_WT Soft-unstrained", "hfLMNA_WT Soft-strained"),
  c("hfLMNA_R377L Rigid-unstrained", "hfLMNA_R377L Soft-strained"),
  c("hfLMNA_R377L Rigid-unstrained", "hfLMNA_R377L Soft-unstrained"),
  c("hfLMNA_R377L Soft-unstrained", "hfLMNA_R377L Soft-strained")
)

sig_main <- get_significance_df(pairwise_results, comparisons_main)
y_max <- max(dat$Intensity, na.rm = TRUE)
sig_main$y.position <- y_max * seq(1.05, by = 0.10, length.out = nrow(sig_main))

p_violin <- ggplot(dat, aes(Condition, Intensity, fill = Condition, color = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.6, linewidth = 0.6) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.3, color = "black", linewidth = 0.5) +
  stat_pvalue_manual(sig_main, label = "stars", xmin = "group1", xmax = "group2",
                     y.position = "y.position", tip.length = 0.01, size = 4) +
  scale_fill_manual(values = condition_colors) +
  scale_color_manual(values = condition_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = NULL, y = "Intensity per nucleus (subtracted by negative control)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11))

print(p_violin)
ggsave(file.path(figure_dir, "main_violin.pdf"), p_violin, width = 7, height = 5, dpi = 300)

# =============================================================================
# 2. Intensity-bin distribution — stacked bar charts (counts and percentages)
# =============================================================================
dat_binned <- dat %>%
  mutate(IntensityBin = cut(
    Intensity,
    breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 100),
    labels = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", ">80"),
    right = FALSE, include.lowest = TRUE
  ))
bin_colors <- c(
  "0-10"  = "#F7F7F7", "10-20" = "#E3E3E3", "20-30" = "#CFCFCF", "30-40" = "#B5B5B5",
  "40-50" = "#9B9B9B", "50-60" = "#818181", "60-70" = "#676767", "70-80" = "#4D4D4D",
  ">80"   = "#6B6B6B"
)

# 2a. Raw counts per bin
p_bins_count <- dat_binned %>%
  count(Condition, IntensityBin) %>%
  ggplot(aes(Condition, n, fill = IntensityBin)) +
  geom_col(position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = bin_colors, name = "Intensity bin") +
  labs(x = NULL, y = "Number of cells") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 11))

print(p_bins_count)
ggsave(file.path(figure_dir, "intensity_bins_count.pdf"), p_bins_count, width = 7, height = 5, dpi = 300)

# 2b. Percentage per bin (better for comparing conditions with unequal n)
bin_order <- dat_binned %>% count(IntensityBin) %>% arrange(desc(n)) %>% pull(IntensityBin)

p_bins_pct <- dat_binned %>%
  count(Condition, IntensityBin) %>%
  group_by(Condition) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(IntensityBin = factor(IntensityBin, levels = bin_order)) %>%
  ggplot(aes(Condition, pct, fill = IntensityBin)) +
  geom_col(position = position_stack(reverse = TRUE), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = bin_colors, name = "Intensity bin") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Percentage of nuclei (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 11))

print(p_bins_pct)
ggsave(file.path(figure_dir, "intensity_bins_pct.pdf"), p_bins_pct, width = 7, height = 5, dpi = 300)

# =============================================================================
# 3. Proportion of cells with intensity > 10 — logistic mixed model + barchart
# =============================================================================
dat_bin2 <- dat_binned %>% mutate(Above10 = as.integer(IntensityBin != "0-10"))

model_above10 <- glmer(Above10 ~ Condition + (1 | ReplicateID/ImageID),
                       data = dat_bin2, family = binomial)
cat("\n=== Logistic mixed model: P(Intensity > 10) ~ Condition ===\n")
print(summary(model_above10))

emm_above10 <- emmeans(model_above10, ~ Condition, type = "response")
pairwise_above10 <- pairs(emm_above10, adjust = "BH")
cat("\n=== Pairwise contrasts, P(Intensity > 10) (BH-adjusted) ===\n")
print(pairwise_above10)

comparisons_within_genotype <- list(
  c("hfLMNA_WT Rigid-unstrained",    "hfLMNA_WT Soft-unstrained"),
  c("hfLMNA_WT Soft-unstrained",     "hfLMNA_WT Soft-strained"),
  c("hfLMNA_WT Rigid-unstrained",    "hfLMNA_WT Soft-strained"),
  c("hfLMNA_R377L Rigid-unstrained", "hfLMNA_R377L Soft-unstrained"),
  c("hfLMNA_R377L Soft-unstrained",  "hfLMNA_R377L Soft-strained"),
  c("hfLMNA_R377L Rigid-unstrained", "hfLMNA_R377L Soft-strained")
)

sig_above10 <- get_significance_df(pairwise_above10, comparisons_within_genotype)

emm_above10_df <- as.data.frame(emm_above10) %>%
  mutate(Condition = factor(Condition, levels = condition_levels))

prob_max <- max(emm_above10_df$asymp.UCL, na.rm = TRUE)
sig_above10$y.position <- prob_max * seq(1.08, by = 0.10, length.out = nrow(sig_above10))

p_above10 <- ggplot(emm_above10_df, aes(Condition, prob, fill = Condition)) +
  geom_col(alpha = 0.8, color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15, linewidth = 0.5) +
  stat_pvalue_manual(sig_above10, label = "stars", xmin = "group1", xmax = "group2",
                     y.position = "y.position", tip.length = 0.01, size = 4) +
  scale_fill_manual(values = condition_colors) +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0),
                     labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "MAPK3 intensity > 10") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 11))

print(p_above10)
ggsave(file.path(figure_dir, "above10_barchart.pdf"), p_above10, width = 7, height = 5, dpi = 300)

# =============================================================================
# Supplement: Correlation 'Area vs Intensity' and 'Perimeter vs Intensity' (6-panel)
# =============================================================================
dat_corr <- dat %>%
  mutate(
    Substrate = case_when(
      grepl("Rigid-unstrained", Condition) ~ "Rigid",
      grepl("Soft-unstrained",  Condition) ~ "Soft",
      grepl("Soft-strained",    Condition) ~ "Stretch"
    ),
    Genotype = case_when(
      grepl("^hfLMNA_WT",    Condition) ~ "WT",
      grepl("^hfLMNA_R377L", Condition) ~ "R377L"
    ),
    Substrate = factor(Substrate, levels = c("Rigid", "Soft", "Stretch"))
  )

genotype_colors <- c("WT" = "chartreuse4", "R377L" = "orangered3")

# Correlation summary (for reporting in text/methods)
dat_corr %>%
  group_by(Genotype, Substrate) %>%
  summarise(
    n_cells = n(),
    r_area  = cor(Area, Intensity, method = "pearson"),
    p_area  = cor.test(Area, Intensity, method = "pearson")$p.value,
    r_peri  = cor(Perimeter, Intensity, method = "pearson"),
    p_peri  = cor.test(Perimeter, Intensity, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  print()

p_supp_area <- ggplot(dat_corr, aes(Area, Intensity)) +
  geom_point(aes(color = Genotype), alpha = 1, size = 1.3) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.6, se = TRUE) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3.2) +
  facet_wrap(~ Genotype + Substrate, nrow = 2) +
  scale_color_manual(values = genotype_colors) +
  labs(x = expression("Area (" * mu * "m"^2 * ")"),
       y = "Intensity per nucleus (subtracted by negative control)",
       color = "Genotype") +
  theme_classic() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 9), axis.title = element_text(size = 10),
        legend.position = "bottom")

print(p_supp_area)
ggsave(file.path(figure_dir, "supplement_area_vs_intensity.pdf"), p_supp_area,
       width = 12, height = 6, dpi = 300)

p_supp_peri <- ggplot(dat_corr, aes(Perimeter, Intensity)) +
  geom_point(aes(color = Genotype), alpha = 1, size = 1.3) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.6, se = TRUE) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3.2) +
  facet_wrap(~ Genotype + Substrate, nrow = 2) +
  scale_color_manual(values = genotype_colors) +
  labs(x = expression("Perimeter (" * mu * "m" * ")"),
       y = "Intensity per nucleus (subtracted by negative control)",
       color = "Genotype") +
  theme_classic() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 9), axis.title = element_text(size = 10),
        legend.position = "bottom")

print(p_supp_peri)
ggsave(file.path(figure_dir, "supplement_perimeter_vs_intensity.pdf"), p_supp_peri,
       width = 12, height = 6, dpi = 300)


sessionInfo()
