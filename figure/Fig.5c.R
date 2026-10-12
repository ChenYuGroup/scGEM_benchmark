# ============================================================
# Before-after comparison of model-size adjustment
#
# Corrections:
#   1. Before adjustment scores are calculated directly from score_long,
#      and are shared across Reactions / Genes / Metabolites rows.
#
#   2. After adjustment scores are constrained within valid metric ranges:
#      Housekeeping genes     : 0-100
#      Metabolic tasks        : 0-100
#      Marker genes           : 0-100
#      Clustering performance : 0-1
#      Essential genes        : 0-observed maximum
#
# Output:
#   PowerPoint only
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(officer)
library(rvg)

# ==============================
# 1. File paths
# ==============================

file <- "../scGEM_benchmark_NCS/result/accuracy_performance.xlsx"
out_dir <- "../scGEM_benchmark/R_figures/"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ==============================
# 2. Sheet names
# ==============================

score_sheet <- "Raw_scores"

size_sheets <- c(
  "Reactions"   = "rxns",
  "Genes"       = "genes",
  "Metabolites" = "mets"
)

# ==============================
# 3. Metric columns
# ==============================

# Raw_scores 工作表
score_metric_cols <- list(
  "Housekeeping genes"     = 14:22,
  "Metabolic tasks"        = 32:40,
  "Marker genes"           = 41:50,
  "Clustering performance" = 51:59,
  "Essential genes"        = 69:73
)

# rxns、genes、mets 工作表
size_metric_cols <- list(
  "Housekeeping genes"     = 5:13,
  "Metabolic tasks"        = 14:22,
  "Marker genes"           = 23:32,
  "Clustering performance" = 33:41,
  "Essential genes"        = 42:46
)

metric_order <- c(
  "Housekeeping genes",
  "Metabolic tasks",
  "Marker genes",
  "Clustering performance",
  "Essential genes"
)

mem_order <- c(
  "ftINIT",
  "GIMME",
  "iMAT",
  "rFASTCORMICS"
)

size_order <- c(
  "Reactions",
  "Genes",
  "Metabolites"
)

# ==============================
# 4. Read score sheet
# ==============================

read_score_sheet <- function() {
  
  raw <- read_excel(
    file,
    sheet = score_sheet,
    col_names = FALSE
  )
  
  data_raw <- raw[-c(1, 2), ]
  
  meta <- data_raw[, 1:4]
  colnames(meta) <- c("Strategy", "MEM", "Data", "Threshold")
  meta <- as.data.frame(meta)
  
  df_list <- lapply(names(score_metric_cols), function(m) {
    
    cols <- score_metric_cols[[m]]
    
    tmp <- bind_cols(
      meta,
      data_raw[, cols]
    )
    
    dataset_names <- as.character(unlist(raw[2, cols]))
    dataset_names <- make.unique(dataset_names)
    
    colnames(tmp) <- c(
      "Strategy", "MEM", "Data", "Threshold",
      dataset_names
    )
    
    tmp %>%
      pivot_longer(
        cols = all_of(dataset_names),
        names_to = "Dataset",
        values_to = "Score"
      ) %>%
      mutate(
        Strategy = as.character(Strategy),
        MEM = as.character(MEM),
        Data = as.character(Data),
        Threshold = as.character(Threshold),
        Dataset = as.character(Dataset),
        Metric = m,
        Score = suppressWarnings(as.numeric(Score))
      )
  })
  
  bind_rows(df_list) %>%
    filter(
      !is.na(Score),
      is.finite(Score),
      !is.na(MEM),
      MEM != ""
    ) %>%
    mutate(
      MEM = case_when(
        grepl("ftINIT", MEM, ignore.case = TRUE) ~ "ftINIT",
        grepl("GIMME", MEM, ignore.case = TRUE) ~ "GIMME",
        grepl("iMAT", MEM, ignore.case = TRUE) ~ "iMAT",
        grepl("rFAST", MEM, ignore.case = TRUE) ~ "rFASTCORMICS",
        TRUE ~ MEM
      ),
      Metric = factor(Metric, levels = metric_order),
      MEM = factor(MEM, levels = mem_order)
    )
}

# ==============================
# 5. Read one model size sheet
# ==============================

read_size_sheet <- function(sheet_name, size_feature) {
  
  raw <- read_excel(
    file,
    sheet = sheet_name,
    col_names = FALSE
  )
  
  data_raw <- raw[-c(1, 2), ]
  
  meta <- data_raw[, 1:4]
  colnames(meta) <- c("Strategy", "MEM", "Data", "Threshold")
  meta <- as.data.frame(meta)
  
  df_list <- lapply(names(size_metric_cols), function(m) {
    
    cols <- size_metric_cols[[m]]
    
    tmp <- bind_cols(
      meta,
      data_raw[, cols]
    )
    
    dataset_names <- as.character(unlist(raw[2, cols]))
    dataset_names <- make.unique(dataset_names)
    
    colnames(tmp) <- c(
      "Strategy", "MEM", "Data", "Threshold",
      dataset_names
    )
    
    tmp %>%
      pivot_longer(
        cols = all_of(dataset_names),
        names_to = "Dataset",
        values_to = "SizeValue"
      ) %>%
      mutate(
        Strategy = as.character(Strategy),
        MEM = as.character(MEM),
        Data = as.character(Data),
        Threshold = as.character(Threshold),
        Dataset = as.character(Dataset),
        Metric = m,
        SizeFeature = size_feature,
        SizeValue = suppressWarnings(as.numeric(SizeValue))
      )
  })
  
  bind_rows(df_list) %>%
    filter(
      !is.na(SizeValue),
      is.finite(SizeValue),
      !is.na(MEM),
      MEM != ""
    ) %>%
    mutate(
      MEM = case_when(
        grepl("ftINIT", MEM, ignore.case = TRUE) ~ "ftINIT",
        grepl("GIMME", MEM, ignore.case = TRUE) ~ "GIMME",
        grepl("iMAT", MEM, ignore.case = TRUE) ~ "iMAT",
        grepl("rFAST", MEM, ignore.case = TRUE) ~ "rFASTCORMICS",
        TRUE ~ MEM
      ),
      Metric = factor(Metric, levels = metric_order),
      MEM = factor(MEM, levels = mem_order),
      SizeFeature = factor(SizeFeature, levels = size_order)
    )
}

# ==============================
# 6. Read and merge data
# ==============================

score_long <- read_score_sheet()

size_long <- bind_rows(
  read_size_sheet(size_sheets["Reactions"],   "Reactions"),
  read_size_sheet(size_sheets["Genes"],       "Genes"),
  read_size_sheet(size_sheets["Metabolites"], "Metabolites")
)


df <- score_long %>%
  inner_join(
    size_long,
    by = c(
      "Strategy",
      "MEM",
      "Data",
      "Threshold",
      "Dataset",
      "Metric"
    )
  ) %>%
  filter(
    !is.na(Score),
    !is.na(SizeValue),
    is.finite(Score),
    is.finite(SizeValue),
    MEM %in% mem_order
  ) %>%
  mutate(
    Metric = factor(Metric, levels = metric_order),
    MEM = factor(MEM, levels = mem_order),
    SizeFeature = factor(SizeFeature, levels = size_order)
  )

cat("Merged data rows:", nrow(df), "\n")

# ==============================
# 7. Metric bounds
# ==============================

metric_bounds <- score_long %>%
  group_by(Metric) %>%
  summarise(
    Observed_max = max(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Bound_low = 0,
    Bound_high = case_when(
      as.character(Metric) %in% c(
        "Housekeeping genes",
        "Metabolic tasks",
        "Marker genes"
      ) ~ 100,
      as.character(Metric) == "Clustering performance" ~ 1,
      as.character(Metric) == "Essential genes" ~ Observed_max,
      TRUE ~ Observed_max
    )
  ) %>%
  select(
    Metric,
    Bound_low,
    Bound_high
  )

print(metric_bounds)

# ==============================
# 8. Before adjustment
# Important correction:
# Calculate directly from score_long, not from df.
# Then copy the same before-adjustment values to all SizeFeature rows.
# ==============================

before_base <- score_long %>%
  filter(
    !is.na(Score),
    is.finite(Score),
    MEM %in% mem_order
  ) %>%
  group_by(Metric, MEM) %>%
  summarise(
    Score_value = mean(Score, na.rm = TRUE),
    SE = sd(Score, na.rm = TRUE) / sqrt(n()),
    N = n(),
    .groups = "drop"
  ) %>%
  mutate(
    CI_low = Score_value - 1.96 * SE,
    CI_high = Score_value + 1.96 * SE,
    Stage = "Before adjustment"
  )

before_results <- before_base %>%
  crossing(
    SizeFeature = size_order
  ) %>%
  mutate(
    SizeFeature = factor(SizeFeature, levels = size_order)
  ) %>%
  select(
    Metric,
    SizeFeature,
    MEM,
    Stage,
    Score_value,
    SE,
    CI_low,
    CI_high,
    N
  )

# ==============================
# 9. After adjustment
# Linear model: Score ~ MEM + SizeValue
# ==============================

get_adjusted_MEM <- function(dat) {
  
  dat <- dat %>%
    filter(
      !is.na(MEM),
      !is.na(Score),
      !is.na(SizeValue),
      is.finite(Score),
      is.finite(SizeValue)
    ) %>%
    mutate(
      MEM = factor(as.character(MEM), levels = mem_order)
    )
  
  metric_name <- unique(as.character(dat$Metric))[1]
  size_name   <- unique(as.character(dat$SizeFeature))[1]
  mean_size   <- mean(dat$SizeValue, na.rm = TRUE)
  
  present_mems <- mem_order[mem_order %in% unique(as.character(dat$MEM))]
  
  fit <- lm(
    Score ~ MEM + SizeValue,
    data = dat
  )
  
  newdat <- data.frame(
    MEM = factor(present_mems, levels = mem_order),
    SizeValue = mean_size
  )
  
  pred <- predict(
    fit,
    newdata = newdat,
    se.fit = TRUE
  )
  
  data.frame(
    Metric = metric_name,
    SizeFeature = size_name,
    MEM = present_mems,
    Stage = "After adjustment",
    Score_value = as.numeric(pred$fit),
    SE = as.numeric(pred$se.fit),
    CI_low = as.numeric(pred$fit - 1.96 * pred$se.fit),
    CI_high = as.numeric(pred$fit + 1.96 * pred$se.fit),
    N = nrow(dat),
    stringsAsFactors = FALSE
  )
}

after_results <- df %>%
  group_by(Metric, SizeFeature) %>%
  group_split() %>%
  lapply(get_adjusted_MEM) %>%
  bind_rows()

# ==============================
# 10. Bound after-adjustment scores and confidence intervals
# Important correction:
# Prevent adjusted scores from exceeding valid metric ranges.
# ==============================

after_results <- after_results %>%
  left_join(
    metric_bounds,
    by = "Metric"
  ) %>%
  mutate(
    Score_value = pmin(
      pmax(Score_value, Bound_low),
      Bound_high
    ),
    CI_low = pmin(
      pmax(CI_low, Bound_low),
      Bound_high
    ),
    CI_high = pmin(
      pmax(CI_high, Bound_low),
      Bound_high
    )
  ) %>%
  select(
    -Bound_low,
    -Bound_high
  )

# Also bound before-results for safety
before_results <- before_results %>%
  left_join(
    metric_bounds,
    by = "Metric"
  ) %>%
  mutate(
    Score_value = pmin(
      pmax(Score_value, Bound_low),
      Bound_high
    ),
    CI_low = pmin(
      pmax(CI_low, Bound_low),
      Bound_high
    ),
    CI_high = pmin(
      pmax(CI_high, Bound_low),
      Bound_high
    )
  ) %>%
  select(
    -Bound_low,
    -Bound_high
  )

# ==============================
# 11. Combine before and after
# ==============================

compare_results <- bind_rows(
  before_results,
  after_results
) %>%
  mutate(
    Metric = factor(Metric, levels = metric_order),
    MEM = factor(MEM, levels = mem_order),
    SizeFeature = factor(SizeFeature, levels = size_order),
    Stage = factor(
      Stage,
      levels = c(
        "Before adjustment",
        "After adjustment"
      )
    )
  )

print(compare_results)

# ==============================
# 12. Plot
# Independent y-axis for each panel
# ==============================

if (!requireNamespace("ggh4x", quietly = TRUE)) {
  install.packages("ggh4x")
}

library(ggh4x)

p_compare <- ggplot(
  compare_results,
  aes(
    x = MEM,
    y = Score_value,
    fill = Stage
  )
) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62,
    color = "grey35",
    linewidth = 0.25
  ) +
  geom_errorbar(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    position = position_dodge(width = 0.72),
    width = 0.18,
    linewidth = 0.30
  ) +
  ggh4x::facet_grid2(
    rows = vars(SizeFeature),
    cols = vars(Metric),
    scales = "free_y",
    independent = "y"
  ) +
  scale_fill_manual(
    values = c(
      "Before adjustment" = "#C9C9C9",
      "After adjustment"  = "#8E7CC3"
    )
  ) +
  labs(
    title = "Before and after model-size adjustment",
    subtitle = "Before = observed mean score; after = size-adjusted score from Score ~ MEM + model size",
    x = NULL,
    y = "Score",
    fill = NULL
  ) +
  theme_bw(
    base_size = 6,
    base_family = "Helvetica"
  ) +
  theme(
    text = element_text(
      family = "Helvetica",
      size = 6
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 8,
      family = "Helvetica"
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 6,
      family = "Helvetica"
    ),
    strip.background = element_rect(
      fill = "#D9D5E5",
      color = "grey40",
      linewidth = 0.35
    ),
    strip.text.x = element_text(
      face = "bold",
      size = 6,
      family = "Helvetica"
    ),
    strip.text.y = element_text(
      face = "bold",
      size = 6,
      family = "Helvetica",
      angle = 90
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey88",
      linewidth = 0.20
    ),
    panel.border = element_rect(
      color = "grey45",
      linewidth = 0.35
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 5.5,
      family = "Helvetica",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 5.5,
      family = "Helvetica",
      color = "black"
    ),
    axis.title.y = element_text(
      size = 6,
      family = "Helvetica",
      face = "bold"
    ),
    axis.ticks = element_line(
      linewidth = 0.25,
      color = "black"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(
      size = 6,
      family = "Helvetica"
    ),
    legend.key.size = unit(0.35, "cm"),
    panel.spacing.x = unit(0.15, "cm"),
    panel.spacing.y = unit(0.18, "cm"),
    plot.margin = margin(4, 8, 4, 4)
  )

print(p_compare)

# ==============================
# 13. Save three-page PowerPoint
# Page 1: Reactions / Genes / Metabolites
# Page 2: Reactions / Genes only
# Page 3: Metabolites only
# ==============================

make_compare_plot <- function(size_features_keep, plot_title_suffix = NULL) {
  
  plot_data <- compare_results %>%
    filter(SizeFeature %in% size_features_keep) %>%
    mutate(
      SizeFeature = factor(SizeFeature, levels = size_features_keep)
    )
  
  ggplot(
    plot_data,
    aes(x = MEM, y = Score_value, fill = Stage)
  ) +
    geom_col(
      position = position_dodge(width = 0.72),
      width = 0.62,
      color = "grey35",
      linewidth = 0.25
    ) +
    geom_errorbar(
      aes(ymin = CI_low, ymax = CI_high),
      position = position_dodge(width = 0.72),
      width = 0.18,
      linewidth = 0.30
    ) +
    ggh4x::facet_grid2(
      rows = vars(SizeFeature),
      cols = vars(Metric),
      scales = "free_y",
      independent = "y"
    ) +
    scale_fill_manual(
      values = c(
        "Before adjustment" = "#C9C9C9",
        "After adjustment"  = "#8E7CC3"
      )
    ) +
    labs(
      title = paste0("Before and after model-size adjustment", plot_title_suffix),
      x = NULL,
      y = "Score",
      fill = NULL
    ) +
    theme_bw(base_size = 6, base_family = "Helvetica") +
    theme(
      text = element_text(family = "Helvetica", size = 6),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
      strip.background = element_rect(
        fill = "#D9D5E5",
        color = "grey40",
        linewidth = 0.35
      ),
      strip.text.x = element_text(face = "bold", size = 6),
      strip.text.y = element_text(face = "bold", size = 6, angle = 90),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.20),
      panel.border = element_rect(color = "grey45", linewidth = 0.35),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 5.5,
        color = "black"
      ),
      axis.text.y = element_text(size = 5.5, color = "black"),
      axis.title.y = element_text(size = 6, face = "bold"),
      axis.ticks = element_line(linewidth = 0.25, color = "black"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.35, "cm"),
      panel.spacing.x = unit(0.15, "cm"),
      panel.spacing.y = unit(0.18, "cm"),
      plot.margin = margin(4, 8, 4, 4)
    )
}

p_all <- make_compare_plot(
  c("Reactions", "Genes", "Metabolites"),
  ""
)

p_rg <- make_compare_plot(
  c("Reactions", "Genes"),
  " (reactions and genes)"
)

p_met <- make_compare_plot(
  c("Metabolites"),
  " (metabolites)"
)

ppt_file <- file.path(
  out_dir,
  "Figure_before_after_model_size_adjustment_three_pages.pptx"
)

fig_width <- 19 / 2.54
fig_height <- 8 / 2.54

ppt <- read_pptx() %>%
  add_slide(layout = "Blank", master = "Office Theme") %>%
  ph_with(
    dml(ggobj = p_all),
    location = ph_location(
      left = 0, top = 0,
      width = fig_width, height = fig_height
    )
  ) %>%
  add_slide(layout = "Blank", master = "Office Theme") %>%
  ph_with(
    dml(ggobj = p_rg),
    location = ph_location(
      left = 0, top = 0,
      width = fig_width, height = fig_height
    )
  ) %>%
  add_slide(layout = "Blank", master = "Office Theme") %>%
  ph_with(
    dml(ggobj = p_met),
    location = ph_location(
      left = 0, top = 0,
      width = fig_width, height = fig_height
    )
  )

print(ppt, target = ppt_file)

cat("PowerPoint saved to:\n", ppt_file, "\n")