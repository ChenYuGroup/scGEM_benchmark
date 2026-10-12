# ============================================================
# Combined factor-level summary plots
# A: Model extraction method
# B: Gene expression threshold
# C: Preprocessing method
#
# Y-axis levels are ordered by Overall score from high to low
# No A/B/C panel labels
# Save PowerPoint only
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
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
# 2. Read Sheet1
# ==============================

raw <- read_excel(
  file,
  sheet = "Norm_scores",
  col_names = FALSE
)

data_raw <- raw[-c(1, 2), ]

meta <- data_raw[, 1:4]
colnames(meta) <- c("Strategy", "MEM", "Data", "Threshold")
meta <- as.data.frame(meta)

# ==============================
# 3. Define metric columns
# ==============================

metric_cols <- list(
  "Housekeeping genes"     = 5:13,
  "Metabolic tasks"        = 14:22,
  "Marker genes"           = 23:32,
  "Clustering performance" = 33:41,
  "Essential genes"        = 42:46
)

metric_order <- c(
  "Overall",
  "Housekeeping genes",
  "Metabolic tasks",
  "Marker genes",
  "Clustering performance",
  "Essential genes"
)

# ==============================
# 4. Convert Sheet1 to long format
# ==============================

df_list <- lapply(names(metric_cols), function(m) {
  
  cols <- metric_cols[[m]]
  
  tmp <- bind_cols(
    meta,
    data_raw[, cols]
  )
  
  dataset_names <- as.character(unlist(raw[2, cols]))
  dataset_names <- make.unique(dataset_names)
  
  colnames(tmp) <- c(
    "Strategy",
    "MEM",
    "Data",
    "Threshold",
    dataset_names
  )
  
  tmp %>%
    pivot_longer(
      cols = all_of(dataset_names),
      names_to = "Dataset",
      values_to = "NormScore"
    ) %>%
    mutate(
      Metric = m,
      NormScore = as.numeric(NormScore)
    )
})

df_long <- bind_rows(df_list) %>%
  filter(
    !is.na(NormScore),
    is.finite(NormScore)
  )

# ==============================
# 5. Strategy-level metric scores
# ==============================

strategy_metric <- df_long %>%
  group_by(Strategy, MEM, Data, Threshold, Metric) %>%
  summarise(
    Score = mean(NormScore, na.rm = TRUE),
    .groups = "drop"
  )

strategy_overall <- strategy_metric %>%
  group_by(Strategy, MEM, Data, Threshold) %>%
  summarise(
    Score = mean(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Metric = "Overall"
  )

strategy_score <- bind_rows(
  strategy_overall,
  strategy_metric
) %>%
  mutate(
    Metric = factor(Metric, levels = metric_order),
    MEM = factor(MEM, levels = c("ftINIT", "GIMME", "iMAT", "rFASTCORMICS")),
    Data = factor(Data, levels = c("R", "LS")),
    Threshold = factor(
      Threshold,
      levels = c("25th", "50th", "75th", "90th", "No threshold")
    )
  )

# ==============================
# 6. General plotting function
# Levels are ordered by Overall score from high to low
# ==============================

make_factor_plot <- function(data, factor_col, title_text) {
  
  plot_points <- data %>%
    mutate(
      Level = as.character(.data[[factor_col]])
    ) %>%
    filter(
      !is.na(Level),
      Level != ""
    ) %>%
    mutate(
      Metric = factor(Metric, levels = metric_order)
    )
  
  mean_df <- plot_points %>%
    group_by(Level, Metric) %>%
    summarise(
      MeanScore = mean(Score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(Metric) %>%
    arrange(desc(MeanScore), .by_group = TRUE) %>%
    mutate(
      Rank = row_number()
    ) %>%
    ungroup()
  
  level_order <- mean_df %>%
    filter(Metric == "Overall") %>%
    arrange(desc(MeanScore)) %>%
    pull(Level)
  
  plot_points <- plot_points %>%
    mutate(
      Level = factor(Level, levels = rev(level_order))
    )
  
  mean_df <- mean_df %>%
    mutate(
      Level = factor(Level, levels = rev(level_order)),
      RankX = 1.10
    )
  
  p <- ggplot() +
    
    geom_col(
      data = mean_df,
      aes(x = MeanScore, y = Level),
      width = 0.38,
      fill = "#8E7CC3",
      color = "grey35",
      linewidth = 0.25
    ) +
    
    geom_point(
      data = plot_points,
      aes(x = Score, y = Level),
      size = 0.9,
      color = "black",
      alpha = 0.75,
      position = position_jitter(height = 0.10, width = 0)
    ) +
    
    geom_text(
      data = mean_df,
      aes(x = RankX, y = Level, label = Rank),
      color = "#C44536",
      size = 2.1,
      family = "Helvetica",
      fontface = "bold",
      hjust = 0.5
    ) +
    
    facet_grid(
      . ~ Metric,
      scales = "fixed",
      space = "fixed"
    ) +
    
    scale_x_continuous(
      limits = c(0, 1.16),
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      expand = c(0, 0)
    ) +
    
    coord_cartesian(
      clip = "off"
    ) +
    
    labs(
      title = title_text,
      x = "Normalized score",
      y = NULL
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
        size = 6,
        family = "Helvetica"
      ),
      strip.background = element_rect(
        fill = "#D9D5E5",
        color = "grey40",
        linewidth = 0.4
      ),
      strip.text = element_text(
        face = "bold",
        size = 6,
        family = "Helvetica"
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(
        color = "grey88",
        linewidth = 0.2
      ),
      panel.border = element_rect(
        color = "grey45",
        linewidth = 0.4
      ),
      axis.text.y = element_text(
        size = 6,
        family = "Helvetica"
      ),
      axis.text.x = element_text(
        size = 6,
        family = "Helvetica"
      ),
      axis.ticks.x = element_line(
        color = "black",
        linewidth = 0.25
      ),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(
        size = 6,
        family = "Helvetica",
        face = "bold"
      ),
      legend.position = "none",
      panel.spacing.x = unit(0.15, "cm"),
      plot.margin = margin(4, 18, 4, 4)
    )
  
  return(p)
}

# ==============================
# 7. Generate three panels
# ==============================

p_mem <- make_factor_plot(
  data = strategy_score,
  factor_col = "MEM",
  title_text = "Model extraction method"
)

threshold_data <- strategy_score %>%
  filter(
    MEM != "rFASTCORMICS",
    Threshold %in% c("25th", "50th", "75th", "90th")
  )

p_threshold <- make_factor_plot(
  data = threshold_data,
  factor_col = "Threshold",
  title_text = "Gene expression threshold"
)

p_data <- make_factor_plot(
  data = strategy_score,
  factor_col = "Data",
  title_text = "Preprocessing method"
)

# ==============================
# 8. Combine as Figure 3
# No A/B/C labels
# ==============================

figure3 <- p_mem / p_threshold / p_data +
  plot_layout(
    ncol = 1,
    heights = c(0.75, 0.75, 0.50)
  )

print(figure3)

# ==============================
# 9. Save as PowerPoint only
# ==============================

# A4 portrait usable width: 19.05 cm = 7.5 inch
fig_width <- 19.05 / 2.54

# Keep similar aspect ratio to previous 13.5 x 11.5 figure
fig_height <- 7.5 * fig_width / 13.5

ppt <- read_pptx() %>%
  add_slide(layout = "Blank", master = "Office Theme") %>%
  ph_with(
    dml(ggobj = figure3),
    location = ph_location(
      left = 0,
      top = 0,
      width = fig_width,
      height = fig_height
    )
  )

print(
  ppt,
  target = file.path(out_dir, "factor_summary_combined_no_ABC.pptx")
)