# ============================================================
# Generate five metric figures:
# Sheet1: normalized scores -> Overall ranking
# Sheet1: normalized scores -> dataset-level score bars and ranks
#
# Ranking rule:
# Use min_rank(desc(score)); tied scores have the same rank
# Example: 1, 1, 3, 4
#
# Rank labels:
# Red rank labels are aligned in one fixed column near the right side
#
# Save PowerPoint only
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

file <- "../scGEM_benchmark/result/accuracy_performance.xlsx"
out_dir <- "../scGEM_benchmark/R_figures/"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ==============================
# 2. Read Sheet1 only
# ==============================

norm_raw <- read_excel(
  file,
  sheet = "Norm_scores",
  col_names = FALSE
)

# Row 1: metric names
# Row 2: dataset names
# Row 3+: normalized scores
norm_data <- norm_raw[-c(1, 2), ]

meta <- norm_data[, 1:4]
colnames(meta) <- c("Strategy", "MEM", "Data", "Threshold")
meta <- as.data.frame(meta)

# ==============================
# 3. Define metric column ranges
# ==============================

metric_cols <- list(
  "Housekeeping genes"     = 5:13,
  "Metabolic tasks"        = 14:22,
  "Marker genes"           = 23:32,
  "Clustering performance" = 33:41,
  "Essential genes"        = 42:46
)

metric_file_names <- c(
  "Housekeeping genes"     = "Housekeeping_genes",
  "Metabolic tasks"        = "Metabolic_tasks",
  "Marker genes"           = "Marker_genes",
  "Clustering performance" = "Clustering_performance",
  "Essential genes"        = "Essential_genes"
)

# ==============================
# 4. Function to make one metric figure
# ==============================

make_metric_plot <- function(metric_name, cols) {
  
  # ------------------------------
  # 4.1 Extract normalized scores from Sheet1
  # ------------------------------
  
  norm_tmp <- bind_cols(
    meta,
    norm_data[, cols]
  )
  
  norm_dataset_names <- as.character(unlist(norm_raw[2, cols]))
  norm_dataset_names <- make.unique(norm_dataset_names)
  
  colnames(norm_tmp) <- c(
    "Strategy",
    "MEM",
    "Data",
    "Threshold",
    norm_dataset_names
  )
  
  norm_df <- norm_tmp %>%
    pivot_longer(
      cols = all_of(norm_dataset_names),
      names_to = "Dataset",
      values_to = "NormScore"
    ) %>%
    mutate(
      Strategy = as.character(Strategy),
      MEM = as.character(MEM),
      Data = as.character(Data),
      Threshold = as.character(Threshold),
      Dataset = as.character(Dataset),
      NormScore = suppressWarnings(as.numeric(NormScore))
    ) %>%
    filter(
      !is.na(NormScore),
      is.finite(NormScore)
    )
  
  # ------------------------------
  # 4.2 Overall ranking from Sheet1 normalized scores
  # Tied scores share the same rank
  # ------------------------------
  
  strategy_rank <- norm_df %>%
    group_by(
      Strategy,
      MEM,
      Data,
      Threshold
    ) %>%
    summarise(
      MeanNormScore = mean(NormScore, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(
      !is.na(MeanNormScore),
      is.finite(MeanNormScore)
    ) %>%
    arrange(
      desc(MeanNormScore),
      Strategy
    ) %>%
    mutate(
      OverallRank = min_rank(desc(MeanNormScore))
    )
  
  # ------------------------------
  # 4.3 Dataset-level rank using normalized scores
  # Tied scores share the same rank
  # ------------------------------
  
  dataset_rank <- norm_df %>%
    group_by(Dataset) %>%
    arrange(
      desc(NormScore),
      Strategy,
      .by_group = TRUE
    ) %>%
    mutate(
      DatasetRank = min_rank(desc(NormScore))
    ) %>%
    ungroup()
  
  # ------------------------------
  # 4.4 Order strategies by Overall ranking
  # If ranks are tied, use score and strategy name for stable display order
  # ------------------------------
  
  strategy_order <- strategy_rank %>%
    arrange(
      OverallRank,
      desc(MeanNormScore),
      Strategy
    ) %>%
    pull(Strategy)
  
  # ------------------------------
  # 4.5 Dataset panels: Sheet1 normalized scores
  # ------------------------------
  
  dataset_plot <- dataset_rank %>%
    left_join(
      strategy_rank %>%
        select(
          Strategy,
          MeanNormScore,
          OverallRank
        ),
      by = "Strategy"
    ) %>%
    mutate(
      Panel = Dataset,
      PlotScore = NormScore,
      PlotRank = DatasetRank,
      OriginalScore = NormScore,
      PanelType = "Dataset",
      Strategy = factor(
        Strategy,
        levels = rev(strategy_order)
      )
    ) %>%
    select(
      Strategy,
      MEM,
      Data,
      Threshold,
      Panel,
      Dataset,
      PlotScore,
      PlotRank,
      OriginalScore,
      MeanNormScore,
      OverallRank,
      PanelType
    )
  
  # ------------------------------
  # 4.6 Overall panel: Sheet1 mean normalized score
  # ------------------------------
  
  overall_plot <- strategy_rank %>%
    mutate(
      Panel = "Overall",
      Dataset = "Overall",
      PlotScore = MeanNormScore,
      PlotRank = OverallRank,
      OriginalScore = MeanNormScore,
      PanelType = "Overall",
      Strategy = factor(
        Strategy,
        levels = rev(strategy_order)
      )
    ) %>%
    select(
      Strategy,
      MEM,
      Data,
      Threshold,
      Panel,
      Dataset,
      PlotScore,
      PlotRank,
      OriginalScore,
      MeanNormScore,
      OverallRank,
      PanelType
    )
  
  # ------------------------------
  # 4.7 Combine Overall + dataset panels
  # ------------------------------
  
  panel_order <- c(
    "Overall",
    norm_dataset_names
  )
  
  plot_df <- bind_rows(
    overall_plot,
    dataset_plot
  ) %>%
    mutate(
      Panel = factor(
        Panel,
        levels = panel_order
      )
    )
  
  # ------------------------------
  # 4.8 Rank label positions
  # Fixed rank column aligned near the right side
  # ------------------------------
  
  rank_x <- 1.055
  
  label_df <- plot_df %>%
    mutate(
      LabelX = rank_x
    )
  
  # ------------------------------
  # 4.9 Plot
  # ------------------------------
  
  p <- ggplot(
    plot_df,
    aes(
      x = PlotScore,
      y = Strategy
    )
  ) +
    geom_col(
      width = 0.62,
      fill = "#8E7CC3",
      color = "grey35",
      linewidth = 0.30
    ) +
    geom_text(
      data = label_df,
      aes(
        x = LabelX,
        y = Strategy,
        label = PlotRank
      ),
      color = "#C44536",
      size = 2.1,
      family = "Helvetica",
      fontface = "bold",
      hjust = 0.5
    ) +
    facet_grid(
      . ~ Panel,
      scales = "fixed",
      space = "fixed"
    ) +
    scale_x_continuous(
      limits = c(0, 1.12),
      breaks = c(0, 0.5, 1),
      labels = c("", "", ""),
      expand = c(0, 0)
    ) +
    coord_cartesian(
      xlim = c(0, 1.12),
      clip = "off"
    ) +
    labs(
      title = metric_name,
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
        family = "Helvetica",
        lineheight = 0.95
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
      axis.text.x = element_blank(),
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
      panel.spacing.x = unit(0.25, "cm"),
      plot.margin = margin(6, 18, 6, 6)
    )
  
  print(p)
  
  # ------------------------------
  # 4.10 Save PowerPoint only
  # Figure width fitted to A4 portrait printable width: 19.05 cm
  # ------------------------------
  
  out_file <- file.path(
    out_dir,
    paste0(
      "Sheet1_normalized_score_with_Overall_tied_rank_aligned_",
      metric_file_names[metric_name],
      ".pptx"
    )
  )
  
  # A4 portrait usable width = 19.05 cm = 7.5 inch
  slide_width <- 19.05 / 2.54
  
  # Keep original aspect ratio
  original_width <- ifelse(
    metric_name == "Marker genes",
    21,
    17
  )
  
  original_height <- 8.5
  
  slide_height <- original_height * slide_width / original_width
  
  ppt <- read_pptx() %>%
    add_slide(
      layout = "Blank",
      master = "Office Theme"
    ) %>%
    ph_with(
      dml(ggobj = p),
      location = ph_location(
        left = 0,
        top = 0,
        width = slide_width,
        height = slide_height
      )
    )
  
  print(
    ppt,
    target = out_file
  )
  
  return(p)
}

# ==============================
# 5. Generate five figures
# ==============================

plots <- lapply(
  names(metric_cols),
  function(m) {
    make_metric_plot(
      metric_name = m,
      cols = metric_cols[[m]]
    )
  }
)

names(plots) <- names(metric_cols)