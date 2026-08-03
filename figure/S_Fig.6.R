# ============================================================
# Generate five metric figures:
# Sheet1: normalized scores -> Overall ranking
# Sheet2: raw scores -> dataset-level raw score bars and raw ranks
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

out_dir <- "../scGEM_benchmark/result"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ==============================
# 2. Read 
# ==============================

norm_raw <- read_excel(
  file,
  sheet = "Norm_scores",
  col_names = FALSE
)

raw_raw <- read_excel(
  file,
  sheet = "Raw_scores",
  col_names = FALSE
)

norm_data <- norm_raw[-c(1, 2), ]
raw_data  <- raw_raw[-c(1, 2), ]

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
# 4. Axis settings
# ==============================

get_dataset_axis_setting <- function(metric_name, raw_scores) {
  
  if (metric_name %in% c("Housekeeping genes", "Metabolic tasks", "Marker genes")) {
    return(list(xlim = c(0, 100)))
  }
  
  if (metric_name == "Clustering performance") {
    return(list(xlim = c(0, 1)))
  }
  
  if (metric_name == "Essential genes") {
    xmax <- max(raw_scores, na.rm = TRUE)
    
    if (!is.finite(xmax) || xmax <= 0) {
      xmax <- 1
    }
    
    return(list(xlim = c(0, xmax)))
  }
}

# ==============================
# 5. Function to make one metric figure
# ==============================

make_metric_plot <- function(metric_name, cols) {
  
  # ------------------------------
  # 5.1 Extract normalized scores from Sheet1
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
      NormScore = suppressWarnings(as.numeric(NormScore))
    ) %>%
    filter(
      !is.na(NormScore),
      is.finite(NormScore)
    )
  
  # ------------------------------
  # 5.2 Overall ranking from Sheet1 normalized scores
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
  # 5.3 Extract raw scores from Sheet2
  # ------------------------------
  
  raw_meta <- raw_data[, 1:4]
  colnames(raw_meta) <- c("Strategy", "MEM", "Data", "Threshold")
  raw_meta <- as.data.frame(raw_meta)
  
  raw_tmp <- bind_cols(
    raw_meta,
    raw_data[, cols]
  )
  
  raw_dataset_names <- as.character(unlist(raw_raw[2, cols]))
  raw_dataset_names <- make.unique(raw_dataset_names)
  
  colnames(raw_tmp) <- c(
    "Strategy",
    "MEM",
    "Data",
    "Threshold",
    raw_dataset_names
  )
  
  raw_df <- raw_tmp %>%
    pivot_longer(
      cols = all_of(raw_dataset_names),
      names_to = "Dataset",
      values_to = "RawScore"
    ) %>%
    mutate(
      RawScore = suppressWarnings(as.numeric(RawScore))
    ) %>%
    filter(
      !is.na(RawScore),
      is.finite(RawScore)
    )
  
  # ------------------------------
  # 5.4 Raw-score rank within each dataset
  # Tied raw scores share the same rank
  # ------------------------------
  
  raw_df <- raw_df %>%
    group_by(Dataset) %>%
    arrange(
      desc(RawScore),
      Strategy,
      .by_group = TRUE
    ) %>%
    mutate(
      RawRank = min_rank(desc(RawScore))
    ) %>%
    ungroup()
  
  # ------------------------------
  # 5.5 Order strategies by Overall ranking from Sheet1
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
  # 5.6 Axis transformation
  # ------------------------------
  
  dataset_axis <- get_dataset_axis_setting(
    metric_name,
    raw_df$RawScore
  )
  
  dataset_x_min <- dataset_axis$xlim[1]
  dataset_x_max <- dataset_axis$xlim[2]
  
  rescale_to01 <- function(x, xmin, xmax) {
    if (xmax == xmin) {
      return(rep(0, length(x)))
    }
    (x - xmin) / (xmax - xmin)
  }
  
  # ------------------------------
  # 5.7 Dataset panels: Sheet2 raw scores
  # ------------------------------
  
  dataset_plot <- raw_df %>%
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
      PlotScore = rescale_to01(
        RawScore,
        dataset_x_min,
        dataset_x_max
      ),
      PlotRank = RawRank,
      OriginalScore = RawScore,
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
  # 5.8 Overall panel: Sheet1 mean normalized score
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
  # 5.9 Combine Overall + dataset panels
  # ------------------------------
  
  panel_order <- c(
    "Overall",
    raw_dataset_names
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
  # 5.10 Rank label positions
  # Fixed rank column for aligned red rank labels
  # ------------------------------
  
  rank_x <- 1.055
  
  label_df <- plot_df %>%
    mutate(
      LabelX = rank_x
    )
  
  # ------------------------------
  # 5.11 Plot
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
      x = "score",
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
  # 5.12 Save PowerPoint only
  # Figure width fitted to A4 portrait printable width: 19.05 cm
  # ------------------------------
  
  out_file <- file.path(
    out_dir,
    paste0(
      "Raw_score_with_Sheet1_Overall_keep_ticks_tied_rank_aligned_",
      metric_file_names[metric_name],
      ".pptx"
    )
  )
  
  slide_width <- 19.05 / 2.54
  
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
# 6. Generate five figures
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