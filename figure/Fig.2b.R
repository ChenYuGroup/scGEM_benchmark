# ============================================================
# Raw-score plots: 8 metrics + slide 9 selected 15 panels
#
# Slide 1-8:
#   Each slide = one metric
#
# Slide 9:
#   5 rows × 3 columns
#   Rows = selected metrics
#   Columns = MEM / Data / Threshold
#   Metric names are vertical
#   No overall title
#
# Output:
#   PPTX only
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
if (!requireNamespace("rvg", quietly = TRUE)) install.packages("rvg")

library(officer)
library(rvg)

# ==============================
# 1. File paths
# ==============================

file <- "../scGEM_benchmark/result/accuracy_performance.xlsx"


dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================
# 2. Read Raw_scores
# ==============================

raw <- read_excel(
  file,
  sheet = "Raw_scores",
  col_names = FALSE
)

data_raw <- raw[-c(1, 2), ]

meta <- data_raw[, 1:4]
colnames(meta) <- c("Strategy", "MEM", "Data", "Threshold")

meta <- as.data.frame(meta) %>%
  mutate(
    Strategy  = as.character(Strategy),
    MEM       = as.character(MEM),
    Data      = as.character(Data),
    Threshold = as.character(Threshold)
  )

# ==============================
# 3. Header information
# ==============================

metric_header <- as.character(unlist(raw[1, 5:ncol(raw)]))
dataset_header <- as.character(unlist(raw[2, 5:ncol(raw)]))

for (i in seq_along(metric_header)) {
  if ((is.na(metric_header[i]) || metric_header[i] == "") && i > 1) {
    metric_header[i] <- metric_header[i - 1]
  }
}

header_df <- data.frame(
  Col     = 5:ncol(raw),
  Metric  = metric_header,
  Dataset = dataset_header,
  stringsAsFactors = FALSE
)

# ==============================
# 4. Metric order
# ==============================

metric_order <- c(
  "Housekeeping gene (2019)",
  "Housekeeping gene (2022)",
  "Housekeeping gene (2025)",
  "Metabolic task",
  "Marker gene",
  "Cluster performance (ARI)",
  "Cluster performance (NMI)",
  "Essential gene"
)

# Keep original metric names
metric_short <- metric_order
names(metric_short) <- metric_order

# ==============================
# 5. Convert to long format
# ==============================

all_long <- lapply(metric_order, function(met) {
  
  cols_use <- header_df %>%
    filter(Metric == met) %>%
    pull(Col)
  
  dataset_names <- header_df %>%
    filter(Metric == met) %>%
    pull(Dataset)
  
  dataset_names <- as.character(dataset_names)
  dataset_names <- make.unique(dataset_names)
  
  tmp <- bind_cols(
    meta,
    data_raw[, cols_use]
  )
  
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
      Metric = met,
      Metric_short = metric_short[met],
      Score = suppressWarnings(as.numeric(Score))
    ) %>%
    filter(
      !is.na(Score),
      is.finite(Score)
    )
}) %>%
  bind_rows()

# ==============================
# 6. Orders and colors
# ==============================

# Default dataset order for metrics other than marker-gene coverage
dataset_order_all <- c(
  "LUAD", "PAAD", "AML", "T2D", "BRCA",
  "MESO", "Cell Lines", "CellLines", "ESC", "SKCM",
  "HCC827", "MCF7", "A549", "HCT116", "K562"
)

# Dataset labels and order specifically for marker-gene coverage
marker_dataset_order <- c(
  "LUAD-1",
  "BRCA",
  "MESO",
  "LUAD-2",
  "HCT116",
  "IMR90",
  "GM12878",
  "H1-1",
  "H1-2",
  "DEC"
)

mem_order <- c("ftINIT", "GIMME", "iMAT", "rFASTCORMICS")
data_order <- c("R", "LS")
threshold_order <- c("25th", "50th", "75th", "90th")

group_colors <- c(
  "ftINIT"        = "#6A51A3",
  "GIMME"         = "#9E9AC8",
  "iMAT"          = "#74C476",
  "rFASTCORMICS"  = "#FDAE6B",
  "R"             = "#BDBDBD",
  "LS"            = "#F46D43",
  "25th"          = "#3182BD",
  "50th"          = "#6BAED6",
  "75th"          = "#9ECAE1",
  "90th"          = "#C6DBEF"
)

# ==============================
# 7. Plot function for one factor
# ==============================

make_one_factor_plot <- function(df_metric, factor_name) {
  
  if (factor_name == "MEM") {
    plot_df <- df_metric %>%
      mutate(Group = MEM) %>%
      filter(Group %in% mem_order) %>%
      mutate(Group = factor(Group, levels = mem_order))
    colors_use <- group_colors[mem_order]
    
  } else if (factor_name == "Data") {
    plot_df <- df_metric %>%
      mutate(Group = Data) %>%
      filter(Group %in% data_order) %>%
      mutate(Group = factor(Group, levels = data_order))
    colors_use <- group_colors[data_order]
    
  } else if (factor_name == "Threshold") {
    plot_df <- df_metric %>%
      mutate(Group = Threshold) %>%
      filter(Group %in% threshold_order) %>%
      mutate(Group = factor(Group, levels = threshold_order))
    colors_use <- group_colors[threshold_order]
  }
  
  metric_name <- unique(as.character(plot_df$Metric))
  available_datasets <- unique(as.character(plot_df$Dataset))
  
  if (length(metric_name) == 1 && metric_name == "Marker gene") {
    dataset_order <- marker_dataset_order[
      marker_dataset_order %in% available_datasets
    ]
  } else {
    dataset_order <- dataset_order_all[
      dataset_order_all %in% available_datasets
    ]
  }
  
  # Append unexpected datasets to prevent accidental omission
  dataset_order <- unique(c(
    dataset_order,
    setdiff(available_datasets, dataset_order)
  ))
  
  plot_df <- plot_df %>%
    mutate(Dataset = factor(Dataset, levels = dataset_order))
  
  plot_sum <- plot_df %>%
    group_by(Dataset, Group) %>%
    summarise(
      Mean = mean(Score, na.rm = TRUE),
      SEM  = sd(Score, na.rm = TRUE) / sqrt(n()),
      N    = n(),
      .groups = "drop"
    )
  
  pd <- position_dodge(width = 0.78)
  
  ggplot() +
    geom_col(
      data = plot_sum,
      aes(x = Dataset, y = Mean, fill = Group),
      position = pd,
      width = 0.68,
      color = "grey35",
      linewidth = 0.25
    ) +
    geom_errorbar(
      data = plot_sum,
      aes(
        x = Dataset,
        ymin = Mean - SEM,
        ymax = Mean + SEM,
        group = Group
      ),
      position = pd,
      width = 0.18,
      linewidth = 0.30
    ) +
    geom_point(
      data = plot_df,
      aes(x = Dataset, y = Score, color = Group),
      position = position_jitterdodge(
        jitter.width = 0.12,
        dodge.width = 0.78
      ),
      size = 0.75,
      alpha = 0.55
    ) +
    scale_fill_manual(
      name = factor_name,
      values = colors_use,
      drop = FALSE
    ) +
    scale_color_manual(
      name = factor_name,
      values = colors_use,
      drop = FALSE
    ) +
    guides(
      color = "none",
      fill = guide_legend(
        title = factor_name,
        nrow = 1,
        override.aes = list(alpha = 0.9, color = "grey35")
      )
    ) +
    labs(
      title = factor_name,
      x = NULL,
      y = "Raw score"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
      panel.border = element_rect(color = "grey45", linewidth = 0.55),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(face = "bold", size = 10),
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.40, "cm"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

# ==============================
# 8. Plot function for slides 1-8
# ==============================

make_metric_slide_plot <- function(metric_name) {
  
  df_metric <- all_long %>%
    filter(Metric == metric_name)
  
  metric_title <- unique(df_metric$Metric_short)
  
  p_mem <- make_one_factor_plot(df_metric, "MEM")
  p_data <- make_one_factor_plot(df_metric, "Data")
  p_threshold <- make_one_factor_plot(df_metric, "Threshold")
  
  p_mem / p_data / p_threshold +
    plot_layout(heights = c(1, 1, 1)) +
    plot_annotation(title = metric_title) &
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
}

# ==============================
# 9. Helper functions for slide 9
# ==============================

make_text_panel <- function(label = "", size = 10, face = "bold", angle = 0) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = label,
      size = size / 2.8,
      fontface = face,
      angle = angle,
      hjust = 0.5,
      vjust = 0.5
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}

make_one_factor_plot_small <- function(df_metric, factor_name, show_y_axis = FALSE) {
  
  if (factor_name == "MEM") {
    plot_df <- df_metric %>%
      mutate(Group = MEM) %>%
      filter(Group %in% mem_order) %>%
      mutate(Group = factor(Group, levels = mem_order))
    colors_use <- group_colors[mem_order]
    
  } else if (factor_name == "Data") {
    plot_df <- df_metric %>%
      mutate(Group = Data) %>%
      filter(Group %in% data_order) %>%
      mutate(Group = factor(Group, levels = data_order))
    colors_use <- group_colors[data_order]
    
  } else if (factor_name == "Threshold") {
    plot_df <- df_metric %>%
      mutate(Group = Threshold) %>%
      filter(Group %in% threshold_order) %>%
      mutate(Group = factor(Group, levels = threshold_order))
    colors_use <- group_colors[threshold_order]
  }
  
  metric_name <- unique(as.character(plot_df$Metric))
  available_datasets <- unique(as.character(plot_df$Dataset))
  
  if (length(metric_name) == 1 && metric_name == "Marker gene") {
    dataset_order <- marker_dataset_order[
      marker_dataset_order %in% available_datasets
    ]
  } else {
    dataset_order <- dataset_order_all[
      dataset_order_all %in% available_datasets
    ]
  }
  
  dataset_order <- unique(c(
    dataset_order,
    setdiff(available_datasets, dataset_order)
  ))
  
  plot_df <- plot_df %>%
    mutate(Dataset = factor(Dataset, levels = dataset_order))
  
  plot_sum <- plot_df %>%
    group_by(Dataset, Group) %>%
    summarise(
      Mean = mean(Score, na.rm = TRUE),
      SEM  = sd(Score, na.rm = TRUE) / sqrt(n()),
      N    = n(),
      .groups = "drop"
    )
  
  pd <- position_dodge(width = 0.72)
  
  ggplot() +
    geom_col(
      data = plot_sum,
      aes(x = Dataset, y = Mean, fill = Group),
      position = pd,
      width = 0.62,
      color = "grey35",
      linewidth = 0.20
    ) +
    geom_errorbar(
      data = plot_sum,
      aes(
        x = Dataset,
        ymin = Mean - SEM,
        ymax = Mean + SEM,
        group = Group
      ),
      position = pd,
      width = 0.14,
      linewidth = 0.25
    ) +
    geom_point(
      data = plot_df,
      aes(x = Dataset, y = Score, color = Group),
      position = position_jitterdodge(
        jitter.width = 0.10,
        dodge.width = 0.72
      ),
      size = 0.45,
      alpha = 0.50
    ) +
    scale_fill_manual(values = colors_use, drop = FALSE) +
    scale_color_manual(values = colors_use, drop = FALSE) +
    guides(fill = "none", color = "none") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 7) +
    theme(
      plot.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.25),
      panel.border = element_rect(color = "grey45", linewidth = 0.45),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 4.8),
      axis.text.y = if (show_y_axis) element_text(size = 5.2) else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line() else element_blank(),
      axis.ticks.x = element_line(linewidth = 0.25),
      axis.title = element_blank(),
      plot.margin = margin(1.5, 1.5, 1.5, 1.5)
    )
}

# ==============================
# 10. Create PPTX with 9 slides
# ==============================

ppt_file <- file.path(
  out_dir,
  "Raw_scores_8_metrics_plus_slide9_vertical_labels.pptx"
)

doc <- read_pptx()
slide_size <- slide_size(doc)

# ---------- Slides 1-8 ----------
for (met in metric_order) {
  
  p_slide <- make_metric_slide_plot(met)
  
  doc <- add_slide(
    doc,
    layout = "Blank",
    master = "Office Theme"
  )
  
  doc <- ph_with(
    doc,
    value = rvg::dml(ggobj = p_slide),
    location = ph_location(
      left = 0.15,
      top = 0.15,
      width = slide_size$width - 0.30,
      height = slide_size$height - 0.30
    )
  )
}

# ---------- Slide 9 ----------
doc <- add_slide(
  doc,
  layout = "Blank",
  master = "Office Theme"
)

# Figure size:
# 19.05 cm × (27.517 cm × 2/3)
use_w <- 19.05 / 2.54
use_h <- (25.517 * 2 / 3) / 2.54

# Center on the slide
left_pos <- (slide_size$width - use_w) / 2
top_pos  <- (slide_size$height - use_h) / 2

doc <- ph_with(
  doc,
  value = rvg::dml(ggobj = p_slide9),
  location = ph_location(
    left   = left_pos,
    top    = top_pos,
    width  = use_w,
    height = use_h
  )
)

print(doc, target = ppt_file)

message("PPTX saved to: ", ppt_file)
