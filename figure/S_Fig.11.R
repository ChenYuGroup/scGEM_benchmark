# ============================================================
# Pairwise correlations of model size
# Rows    : grouping by MEM / Data / Threshold
# Columns : pairwise model size comparisons
#
# Each subplot has its own x/y axis labels:
# Reaction number, Gene number, Metabolite number
#
# Font settings are consistent with previous figures:
# Helvetica, base size = 6
#
# Save PowerPoint only
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(officer)
library(rvg)
library(patchwork)

# ==============================
# 1. File paths
# ==============================

file <- "../scGEM_benchmark/result/accuracy_performance.xlsx"
out_dir <- "../scGEM_benchmark/R_figures/"


if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ==============================
# 2. Colors
# ==============================

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
# 3. Read 
# ==============================

raw <- read_excel(
  file,
  sheet = "Size",
  col_names = FALSE
)

data_raw <- raw[-c(1, 2), ]

meta <- data_raw[, 1:4]
colnames(meta) <- c("Strategy", "MEM", "Data", "Threshold")
meta <- as.data.frame(meta)

# ==============================
# 4. Define model size columns
# ==============================

metric_cols <- list(
  Reactions   = 5:13,
  Genes       = 14:22,
  Metabolites = 23:31
)

# ==============================
# 5. Convert each model size metric to long format
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
      values_to = "Value"
    ) %>%
    mutate(
      Metric = m,
      Value = suppressWarnings(as.numeric(Value))
    )
})

df_long <- bind_rows(df_list) %>%
  filter(
    !is.na(Value),
    is.finite(Value)
  ) %>%
  mutate(
    MEM = case_when(
      grepl("ftINIT", MEM, ignore.case = TRUE) ~ "ftINIT",
      grepl("GIMME", MEM, ignore.case = TRUE) ~ "GIMME",
      grepl("iMAT", MEM, ignore.case = TRUE) ~ "iMAT",
      grepl("rFAST", MEM, ignore.case = TRUE) ~ "rFASTCORMICS",
      TRUE ~ as.character(MEM)
    ),
    Data = case_when(
      Data %in% c("R", "Raw", "raw") ~ "R",
      Data %in% c("LS", "Linnorm+SAVER", "linnorm+saver") ~ "LS",
      TRUE ~ as.character(Data)
    ),
    Threshold = trimws(as.character(Threshold)),
    Threshold = gsub("\\s+", "", Threshold),
    Threshold = gsub("%", "", Threshold),
    Threshold = case_when(
      Threshold %in% c("25th", "25", "0.25", "0.250", "0.2500") ~ "25th",
      Threshold %in% c("50th", "50", "0.5", "0.50", "0.500", "0.5000") ~ "50th",
      Threshold %in% c("75th", "75", "0.75", "0.750", "0.7500") ~ "75th",
      Threshold %in% c("90th", "90", "0.9", "0.90", "0.900", "0.9000") ~ "90th",
      TRUE ~ Threshold
    )
  )

# ==============================
# 6. Convert to wide format
# ==============================

df_wide <- df_long %>%
  select(
    Strategy,
    MEM,
    Data,
    Threshold,
    Dataset,
    Metric,
    Value
  ) %>%
  pivot_wider(
    names_from = Metric,
    values_from = Value,
    values_fn = mean
  ) %>%
  filter(
    !is.na(Reactions),
    !is.na(Genes),
    !is.na(Metabolites)
  ) %>%
  mutate(
    MEM = factor(
      MEM,
      levels = c("ftINIT", "GIMME", "iMAT", "rFASTCORMICS")
    ),
    Data = factor(
      Data,
      levels = c("R", "LS")
    ),
    Threshold = factor(
      Threshold,
      levels = c("25th", "50th", "75th", "90th")
    )
  )

# ==============================
# 7. Pairwise comparison settings
# ==============================

pair_settings <- data.frame(
  Pair = c(
    "Reactions vs Genes",
    "Reactions vs Metabolites",
    "Genes vs Metabolites"
  ),
  XVar = c(
    "Reactions",
    "Reactions",
    "Genes"
  ),
  YVar = c(
    "Genes",
    "Metabolites",
    "Metabolites"
  ),
  XLabel = c(
    "Number of reactions",
    "Number of reactions",
    "Number of genes"
  ),
  YLabel = c(
    "Number of genes",
    "Number of metabolites",
    "Number of metabolites"
  ),
  stringsAsFactors = FALSE
)

# ==============================
# 8. Function to make one subplot
# ==============================

make_one_plot <- function(group_type, pair_name, x_var, y_var, x_lab, y_lab) {
  
  if (group_type == "MEM") {
    
    plot_df <- df_wide %>%
      mutate(Group = as.character(MEM)) %>%
      filter(!is.na(Group))
    
  } else if (group_type == "Data") {
    
    plot_df <- df_wide %>%
      mutate(Group = as.character(Data)) %>%
      filter(!is.na(Group))
    
  } else if (group_type == "Threshold") {
    
    plot_df <- df_wide %>%
      filter(
        !is.na(Threshold),
        Threshold %in% c("25th", "50th", "75th", "90th")
      ) %>%
      mutate(Group = as.character(Threshold))
  }
  
  r_value <- cor(
    plot_df[[x_var]],
    plot_df[[y_var]],
    method = "spearman",
    use = "complete.obs"
  )
  
  r_label <- paste0("r = ", round(r_value, 3))
  
  p <- ggplot(
    plot_df,
    aes(
      x = .data[[x_var]],
      y = .data[[y_var]]
    )
  ) +
    geom_point(
      aes(color = Group),
      size = 0.80,
      shape = 16,
      alpha = 0.78
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.30,
      color = "grey35",
      fill = "grey80"
    ) +
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = r_label,
      hjust = 1.05,
      vjust = -0.8,
      size = 1.8,
      family = "Helvetica",
      fontface = "bold"
    ) +
    scale_color_manual(
      values = group_colors,
      drop = FALSE
    ) +
    labs(
      title = pair_name,
      x = x_lab,
      y = y_lab,
      color = NULL
    ) +
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 1.5)
      )
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
        size = 5.5,
        family = "Helvetica"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        color = "grey88",
        linewidth = 0.18
      ),
      panel.border = element_rect(
        color = "grey45",
        linewidth = 0.35
      ),
      axis.text = element_text(
        size = 4.8,
        family = "Helvetica",
        color = "black"
      ),
      axis.title.x = element_text(
        size = 5.2,
        family = "Helvetica",
        face = "bold",
        color = "black"
      ),
      axis.title.y = element_text(
        size = 5.2,
        family = "Helvetica",
        face = "bold",
        color = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.25,
        color = "black"
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(
        size = 5.5,
        family = "Helvetica"
      ),
      legend.key.size = unit(0.30, "cm"),
      plot.margin = margin(2, 2, 2, 2)
    )
  
  return(p)
}

# ==============================
# 9. Generate plots by row
# ==============================

make_row <- function(group_type) {
  
  plot_list <- lapply(seq_len(nrow(pair_settings)), function(i) {
    
    make_one_plot(
      group_type = group_type,
      pair_name = pair_settings$Pair[i],
      x_var = pair_settings$XVar[i],
      y_var = pair_settings$YVar[i],
      x_lab = pair_settings$XLabel[i],
      y_lab = pair_settings$YLabel[i]
    )
  })
  
  row_plot <- wrap_plots(
    plot_list,
    nrow = 1,
    guides = "collect"
  ) &
    theme(
      legend.position = "bottom",
      legend.text = element_text(
        size = 5.5,
        family = "Helvetica"
      ),
      legend.key.size = unit(0.30, "cm")
    )
  
  return(row_plot)
}

p_mem <- make_row("MEM")
p_data <- make_row("Data")
p_threshold <- make_row("Threshold")

# ==============================
# 10. Combine 3 rows
# ==============================

p_final <- p_mem / p_data / p_threshold +
  plot_annotation(
    title = "Pairwise correlations of model size",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 7,
        family = "Helvetica"
      )
    )
  )

print(p_final)

# ==============================
# 11. Save as PowerPoint only
# ==============================

fig_width <- 19.05 / 2.54
fig_height <- 14.85 / 2.54

ppt_file <- file.path(
  out_dir,
  "Figure_model_size_pairwise_3x3_each_axis_labeled_same_font.pptx"
)

ppt <- read_pptx() %>%
  add_slide(
    layout = "Blank",
    master = "Office Theme"
  ) %>%
  ph_with(
    dml(ggobj = p_final),
    location = ph_location(
      left = 0,
      top = 0,
      width = fig_width,
      height = fig_height
    )
  )

print(
  ppt,
  target = ppt_file
)