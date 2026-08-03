# ============================================================
# Model size panels with raw data points + mean bars
#
# Page 1: Reactions / Genes / Metabolites
# Page 2: Reactions / Genes only
#
# No Overall panel
# No rank labels
#
# Strategy order:
# Ordered by mean reaction number from high to low
#
# Save PowerPoint only
# Figure width = half of A4 paper width
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
# 2. Read Excel 
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
# 3. Define model size columns
# ==============================

metric_cols <- list(
  "Reactions"   = 5:13,
  "Genes"       = 14:22,
  "Metabolites" = 23:31
)

panel_order_all <- c("Reactions", "Genes", "Metabolites")
panel_order_rg  <- c("Reactions", "Genes")

# ==============================
# 4. Convert wide table to long table
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

df_all <- bind_rows(df_list)

# ==============================
# 5. Filter empty / invalid values
# ==============================

df <- df_all %>%
  filter(
    !is.na(Value),
    is.finite(Value)
  )

# ==============================
# 6. Calculate mean value
# ==============================

metric_value <- df %>%
  group_by(
    Strategy,
    MEM,
    Data,
    Threshold,
    Metric
  ) %>%
  summarise(
    MeanValue = mean(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    !is.na(MeanValue),
    is.finite(MeanValue)
  )

# ==============================
# 7. Order strategies by reaction number from high to low
# ==============================

strategy_order <- metric_value %>%
  filter(Metric == "Reactions") %>%
  arrange(desc(MeanValue)) %>%
  pull(Strategy)

# ==============================
# 8. Function to make model size plot
# ==============================

make_model_size_plot <- function(
    panel_order,
    plot_title = "Model size"
) {
  
  df_sub <- df %>%
    filter(Metric %in% panel_order)
  
  metric_value_sub <- metric_value %>%
    filter(Metric %in% panel_order)
  
  raw_points <- df_sub %>%
    mutate(
      Panel = Metric,
      Strategy = factor(
        Strategy,
        levels = rev(strategy_order)
      ),
      Panel = factor(
        Panel,
        levels = panel_order
      )
    )
  
  mean_plot <- metric_value_sub %>%
    mutate(
      Panel = Metric,
      Strategy = factor(
        Strategy,
        levels = rev(strategy_order)
      ),
      Panel = factor(
        Panel,
        levels = panel_order
      )
    )
  
  range_df <- df_sub %>%
    group_by(
      Strategy,
      Metric
    ) %>%
    summarise(
      xmin = min(Value, na.rm = TRUE),
      xmax = max(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(
      is.finite(xmin),
      is.finite(xmax)
    ) %>%
    mutate(
      Panel = Metric,
      Strategy = factor(
        Strategy,
        levels = rev(strategy_order)
      ),
      Panel = factor(
        Panel,
        levels = panel_order
      )
    )
  
  p <- ggplot() +
    
    geom_segment(
      data = range_df,
      aes(
        x = xmin,
        xend = xmax,
        y = Strategy,
        yend = Strategy
      ),
      linewidth = 0.25,
      color = "grey55",
      alpha = 0.8
    ) +
    
    geom_col(
      data = mean_plot,
      aes(
        x = MeanValue,
        y = Strategy
      ),
      width = 0.58,
      fill = "#8E7CC3",
      color = "grey35",
      linewidth = 0.25
    ) +
    
    geom_point(
      data = raw_points,
      aes(
        x = Value,
        y = Strategy
      ),
      size = 0.65,
      shape = 16,
      color = "black",
      alpha = 0.75,
      position = position_jitter(
        height = 0.12,
        width = 0
      )
    ) +
    
    facet_grid(
      . ~ Panel,
      scales = "free_x",
      space = "fixed"
    ) +
    
    labs(
      title = plot_title,
      x = "Count",
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
        size = 7,
        family = "Helvetica"
      ),
      strip.background = element_rect(
        fill = "#D9D5E5",
        color = "grey40",
        linewidth = 0.35
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
        linewidth = 0.20
      ),
      panel.border = element_rect(
        color = "grey45",
        linewidth = 0.35
      ),
      axis.text.y = element_text(
        size = 5.5,
        family = "Helvetica",
        color = "black"
      ),
      axis.text.x = element_text(
        size = 5,
        family = "Helvetica",
        angle = 45,
        hjust = 1,
        color = "black"
      ),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(
        size = 6,
        family = "Helvetica",
        face = "bold"
      ),
      legend.position = "none",
      panel.spacing.x = unit(0.10, "cm"),
      plot.margin = margin(3, 3, 3, 3)
    )
  
  return(p)
}

# ==============================
# 9. Generate two plots
# ==============================

p_all <- make_model_size_plot(
  panel_order = panel_order_all,
  plot_title = "Model size"
)

p_rg <- make_model_size_plot(
  panel_order = panel_order_rg,
  plot_title = "Model size"
)

print(p_all)
print(p_rg)

# ==============================
# 10. Save two-page PowerPoint only
# Figure width = half of A4 paper width
# Figure height = half of A4 paper height
# ==============================

# A4 width = 21 cm
# Half A4 width = 10.5 cm
fig_width <- 9 / 2.54

# A4 height = 29.7 cm
# Half A4 height = 14.85 cm
fig_height <- 14.85 / 2.54

ppt_file <- file.path(
  out_dir,
  "Figure_model_size_two_pages_all_and_reactions_genes_half_A4_size.pptx"
)

ppt <- read_pptx()

# Page 1: Reactions / Genes / Metabolites
ppt <- ppt %>%
  add_slide(
    layout = "Blank",
    master = "Office Theme"
  ) %>%
  ph_with(
    dml(ggobj = p_all),
    location = ph_location(
      left = 0,
      top = 0,
      width = fig_width,
      height = fig_height
    )
  )

# Page 2: Reactions / Genes only
ppt <- ppt %>%
  add_slide(
    layout = "Blank",
    master = "Office Theme"
  ) %>%
  ph_with(
    dml(ggobj = p_rg),
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

cat("PowerPoint saved to:\n", ppt_file, "\n")