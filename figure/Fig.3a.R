# ============================================================
# Sensitivity analysis line plots
#
# Sheet1:
# Row 1 = metric names
# Row 2 = correlation levels
# Row 3+ = strategy scores
#
# Outputs:
# Pages 1-9  : one metric per slide
# Page 10    : 9 plots in 3 x 3 layout, one legend
# Page 11    : 6 selected plots in 3 x 2 layout, one legend
#
# Legend:
# 4 columns; each column corresponds to one MEM
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
file <- "../scGEM_benchmark/result/sensitivity_performance.xlsx"

out_dir <- "../scGEM_benchmark/result"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ==============================
# 2. Read Excel
# Read all columns as text to avoid type conflicts
# ==============================

raw <- read_excel(
  file,
  sheet = "Sheet1",
  col_names = FALSE,
  col_types = "text"
)

metric_header <- as.character(unlist(raw[1, ]))
rho_header    <- as.character(unlist(raw[2, ]))

# Fill merged metric names to the right
for (i in seq_along(metric_header)) {
  if (i == 1) next
  
  if (is.na(metric_header[i]) || metric_header[i] == "") {
    metric_header[i] <- metric_header[i - 1]
  }
}

data_raw <- raw[-c(1, 2), ]

# ==============================
# 3. Clean empty columns and build column names
# ==============================

metric_header[1] <- "Strategy"
rho_header[1] <- "Strategy"

# Remove completely empty columns
keep_cols <- !(
  (is.na(metric_header) | metric_header == "") &
    (is.na(rho_header) | rho_header == "")
)

metric_header <- metric_header[keep_cols]
rho_header    <- rho_header[keep_cols]
data_raw      <- data_raw[, keep_cols]

# Remove columns with empty correlation except Strategy column
keep_cols2 <- seq_along(metric_header) == 1 |
  !(is.na(rho_header) | rho_header == "")

metric_header <- metric_header[keep_cols2]
rho_header    <- rho_header[keep_cols2]
data_raw      <- data_raw[, keep_cols2]

new_colnames <- rep(NA, length(metric_header))
new_colnames[1] <- "Strategy"

for (i in 2:length(metric_header)) {
  new_colnames[i] <- paste0(
    metric_header[i],
    "__rho_",
    rho_header[i]
  )
}

colnames(data_raw) <- new_colnames

valid_cols <- !is.na(colnames(data_raw)) & colnames(data_raw) != ""
data_raw <- data_raw[, valid_cols]

data_raw <- data_raw %>%
  mutate(
    across(
      everything(),
      as.character
    )
  )

# ==============================
# 4. Convert to long format
# ==============================

df_long <- data_raw %>%
  mutate(
    Strategy = as.character(Strategy)
  ) %>%
  filter(
    !is.na(Strategy),
    Strategy != ""
  ) %>%
  pivot_longer(
    cols = -Strategy,
    names_to = "Metric_Rho",
    values_to = "Score",
    values_transform = list(Score = as.character)
  ) %>%
  separate(
    Metric_Rho,
    into = c("Metric", "Correlation"),
    sep = "__rho_",
    remove = TRUE
  ) %>%
  mutate(
    Strategy = as.character(Strategy),
    Metric = as.character(Metric),
    Correlation = as.numeric(Correlation),
    Score = as.numeric(Score)
  ) %>%
  filter(
    !is.na(Strategy),
    Strategy != "",
    !is.na(Metric),
    Metric != "",
    !is.na(Correlation),
    !is.na(Score),
    is.finite(Score)
  )

# ==============================
# 5. Order settings
# ==============================

rho_levels <- c(1, 0.9, 0.7, 0.5, 0.3, 0.1, 0)

metric_order <- c(
  "Model similarity",
  "Housekeeping gene (2019)",
  "Housekeeping gene (2022)",
  "Housekeeping gene (2025)",
  "Metabolic task",
  "Marker gene",
  "Clustering performance (ARI)",
  "Clustering performance (NMI)",
  "Essential gene"
)

metric_order <- metric_order[metric_order %in% unique(df_long$Metric)]

# Important:
# This order is designed for guide_legend(ncol = 4, byrow = FALSE)
# so each legend column corresponds to one MEM.
preferred_strategy_order <- c(
  "ftINIT-25th", "ftINIT-50th", "ftINIT-75th", "ftINIT-90th",
  "GIMME-25th", "GIMME-50th", "GIMME-75th", "GIMME-90th",
  "iMAT-25th", "iMAT-50th", "iMAT-75th", "iMAT-90th",
  "rFASTCORMICS-LS"
)

actual_strategies <- unique(df_long$Strategy)

strategy_order <- c(
  preferred_strategy_order[preferred_strategy_order %in% actual_strategies],
  actual_strategies[!actual_strategies %in% preferred_strategy_order]
)

df_long <- df_long %>%
  mutate(
    Strategy = factor(Strategy, levels = strategy_order),
    Metric = factor(Metric, levels = metric_order)
  )

# ==============================
# 6. Strategy colors
# ==============================

strategy_cols <- c(
  "ftINIT-25th" = "#6A51A3",
  "ftINIT-50th" = "#807DBA",
  "ftINIT-75th" = "#9E9AC8",
  "ftINIT-90th" = "#BCBDDC",
  
  "GIMME-25th" = "#238B45",
  "GIMME-50th" = "#41AB5D",
  "GIMME-75th" = "#74C476",
  "GIMME-90th" = "#A1D99B",
  
  "iMAT-25th" = "#CB181D",
  "iMAT-50th" = "#EF3B2C",
  "iMAT-75th" = "#FB6A4A",
  "iMAT-90th" = "#FC9272",
  
  "rFASTCORMICS-LS" = "#FDAE6B"
)

missing_strategies <- setdiff(strategy_order, names(strategy_cols))

if (length(missing_strategies) > 0) {
  extra_cols <- grDevices::hcl.colors(
    length(missing_strategies),
    palette = "Dark 3"
  )
  names(extra_cols) <- missing_strategies
  strategy_cols <- c(strategy_cols, extra_cols)
}

# ==============================
# 7. Y-axis settings by metric
# ==============================

get_y_scale <- function(metric_name, scores) {
  
  metric_name <- as.character(metric_name)
  
  if (metric_name %in% c(
    "Housekeeping gene (2019)",
    "Housekeeping gene (2022)",
    "Housekeeping gene (2025)",
    "Metabolic task",
    "Marker gene"
  )) {
    return(
      scale_y_continuous(
        limits = c(0, 100),
        breaks = c(0, 25, 50, 75, 100),
        expand = expansion(mult = c(0, 0.02))
      )
    )
  }
  
  if (metric_name %in% c(
    "Model similarity",
    "Clustering performance (ARI)",
    "Clustering performance (NMI)"
  )) {
    return(
      scale_y_continuous(
        limits = c(0, 1),
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        expand = expansion(mult = c(0, 0.02))
      )
    )
  }
  
  if (metric_name == "Essential gene") {
    
    ymin <- min(scores, na.rm = TRUE)
    ymax <- max(scores, na.rm = TRUE)
    
    if (!is.finite(ymin)) ymin <- 0
    if (!is.finite(ymax)) ymax <- 1
    
    if (ymin == ymax) {
      ymin <- ymin - 0.05
      ymax <- ymax + 0.05
    }
    
    return(
      scale_y_continuous(
        limits = c(ymin, ymax),
        breaks = pretty(c(ymin, ymax), n = 5),
        expand = expansion(mult = c(0.03, 0.05))
      )
    )
  }
  
  return(
    scale_y_continuous(
      breaks = pretty(scores, n = 5),
      expand = expansion(mult = c(0.03, 0.05))
    )
  )
}

# ==============================
# 8. Function to make one line plot
# ==============================

make_sensitivity_plot <- function(metric_name) {
  
  plot_df <- df_long %>%
    filter(Metric == metric_name) %>%
    droplevels()
  
  y_scale_obj <- get_y_scale(
    metric_name = metric_name,
    scores = plot_df$Score
  )
  
  p <- ggplot(
    plot_df,
    aes(
      x = Correlation,
      y = Score,
      group = Strategy,
      color = Strategy
    )
  ) +
    geom_line(
      linewidth = 0.45,
      alpha = 0.95
    ) +
    geom_point(
      size = 1.0,
      alpha = 0.95
    ) +
    scale_color_manual(
      values = strategy_cols,
      breaks = strategy_order,
      drop = FALSE,
      guide = guide_legend(
        ncol = 4,
        byrow = FALSE,
        override.aes = list(
          linewidth = 0.6,
          size = 1.2
        )
      )
    ) +
    scale_x_reverse(
      breaks = rho_levels,
      limits = c(1, 0)
    ) +
    y_scale_obj +
    labs(
      title = as.character(metric_name),
      x = "Spearman R",
      y = "Score",
      color = NULL
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
      axis.title = element_text(
        size = 6,
        family = "Helvetica"
      ),
      axis.text = element_text(
        size = 6,
        family = "Helvetica",
        color = "black"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        color = "grey88",
        linewidth = 0.2
      ),
      panel.border = element_rect(
        color = "grey45",
        linewidth = 0.4
      ),
      axis.ticks = element_line(
        linewidth = 0.25
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(
        size = 5,
        family = "Helvetica"
      ),
      legend.key.height = unit(0.22, "cm"),
      legend.key.width = unit(0.45, "cm"),
      legend.spacing.x = unit(0.08, "cm"),
      legend.spacing.y = unit(0.02, "cm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      plot.margin = margin(6, 6, 6, 6)
    )
  
  return(p)
}

# ==============================
# 9. Generate all plots
# ==============================

plots <- lapply(
  metric_order,
  make_sensitivity_plot
)

names(plots) <- metric_order

# ==============================
# 10. Save PowerPoint
# ==============================

a4_width_cm  <- 21
a4_height_cm <- 29.7
a4_usable_width_cm <- 19.05

fig_width_single  <- a4_usable_width_cm / 2.54
fig_height_single <- 10 / 2.54

fig_width_9  <- a4_usable_width_cm / 2.54
fig_height_9 <- 24 / 2.54

fig_width_6  <- (a4_width_cm * 1 / 2) / 2.54
fig_height_6 <- (a4_height_cm * 1 / 2) / 2.54

ppt <- read_pptx()

# ==============================
# Pages 1-9: one plot per slide
# ==============================

for (m in metric_order) {
  
  p <- plots[[m]]
  
  ppt <- ppt %>%
    add_slide(layout = "Blank", master = "Office Theme") %>%
    ph_with(
      dml(ggobj = p),
      location = ph_location(
        left = 0,
        top = 0,
        width = fig_width_single,
        height = fig_height_single
      )
    )
}

# ==============================
# Page 10: all 9 plots in 3 x 3 layout
# One legend, 4 columns by MEM
# ==============================

combined_plot_9 <- wrap_plots(
  plots,
  ncol = 3,
  nrow = 3,
  guides = "collect"
) &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 5,
      family = "Helvetica"
    ),
    legend.key.height = unit(0.22, "cm"),
    legend.key.width = unit(0.45, "cm"),
    legend.spacing.x = unit(0.08, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ppt <- ppt %>%
  add_slide(layout = "Blank", master = "Office Theme") %>%
  ph_with(
    dml(ggobj = combined_plot_9),
    location = ph_location(
      left = 0,
      top = 0,
      width = fig_width_9,
      height = fig_height_9
    )
  )

# ==============================
# Page 11: selected 6 plots in 3 x 2 layout
# Row 1: Model similarity + Clustering performance (ARI)
# Row 2: Housekeeping gene (2022) + Marker gene
# Row 3: Metabolic task + Essential gene
# One legend, 4 columns by MEM
# ==============================

selected_metrics_6 <- c(
  "Model similarity",
  "Clustering performance (ARI)",
  "Housekeeping gene (2022)",
  "Marker gene",
  "Metabolic task",
  "Essential gene"
)

selected_metrics_6 <- selected_metrics_6[selected_metrics_6 %in% names(plots)]

plots_6 <- plots[selected_metrics_6]

combined_plot_6 <- wrap_plots(
  plots_6,
  ncol = 2,
  nrow = 3,
  guides = "collect"
) &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 5,
      family = "Helvetica"
    ),
    legend.key.height = unit(0.22, "cm"),
    legend.key.width = unit(0.45, "cm"),
    legend.spacing.x = unit(0.08, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ppt <- ppt %>%
  add_slide(layout = "Blank", master = "Office Theme") %>%
  ph_with(
    dml(ggobj = combined_plot_6),
    location = ph_location(
      left = 0,
      top = 0,
      width = fig_width_6,
      height = fig_height_6
    )
  )

# ==============================
# 11. Save file
# ==============================

print(
  ppt,
  target = file.path(
    out_dir,
    "Sensitivity_line_plots_by_metric_9panel_6panel_MEMlegend.pptx"
  )
)