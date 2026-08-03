# ============================================================
# Variance decomposition across five accuracy metrics
#
# Input:
# variance_decomposition.xlsx
# Sheet: Accuracy_variance
#
# Output:
# Editable PowerPoint
# ============================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(officer)
library(rvg)

# ==============================
# 1. File paths
# ==============================

source_file <- paste0(
  "../scGEM_benchmark/result/",
  "variance_decomposition.xlsx"
)

out_dir <- "../scGEM_benchmark/result"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ppt_file <- file.path(
  out_dir,
  "variance_explained_stacked_bar_half_A4_width.pptx"
)

if (!file.exists(source_file)) {
  stop(
    paste(
      "Input file was not found:",
      source_file
    )
  )
}

# ==============================
# 2. Read variance data
# ==============================

df <- read_excel(
  source_file,
  sheet = "Accuracy_variance"
)

required_columns <- c(
  "Metric",
  "MEM",
  "Data",
  "Threshold",
  "Residual"
)

missing_columns <- setdiff(
  required_columns,
  colnames(df)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}

df <- df %>%
  select(
    Metric,
    MEM,
    Data,
    Threshold,
    Residual
  ) %>%
  filter(
    !is.na(Metric)
  ) %>%
  mutate(
    Metric = as.character(Metric),
    across(
      c(MEM, Data, Threshold, Residual),
      as.numeric
    )
  )

# Check whether each row sums to 100%
print(
  df %>%
    mutate(
      CheckSum = MEM + Data + Threshold + Residual
    ) %>%
    select(
      Metric,
      CheckSum
    )
)

# ==============================
# 3. Convert to long format
# ==============================

df_long <- df %>%
  pivot_longer(
    cols = c(
      MEM,
      Data,
      Threshold,
      Residual
    ),
    names_to = "Factor",
    values_to = "Variance"
  )

# ==============================
# 4. Set plotting order
# ==============================

metric_order <- c(
  "Housekeeping genes",
  "Metabolic tasks",
  "Marker genes",
  "Clustering performance",
  "Essential genes"
)

df_long$Metric <- factor(
  df_long$Metric,
  levels = metric_order
)

df_long$Factor <- factor(
  df_long$Factor,
  levels = c(
    "MEM",
    "Data",
    "Threshold",
    "Residual"
  )
)

# ==============================
# 5. Set colors
# ==============================

factor_colors <- c(
  "MEM"       = "#8E7CC3",
  "Data"      = "#F08A65",
  "Threshold" = "#6BAED6",
  "Residual"  = "#C7C7C7"
)

# ==============================
# 6. Plot
# ==============================

p <- ggplot(
  df_long,
  aes(
    x = Metric,
    y = Variance,
    fill = Factor
  )
) +
  
  geom_col(
    width = 0.70,
    color = "white",
    linewidth = 0.25
  ) +
  
  scale_fill_manual(
    values = factor_colors
  ) +
  
  scale_y_continuous(
    limits = c(
      0,
      100
    ),
    breaks = seq(
      0,
      100,
      by = 20
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  labs(
    x = NULL,
    y = "Proportion of variance explained (%)",
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 6,
    base_family = "Helvetica"
  ) +
  
  theme(
    text = element_text(
      family = "Helvetica",
      size = 6
    ),
    
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1,
      size = 6,
      family = "Helvetica",
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 6,
      family = "Helvetica",
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 6,
      family = "Helvetica",
      face = "bold",
      color = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.30,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.25,
      color = "black"
    ),
    
    legend.position = "top",
    legend.direction = "horizontal",
    
    legend.text = element_text(
      size = 6,
      family = "Helvetica",
      color = "black"
    ),
    
    legend.key.size = unit(
      0.30,
      "cm"
    ),
    
    legend.spacing.x = unit(
      0.12,
      "cm"
    ),
    
    plot.margin = margin(
      3,
      3,
      3,
      3
    )
  )

print(p)

# ==============================
# 7. Save editable PowerPoint
# ==============================

figure_width_cm  <- 10.5
figure_height_cm <- 7.5

figure_width_in  <- figure_width_cm / 2.54
figure_height_in <- figure_height_cm / 2.54

ppt <- read_pptx() %>%
  add_slide(
    layout = "Blank",
    master = "Office Theme"
  ) %>%
  ph_with(
    dml(
      ggobj = p
    ),
    location = ph_location(
      left = 0,
      top = 0,
      width = figure_width_in,
      height = figure_height_in
    )
  )

print(
  ppt,
  target = ppt_file
)

cat(
  "PowerPoint saved to:\n",
  ppt_file,
  "\n"
)