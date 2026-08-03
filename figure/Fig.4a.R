# ============================================================
# Plot Kendall's W for five accuracy metrics
#
# Interpretation:
# Higher Kendall's W indicates greater consistency in strategy
# rankings across datasets and therefore weaker dataset dependence.
#
# Input:
# Kendall_W_ranking_agreement.xlsx
#
# Output:
# Editable PowerPoint only
# Font: Helvetica
# Font size: 6 pt
# ============================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(officer)
library(rvg)

# ==============================
# 1. File paths
# ==============================

source_file <- paste0(
  "../scGEM_benchmark/result/",
  "Kendall_W_ranking_agreement.xlsx"
)

out_dir <- "../scGEM_benchmark/result"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

out_file <- file.path(
  out_dir,
  "Kendall_W_dataset_concordance.pptx"
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
# 2. Read Kendall's W results
# ==============================

df <- read_excel(
  source_file,
  sheet = 1
)

# Check required columns
required_columns <- c(
  "Metric",
  "Kendall_W",
  "N_datasets"
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

# Retain columns required for plotting
df <- df %>%
  select(
    Metric,
    Kendall_W,
    N_datasets
  ) %>%
  filter(
    !is.na(Metric),
    !is.na(Kendall_W),
    !is.na(N_datasets)
  ) %>%
  mutate(
    Metric = as.character(Metric),
    Kendall_W = as.numeric(Kendall_W),
    N_datasets = as.integer(N_datasets)
  )

# ==============================
# 3. Order metrics
# ==============================

df <- df %>%
  arrange(Kendall_W) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = Metric
    ),
    W_label = sprintf(
      "%.2f",
      Kendall_W
    ),
    Dataset_label = paste0(
      "n = ",
      N_datasets
    )
  )

# ==============================
# 4. Plot
# ==============================

p <- ggplot(
  df,
  aes(
    x = Kendall_W,
    y = Metric
  )
) +
  
  # Reference grid lines
  geom_vline(
    xintercept = c(
      0.25,
      0.50,
      0.75
    ),
    color = "grey85",
    linewidth = 0.3
  ) +
  
  # Bars
  geom_col(
    width = 0.58,
    fill = "#8E7CC3",
    color = "grey35",
    linewidth = 0.35
  ) +
  
  # Kendall's W labels
  geom_text(
    aes(
      x = Kendall_W + 0.025,
      label = W_label
    ),
    hjust = 0,
    family = "Helvetica",
    size = 2.1,
    fontface = "bold"
  ) +
  
  # Number of datasets
  geom_text(
    aes(
      x = 1.08,
      label = Dataset_label
    ),
    hjust = 1,
    family = "Helvetica",
    size = 2.0,
    color = "grey35"
  ) +
  
  scale_x_continuous(
    limits = c(
      0,
      1.10
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    ),
    labels = sprintf(
      "%.2f",
      seq(
        0,
        1,
        by = 0.25
      )
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  labs(
    x = "Kendall's W",
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
    
    axis.text.x = element_text(
      size = 6,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 6,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 6,
      face = "bold",
      margin = margin(
        t = 4
      )
    ),
    
    axis.ticks.y = element_blank(),
    
    axis.ticks.x = element_line(
      color = "black",
      linewidth = 0.25
    ),
    
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      linewidth = 0.4
    ),
    
    plot.margin = margin(
      t = 5,
      r = 6,
      b = 5,
      l = 5
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  )

print(p)

# ==============================
# 5. Save editable PowerPoint
# ==============================

figure_width_cm  <- 8
figure_height_cm <- 5

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
  target = out_file
)

cat(
  "PowerPoint saved to:\n",
  out_file,
  "\n"
)