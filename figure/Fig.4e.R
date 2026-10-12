# ============================================================
# Variance decomposition of ARI-based clustering performance
#
# Each stacked bar is normalized to 100%
# Output: Editable PowerPoint
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

out_dir <- "C:/Users/jjyan/Documents/GitHub/scGEM_benchmark_NCS/R_figures"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

out_file <- file.path(
  out_dir,
  "ARI_variance_explained_stacked_bar_half_A4_width.pptx"
)

if (!file.exists(source_file)) {
  stop("Input file was not found: ", source_file)
}

# ==============================
# 2. Read and normalize data
# ==============================

df <- read_excel(
  source_file,
  sheet = "ARI_variance"
)

required_columns <- c(
  "Dataset",
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
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

factor_columns <- c(
  "MEM",
  "Data",
  "Threshold",
  "Residual"
)

df <- df %>%
  select(all_of(required_columns)) %>%
  filter(!is.na(Dataset)) %>%
  mutate(
    Dataset = trimws(as.character(Dataset)),
    across(all_of(factor_columns), as.numeric)
  ) %>%
  filter(Dataset != "")

# Check values before normalization
value_matrix <- as.matrix(df[, factor_columns])

invalid_rows <- apply(
  value_matrix,
  1,
  function(x) any(!is.finite(x) | x < 0)
)

if (any(invalid_rows)) {
  stop(
    "Missing, non-finite, or negative values in datasets: ",
    paste(df$Dataset[invalid_rows], collapse = ", ")
  )
}

row_totals <- rowSums(value_matrix)

if (any(row_totals <= 0)) {
  stop(
    "The total variance contribution must be positive: ",
    paste(df$Dataset[row_totals <= 0], collapse = ", ")
  )
}

# Normalize each dataset to 100%
df <- df %>%
  mutate(
    OriginalSum = row_totals,
    across(
      all_of(factor_columns),
      ~ .x / OriginalSum * 100
    )
  )

# Print totals before and after normalization
print(
  df %>%
    mutate(
      CheckSum = rowSums(
        across(all_of(factor_columns))
      )
    ) %>%
    select(
      Dataset,
      OriginalSum,
      CheckSum
    )
)

# ==============================
# 3. Convert to long format
# ==============================

df_long <- df %>%
  pivot_longer(
    cols = all_of(factor_columns),
    names_to = "Factor",
    values_to = "Variance"
  )

# ==============================
# 4. Set factor order
# ==============================

dataset_order <- c(
  "LUAD",
  "PAAD",
  "AML",
  "T2D",
  "BRCA",
  "MESO",
  "CellLines",
  "ESC",
  "SKCM"
)

# Retain any additional datasets
dataset_order <- c(
  dataset_order,
  setdiff(unique(df$Dataset), dataset_order)
)

df_long$Dataset <- factor(
  df_long$Dataset,
  levels = dataset_order
)

df_long$Factor <- factor(
  df_long$Factor,
  levels = factor_columns
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
    x = Dataset,
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
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0)
  ) +
  coord_cartesian(
    ylim = c(0, 100)
  ) +
  labs(
    title = "Clustering performance (ARI)",
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
    plot.title = element_text(
      size = 7,
      face = "bold",
      hjust = 0.5,
      family = "Helvetica",
      color = "black"
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
    legend.key.size = grid::unit(
      0.30,
      "cm"
    ),
    legend.spacing.x = grid::unit(
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

figure_width_cm <- 9
figure_height_cm <- 7

figure_width_in <- figure_width_cm / 2.54
figure_height_in <- figure_height_cm / 2.54

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