# ============================================================
# Score panels with raw data points + mean bars + separated rank labels
#
# Ranking rule:
# If scores are tied, they share the same rank.

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


if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ==============================
# 2. Read Excel
# ==============================

raw <- read_excel(
  file,
  sheet = "Norm_scores",
  col_names = FALSE
)

# Row 1: metric names
# Row 2: dataset names
# Row 3+: scores
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

panel_order <- c(
  "Overall",
  "Housekeeping genes",
  "Metabolic tasks",
  "Marker genes",
  "Clustering performance",
  "Essential genes"
)

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
      values_to = "Score"
    ) %>%
    mutate(
      Metric = m,
      Score = suppressWarnings(as.numeric(Score))
    )
})

df_all <- bind_rows(df_list)

# ==============================
# 5. Filter empty / invalid values
# ==============================

df <- df_all %>%
  mutate(
    Score = suppressWarnings(as.numeric(Score))
  ) %>%
  filter(
    !is.na(Score),
    is.finite(Score)
  )

# ==============================
# 6. Calculate mean score for each strategy within each metric
# ==============================

metric_score <- df %>%
  group_by(
    Strategy,
    MEM,
    Data,
    Threshold,
    Metric
  ) %>%
  summarise(
    MeanScore = mean(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    !is.na(MeanScore),
    is.finite(MeanScore)
  )

# ==============================
# 7. Metric-level ranking
# Higher score = better rank
# Tied scores share the same rank
# ==============================

metric_rank <- metric_score %>%
  group_by(Metric) %>%
  arrange(
    desc(MeanScore),
    Strategy,
    .by_group = TRUE
  ) %>%
  mutate(
    Rank = min_rank(desc(MeanScore))
  ) %>%
  ungroup()

# ==============================
# 8. Overall score and overall ranking
# First average across datasets within each metric,
# then average across the five metrics
# Tied scores share the same rank
# ==============================

overall_score <- metric_score %>%
  group_by(
    Strategy,
    MEM,
    Data,
    Threshold
  ) %>%
  summarise(
    MeanScore = mean(MeanScore, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    !is.na(MeanScore),
    is.finite(MeanScore)
  ) %>%
  arrange(
    desc(MeanScore),
    Strategy
  ) %>%
  mutate(
    Metric = "Overall",
    Rank = min_rank(desc(MeanScore))
  )

# ==============================
# 9. Order strategies by overall ranking
# If overall rank is tied, use MeanScore and Strategy for stable display order
# ==============================

strategy_order <- overall_score %>%
  arrange(
    Rank,
    desc(MeanScore),
    Strategy
  ) %>%
  pull(Strategy)

# ==============================
# 10. Prepare plotting data
# ==============================

raw_points <- df %>%
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

metric_plot <- metric_rank %>%
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

overall_plot <- overall_score %>%
  mutate(
    Panel = "Overall",
    Strategy = factor(
      Strategy,
      levels = rev(strategy_order)
    ),
    Panel = factor(
      Panel,
      levels = panel_order
    )
  )

mean_plot <- bind_rows(
  overall_plot %>%
    select(
      Strategy,
      MEM,
      Data,
      Threshold,
      Panel,
      MeanScore,
      Rank
    ),
  metric_plot %>%
    select(
      Strategy,
      MEM,
      Data,
      Threshold,
      Panel,
      MeanScore,
      Rank
    )
) %>%
  filter(
    !is.na(MeanScore),
    is.finite(MeanScore)
  )

range_df <- df %>%
  group_by(
    Strategy,
    Metric
  ) %>%
  summarise(
    xmin = min(Score, na.rm = TRUE),
    xmax = max(Score, na.rm = TRUE),
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

# ==============================
# 11. Add separated rank column
# ==============================

rank_x <- 1.13

rank_plot <- mean_plot %>%
  mutate(
    RankX = rank_x
  )

# ==============================
# 12. Plot
# ==============================

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
      x = MeanScore,
      y = Strategy
    ),
    width = 0.62,
    fill = "#8E7CC3",
    color = "grey35",
    linewidth = 0.25
  ) +
  
  geom_point(
    data = raw_points,
    aes(
      x = Score,
      y = Strategy
    ),
    size = 0.75,
    shape = 16,
    color = "black",
    alpha = 0.75,
    position = position_jitter(
      height = 0.14,
      width = 0
    )
  ) +
  
  geom_text(
    data = rank_plot,
    aes(
      x = RankX,
      y = Strategy,
      label = Rank
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
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0)
  ) +
  
  coord_cartesian(
    xlim = c(0, 1.20),
    clip = "off"
  ) +
  
  labs(
    title = "Accuracy",
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
    axis.title.x = element_text(
      size = 6,
      family = "Helvetica"
    ),
    axis.text.y = element_text(
      size = 6,
      family = "Helvetica"
    ),
    axis.text.x = element_text(
      size = 6,
      angle = 45,
      hjust = 1,
      family = "Helvetica"
    ),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(
      linewidth = 0.25
    ),
    legend.position = "none",
    panel.spacing.x = unit(0.10, "cm"),
    plot.margin = margin(6, 25, 6, 6)
  )

print(p)

# ==============================
# 13. Save as PowerPoint only
# Figure width fitted to A4 portrait width: 19.05 cm
# ==============================

fig_width <- 19.05 / 2.54

# Keep original aspect ratio of previous figure: 15 x 8.5
fig_height <- 8.5 * fig_width / 15

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
      width = fig_width,
      height = fig_height
    )
  )

print(
  ppt,
  target = file.path(
    out_dir,
    "Figure_score_panels_points_rank_tied_rank.pptx"
  )
)