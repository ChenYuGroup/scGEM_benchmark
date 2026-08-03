# ============================================================
# Dataset-level grouped ranking plot
# Rows grouped by MEM and Threshold
# Data source: Sheet8
# Sheet8 structure:
# Strategy | MEM | Threshold | datasets...
#
# MEM score:
# For each MEM, use the maximum score across thresholds
# rFASTCORMICS has no threshold, so use its own score
#
# Threshold score:
# Keep unchanged; average across MEMs under each threshold
#
# Automatically add Overall if missing
# Ranking rule: min_rank(desc(score))
#
# Rank labels:
# Red rank labels are aligned in one fixed column on the right side
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

out_dir <- "../scGEM_benchmark/result/"


if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ==============================
# 2. Read Clustering
# Row 1 = metric title
# Row 2 = dataset names
# Row 3+ = scores
# ==============================

raw <- read_excel(
  file,
  sheet = "Clustering",
  col_names = FALSE
)

header1 <- trimws(as.character(unlist(raw[1, ])))
header2 <- trimws(as.character(unlist(raw[2, ])))

title_candidates <- header1[!is.na(header1) & header1 != ""]

if (length(title_candidates) > 0) {
  metric_title <- title_candidates[length(title_candidates)]
} else {
  metric_title <- "Clustering performance"
}

data_raw <- raw[-c(1, 2), ]

col_names <- header2
col_names[1:3] <- c("Strategy", "MEM", "Threshold")

empty_cols <- is.na(col_names) | col_names == ""
col_names[empty_cols] <- paste0("Empty_", which(empty_cols))

col_names <- make.unique(col_names)

colnames(data_raw) <- col_names
data_raw <- as.data.frame(data_raw)

data_raw <- data_raw[, colSums(!is.na(data_raw)) > 0]
colnames(data_raw) <- make.unique(colnames(data_raw))

# ==============================
# 3. Clean metadata
# ==============================

data_raw <- data_raw %>%
  mutate(
    Strategy = trimws(as.character(Strategy)),
    MEM = trimws(as.character(MEM)),
    Threshold = trimws(as.character(Threshold))
  ) %>%
  filter(
    !is.na(Strategy),
    Strategy != ""
  )

data_raw <- data_raw %>%
  mutate(
    MEM = case_when(
      grepl("ftINIT", MEM, ignore.case = TRUE) ~ "ftINIT",
      grepl("GIMME", MEM, ignore.case = TRUE) ~ "GIMME",
      grepl("iMAT", MEM, ignore.case = TRUE) ~ "iMAT",
      grepl("rFAST", MEM, ignore.case = TRUE) ~ "rFASTCORMICS",
      TRUE ~ MEM
    )
  )

data_raw <- data_raw %>%
  mutate(
    Threshold_raw = Threshold,
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

print(unique(data_raw$Threshold))

dataset_cols <- colnames(data_raw)[4:ncol(data_raw)]

dataset_cols <- dataset_cols[
  !is.na(dataset_cols) &
    dataset_cols != "" &
    dataset_cols != "NA" &
    !grepl("^Empty_", dataset_cols) &
    dataset_cols != "Threshold_raw"
]

print(dataset_cols)

# ==============================
# 4. Convert to long format
# ==============================

df_long <- data_raw %>%
  pivot_longer(
    cols = all_of(dataset_cols),
    names_to = "Dataset",
    values_to = "Score_raw"
  ) %>%
  mutate(
    Strategy = trimws(as.character(Strategy)),
    MEM = trimws(as.character(MEM)),
    Threshold = trimws(as.character(Threshold)),
    Dataset = trimws(as.character(Dataset)),
    Dataset = sub("\\.\\.\\.[0-9]+$", "", Dataset),
    Score_raw = trimws(as.character(Score_raw)),
    Score = suppressWarnings(as.numeric(Score_raw))
  ) %>%
  filter(
    !is.na(Strategy),
    Strategy != "",
    !is.na(Dataset),
    Dataset != "",
    !is.na(Score),
    is.finite(Score)
  )

# ==============================
# 5. Add Overall if missing
# ==============================

if (!"Overall" %in% unique(df_long$Dataset)) {
  
  overall_df <- df_long %>%
    group_by(Strategy, MEM, Threshold) %>%
    summarise(
      Score = mean(Score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Dataset = "Overall"
    )
  
  df_long <- bind_rows(
    overall_df,
    df_long %>%
      select(Strategy, MEM, Threshold, Dataset, Score)
  )
  
} else {
  
  df_long <- df_long %>%
    select(Strategy, MEM, Threshold, Dataset, Score)
}

# ==============================
# 6. Average strategies within each factor group
# ==============================

mem_threshold_score <- df_long %>%
  group_by(Dataset, MEM, Threshold) %>%
  summarise(
    ThresholdScore = mean(Score, na.rm = TRUE),
    .groups = "drop"
  )

mem_df <- mem_threshold_score %>%
  group_by(Dataset, MEM) %>%
  summarise(
    MeanScore = max(ThresholdScore, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Factor = "MEM",
    Level = MEM
  ) %>%
  select(Dataset, Factor, Level, MeanScore)

threshold_df <- df_long %>%
  filter(
    MEM != "rFASTCORMICS",
    Threshold %in% c("25th", "50th", "75th", "90th")
  ) %>%
  group_by(Dataset, Threshold) %>%
  summarise(
    MeanScore = mean(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Factor = "Threshold",
    Level = Threshold
  ) %>%
  select(Dataset, Factor, Level, MeanScore)

plot_df <- bind_rows(
  mem_df,
  threshold_df
)

# ==============================
# 7. Rank within each dataset and factor
# ==============================

plot_df <- plot_df %>%
  group_by(Factor, Dataset) %>%
  arrange(
    desc(MeanScore),
    .by_group = TRUE
  ) %>%
  mutate(
    Rank = min_rank(desc(MeanScore))
  ) %>%
  ungroup()

# ==============================
# 8. Order settings
# ==============================

dataset_order <- c(
  "Overall",
  "LUAD",
  "PAAD",
  "AML",
  "T2D",
  "BRCA",
  "MESO",
  "Cell Lines",
  "ESC",
  "SKCM"
)

dataset_order <- dataset_order[
  dataset_order %in% unique(plot_df$Dataset)
]

factor_order <- c("MEM", "Threshold")

overall_order <- plot_df %>%
  filter(Dataset == "Overall") %>%
  group_by(Factor) %>%
  arrange(
    desc(MeanScore),
    Level,
    .by_group = TRUE
  ) %>%
  mutate(
    OverallOrder = row_number()
  ) %>%
  ungroup() %>%
  select(
    Factor,
    Level,
    OverallOrder
  )

print(overall_order)

plot_df <- plot_df %>%
  left_join(
    overall_order,
    by = c("Factor", "Level")
  ) %>%
  mutate(
    Dataset = factor(Dataset, levels = dataset_order),
    Factor = factor(Factor, levels = factor_order),
    LevelAxis = paste0(
      Factor,
      "__",
      sprintf("%02d", OverallOrder),
      "__",
      Level
    ),
    
    # Fixed rank position:
    # all red rank labels are aligned in one column near the right side
    RankX = 1.055
  )

level_axis_order <- plot_df %>%
  distinct(
    Factor,
    OverallOrder,
    Level,
    LevelAxis
  ) %>%
  mutate(
    Factor = factor(Factor, levels = factor_order)
  ) %>%
  arrange(
    Factor,
    desc(OverallOrder)
  ) %>%
  pull(LevelAxis)

plot_df <- plot_df %>%
  mutate(
    LevelAxis = factor(
      LevelAxis,
      levels = level_axis_order
    )
  )

# ==============================
# 9. Plot
# ==============================

p <- ggplot(
  plot_df,
  aes(
    x = MeanScore,
    y = LevelAxis
  )
) +
  geom_col(
    width = 0.55,
    fill = "#8E7CC3",
    color = "grey35",
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      x = RankX,
      y = LevelAxis,
      label = Rank
    ),
    color = "#C44536",
    size = 2.0,
    family = "Helvetica",
    fontface = "bold",
    hjust = 0.5
  ) +
  facet_grid(
    Factor ~ Dataset,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_x_continuous(
    limits = c(0, 1.12),
    breaks = c(0, 0.5, 1),
    labels = c("", "", ""),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    labels = function(x) sub("^.*__[0-9]+__", "", x)
  ) +
  coord_cartesian(
    xlim = c(0, 1.12),
    clip = "off"
  ) +
  labs(
    title = metric_title,
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
      size = 7,
      family = "Helvetica"
    ),
    strip.background.x = element_rect(
      fill = "#D9D5E5",
      color = "grey40",
      linewidth = 0.4
    ),
    strip.text.x = element_text(
      face = "bold",
      size = 6,
      family = "Helvetica"
    ),
    strip.text.y.left = element_blank(),
    strip.background.y = element_blank(),
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
      family = "Helvetica",
      color = "black"
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
    legend.position = "none",
    panel.spacing.x = unit(0.12, "cm"),
    panel.spacing.y = unit(0.35, "cm"),
    plot.margin = margin(4, 18, 4, 4)
  )

print(p)

# ==============================
# 10. Save as PowerPoint only
# ==============================

fig_width <- 19 / 2.54
fig_height <- 7.0 / 2.54

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
    "Sheet8_clustering_performance_MEM_max_Threshold_mean_with_Overall_rank_aligned.pptx"
  )
)