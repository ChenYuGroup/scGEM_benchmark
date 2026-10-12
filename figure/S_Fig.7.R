# ============================================================
# Dataset-level ranking plot: MEM + Data + Threshold
#
# Sheet: Clustering_all
# Row 1: metric title
# Row 2: column names
# Row 3+: scores
#
# First four columns:
# Strategy | MEM | Data | Threshold
#
# MEM: maximum across threshold-level mean scores
# Data: mean across all strategies within LS / R
# Threshold: mean across strategies, excluding rFASTCORMICS
#
# Output: editable PowerPoint
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(officer)
library(rvg)

# ==============================
# 1. File paths and column positions
# ==============================

file <- paste0(
  "../scGEM_benchmark/result/",
  "accuracy_performance.xlsx"
)

out_dir <- paste0(
  "../scGEM_benchmark/",
  "R_figures/"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Adjust these positions if the Excel column order differs
metadata_cols <- c(
  Strategy  = 1L,
  MEM       = 2L,
  Data      = 3L,
  Threshold = 4L
)

if (!file.exists(file)) {
  stop("Input file was not found: ", file)
}

# ==============================
# 2. Read worksheet
# ==============================

raw <- read_excel(
  file,
  sheet = "Clustering_all",
  col_names = FALSE
)

if (nrow(raw) < 3 || ncol(raw) < 5) {
  stop("The worksheet must contain two header rows and at least five columns.")
}

header1 <- trimws(as.character(unlist(raw[1, ])))
header2 <- trimws(as.character(unlist(raw[2, ])))

title_candidates <- header1[
  !is.na(header1) & header1 != ""
]

metric_title <- if (length(title_candidates) > 0) {
  tail(title_candidates, 1)
} else {
  "Clustering performance"
}

col_names <- header2
col_names[unname(metadata_cols)] <- names(metadata_cols)

empty_cols <- is.na(col_names) | col_names == ""
col_names[empty_cols] <- paste0("Empty_", which(empty_cols))
col_names <- make.unique(col_names)

data_raw <- as.data.frame(raw[-c(1, 2), ])
colnames(data_raw) <- col_names

dataset_cols <- col_names[
  setdiff(seq_along(col_names), unname(metadata_cols))
]

dataset_cols <- dataset_cols[
  !grepl("^Empty_", dataset_cols) &
    dataset_cols != "NA"
]

if (length(dataset_cols) == 0) {
  stop("No dataset columns were found.")
}

# ==============================
# 3. Clean metadata
# ==============================

data_raw <- data_raw %>%
  mutate(
    across(
      all_of(names(metadata_cols)),
      ~ trimws(as.character(.x))
    )
  ) %>%
  filter(!is.na(Strategy), Strategy != "") %>%
  mutate(
    MEM = case_when(
      grepl("ftINIT", MEM, ignore.case = TRUE) ~ "ftINIT",
      grepl("GIMME", MEM, ignore.case = TRUE) ~ "GIMME",
      grepl("iMAT", MEM, ignore.case = TRUE) ~ "iMAT",
      grepl("rFAST", MEM, ignore.case = TRUE) ~ "rFASTCORMICS",
      TRUE ~ MEM
    ),
    
    Data = toupper(gsub("\\s+", "", Data)),
    
    Threshold = tolower(
      gsub("[[:space:]%]+", "", Threshold)
    ),
    
    Threshold = case_when(
      Threshold %in% c(
        "25th", "25", "0.25", "0.250", "0.2500"
      ) ~ "25th",
      
      Threshold %in% c(
        "50th", "50", "0.5", "0.50", "0.500", "0.5000"
      ) ~ "50th",
      
      Threshold %in% c(
        "75th", "75", "0.75", "0.750", "0.7500"
      ) ~ "75th",
      
      Threshold %in% c(
        "90th", "90", "0.9", "0.90", "0.900", "0.9000"
      ) ~ "90th",
      
      TRUE ~ Threshold
    )
  )

# Print metadata to check column positions
print(head(data_raw[c("Strategy", "MEM", "Data", "Threshold")], 12))

# ==============================
# 4. Convert to long format
# ==============================

df_long <- data_raw %>%
  mutate(
    across(all_of(dataset_cols), as.character)
  ) %>%
  pivot_longer(
    cols = all_of(dataset_cols),
    names_to = "Dataset",
    values_to = "Score_raw"
  ) %>%
  mutate(
    Dataset = trimws(Dataset),
    Dataset = sub("\\.\\.\\.[0-9]+$", "", Dataset),
    Score_raw = trimws(Score_raw),
    
    # Convert values such as (0.12) to -0.12
    Score_raw = sub(
      "^\\((.*)\\)$",
      "-\\1",
      Score_raw
    ),
    
    Score = suppressWarnings(as.numeric(Score_raw))
  ) %>%
  filter(
    !is.na(MEM),
    MEM != "",
    !is.na(Dataset),
    Dataset != "",
    is.finite(Score)
  ) %>%
  select(
    Strategy, MEM, Data, Threshold, Dataset, Score
  )

if (nrow(df_long) == 0) {
  stop("No valid numeric scores were found.")
}

# ==============================
# 5. Add Overall if missing
# ==============================

if (!"Overall" %in% df_long$Dataset) {
  
  overall_df <- df_long %>%
    group_by(Strategy, MEM, Data, Threshold) %>%
    summarise(
      Score = mean(Score),
      .groups = "drop"
    ) %>%
    mutate(Dataset = "Overall")
  
  df_long <- bind_rows(overall_df, df_long)
}

# ==============================
# 6. Calculate MEM comparison
# ==============================

# Average LS / R within each MEM and threshold
mem_threshold_score <- df_long %>%
  group_by(Dataset, MEM, Threshold) %>%
  summarise(
    ThresholdScore = mean(Score),
    .groups = "drop"
  )

# Select the maximum threshold-level mean for each MEM
mem_df <- mem_threshold_score %>%
  group_by(Dataset, MEM) %>%
  summarise(
    MeanScore = max(ThresholdScore),
    .groups = "drop"
  ) %>%
  transmute(
    Dataset,
    Factor = "MEM",
    Level = MEM,
    MeanScore
  )

# ==============================
# 7. Calculate Data comparison
# ==============================

data_df <- df_long %>%
  filter(Data %in% c("LS", "R")) %>%
  group_by(Dataset, Data) %>%
  summarise(
    MeanScore = mean(Score),
    .groups = "drop"
  ) %>%
  transmute(
    Dataset,
    Factor = "Data",
    Level = Data,
    MeanScore
  )

# ==============================
# 8. Calculate Threshold comparison
# ==============================

threshold_df <- df_long %>%
  filter(
    MEM != "rFASTCORMICS",
    Threshold %in% c("25th", "50th", "75th", "90th")
  ) %>%
  group_by(Dataset, Threshold) %>%
  summarise(
    MeanScore = mean(Score),
    .groups = "drop"
  ) %>%
  transmute(
    Dataset,
    Factor = "Threshold",
    Level = Threshold,
    MeanScore
  )

# Stop instead of silently omitting comparison groups
if (nrow(data_df) == 0) {
  stop(
    "Data comparison is empty. Check the Data column position ",
    "and whether its values are LS / R."
  )
}

if (nrow(threshold_df) == 0) {
  stop(
    "Threshold comparison is empty. Check the Threshold column ",
    "position and its values."
  )
}

plot_df <- bind_rows(mem_df, data_df, threshold_df)

print(plot_df %>% count(Factor))

# ==============================
# 9. Rank within each dataset and factor
# ==============================

plot_df <- plot_df %>%
  group_by(Factor, Dataset) %>%
  mutate(
    Rank = min_rank(desc(MeanScore))
  ) %>%
  ungroup()

# ==============================
# 10. Order settings
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

# Preserve additional dataset names
dataset_order <- c(
  dataset_order[dataset_order %in% plot_df$Dataset],
  setdiff(unique(plot_df$Dataset), dataset_order)
)

factor_order <- c("MEM", "Threshold", "Data")

overall_order <- plot_df %>%
  filter(Dataset == "Overall") %>%
  group_by(Factor) %>%
  arrange(desc(MeanScore), Level, .by_group = TRUE) %>%
  mutate(OverallOrder = row_number()) %>%
  ungroup() %>%
  select(Factor, Level, OverallOrder)

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
    
    RankX = 1.055
  )

level_axis_order <- plot_df %>%
  distinct(Factor, OverallOrder, Level, LevelAxis) %>%
  arrange(Factor, desc(OverallOrder)) %>%
  pull(LevelAxis)

plot_df <- plot_df %>%
  mutate(
    LevelAxis = factor(
      LevelAxis,
      levels = level_axis_order
    )
  )

# ==============================
# 11. Plot
# ==============================

p <- ggplot(
  plot_df,
  aes(x = MeanScore, y = LevelAxis)
) +
  geom_col(
    width = 0.55,
    fill = "#8E7CC3",
    color = "grey35",
    linewidth = 0.25
  ) +
  geom_text(
    aes(x = RankX, label = Rank),
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
    breaks = c(0, 0.5, 1),
    labels = c("", "", ""),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    labels = function(x) {
      sub("^.*__[0-9]+__", "", x)
    }
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
      size = 7
    ),
    
    strip.background.x = element_rect(
      fill = "#D9D5E5",
      color = "grey40",
      linewidth = 0.4
    ),
    
    strip.text.x = element_text(
      face = "bold",
      size = 6
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
      face = "bold"
    ),
    
    legend.position = "none",
    
    panel.spacing.x = grid::unit(0.12, "cm"),
    panel.spacing.y = grid::unit(0.35, "cm"),
    
    plot.margin = margin(4, 18, 4, 4)
  )

print(p)

# ==============================
# 12. Save editable PowerPoint
# ==============================

fig_width <- 19 / 2.54
fig_height <- 9 / 2.54

ppt_file <- file.path(
  out_dir,
  paste0(
    "Clustering_performance_MEM_max_Data_mean_",
    "Threshold_mean_with_Overall_rank_aligned.pptx"
  )
)

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

print(ppt, target = ppt_file)

cat("PowerPoint saved to:\n", ppt_file, "\n")