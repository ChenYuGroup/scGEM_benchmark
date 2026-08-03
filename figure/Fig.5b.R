# ============================================================
# Correlation between five metric scores and model size
#
# Raw_scores    : raw scores of five metrics
# rxns : reaction numbers
# genes : gene numbers
# mets : metabolite numbers
#
# Each point = one strategy in one dataset
#
# Page 1 color = MEM
# Page 2 color = Data
# Page 3 color = Threshold
# Page 4 color = MEM, only Reactions and Genes
#
# Save PowerPoint only
# Figure width  = half of A4 paper width
# Figure height = half of A4 paper height
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

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================
# 2. Sheet names
# ==============================

score_sheet <- "Raw_scores"

size_sheets <- c(
  "Reactions"   = "rxns",
  "Genes"       = "genes",
  "Metabolites" = "mets"
)

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

metric_order <- c(
  "Housekeeping genes",
  "Metabolic tasks",
  "Marker genes",
  "Clustering performance",
  "Essential genes"
)

size_order <- c(
  "Reactions",
  "Genes",
  "Metabolites"
)

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
# 4. Function: clean metadata
# ==============================

clean_meta <- function(df) {
  
  df %>%
    mutate(
      Strategy = as.character(Strategy),
      MEM = trimws(as.character(MEM)),
      Data = trimws(as.character(Data)),
      Threshold = trimws(as.character(Threshold)),
      
      MEM = case_when(
        grepl("ftINIT", MEM, ignore.case = TRUE) ~ "ftINIT",
        grepl("GIMME", MEM, ignore.case = TRUE) ~ "GIMME",
        grepl("iMAT", MEM, ignore.case = TRUE) ~ "iMAT",
        grepl("rFAST", MEM, ignore.case = TRUE) ~ "rFASTCORMICS",
        TRUE ~ MEM
      ),
      
      Data = case_when(
        Data %in% c("R", "Raw", "raw") ~ "R",
        Data %in% c("LS", "Linnorm+SAVER", "linnorm+saver") ~ "LS",
        TRUE ~ Data
      ),
      
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
}

# ==============================
# 5. Function: read score sheet
# ==============================

read_score_sheet <- function() {
  
  raw <- read_excel(
    file,
    sheet = score_sheet,
    col_names = FALSE
  )
  
  data_raw <- raw[-c(1, 2), ]
  
  meta <- data_raw[, 1:4]
  colnames(meta) <- c("Strategy", "MEM", "Data", "Threshold")
  meta <- as.data.frame(meta)
  
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
        Dataset = as.character(Dataset),
        Metric = m,
        Score = suppressWarnings(as.numeric(Score))
      )
  })
  
  bind_rows(df_list) %>%
    clean_meta() %>%
    filter(
      !is.na(Score),
      is.finite(Score)
    )
}

# ==============================
# 6. Function: read one model size sheet
# ==============================

read_size_sheet <- function(sheet_name, size_feature) {
  
  raw <- read_excel(
    file,
    sheet = sheet_name,
    col_names = FALSE
  )
  
  data_raw <- raw[-c(1, 2), ]
  
  meta <- data_raw[, 1:4]
  colnames(meta) <- c("Strategy", "MEM", "Data", "Threshold")
  meta <- as.data.frame(meta)
  
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
        values_to = "SizeValue"
      ) %>%
      mutate(
        Dataset = as.character(Dataset),
        Metric = m,
        SizeFeature = size_feature,
        SizeValue = suppressWarnings(as.numeric(SizeValue))
      )
  })
  
  bind_rows(df_list) %>%
    clean_meta() %>%
    filter(
      !is.na(SizeValue),
      is.finite(SizeValue)
    )
}

# ==============================
# 7. Read and merge data
# ==============================

score_long <- read_score_sheet()

size_long <- bind_rows(
  read_size_sheet(size_sheets["Reactions"],   "Reactions"),
  read_size_sheet(size_sheets["Genes"],       "Genes"),
  read_size_sheet(size_sheets["Metabolites"], "Metabolites")
)

plot_df <- score_long %>%
  inner_join(
    size_long,
    by = c(
      "Strategy",
      "MEM",
      "Data",
      "Threshold",
      "Dataset",
      "Metric"
    )
  ) %>%
  filter(
    !is.na(Score),
    !is.na(SizeValue),
    is.finite(Score),
    is.finite(SizeValue)
  ) %>%
  mutate(
    Metric = factor(Metric, levels = metric_order),
    SizeFeature = factor(SizeFeature, levels = size_order),
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
# 8. Function: calculate Spearman correlations
# ==============================

get_cor_results <- function(df) {
  
  df %>%
    group_by(Metric, SizeFeature) %>%
    summarise(
      Spearman_r = cor(
        SizeValue,
        Score,
        method = "spearman",
        use = "complete.obs"
      ),
      N = n(),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0("r = ", round(Spearman_r, 3))
    )
}

# ==============================
# 9. Function: make correlation plot
# ==============================

make_correlation_plot <- function(
    df,
    color_var,
    plot_title,
    size_features_to_keep = c("Reactions", "Genes", "Metabolites")
) {
  
  color_var <- rlang::ensym(color_var)
  color_name <- rlang::as_string(color_var)
  
  df_plot <- df %>%
    filter(
      SizeFeature %in% size_features_to_keep
    ) %>%
    mutate(
      SizeFeature = factor(
        as.character(SizeFeature),
        levels = size_features_to_keep
      )
    )
  
  if (color_name == "Threshold") {
    df_plot <- df_plot %>%
      filter(
        !is.na(Threshold),
        Threshold %in% c("25th", "50th", "75th", "90th")
      )
  }
  
  cor_results <- get_cor_results(df_plot)
  
  p <- ggplot(
    df_plot,
    aes(
      x = SizeValue,
      y = Score
    )
  ) +
    geom_point(
      aes(color = !!color_var),
      size = 0.75,
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
    geom_text(
      data = cor_results,
      aes(
        x = Inf,
        y = -Inf,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = -0.8,
      size = 1.8,
      family = "Helvetica",
      fontface = "bold"
    ) +
    facet_grid(
      Metric ~ SizeFeature,
      scales = "free"
    ) +
    scale_color_manual(
      values = group_colors,
      drop = FALSE
    ) +
    labs(
      title = plot_title,
      x = "Model size",
      y = "Raw score",
      color = color_name
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
      axis.title = element_text(
        size = 6,
        family = "Helvetica",
        face = "bold"
      ),
      axis.ticks = element_line(
        linewidth = 0.25,
        color = "black"
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(
        face = "bold",
        size = 6,
        family = "Helvetica"
      ),
      legend.text = element_text(
        size = 5.5,
        family = "Helvetica"
      ),
      legend.key.size = unit(0.30, "cm"),
      panel.spacing.x = unit(0.08, "cm"),
      panel.spacing.y = unit(0.10, "cm"),
      plot.margin = margin(3, 3, 3, 3)
    )
  
  return(p)
}

# ==============================
# 10. Generate four plots
# ==============================

p_mem <- make_correlation_plot(
  df = plot_df,
  color_var = MEM,
  plot_title = "Correlations between model size and metric performance by MEM",
  size_features_to_keep = c("Reactions", "Genes", "Metabolites")
)

p_data <- make_correlation_plot(
  df = plot_df,
  color_var = Data,
  plot_title = "Correlations between model size and metric performance by data preprocessing",
  size_features_to_keep = c("Reactions", "Genes", "Metabolites")
)

p_threshold <- make_correlation_plot(
  df = plot_df,
  color_var = Threshold,
  plot_title = "Correlations between model size and metric performance by threshold",
  size_features_to_keep = c("Reactions", "Genes", "Metabolites")
)

p_mem_reactions_genes <- make_correlation_plot(
  df = plot_df,
  color_var = MEM,
  plot_title = "Correlations with reactions and genes by MEM",
  size_features_to_keep = c("Reactions", "Genes")
)

print(p_mem)
print(p_data)
print(p_threshold)
print(p_mem_reactions_genes)

# ==============================
# 11. Save four-page PowerPoint only
# Figure width  = half of A4 paper width
# Figure height = half of A4 paper height
# ==============================

# A4 paper size:
# width  = 21 cm
# height = 29.7 cm
# Half width  = 10.5 cm
# Half height = 14.85 cm

fig_width <- 10.5 / 2.54
fig_height <- 14.85 / 2.54

ppt_file <- file.path(
  out_dir,
  "Figure_model_size_vs_five_metrics_spearman_four_pages_half_A4_size.pptx"
)

ppt <- read_pptx()

ppt <- ppt %>%
  add_slide(
    layout = "Blank",
    master = "Office Theme"
  ) %>%
  ph_with(
    dml(ggobj = p_mem),
    location = ph_location(
      left = 0,
      top = 0,
      width = fig_width,
      height = fig_height
    )
  )

ppt <- ppt %>%
  add_slide(
    layout = "Blank",
    master = "Office Theme"
  ) %>%
  ph_with(
    dml(ggobj = p_data),
    location = ph_location(
      left = 0,
      top = 0,
      width = fig_width,
      height = fig_height
    )
  )

ppt <- ppt %>%
  add_slide(
    layout = "Blank",
    master = "Office Theme"
  ) %>%
  ph_with(
    dml(ggobj = p_threshold),
    location = ph_location(
      left = 0,
      top = 0,
      width = fig_width,
      height = fig_height
    )
  )

ppt <- ppt %>%
  add_slide(
    layout = "Blank",
    master = "Office Theme"
  ) %>%
  ph_with(
    dml(ggobj = p_mem_reactions_genes),
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