# ============================================================
# Model size versus Metabolic tasks / Essential genes
# Exclude ftINIT before correlation calculation and plotting
#
# Page 1: color = MEM
# Page 2: color = Data
# Page 3: color = Threshold (valid thresholds only)
# Page 4: color = MEM, Reactions and Genes only
#
# Output: editable four-page PowerPoint
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(officer)
library(rvg)

if (!requireNamespace("ggh4x", quietly = TRUE)) {
  stop("Please install ggh4x: install.packages('ggh4x')")
}

# ==============================
# 1. Paths and settings
# ==============================

file <- paste0(
  "../scGEM_benchmark/result/",
  "accuracy_performance.xlsx"
)

out_dir <- paste0(
  "C:/Users/jjyan/Documents/GitHub/scGEM_benchmark/",
  "R_figures/"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

score_metric_cols <- list(
  "Metabolic tasks" = 32:40,
  "Essential genes" = 69:73
)

size_metric_cols <- list(
  "Metabolic tasks" = 14:22,
  "Essential genes" = 42:46
)

size_sheets <- c(
  Reactions = "rxns",
  Genes = "genes",
  Metabolites = "mets"
)

metric_order <- c("Metabolic tasks", "Essential genes")
size_order <- names(size_sheets)
mem_order <- c("GIMME", "iMAT", "rFASTCORMICS")
threshold_order <- c("25th", "50th", "75th", "90th")

group_colors <- c(
  GIMME        = "#9E9AC8",
  iMAT         = "#74C476",
  rFASTCORMICS = "#FDAE6B",
  R            = "#BDBDBD",
  LS           = "#F46D43",
  "25th"       = "#3182BD",
  "50th"       = "#6BAED6",
  "75th"       = "#9ECAE1",
  "90th"       = "#C6DBEF"
)

# ==============================
# 2. Clean metadata
# ==============================

clean_meta <- function(df) {
  df %>%
    mutate(
      across(
        c(Strategy, MEM, Data, Threshold, Dataset),
        ~ trimws(as.character(.x))
      ),
      
      Dataset = recode(Dataset, "Cell Lines" = "CellLines"),
      
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
      
      Threshold = gsub("[[:space:]%]+", "", Threshold),
      
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
}

# ==============================
# 3. Read scores or model sizes
# ==============================

read_metric_sheet <- function(sheet, metric_cols, value_name) {
  raw <- read_excel(file, sheet = sheet, col_names = FALSE)
  dat <- raw[-c(1, 2), ]
  
  meta <- dat[, 1:4]
  names(meta) <- c("Strategy", "MEM", "Data", "Threshold")
  
  bind_rows(lapply(names(metric_cols), function(m) {
    cols <- metric_cols[[m]]
    dataset_names <- make.unique(
      as.character(unlist(raw[2, cols]))
    )
    
    tmp <- bind_cols(meta, dat[, cols])
    names(tmp) <- c(names(meta), dataset_names)
    
    tmp %>%
      mutate(across(all_of(dataset_names), as.character)) %>%
      pivot_longer(
        cols = all_of(dataset_names),
        names_to = "Dataset",
        values_to = value_name
      ) %>%
      mutate(
        Metric = m,
        across(
          all_of(value_name),
          ~ suppressWarnings(as.numeric(.x))
        )
      )
  })) %>%
    clean_meta()
}

score_long <- read_metric_sheet(
  "Raw_scores", score_metric_cols, "Score"
) %>%
  filter(is.finite(Score), MEM %in% mem_order)

size_long <- bind_rows(lapply(size_order, function(s) {
  read_metric_sheet(
    size_sheets[[s]], size_metric_cols, "SizeValue"
  ) %>%
    mutate(SizeFeature = s)
})) %>%
  filter(is.finite(SizeValue), MEM %in% mem_order)

# ==============================
# 4. Merge and exclude ftINIT
# ==============================

plot_df <- inner_join(
  score_long,
  size_long,
  by = c(
    "Strategy", "MEM", "Data", "Threshold",
    "Dataset", "Metric"
  )
) %>%
  filter(
    MEM %in% mem_order,
    is.finite(Score),
    is.finite(SizeValue)
  ) %>%
  mutate(
    Metric = factor(Metric, levels = metric_order),
    SizeFeature = factor(SizeFeature, levels = size_order),
    MEM = factor(MEM, levels = mem_order),
    Data = factor(Data, levels = c("R", "LS")),
    Threshold = factor(Threshold, levels = threshold_order)
  )

if (nrow(plot_df) == 0) {
  stop("No matched data. Check worksheet headers and matching fields.")
}

print(plot_df %>% count(Metric, SizeFeature, MEM), n = Inf)

# ==============================
# 5. Calculate Spearman correlations
# ==============================

get_cor_results <- function(df) {
  df %>%
    group_by(Metric, SizeFeature) %>%
    summarise(
      Spearman_r = if (
        n() >= 2 &&
        n_distinct(SizeValue) > 1 &&
        n_distinct(Score) > 1
      ) {
        cor(SizeValue, Score, method = "spearman")
      } else {
        NA_real_
      },
      N = n(),
      .groups = "drop"
    ) %>%
    mutate(label = paste0("r = ", round(Spearman_r, 3)))
}

cat("\nOverall correlations excluding ftINIT:\n")
print(get_cor_results(plot_df))

cat("\nCorrelations using valid thresholds only:\n")
print(get_cor_results(filter(plot_df, !is.na(Threshold))))

# ==============================
# 6. Plot function
# ==============================

make_correlation_plot <- function(
    df, color_var, plot_title,
    size_features_to_keep = size_order
) {
  color_var <- rlang::ensym(color_var)
  color_name <- rlang::as_string(color_var)
  
  df_plot <- df %>%
    filter(SizeFeature %in% size_features_to_keep) %>%
    mutate(
      SizeFeature = factor(
        as.character(SizeFeature),
        levels = size_features_to_keep
      )
    )
  
  # Preserve the original threshold-page filtering rule
  if (color_name == "Threshold") {
    df_plot <- df_plot %>%
      filter(!is.na(Threshold))
  }
  
  cor_results <- get_cor_results(df_plot)
  
  ggplot(df_plot, aes(x = SizeValue, y = Score)) +
    geom_point(
      aes(color = !!color_var),
      size = 0.75, shape = 16, alpha = 0.78
    ) +
    geom_smooth(
      method = "lm", se = TRUE,
      linewidth = 0.30,
      color = "grey35", fill = "grey80"
    ) +
    geom_text(
      data = cor_results,
      aes(x = Inf, y = -Inf, label = label),
      inherit.aes = FALSE,
      hjust = 1.05, vjust = -0.8,
      size = 1.8, family = "Helvetica",
      fontface = "bold"
    ) +
    facet_grid(Metric ~ SizeFeature, scales = "free") +
    ggh4x::facetted_pos_scales(
      y = list(
        Metric == "Metabolic tasks" ~ scale_y_continuous(
          limits = c(0, 100),
          breaks = c(0, 25, 50, 75, 100),
          labels = function(x) sprintf("%.0f", x),
          expand = expansion(mult = c(0.04, 0.04))
        ),
        Metric == "Essential genes" ~ scale_y_continuous(
          labels = function(x) sprintf("%.1f", x),
          expand = expansion(mult = c(0.05, 0.05))
        )
      )
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
    theme_bw(base_size = 6, base_family = "Helvetica") +
    theme(
      text = element_text(family = "Helvetica", size = 6),
      plot.title = element_text(
        hjust = 0.5, face = "bold", size = 7
      ),
      strip.background = element_rect(
        fill = "#D9D5E5", color = "grey40", linewidth = 0.35
      ),
      strip.text = element_text(face = "bold", size = 5.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        color = "grey88", linewidth = 0.18
      ),
      panel.border = element_rect(
        color = "grey45", linewidth = 0.35
      ),
      axis.text = element_text(size = 4.8, color = "black"),
      axis.title = element_text(size = 6, face = "bold"),
      axis.ticks = element_line(
        linewidth = 0.25, color = "black"
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold", size = 6),
      legend.text = element_text(size = 5.5),
      legend.key.size = grid::unit(0.30, "cm"),
      panel.spacing.x = grid::unit(0.08, "cm"),
      panel.spacing.y = grid::unit(0.10, "cm"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

# ==============================
# 7. Generate four plots
# ==============================

plots <- list(
  make_correlation_plot(
    plot_df, MEM,
    "Model size and performance by MEM (excluding ftINIT)"
  ),
  make_correlation_plot(
    plot_df, Data,
    "Model size and performance by preprocessing (excluding ftINIT)"
  ),
  make_correlation_plot(
    plot_df, Threshold,
    "Model size and performance by threshold (excluding ftINIT)"
  ),
  make_correlation_plot(
    plot_df, MEM,
    "Reactions and genes versus performance (excluding ftINIT)",
    c("Reactions", "Genes")
  )
)

for (p in plots) print(p)

# ==============================
# 8. Save editable PowerPoint
# ==============================

fig_width <- 10.5 / 2.54
fig_height <- 7 / 2.54

ppt_file <- file.path(
  out_dir,
  "Figure_tasks_essential_genes_no_ftINIT_spearman_four_pages.pptx"
)

ppt <- read_pptx()

for (p in plots) {
  ppt <- ppt %>%
    add_slide(layout = "Blank", master = "Office Theme") %>%
    ph_with(
      dml(ggobj = p),
      location = ph_location(
        left = 0, top = 0,
        width = fig_width, height = fig_height
      )
    )
}

print(ppt, target = ppt_file)
cat("PowerPoint saved to:\n", ppt_file, "\n")