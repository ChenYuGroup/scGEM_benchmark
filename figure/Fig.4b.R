# 在原代码基础上：
# 1. 不再保存 PNG/CSV
# 2. 仅保存 5 张相关性热图到一个 PPT
# 3. 字体 Helvetica，字号 6

library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(officer)
library(rvg)
library(grid)


file <- "../scGEM_benchmark/result/accuracy_performance.xlsx"
out_dir <- "../scGEM_benchmark/result/dataset_correlation_heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read_excel(file, sheet = "Raw_scores", col_names = FALSE)
raw_data <- raw[-c(1, 2), ]

meta <- raw_data[, 1:4]
colnames(meta) <- c("Strategy", "MEM", "Data", "Threshold")
meta <- as.data.frame(meta) %>%
  mutate(across(everything(), as.character)) %>%
  filter(!is.na(Strategy), Strategy != "")

metric_cols <- list(
  "Housekeeping genes"     = 5:13,
  "Metabolic tasks"        = 14:22,
  "Marker genes"           = 23:32,
  "Clustering performance" = 33:41,
  "Essential genes"        = 42:46
)

extract_metric_matrix <- function(metric_name, cols) {
  tmp <- bind_cols(meta, raw_data[, cols])
  dataset_names <- as.character(unlist(raw[2, cols]))
  dataset_names <- make.unique(dataset_names)
  
  colnames(tmp) <- c("Strategy", "MEM", "Data", "Threshold", dataset_names)
  
  tmp <- tmp %>%
    mutate(across(all_of(dataset_names), ~ suppressWarnings(as.numeric(.)))) %>%
    filter(!is.na(Strategy), Strategy != "")
  
  score_matrix <- tmp %>%
    select(Strategy, all_of(dataset_names)) %>%
    column_to_rownames("Strategy") %>%
    as.matrix()
  
  storage.mode(score_matrix) <- "numeric"
  
  score_matrix
}

get_cluster_order <- function(cor_matrix) {
  cor_cluster <- cor_matrix
  cor_cluster[!is.finite(cor_cluster)] <- 0
  diag(cor_cluster) <- 1
  
  hc <- hclust(as.dist(1 - cor_cluster), method = "average")
  colnames(cor_matrix)[hc$order]
}

correlation_to_long <- function(cor_matrix, dataset_order) {
  cor_long <- as.data.frame(as.table(cor_matrix))
  colnames(cor_long) <- c("Dataset1", "Dataset2", "SpearmanR")
  
  cor_long %>%
    mutate(
      Dataset1 = factor(Dataset1, levels = dataset_order),
      Dataset2 = factor(Dataset2, levels = rev(dataset_order)),
      CorLabel = ifelse(is.finite(SpearmanR), sprintf("%.2f", SpearmanR), "NA")
    )
}

theme_heatmap <- theme_bw(base_size = 6, base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 6),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey35", linewidth = 0.35),
    legend.position = "right",
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6),
    plot.margin = margin(4, 4, 4, 4)
  )

plots <- list()

for (metric_name in names(metric_cols)) {
  
  score_matrix <- extract_metric_matrix(metric_name, metric_cols[[metric_name]])
  
  valid_dataset <- apply(score_matrix, 2, function(x) sum(is.finite(x)) >= 3)
  score_matrix <- score_matrix[, valid_dataset, drop = FALSE]
  
  cor_matrix <- cor(
    score_matrix,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  diag(cor_matrix) <- 1
  
  dataset_order <- get_cluster_order(cor_matrix)
  cor_long <- correlation_to_long(cor_matrix, dataset_order)
  
  p <- ggplot(
    cor_long,
    aes(
      x = Dataset1,
      y = Dataset2,
      fill = SpearmanR
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.35
    ) +
    scale_fill_gradientn(
      colours = c(
        "#4575B4",
        "#D9EAF4",
        "white",
        "#FDD0C0",
        "#F46D43",
        "#B2182B"
      ),
      values = scales::rescale(
        c(-1, -0.25, 0, 0.5, 0.8, 1)
      ),
      limits = c(-1, 1),
      breaks = c(-1, -0.5, 0, 0.5, 0.8, 1),
      name = "Spearman\ncorrelation",
      na.value = "grey85"
    ) +
    guides(
      fill = guide_colourbar(
        barwidth = unit(0.28, "cm"),
        barheight = unit(2.8, "cm"),
        title.position = "top",
        title.hjust = 0.5,
        ticks = TRUE,
        frame.colour = "grey35",
        frame.linewidth = 0.3
      )
    ) +
    coord_fixed() +
    labs(
      title = metric_name
    ) +
    theme_heatmap
  
  plots[[metric_name]] <- p
}

ppt <- read_pptx()

slide_width <- 8 / 2.54
slide_height <- 5 / 2.54

for (metric_name in names(plots)) {
  ppt <- ppt %>%
    add_slide(layout = "Blank", master = "Office Theme") %>%
    ph_with(
      dml(ggobj = plots[[metric_name]]),
      location = ph_location(
        left = 0,
        top = 0,
        width = slide_width,
        height = slide_height
      )
    )
}

out_file <- file.path(out_dir, "Dataset_correlation_heatmaps_5_metrics.pptx")
print(ppt, target = out_file)

cat("PPT saved to:\n", out_file, "\n")