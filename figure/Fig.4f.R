# ============================================================
# Compare scGEM- and expression-based clustering performance
#
# Input:
# clustering_performance_scRNA.xlsx
# Sheets: ARI and NMI
#
# Output:
# Editable PowerPoint with ARI and NMI plots
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
  "clustering_performance_scRNA.xlsx"
)

out_dir <- "../scGEM_benchmark/result"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ppt_file <- file.path(
  out_dir,
  "ARI_NMI_barplot_two_pages_legend_under_title.pptx"
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
# 2. Read ARI and NMI data
# ==============================

ari_data <- read_excel(
  source_file,
  sheet = "ARI"
)

nmi_data <- read_excel(
  source_file,
  sheet = "NMI"
)

# Check required columns
required_ari_columns <- c(
  "Dataset",
  "Group",
  "ARI"
)

required_nmi_columns <- c(
  "Dataset",
  "Group",
  "NMI"
)

missing_ari_columns <- setdiff(
  required_ari_columns,
  colnames(ari_data)
)

missing_nmi_columns <- setdiff(
  required_nmi_columns,
  colnames(nmi_data)
)

if (length(missing_ari_columns) > 0) {
  stop(
    paste(
      "Missing columns in the ARI sheet:",
      paste(
        missing_ari_columns,
        collapse = ", "
      )
    )
  )
}

if (length(missing_nmi_columns) > 0) {
  stop(
    paste(
      "Missing columns in the NMI sheet:",
      paste(
        missing_nmi_columns,
        collapse = ", "
      )
    )
  )
}

# Format and combine the two sheets
ari_data <- ari_data %>%
  select(
    Dataset,
    Group,
    ARI
  ) %>%
  filter(
    !is.na(Dataset),
    !is.na(Group),
    !is.na(ARI)
  ) %>%
  mutate(
    Dataset = as.character(Dataset),
    Group = as.character(Group),
    ARI = as.numeric(ARI)
  )

nmi_data <- nmi_data %>%
  select(
    Dataset,
    Group,
    NMI
  ) %>%
  filter(
    !is.na(Dataset),
    !is.na(Group),
    !is.na(NMI)
  ) %>%
  mutate(
    Dataset = as.character(Dataset),
    Group = as.character(Group),
    NMI = as.numeric(NMI)
  )

df <- full_join(
  ari_data,
  nmi_data,
  by = c(
    "Dataset",
    "Group"
  )
)

# Check for unmatched or missing observations
if (anyNA(df$ARI) || anyNA(df$NMI)) {
  warning(
    paste(
      "Some Dataset–Group combinations were not matched",
      "between the ARI and NMI sheets."
    )
  )
}

# ==============================
# 3. Set plotting order and colors
# ==============================

dataset_order <- c(
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

group_order <- c(
  "scGEM Best Performance",
  "scRNA MeGenes",
  "scRNA AllGenes"
)

group_colors <- c(
  "scGEM Best Performance" = "#A1D99B",
  "scRNA MeGenes" = "#9ECAE1",
  "scRNA AllGenes" = "#FDD0A2"
)

# ==============================
# 4. Plot function
# ==============================

make_plot <- function(metric_name) {
  
  plot_df <- df %>%
    select(
      Dataset,
      Group,
      all_of(metric_name)
    ) %>%
    rename(
      Score = all_of(metric_name)
    ) %>%
    mutate(
      Dataset = factor(
        Dataset,
        levels = dataset_order
      ),
      Group = factor(
        Group,
        levels = group_order
      )
    )
  
  p <- ggplot(
    plot_df,
    aes(
      x = Dataset,
      y = Score,
      fill = Group
    )
  ) +
    
    geom_col(
      position = position_dodge(
        width = 0.75
      ),
      width = 0.65,
      color = "grey35",
      linewidth = 0.25
    ) +
    
    scale_fill_manual(
      values = group_colors
    ) +
    
    scale_y_continuous(
      limits = c(
        0,
        1.05
      ),
      breaks = seq(
        0,
        1,
        by = 0.2
      ),
      expand = c(
        0,
        0
      )
    ) +
    
    labs(
      title = metric_name,
      x = NULL,
      y = "Score",
      fill = NULL
    ) +
    
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE
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
        family = "Helvetica",
        margin = margin(
          b = 2
        )
      ),
      
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      
      panel.grid.major.y = element_line(
        color = "grey88",
        linewidth = 0.25
      ),
      
      panel.border = element_rect(
        color = "grey45",
        linewidth = 0.35
      ),
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
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
        face = "bold",
        size = 6,
        family = "Helvetica"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.25,
        color = "black"
      ),
      
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "center",
      
      legend.text = element_text(
        size = 6,
        family = "Helvetica"
      ),
      
      legend.key.size = unit(
        0.30,
        "cm"
      ),
      
      legend.margin = margin(
        t = 0,
        r = 0,
        b = 2,
        l = 0
      ),
      
      legend.box.margin = margin(
        t = 0,
        r = 0,
        b = 2,
        l = 0
      ),
      
      plot.margin = margin(
        3,
        3,
        3,
        3
      )
    )
  
  return(p)
}

# ==============================
# 5. Generate plots
# ==============================

p_ari <- make_plot("ARI")
p_nmi <- make_plot("NMI")

print(p_ari)
print(p_nmi)

# ==============================
# 6. Save editable PowerPoint
# ==============================

figure_width_cm  <- 9
figure_height_cm <- 7

figure_width_in  <- figure_width_cm / 2.54
figure_height_in <- figure_height_cm / 2.54

ppt <- read_pptx()

ppt <- ppt %>%
  add_slide(
    layout = "Blank",
    master = "Office Theme"
  ) %>%
  ph_with(
    dml(
      ggobj = p_ari
    ),
    location = ph_location(
      left = 0,
      top = 0,
      width = figure_width_in,
      height = figure_height_in
    )
  )

ppt <- ppt %>%
  add_slide(
    layout = "Blank",
    master = "Office Theme"
  ) %>%
  ph_with(
    dml(
      ggobj = p_nmi
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