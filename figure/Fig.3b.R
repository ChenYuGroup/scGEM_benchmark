# ==============================
#
# Input:
# feasibility_performance.xlsx
# Sheet: feasibility_performance
#
# Layout:
# Strategy | Status columns | Mean construction time
#
# Save as editable PowerPoint.
# ==============================

library(readxl)
library(ggplot2)
library(dplyr)
library(officer)
library(rvg)
library(grid)

# ------------------------------
# 1. File paths
# ------------------------------
file <- "../scGEM_benchmark/result/sensitivity_performance.xlsx"
out_dir <- "../scGEM_benchmark/result"
source_file <- file


if (!file.exists(source_file)) {
  stop(
    paste(
      "Source-data file was not found:",
      source_file
    )
  )
}

if (!dir.exists(out_dir)) {
  dir.create(
    out_dir,
    recursive = TRUE
  )
}

# ------------------------------
# 2. Read feasibility data
# ------------------------------
df_plot <- read_excel(
  source_file,
  sheet = "feasibility_performance"
) %>%
  transmute(
    Method = as.character(Strategy),
    Time = as.numeric(Mean_construction_time_min),
    Timeout = as.integer(Timeout),
    Failure = as.integer(Failure)
  )

# ------------------------------
# 3. Sort by construction time
# Shorter construction time at top
# ------------------------------
df_plot <- df_plot %>%
  arrange(Time) %>%
  mutate(
    Method = factor(
      Method,
      levels = rev(Method)
    ),
    Time_label = sprintf("%.2f", Time)
  )

# ------------------------------
# 4. Convert status data to long format
# ------------------------------
status_data <- bind_rows(
  df_plot %>%
    filter(Timeout == 1) %>%
    transmute(
      Method,
      Status = "Timeout"
    ),
  
  df_plot %>%
    filter(Failure == 1) %>%
    transmute(
      Method,
      Status = "Failure"
    )
) %>%
  mutate(
    Status = factor(
      Status,
      levels = c("Timeout", "Failure")
    ),
    Status_x = case_when(
      Status == "Timeout" ~ -5.5,
      Status == "Failure" ~ -2.4
    )
  )

# ------------------------------
# 5. Plot positions
# ------------------------------
x_min <- -7.0
x_max <- 21.5

status_split_x <- -3.95

# ------------------------------
# 6. Generate the plot
# ------------------------------
fig2C <- ggplot(
  df_plot,
  aes(
    x = Time,
    y = Method
  )
) +
  
  # Construction-time bars
  geom_col(
    width = 0.68,
    fill = "#8E7CC3",
    color = "black",
    linewidth = 0.25
  ) +
  
  # Construction-time values
  geom_text(
    aes(
      x = Time + 0.28,
      label = Time_label
    ),
    hjust = 0,
    size = 2.3,
    family = "Helvetica",
    color = "black"
  ) +
  
  # Divider between the two status columns
  geom_vline(
    xintercept = status_split_x,
    linetype = "22",
    linewidth = 0.25,
    color = "grey70"
  ) +
  
  # Divider between the status and time regions
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    linewidth = 0.45,
    color = "black"
  ) +
  
  # Timeout and failure symbols
  geom_point(
    data = status_data,
    aes(
      x = Status_x,
      y = Method,
      shape = Status
    ),
    inherit.aes = FALSE,
    size = 2.0,
    stroke = 0.80,
    color = "#D94B3D"
  ) +
  
  scale_shape_manual(
    name = NULL,
    values = c(
      "Timeout" = 16,
      "Failure" = 4
    )
  ) +
  
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = c(0, 5, 10, 15, 20),
    labels = c("0", "5", "10", "15", "20"),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  
  scale_y_discrete(
    expand = expansion(
      add = c(0.35, 0.35)
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  labs(
    x = "Mean construction time (min)",
    y = NULL
  ) +
  
  guides(
    shape = guide_legend(
      direction = "horizontal",
      override.aes = list(
        color = "#D94B3D",
        size = 2.1,
        stroke = 0.8
      )
    )
  ) +
  
  theme_classic(
    base_size = 8,
    base_family = "Helvetica"
  ) +
  
  theme(
    plot.title = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.50
    ),
    
    axis.line = element_blank(),
    
    axis.text.x = element_text(
      size = 7,
      color = "black",
      family = "Helvetica"
    ),
    
    axis.text.y = element_text(
      size = 6.2,
      color = "black",
      family = "Helvetica",
      margin = margin(r = 3)
    ),
    
    axis.title.x = element_text(
      size = 7.5,
      color = "black",
      family = "Helvetica",
      margin = margin(t = 4)
    ),
    
    axis.title.y = element_blank(),
    
    axis.ticks.x = element_line(
      color = "black",
      linewidth = 0.25
    ),
    
    axis.ticks.y = element_line(
      color = "black",
      linewidth = 0.25
    ),
    
    axis.ticks.length = unit(
      0.10,
      "cm"
    ),
    
    # Status legend placed inside the upper-right region
    legend.position = c(0.78, 0.94),
    
    legend.justification = c(
      0.5,
      1
    ),
    
    legend.direction = "horizontal",
    
    legend.text = element_text(
      size = 6.2,
      family = "Helvetica",
      color = "black"
    ),
    
    legend.key.size = unit(
      0.30,
      "cm"
    ),
    
    legend.spacing.x = unit(
      0.10,
      "cm"
    ),
    
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.25
    ),
    
    legend.margin = margin(
      t = 2,
      r = 3,
      b = 2,
      l = 3
    ),
    
    plot.margin = margin(
      t = 2,
      r = 4,
      b = 3,
      l = 3
    )
  )

fig2C

# ------------------------------
# 7. Save editable PowerPoint
# ------------------------------
fig_width  <- 7 / 2.54
fig_height <- 14 / 2.54

ppt <- read_pptx() %>%
  add_slide(
    layout = "Blank",
    master = "Office Theme"
  ) %>%
  ph_with(
    dml(
      ggobj = fig2C
    ),
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
    "Figure_3b_computational_feasibility_status_legend.pptx"
  )
)