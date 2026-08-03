# ============================================================
# Kendall's W ranking agreement across datasets
#
# Input:
# Raw strategy scores for five accuracy metrics
#
# Output:
# Kendall_W_ranking_agreement.xlsx
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(writexl)

# ------------------------------
# 1. File paths
# ------------------------------
file <- "../scGEM_benchmark/result/accuracy_performance.xlsx"

out_dir <- "../scGEM_benchmark/result"

if (!dir.exists(out_dir)) {
  dir.create(
    out_dir,
    recursive = TRUE
  )
}

output_file <- file.path(
  out_dir,
  "Kendall_W_ranking_agreement.xlsx"
)

# ------------------------------
# 2. Read the input data
# ------------------------------
raw <- read_excel(
  file,
  sheet = "Raw_scores",
  col_names = FALSE
)

# Remove completely empty rows and columns
raw <- raw %>%
  select(
    where(~ !all(is.na(.)))
  ) %>%
  filter(
    if_any(
      everything(),
      ~ !is.na(.)
    )
  )

# Row 1: metric names
# Row 2: dataset names
# Row 3 onward: strategy scores
raw_data <- raw[-c(1, 2), ]

# Remove rows without a strategy name
valid_rows <- !is.na(raw_data[[1]]) &
  as.character(raw_data[[1]]) != ""

raw_data <- raw_data[valid_rows, ]

strategy_names <- as.character(
  raw_data[[1]]
)

# ------------------------------
# 3. Metric column ranges
# ------------------------------
metric_cols <- list(
  "Housekeeping genes"     = 5:13,
  "Metabolic tasks"        = 14:22,
  "Marker genes"           = 23:32,
  "Clustering performance" = 33:41,
  "Essential genes"        = 42:46
)

# ------------------------------
# 4. Extract one accuracy metric
# ------------------------------
extract_metric_data <- function(metric_name, cols) {
  
  dataset_names <- as.character(
    unlist(
      raw[2, cols]
    )
  )
  
  dataset_names <- make.unique(
    dataset_names
  )
  
  metric_df <- bind_cols(
    data.frame(
      Strategy = strategy_names
    ),
    as.data.frame(
      raw_data[, cols]
    )
  )
  
  colnames(metric_df) <- c(
    "Strategy",
    dataset_names
  )
  
  metric_df <- metric_df %>%
    pivot_longer(
      cols = -Strategy,
      names_to = "Dataset",
      values_to = "RawScore"
    ) %>%
    mutate(
      RawScore = suppressWarnings(
        as.numeric(RawScore)
      ),
      Metric = metric_name
    ) %>%
    filter(
      !is.na(RawScore),
      is.finite(RawScore)
    )
  
  return(metric_df)
}

# ------------------------------
# 5. Combine the five metrics
# ------------------------------
score_long <- bind_rows(
  lapply(
    names(metric_cols),
    function(metric_name) {
      extract_metric_data(
        metric_name,
        metric_cols[[metric_name]]
      )
    }
  )
)

# ------------------------------
# 6. Kendall's W function
# ------------------------------
calculate_kendall_w <- function(score_matrix) {
  
  score_matrix <- as.matrix(
    score_matrix
  )
  
  storage.mode(score_matrix) <- "numeric"
  
  n_strategies <- nrow(
    score_matrix
  )
  
  n_datasets <- ncol(
    score_matrix
  )
  
  # Higher raw score receives a better rank
  rank_matrix <- apply(
    score_matrix,
    2,
    function(x) {
      rank(
        -x,
        ties.method = "average"
      )
    }
  )
  
  rank_sums <- rowSums(
    rank_matrix
  )
  
  S <- sum(
    (
      rank_sums -
        mean(rank_sums)
    )^2
  )
  
  # Tie correction
  tie_correction <- 0
  
  for (j in seq_len(n_datasets)) {
    
    tie_sizes <- table(
      score_matrix[, j]
    )
    
    tie_sizes <- tie_sizes[
      tie_sizes > 1
    ]
    
    if (length(tie_sizes) > 0) {
      tie_correction <- tie_correction +
        sum(
          tie_sizes^3 -
            tie_sizes
        )
    }
  }
  
  denominator <- (
    n_datasets^2 *
      (
        n_strategies^3 -
          n_strategies
      )
  ) -
    (
      n_datasets *
        tie_correction
    )
  
  if (
    denominator <= 0 ||
    !is.finite(denominator)
  ) {
    W <- NA_real_
  } else {
    W <- 12 * S / denominator
  }
  
  chi_square <- n_datasets *
    (
      n_strategies - 1
    ) *
    W
  
  degrees_freedom <- n_strategies - 1
  
  p_value <- pchisq(
    chi_square,
    df = degrees_freedom,
    lower.tail = FALSE
  )
  
  data.frame(
    Kendall_W = W,
    Chi_square = chi_square,
    Degrees_freedom = degrees_freedom,
    P_value = p_value,
    N_strategies = n_strategies,
    N_datasets = n_datasets
  )
}

# ------------------------------
# 7. Calculate Kendall's W
# ------------------------------
kendall_results <- bind_rows(
  lapply(
    names(metric_cols),
    function(metric_name) {
      
      metric_df <- score_long %>%
        filter(
          Metric == metric_name
        ) %>%
        select(
          Strategy,
          Dataset,
          RawScore
        ) %>%
        pivot_wider(
          names_from = Dataset,
          values_from = RawScore
        )
      
      strategy_ids <- metric_df$Strategy
      
      score_matrix <- as.matrix(
        metric_df[, -1]
      )
      
      rownames(score_matrix) <- strategy_ids
      
      # Retain strategies with scores in every applicable dataset
      score_matrix <- score_matrix[
        complete.cases(score_matrix),
        ,
        drop = FALSE
      ]
      
      result_i <- calculate_kendall_w(
        score_matrix
      )
      
      result_i$Metric <- metric_name
      
      result_i %>%
        select(
          Metric,
          everything()
        )
    }
  )
)

# ------------------------------
# 8. Save the result table
# ------------------------------
write_xlsx(
  list(
    Kendall_W_ranking_agreement = kendall_results
  ),
  path = output_file
)

disp_results <- kendall_results
print(disp_results)

message(
  "Kendall's W results saved to: ",
  output_file
)