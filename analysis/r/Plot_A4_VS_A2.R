library(kableExtra)
library(dplyr)
library(purrr)
library(ggplot2)
library(readr)
library(stringr)

# ==================================================
# Runtime Scaling Analysis for Algorithm A4
# ==================================================
# Description:
# Aggregates and visualizes runtime scaling results
# for algorithm A4, comparing average runtime,
# performance (GFLOPS), and speedup over matrix
# sizes, decomposition parameters (c, i), and shapes.
#
# Author: Johannes Held
# Last Modified: 2025-10-13
# ==================================================

theme_set(theme_minimal(base_size = 10))

color_values <- c(
  "#F00000", # 1. Baseline
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999",
  "#000000", "#8B4513", # Your original 10 + Baseline
  "#3CB371", "#4682B4", "#FFA07A", "#A0522D"  # 5 additions
)

process_a4_data <- function(base_dir, subfolders, gflops_divisor = 1e9) {
  
  # 1. Load Data
  files <- file.path(base_dir, subfolders) |>
    map(~ list.files(.x, pattern = "\\.csv$", full.names = TRUE)) |>
    unlist(use.names = FALSE)
  
  # Read all files and combine into one data frame
  df <- files |>
    map_dfr(~ read.csv(.x, sep = ";", header = TRUE, stringsAsFactors = FALSE))
  
  # 2. Aggregate and Compute Stats
  avg_data <- df |>
    filter(if_all(everything(), ~ !is.na(.))) |>
    # Group by all unique experimental parameters
    group_by(c, i, algorithm, m, n, nodes) |>
    summarise(
      Average_Time = mean(runtime_ms, na.rm = TRUE),
      StdDev_Time  = sd(runtime_ms, na.rm = TRUE),
      Runs         = n(),
      .groups = "drop"
    ) |>
    
    # 3. Mutate and Calculate Metrics
    mutate(
      # Convert factors *after* summarization for efficiency
      c = c,
      i = i,
      nodes = nodes,
      algorithm = as.factor(algorithm),
      
      # Calculate shape type
      shape_type = case_when(
        # Use tolerance for square matrices for both datasets
        abs(m / n - 1) < 0.10 ~ "Square", 
        m < n  ~ "Fat",
        m > n  ~ "Thin"
      ),
      size_metric = m * n,
      
      # Gflops calculation: Uses parameter for divisor
      # NOTE: If runtime_ms is in milliseconds, divisor should be 1e6. 
      # Using 1e9 as per original code, but this assumes nanoseconds.
      Gflops = (n * m^2) / (gflops_divisor * Average_Time)
    ) %>%
    
    # Sort by shape and size_metric for proper factor ordering later
    arrange(shape_type, size_metric) %>%
    
    # Create combined_size factor based on the correct, sorted order
    mutate(
      combined_size = factor(paste0(m, "x", n), 
                             levels = unique(paste0(m, "x", n)))
    )
  
  return(avg_data)
}


# --- Definitions for Block 1 (Older Data) ---
base_dir_1 <- "/Users/johannes/Downloads/Uni/Bachelor Arbeit/A4_Scaling/2025-10-15"

subfolders_1 <- sprintf(
  "A4_scaling_%s",
  c("02", "1_2", "1_3", "1_4", "1_5", "1_6",
    "2_2", "2_3", "2_4", "2_5", "2_6",
    "3_2", "3_3")
)

# --- Definitions for Block 2 (Newer Data) ---
base_dir_2 <- "/Users/johannes/Downloads/Uni/Bachelor Arbeit/A4_Scaling/2025-12-03"
subfolders_2 <- sprintf(
  "A4_scaling_%s",
  c("fat_2", "fat_3", "square_3", "square_2", "thin_2", "thin_3", "reference")
)

# --- Processing ---

# Process Block 1 (using original Gflops divisor)
avg_data_1 <- process_a4_data(base_dir_1, subfolders_1, gflops_divisor = 1e9)

# Process Block 2 (using original Gflops divisor)
avg_data_2 <- process_a4_data(base_dir_2, subfolders_2, gflops_divisor = 1e9)


# --- Apply Filters (These should ideally be inside the function if consistent) ---

# Apply filters from original script:
#avg_data_1 <- avg_data_1 |> filter(shape_type != "Square")
#avg_data_1 <- avg_data_1 |> filter(c == 3 | algorithm == 2)


# --- Final Merge ---

# Ensure consistent column structure for safe row binding
# Since the function is identical, they should share the structure.
# We skip the unnecessary union/subsetting step.

combined_data <- bind_rows(
  mutate(avg_data_1, source = "old"),
  mutate(avg_data_2, source = "new")
)

# Final cleanup (arranging and re-factoring combined_size)
combined_data <- combined_data |>
  arrange(size_metric) |>
  mutate(
    # Re-factor based on the new, combined, sorted data
    combined_size = factor(combined_size, levels = unique(combined_size)),
    
    # Create Group factor for shape/legend (as intended)
    Group = case_when(
      algorithm == 4 ~ "A4",
      algorithm == 2 & nodes == 1 ~ "Baseline (A2 with 1 Node)",
      TRUE ~ "Other"
    ),
    Group = factor(Group, levels = c("Baseline (A2 with 1 Node)", "A4")),
    
    # Create Combined_Color_Factor based on 'c' value for A4
    Combined_Color_Factor = case_when(
      algorithm == 2 & nodes == 1 ~ "Baseline (A2 with 1 node)",
      algorithm == 4 ~ paste0("A4 (c=", c, ", i=", i, ", nodes=", c * (c + 1) * i, ")"),
      TRUE ~ "Other"
    )
  ) 
###############################################################################

# a. Get the unique node counts for A4 and sort them numerically
a0_labels_ordered <- unique(combined_data$Combined_Color_Factor[combined_data$algorithm == 4]) %>% sort()

# c. Define the final, complete, ordered levels list
ordered_levels <- c(
  "Baseline (A2 with 1 node)",
  a0_labels_ordered
)

# d. Apply the new levels to the factor
combined_data <- combined_data %>%
  mutate(
    Combined_Color_Factor = factor(Combined_Color_Factor, levels = ordered_levels)
  )

# 2. Assign the colors to the factors (assuming 'ordered_levels' is a vector of 15 factor levels)
color_map <- setNames(color_values, ordered_levels)

# --------------------------
# 3. Plot: Average Runtime
# --------------------------
p_runtime <- ggplot(combined_data, aes(
  x = combined_size, 
  y = Average_Time,
  color = Combined_Color_Factor, 
  #shape = as.factor(c), 
  shape = as.factor(i),
  linetype = as.factor(c),
)) +
  geom_line(aes(group = interaction(c, i, algorithm)), linewidth = 0.6) +
  geom_point(size = 2) +
  scale_color_manual(values = color_map) +
  facet_wrap(~shape_type, scales = "free_x", nrow = 1) +
  labs(
    title = "Average Runtime vs. Matrix Size",
    subtitle = "Performance grouped by matrix shape",
    x = "Matrix Size (m × n)",
    y = "Average Runtime (ms)",
    color = "",
    shape = "i",
    linetype = "c"
  ) +
  guides(
    color = guide_legend(order = 1), # Color legend appears first
    linetype = guide_legend(order = 2)  # Shape legend appears second
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1)
  )
p_runtime

# Log-scaled variant
p_runtime_log <- p_runtime + scale_y_log10() +
  labs(subtitle = "Logarithmic Y-Axis")
p_runtime_log

# --------------------------
# 4. Plot: GFLOPS
# --------------------------
p_perf <- ggplot(combined_data, aes(
  x = combined_size, 
  y = Gflops,
  color = Combined_Color_Factor,              # color now maps to iteration factor
  shape = as.factor(i),              # shape now maps to column decomposition
  linetype = as.factor(c)    # different line type for each algorithm (e.g., baseline)
)) +
  geom_line(aes(group = interaction(c, i, algorithm)), linewidth = 0.5) +
  geom_point(size = 2) +
  scale_color_manual(values = color_map) +
  facet_wrap(~shape_type, scales = "free_x", nrow = 1) +
  labs(
    title = "Performance (GFLOPS/s) vs. Matrix Size",
    subtitle = "Computed as n × m² / runtime",
    x = "Matrix Size (m × n)",
    y = "Performance (GFLOPS/s)",
    color = "",
    shape = "i",
    linetype = "c"
  ) +
  guides(
    color = guide_legend(order = 1), # Color legend appears first
    linetype = guide_legend(order = 2)  # Shape legend appears second
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1)
  )
p_perf


# --------------------------
# 5. Speedup Computation
# --------------------------
# --------------------------------------------------
# Compute Speedup (using algorithm 2 with c=1, i=2 as baseline)
# --------------------------------------------------
speedup_data <- combined_data |>
  group_by(m, n) |>
  mutate(
    # Define your baseline explicitly — adjust if needed
    baseline_time = Average_Time[algorithm == 2],
    Speedup = baseline_time / Average_Time
  ) |>
  ungroup()

# --------------------------------------------------
# Speedup Plot (using c and i instead of nodes)
# --------------------------------------------------
p_speedup <- ggplot(speedup_data, aes(
  x = combined_size,
  y = Speedup,
  color = Combined_Color_Factor,              # color for iteration factor
  shape = as.factor(i),              # shape for column decomposition
  linetype = as.factor(c)    # separate algorithms by line type
)) +
  geom_line(aes(group = interaction(c, i, algorithm)), linewidth = 0.5) +
  geom_point(size = 2) +
  scale_color_manual(values = color_map) +
  facet_wrap(~shape_type, scales = "free_x", nrow = 1) +
  labs(
    title = "Speedup vs. Matrix Size",
    subtitle = "Normalized to baseline runtime (A2)",
    x = "Matrix Size (m × n)",
    y = "Speedup",
    color = "",
    shape = "i",
    linetype = "c"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1)
  )

p_speedup


# --------------------------
# 6. Speedup Summary Table
# --------------------------
# Ensure algorithm is numeric
combined_data <- combined_data |>
  mutate(algorithm = as.numeric(as.character(algorithm)))

# Select A3 (the improved algorithm)
a4_data <- filter(combined_data, algorithm == 4)

# --------------------------------------------------
# Baseline definitions (adjust if needed)
# --------------------------------------------------
# Baseline for Absolute Speedup: Algorithm 2 with smallest c and i (e.g. c=1, i=2)
baseline_a2 <- combined_data |>
  filter(algorithm == 2) |>
  select(m, n, A2_baseline = Average_Time)

# Baseline for Relative Speedup: Algorithm 3 with smallest c and i (e.g. c=1, i=2)
baseline_a4 <- a4_data |>
  filter(c == 1, i == 2) |>
  select(m, n, A4_baseline = Average_Time)

# --------------------------------------------------
# Compute speedups
# --------------------------------------------------
speedup_stats <- a4_data |>
  left_join(baseline_a2, by = c("m", "n")) |>
  left_join(baseline_a4, by = c("m", "n")) |>
  mutate(
    Absolute_Speedup = A2_baseline / Average_Time,
    Relative_Speedup = A4_baseline / Average_Time
  )

# --------------------------------------------------
# Summarize by shape_type, c, and i
# --------------------------------------------------
summary_table_by_shape_ci <- speedup_stats |>
  group_by(shape_type, c, i) |>
  summarise(
    across(
      .cols = c(Absolute_Speedup, Relative_Speedup),
      .fns = list(Min = min, Max = max, Mean = mean, SD = sd),
      .names = "{.fn}_{.col}"
    ),
    .groups = "drop"
  ) |>
  rename(
    Shape = shape_type,
    `c (column decomposition)` = c,
    `i (iteration factor)` = i
  )
View(summary_table_by_shape_ci)

# --------------------------------------------------
# Output as LaTeX table
# --------------------------------------------------
summary_table_by_shape_ci |>
  kable(
    format = "latex",
    booktabs = TRUE,
    digits = 3,
    col.names = c(
      "Shape", "c", "i",
      "Min Abs", "Max Abs", "Mean Abs", "SD Abs",
      "Min Rel", "Max Rel", "Mean Rel", "SD Rel"
    ),
    caption = "Statistics of Absolute and Relative Speedup of A3 by Matrix Shape, c, and i"
  ) |>
  kable_styling(
    latex_options = c("striped", "hold_position", "scale_down")
  )
