library(kableExtra)
library(dplyr)
library(purrr)
library(ggplot2)
library(readr)
library(stringr)

# --------------------------
# Description:
#
# Last Modified: 2025-12-08
# Author: Johannes Held
# --------------------------
# 1. Read all CSV files
# --------------------------
base_dir <- "../../experiments/raw"
folders <- c(
  file.path(base_dir, "A0_VS_A2", "2025-08-22"),
  file.path(
    base_dir,
    "A1_VS_A2",
    "2025-08-26",
    paste0("A1_scaling_", 1:3)
  )
)

files <- map(folders, ~ list.files(.x, pattern = "\\.csv$", full.names = TRUE)) %>%
  unlist()

df <- files %>% 
  lapply(read.csv, sep = ";", header = TRUE, stringsAsFactors = FALSE) %>% 
  bind_rows()

df <- df %>% filter(algorithm == 1 | algorithm == 2)

avg_data <- df %>%
  # remove rows with any NA values
  filter(if_all(everything(), ~ !is.na(.))) %>%
  # group by experiment parameters
  group_by(nodes, algorithm, m, n) %>%
  # compute summary statistics
  summarise(
    Average_Time = mean(runtime_ms, na.rm = TRUE),
    StdDev_Time  = sd(runtime_ms, na.rm = TRUE),
    Runs         = n(),   # number of runs
    .groups = 'drop'
  ) %>%
  mutate(
    nodes = as.factor(nodes),
    algorithm = as.factor(algorithm),
    combined_size = paste0(m, "x", n)
  ) %>%
  mutate(
    shape_type = case_when(
      m == n ~ "Square",
      m < n  ~ "Fat",
      m > n  ~ "Thin"
    ),
    size_metric = m * n,
    combined_size = paste0(m, "x", n)
  ) %>%
  arrange(shape_type, size_metric) %>%
  mutate(
    combined_size = factor(combined_size, levels = unique(combined_size)),
    Gflops = (n * m^2) / (1e6 * Average_Time)  # runtime_ms to seconds (1e6)
  )

okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000"
)

avg_data_tmp <- avg_data %>%
  # Create a dedicated factor for the legend
  mutate(
    Group = case_when(
      algorithm == 1 ~ "A1",
      algorithm == 2 & nodes == 1 ~ "Baseline (A2 with 1 Node)",
      TRUE ~ "Other"
    ),
    Group = factor(Group, levels = c("Baseline (A2 with 1 Node)", "A1"))
  ) %>%
  mutate(
    # Combine the node count and algorithm for the color aesthetic
    Combined_Color_Factor = case_when(
      algorithm == 2 & nodes == 1 ~ "Baseline (A2 with 1 node)",
      algorithm == 1 ~ paste("A1 -", nodes, "nodes"),
      TRUE ~ "Other"
    )
  )



# ------------------------------------------------------------------
# 1. ORDERING FIX: Set the Levels for Combined_Color_Factor
# ------------------------------------------------------------------

# a. Get the unique node counts for A0 and sort them numerically
a1_nodes_sorted <- unique(avg_data_tmp$nodes[avg_data_tmp$algorithm == 1]) %>% sort()

# b. Create the ordered labels for the A0 cases
a1_labels_ordered <- paste("A1 -", a1_nodes_sorted, "nodes")

# c. Define the final, complete, ordered levels list
ordered_levels <- c(
  "Baseline (A2 with 1 node)",
  a1_labels_ordered
)

# d. Apply the new levels to the factor
avg_data_tmp <- avg_data_tmp %>%
  mutate(
    Combined_Color_Factor = factor(Combined_Color_Factor, levels = ordered_levels)
  )

# ------------------------------------------------------------------
# 2. Rebuild the color_map with the correctly ordered labels
# ------------------------------------------------------------------
BASELINE_COLOR <- "#F00000"

# Re-extract the factors (now correctly ordered)
ordered_a1_factors <- ordered_levels[-1] # Remove Baseline, keep ordered A0

# 2. Get the colors for the A0 lines (excluding the one used for baseline)
#    We use the first 'n' colors from okabe_ito, where 'n' is the number of unique A0 node counts.
A1_NODE_COLORS <- okabe_ito[1:length(unique(avg_data_tmp$nodes[avg_data_tmp$algorithm == 1]))]

# 3. Create the named vector for scale_color_manual using the ordered A0 factors
color_map <- c(
  "Baseline (A2 with 1 node)" = BASELINE_COLOR,
  # This line ensures the colors in A0_NODE_COLORS (which start with #E69F00) 
  # are matched positionally to the labels in ordered_a0_factors (which start 
  # with "A1 - 1 nodes", then "A1 - 2 nodes", etc.)
  A1_NODE_COLORS %>% setNames(ordered_a1_factors)
)
# Ensure the map doesn't contain the "Other" category if it exists.
color_map <- color_map[!names(color_map) %in% "Other"]

# --------------------------------------------------
# Plot 1:
# --------------------------------------------------
p_runtime <- ggplot(avg_data_tmp, aes(
  x = combined_size,
  y = Average_Time,
  color = Combined_Color_Factor,
  shape = Group
)) +
  #geom_errorbar(
  #  aes(ymin = Average_Time - StdDev_Time, ymax = Average_Time + StdDev_Time),
  #  width = 0.2, alpha = 0.5
  #) +
  geom_line(aes(group = interaction(as.factor(nodes), algorithm)), linewidth = .5) +
  geom_point(size = 2, aes(shape = Group)) +
  labs(
    title = "Runtime Comparison: A1 Scaling vs. Baseline (A2 with 1 Node)",
    subtitle = "Performance grouped by matrix shape",
    x = "Matrix Size (m x n)",
    y = "Average Runtime (seconds)",
    color = "",
    shape = "Algorithm"
  ) +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = color_map) + # Select distinct, colorblind-friendly colors
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_wrap(~shape_type, scales = "free_x", nrow = 1)
p_runtime

# --------------------------------------------------
# Plot 2: y log scaled
# --------------------------------------------------
p_runtime_log <- p_runtime + scale_y_log10() +
  labs(
    subtitle = "Logarithmic Y-Axis"
  )
p_runtime_log

# --------------------------------------------------
# Plot 3: Performance (Gflops/s) vs. Matrix Size
# --------------------------------------------------
p_perf <- ggplot(avg_data_tmp %>% filter(algorithm == 1), aes(
  x = combined_size,
  y = Gflops,
  color = Combined_Color_Factor,
  shape = Group
)) +
  geom_line(aes(group = interaction(nodes, algorithm), linetype = as.factor(nodes)), linewidth = .5) +
  geom_point(size = 2) +
  labs(
    title = "Performance (Gflops/s) vs. Matrix Size",
    subtitle = "Computed as n * m² / runtime",
    x = "Matrix Size (m x n)",
    y = "Performance (Gflops/s)",
    shape = "Algorithm",
    color = "",
  ) +
  guides(
    color = guide_legend(order = 1), # Color legend appears first
    shape = guide_legend(order = 2)  # Shape legend appears second
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = color_map) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1)
  ) +
  facet_wrap(~shape_type, scales = "free_x", nrow = 1) +
  guides(linetype = "none")
p_perf

# --------------------------------------------------
# Compute Speedup = T(1 node) / T(k nodes)
# --------------------------------------------------
speedup_data <- avg_data_tmp %>%
  group_by(m, n) %>%
  mutate(
    baseline_time = Average_Time[nodes == "1" & algorithm == "2"],   # runtime for 1 node
    Speedup = baseline_time / Average_Time
  ) %>%
  ungroup() %>%
  filter(algorithm == 1)

# --------------------------------------------------
# Speedup Plot
# --------------------------------------------------
p_speedup <- ggplot(speedup_data, aes(
  x = combined_size,
  y = Speedup,
  color = Combined_Color_Factor,
  shape = Group
)) +
  geom_line(aes(group = interaction(nodes, algorithm), linetype = nodes), linewidth = .5) +
  geom_point(size = 2) +
  labs(
    title = "Absolute Speedup of A1",
    subtitle = "Computed as Tseq(m*n) / Tpar(p, m*n)",
    x = "Matrix Size (m x n)",
    y = "Speedup",
    color = "",
    shape = "Algorithm"
  ) +
  guides(
    color = guide_legend(order = 1), # Color legend appears first
    shape = guide_legend(order = 2)  # Shape legend appears second
  ) +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = color_map) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1)
  ) +
  facet_wrap(~shape_type, scales = "free_x", nrow = 1) +
  guides(linetype = "none")
p_speedup


# --------------------------------------------------
# Table output
# --------------------------------------------------

# Ensure algorithm is numeric
avg_data_tmp <- avg_data_tmp %>%
  mutate(algorithm = as.numeric(as.character(algorithm)))

# Filter only A0
a1_data <- avg_data_tmp %>% filter(algorithm == 1)

# Baseline: 1-node runtime of A2 (for absolute speedup)
baseline_a2 <- avg_data_tmp %>%
  filter(algorithm == 2, nodes == 1) %>%
  select(m, n, A2_1node = Average_Time)

# Baseline: 1-node runtime of A0 (for relative speedup)
baseline_a1 <- a1_data %>%
  filter(nodes == 1) %>%
  select(m, n, A1_1node = Average_Time)

# Join with A0 and compute speedups
speedup_stats <- a1_data %>%
  left_join(baseline_a2, by = c("m", "n")) %>%
  left_join(baseline_a1, by = c("m", "n")) %>%
  mutate(
    Absolute_Speedup = A2_1node / Average_Time,
    Relative_Speedup = A1_1node / Average_Time
  )

# Summarize per shape type
summary_table_by_shape <- speedup_stats %>%
  group_by(shape_type) %>%
  summarise(
    Min_Absolute = min(Absolute_Speedup),
    Max_Absolute = max(Absolute_Speedup),
    Mean_Absolute = mean(Absolute_Speedup),
    SD_Absolute = sd(Absolute_Speedup),
    Min_Relative = min(Relative_Speedup),
    Max_Relative = max(Relative_Speedup),
    Mean_Relative = mean(Relative_Speedup),
    SD_Relative = sd(Relative_Speedup),
    .groups = "drop"
  )
summary_table_by_shape_renamed <- summary_table_by_shape %>%
  rename(
    Shape = shape_type,
    `Min Abs` = Min_Absolute,
    `Max Abs` = Max_Absolute,
    `Mean Abs` = Mean_Absolute,
    `SD Abs` = SD_Absolute,
    `Min Rel` = Min_Relative,
    `Max Rel` = Max_Relative,
    `Mean Rel` = Mean_Relative,
    `SD Rel` = SD_Relative
  )
View(summary_table_by_shape)

# Produce LaTeX table
summary_table_by_shape %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    digits = 4,
    col.names = c(
      "Shape", "Min Abs", "Max Abs", "Mean Abs", "SD Abs",
      "Min Rel", "Max Rel", "Mean Rel", "SD Rel"
    ),
    caption = "Statistics of Absolute and Relative Speedup of A1 by Matrix Shape"
  ) %>%
  kable_styling(
    latex_options = c("striped", "hold_position", "scale_down")
  )


