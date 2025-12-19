library(kableExtra)
library(dplyr)
library(purrr)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# --------------------------
# Description:
#
# Last Modified: 2025-10-09
# Author: Johannes Held
# --------------------------
# 1. Read all CSV files
# --------------------------

process_exp_data <- function(folder_paths, c_value = NULL) {
  
  # 1. Load Data
  files <- map(folder_paths, ~ list.files(.x, pattern = "\\.csv$", full.names = TRUE)) %>%
    unlist()
  
  df <- files %>%
    # Use map_dfr for combined reading and binding, more idiomatic than map + bind_rows
    map_dfr(~ read.csv(.x, sep = ";", header = TRUE, stringsAsFactors = FALSE)) 
  
  # 2. Summarize Data
  avg_data <- df %>%
    # Filter NA rows (using !is.na(.) is slightly cleaner)
    filter(if_all(everything(), ~ !is.na(.))) %>%
    
    group_by(nodes, algorithm, m, n) %>%
    summarise(
      Average_Time = mean(runtime_ms, na.rm = TRUE),
      StdDev_Time  = sd(runtime_ms, na.rm = TRUE),
      Runs         = n(),
      .groups = 'drop'
    ) %>%
    
    # 3. Mutate (Standard Calculations)
    mutate(
      nodes = as.numeric(as.character(nodes)),
      algorithm = as.factor(as.character(algorithm)),
      
      # Use the more sophisticated shape calculation from your second block
      shape_type = case_when(
        abs(m / n - 1) < 0.10 ~ "Square", # within ±10%
        m < n  ~ "Fat",
        m > n  ~ "Thin"
      ),
      size_metric = m * n,
      combined_size = paste0(m, "x", n),
      Gflops = (n * m^2) / (1e6 * Average_Time) # runtime_ms -> seconds
    )
  
  # 4. Handle 'c' Calculation (Specific to A3)
  if (is.null(c_value)) {
    # Calculation for the first dataset (A3 scaling based on nodes)
    avg_data <- avg_data %>%
      mutate(
        c_calc = case_when(
          algorithm == 3 ~ floor((-1 + sqrt(1 + 4 * nodes)) / 2),
          TRUE ~ NA_real_ # Ensure NA is numeric type for consistency
        )
      )
  } else {
    # Assignment for the second dataset (where c is a fixed value like 31)
    avg_data <- avg_data %>%
      mutate(
        c_calc = if_else(algorithm == 3, as.numeric(c_value), NA_real_)
      )
  }
  
  return(avg_data)
}

base_dir <- "../../experiments/raw"
# Block 1 (A3 Scaling)
folders_old <- file.path(
  base_dir,
  "A3_Scaling",
  "2025-10-12",
  paste0("A3_scaling_", 0:2)
)

# Block 2 (New Experiments/Reference)
folders_new <- file.path(
  base_dir,
  "A3_Scaling",
  "2025-11-28",
  paste0("A3_scaling_", c("square", "thin", "fat", "reference"))
)

# --- Process Data Blocks ---
# Process data block 1 (c is calculated based on nodes)
avg_data_old <- process_exp_data(folders_old, c_value = NULL)

# Process data block 2 (c is fixed at 31 for A3, 0 or NA for others)
# I used the assumption that c=31 is only relevant for A3 in this set.
avg_data_new <- process_exp_data(folders_new, c_value = 31)


# --- Final Merge and Cleanup ---

# Rename c_calc back to 'c' and ensure both tables have identical columns
avg_data_old <- avg_data_old %>% rename(c = c_calc) %>% mutate(source = "old")
avg_data_new <- avg_data_new %>% rename(c = c_calc) %>% mutate(source = "new")

# Combine the data frames. Since the function ensured identical column structure, 
# this is simple and safe.
combined_data <- bind_rows(avg_data_old, avg_data_new)

# Final formatting and factor creation
combined_data <- combined_data %>%
  arrange(size_metric) %>%
  mutate(
    # Set combined_size as an ordered factor based on sorted size_metric
    combined_size = factor(combined_size, levels = unique(combined_size)),
    
    # Create Group factor for shape/legend (as intended)
    Group = case_when(
      algorithm == 3 ~ "A3",
      algorithm == 2 & nodes == 1 ~ "Baseline (A2 with 1 Node)",
      TRUE ~ "Other"
    ),
    Group = factor(Group, levels = c("Baseline (A2 with 1 Node)", "A3")),
    
    # Create Combined_Color_Factor based on 'c' value for A3
    Combined_Color_Factor = case_when(
      algorithm == 2 & nodes == 1 ~ "Baseline (A2 with 1 node)",
      algorithm == 3 ~ paste0("A3 (c=", c, ", nodes=", c * (c + 1), ")"),
      TRUE ~ "Other"
    )
  )

okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000"
)

# --------------------------------------------------
# Plot 1:
# --------------------------------------------------
p_runtime <- ggplot(combined_data, aes(
  x = combined_size,
  y = Average_Time,
  color = Combined_Color_Factor,               # ensure c is treated as a factor
  linetype = Group,     # ensure algorithm is treated as a factor
)) +
  geom_line(aes(group = interaction(c, algorithm)), linewidth = 0.6) +
  geom_point(size = 2) +
  labs(
    title = "Runtime Comparison: A3 Scaling vs. Baseline (A2 with 1 Node)",
    subtitle = "Performance grouped by matrix shape",
    x = "Matrix Size (m × n)",
    y = "Average Runtime (seconds)",
    color = "",
    shape = "Algorithm",
    linetype = "Algorithm"
  ) +
  guides(
    color = guide_legend(order = 1), # Color legend appears first
    shape = guide_legend(order = 2)  # Shape legend appears second
  ) +
  scale_color_manual(values = okabe_ito) +
  theme_minimal(base_size = 11) +
  theme( legend.position = "right", axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1), legend.box = "vertical" ) +
  facet_wrap(~shape_type, scales = "free_x", nrow = 1)
p_runtime

# --------------------------------------------------
# Plot 2: y log scaled
# --------------------------------------------------
p_runtime_log <- p_runtime +
  scale_y_log10() +
  labs(
    subtitle = "Logarithmic Y-Axis"
  )
p_runtime_log

# --------------------------------------------------
# Plot 3: Performance (Gflops/s) vs. Matrix Size
# --------------------------------------------------
p_perf <- ggplot(combined_data, aes(
  x = combined_size,
  y = Gflops,
  color = Combined_Color_Factor,              # use c as color (or keep nodes if you want)
  linetype = Group       # ensure algorithm is treated as a factor
)) +
  geom_line(
    aes(group = interaction(c, algorithm)),  # use c for grouping/linetype
    linewidth = 0.5
  ) +
  geom_point(size = 2) +
  scale_color_manual(values = okabe_ito) +
  facet_wrap(~shape_type, scales = "free_x", nrow = 1) +
  labs(
    title = "Performance (GFLOPS/s) vs. Matrix Size",
    subtitle = "Computed as n × m² / runtime",
    x = "Matrix Size (m × n)",
    y = "Performance (GFLOPS/s)",
    color = "",
    shape = "Algorithm",
    linetype = "Algorithm"
  ) +
  guides(
    color = guide_legend(order = 1), # Color legend appears first
    shape = guide_legend(order = 2)  # Shape legend appears second
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1)
  ) 
p_perf

# --------------------------------------------------
# Compute Speedup = T(1 node) / T(k nodes)
# --------------------------------------------------
speedup_data <- combined_data %>%
  group_by(m, n) %>%
  mutate(
    baseline_time = mean(Average_Time[nodes == "1" & algorithm == "2"]),   # runtime for 1 node
    Speedup = baseline_time / Average_Time
  ) %>%
  ungroup() 

# --------------------------------------------------
# Speedup Plot
# --------------------------------------------------
p_speedup <- ggplot(speedup_data, aes(
  x = combined_size,
  y = Speedup,
  color = Combined_Color_Factor,          # use c instead of nodes (discrete)
  linetype = Group,  # ensure algorithm is treated as a factor
)) +
  geom_line(aes(group = interaction(c, algorithm)), linewidth = 0.6) +
  geom_point(size = 2) +
  labs(
    title = "Absolute Speedup vs. Matrix Size",
    subtitle = "Normalized to baseline runtime (A2)",
    x = "Matrix Size (m × n)",
    y = "Speedup",
    color = "",
    shape = "Algorithm",
    linetype = "Algorithm"
  ) +
  guides(
    color = guide_legend(order = 1), # Color legend appears first
    shape = guide_legend(order = 2)  # Shape legend appears second
  ) +
  theme_minimal(base_size = 11) +
  scale_color_manual(values = okabe_ito) +
  facet_wrap(~shape_type, scales = "free_x", nrow = 1) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),
  )

p_speedup


# --------------------------------------------------
# Table output
# --------------------------------------------------
library(dplyr)
library(kableExtra)

# Ensure algorithm is numeric
avg_data <- combined_data %>%
  mutate(algorithm = as.numeric(as.character(algorithm)))

# Filter only A0
a3_data <- combined_data %>% filter(algorithm == 3) 

# Baseline: 1-node runtime of A2 (for absolute speedup)
baseline_a2 <- combined_data %>%
  filter(algorithm == 2) %>%
  select(m, n, A2_1node = Average_Time)

# Baseline: 1-node runtime of A0 (for relative speedup)
baseline_a3 <- a3_data %>%
  filter(nodes == 2) %>%
  select(m, n, A3_2node = Average_Time)

# Join with A0 and compute speedups
speedup_stats <- a3_data %>%
  left_join(baseline_a2, by = c("m", "n")) %>%
  left_join(baseline_a3, by = c("m", "n")) %>%
  mutate(
    Absolute_Speedup = A2_1node / Average_Time,
    Relative_Speedup = A3_2node / Average_Time
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


