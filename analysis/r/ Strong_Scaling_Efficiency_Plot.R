library(dplyr)
library(tidyverse)
library(readr) # Specifically for read_csv
library(stringr)
library(ggplot2)
library(kableExtra)


# ------------------------------------------------------------
# REQUIREMENTS:
# avg_data must contain at minimum:
#   nodes, algorithm, m, n, Average_Time, shape_type
# nodes must be numeric or coercible to numeric
# shape_type must have values: Square, Thin, Fat
# ------------------------------------------------------------

# ----------------------------------------------------------------------
# PART 1: Data Ingestion and Combination (Creating the 'df' table)
# ----------------------------------------------------------------------

# --- 1. Define all folders and paths ---

# A0 data paths
# maybe need to set the working directory using setwd(...)
base_dir <- "../../experiments/raw"
folders_a0 <- c(
  file.path(base_dir, "A0_VS_A2", "A0_scaling_1")
)

# A1 data paths
folders_a1 <- file.path(
  base_dir,
  "A1_VS_A2",
  "2025-08-26",
  paste0("A1_scaling_", 1:3)
)

# A2 data paths
folders_a2 <- c(
  file.path(base_dir, "A0_VS_A2", "2025-08-22"),
  file.path(
    base_dir,
    "A2_Scaling",
    "2025-10-09",
    paste0("A2_scaling_", 1:3)
  )
)

# A3 data paths
folders_a3 <- file.path(
  base_dir,
  "A3_Scaling",
  "2025-10-12",
  paste0("A3_scaling_", 0:2)
)

# A4 data paths (constructed programmatically)
base_dir_a4 <- file.path(base_dir, "A4_Scaling", "2025-10-15")
subfolders_a4 <- sprintf(
  "A4_scaling_%s",
  c("02", "1_2", "1_3", "1_4", "1_5", "1_6",
    "2_2", "2_3", "2_4", "2_5", "2_6",
    "3_2", "3_3")
)
folders_a4 <- file.path(base_dir_a4, subfolders_a4)

# Combine all unique folder paths into one vector
all_folders <- unique(c(folders_a0, folders_a1, folders_a2, folders_a3, folders_a4))


# --- 2. Find all data files (assuming .csv) ---

# Find all CSV files within the listed folders (recursive=TRUE looks in subfolders too)
all_files <- unlist(sapply(all_folders, function(folder) {
  list.files(
    path = folder,
    pattern = "\\.csv$", # Looks for files ending with .csv
    full.names = TRUE,
    recursive = TRUE
  )
}))

# Check if any files were found
if (length(all_files) == 0) {
  stop("No CSV files found in the specified directories. Please check your paths and file extensions.")
}


# --- 3. Read and combine all data files ---

# Use purrr::map_dfr to iterate through the file paths, read each CSV, 
# and combine them by row into a single data frame (df).
# Since files use ';' as a delimiter, we use read_delim().
# map_dfr automatically handles differing columns (like 'c' and 'i') by 
# filling the missing cells with NA, which is exactly what we need.
df <- map_dfr(all_files, ~read_delim(.x, delim = ";", show_col_types = FALSE), 
              .id = "source_file_path", 
              .progress = TRUE)


# ----------------------------------------------------------------------
# PART 2: Data Analysis Pipeline (Your original code, adapted for 'df')
# ----------------------------------------------------------------------

# Now apply your summary and mutation steps to the combined 'df'
avg_data <- df %>%
  # Optional: Add an experiment ID based on the folder structure for easier filtering
  #mutate(experiment_id = str_extract(source_file_path, "A\\d+_scaling_\\d.*")) %>%
  
  # group by experiment parameters
  #group_by(experiment_id, nodes, algorithm, m, n) %>%
  group_by(nodes, algorithm, m, n) %>%
  
  # compute summary statistics
  summarise(
    Average_Time = mean(runtime_ms, na.rm = TRUE),
    StdDev_Time  = sd(runtime_ms, na.rm = TRUE),
    Runs         = n(),    # number of runs
    .groups = 'drop'
  ) %>%
  
  # Data type conversion and feature engineering
  mutate(
    nodes = as.factor(nodes),
    algorithm = as.factor(algorithm)
  ) %>%
  mutate(
    shape_type = case_when(
      m < n  ~ "Fat",
      m == n ~ "Square",
      m > n  ~ "Thin"
    ),
    size_metric = m * n, # Total matrix size for ordering
    combined_size = paste0(m, "x", n)
  ) %>%
  
  # Ordering and factoring 'combined_size' for plotting
  arrange(shape_type, size_metric) %>%
  mutate(
    combined_size = factor(combined_size, levels = unique(combined_size)),
    
    # --- NEW: Compute Gflops per second ---
    # Assuming the matrix operation is O(n*m^2) for simplicity, 
    # and converting runtime_ms to seconds (1e6)
    Gflops = (n * m^2) / (1e6 * Average_Time) 
  )


# Ensure proper types
avg_data <- avg_data %>%
  mutate(
    nodes = as.numeric(as.character(nodes)),
    algorithm = as.factor(algorithm),
    combined_size = paste0(m, "x", n),
    shape_type = factor(shape_type, levels = c("Fat", "Square", "Thin"))
  )

# ------------------------------------------------------------
# Compute Strong Scaling Efficiency per (algorithm, m, n, nodes)
# ------------------------------------------------------------
sse_data <- avg_data %>%
  group_by(algorithm, m, n) %>%
  mutate(
    base_nodes = min(nodes),                     # smallest available node count
    T_base = mean(Average_Time[nodes == base_nodes]),  # baseline runtime
    Strong_Eff = T_base / ((nodes / base_nodes ) * Average_Time)
  ) %>%
  ungroup()

# ------------------------------------------------------------
# Summary table by shape type & algorithm
# ------------------------------------------------------------
sse_summary <- sse_data %>%
  group_by(shape_type, algorithm) %>%
  summarise(
    Min_Efficiency = min(Strong_Eff, na.rm = TRUE),
    Max_Efficiency = max(Strong_Eff, na.rm = TRUE),
    Mean_Efficiency = mean(Strong_Eff, na.rm = TRUE),
    SD_Efficiency = sd(Strong_Eff, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# LaTeX table output
# ------------------------------------------------------------
sse_summary %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    digits = 4,
    col.names = c(
      "Shape", "Algorithm",
      "Min Eff", "Max Eff", "Mean Eff", "SD Eff"
    ),
    caption = "Strong Scaling Efficiency per Shape Category and Algorithm"
  ) %>%
  kable_styling(latex_options = c("striped", "hold_position", "scale_down"))

# ------------------------------------------------------------
# Plot: Strong Scaling Efficiency vs Matrix Size
# ------------------------------------------------------------
# 1. Determine the correct order based on m and then n
#    We use distinct() to get unique size combinations and then arrange them.
size_order <- sse_data %>%
  distinct(m, n, combined_size) %>%
  arrange(m, n) %>%
  pull(combined_size)

# 2. Convert 'combined_size' to a factor using the determined order
sse_data$combined_size <- factor(
  sse_data$combined_size, 
  levels = size_order
)

sse_data_filtered <- sse_data %>% filter(algorithm != 0)

p_sse <- ggplot(sse_data_filtered, aes(
  x = combined_size,
  y = Strong_Eff,
  color = factor(nodes),
  shape = algorithm
)) +
  geom_line(aes(group = interaction(nodes, algorithm)), linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~shape_type, scales = "free_x", nrow = 1) +
  theme_minimal(base_size = 10) +
  labs(
    title = "Strong Scaling Efficiency vs. Matrix Size",
    subtitle = "Efficiency = T(1) / (k × T(k))",
    x = "Matrix Size (m × n)",
    y = "Strong Scaling Efficiency",
    color = "Nodes",
    shape = "Algorithm"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1)
  )

p_sse



################################################################################

# Define the algorithms you want to highlight
algorithms_to_highlight <- c(2)

sse_data_highlight <- sse_data %>%
  # Apply Factor Ordering to the base data frame
  mutate(combined_size = factor(combined_size, levels = size_order)) %>%
  # Create the highlight group
  mutate(
    highlight_group = ifelse(
      algorithm %in% algorithms_to_highlight, 
      "Highlighted", 
      "Other"
    ),
    highlight_group = factor(highlight_group, levels = c("Highlighted", "Other"))
  ) %>%
  filter(algorithm != 0) # Apply the required filter


p_sse_highlighted <- ggplot(sse_data_highlight, aes(
  x = combined_size,
  y = Strong_Eff,
  shape = algorithm
)) +
  # 1. Plot the "Other" algorithms (gray and faded)
  geom_line(
    # 👇 CORRECTED: Use sse_data_highlight
    data = subset(sse_data_highlight, highlight_group == "Other"), 
    aes(group = interaction(nodes, algorithm)), 
    color = "gray60", 
    alpha = 0.3, 
    linewidth = 0.5
  ) +
  geom_point(
    # 👇 CORRECTED: Use sse_data_highlight
    data = subset(sse_data_highlight, highlight_group == "Other"), 
    color = "gray60", 
    alpha = 0.3, 
    size = 2
  ) +
  
  # 2. Plot the "Highlighted" algorithms (full color and opacity)
  geom_line(
    # 👇 CORRECTED: Use sse_data_highlight
    data = subset(sse_data_highlight, highlight_group == "Highlighted"), 
    aes(color = factor(nodes), group = interaction(nodes, algorithm)), 
    linewidth = 0.5
  ) +
  geom_point(
    # 👇 CORRECTED: Use sse_data_highlight
    data = subset(sse_data_highlight, highlight_group == "Highlighted"), 
    aes(color = factor(nodes)), 
    size = 2
  ) +
  
  # --- Standard Plot Elements ---
  facet_wrap(~shape_type, scales = "free_x", nrow = 1) +
  theme_minimal(base_size = 10) +
  labs(
    title = "Strong Scaling Efficiency vs. Matrix Size (Algorithms Highlighted)",
    subtitle = paste("Highlighted:", paste(algorithms_to_highlight, collapse = ", ")),
    x = "Matrix Size (m × n)",
    y = "Strong Scaling Efficiency",
    color = "Nodes",
    shape = "Algorithm"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1)
  )

p_sse_highlighted

