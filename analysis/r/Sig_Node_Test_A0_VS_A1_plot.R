library(kableExtra)
library(dplyr)
library(purrr)
library(ggplot2)

# ------------------------
# Description:
# This R script analyzes the runtime performance of the SYRK implementations A0 and A1 
# across different node counts and matrix sizes. It performs statistical tests to assess:
#   1. Whether the mean runtimes of A0 and A1 differ significantly (using a two-sample t-test),
#   2. Whether the runtime ratio (A0/A1) is approximately equal to 2 within a tolerance of ±0.1,
#      tested via confidence intervals and t-tests against the bounds 1.9 and 2.1.
#
# The script also evaluates the normality of the runtime distributions using 
# the Shapiro–Wilk test and generates two LaTeX tables summarizing:
#   - Mean runtimes, ratios, and significance per matrix and node size,
#   - Test results for whether the runtime ratio differs significantly from 2.
#
#
# Last Modified: 2025-10-07
# --------------------------
# 1. Read all CSV files
# --------------------------
folder <- "/Users/johannes/Downloads/Uni/Bachelor Arbeit/2025-08-12/multi_node"
files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)

df <- files %>% 
  lapply(read.csv, sep = ";", header = TRUE, stringsAsFactors = FALSE) %>% 
  bind_rows()

# --------------------------
# 2. Compute stats
# --------------------------
safe_shapiro_p <- function(x) {
  if (length(x) >= 3) shapiro.test(x)$p.value else NA
}

safe_t_test_p <- function(x, y) {
  if (length(x) >= 2 && length(y) >= 2) t.test(x, y, var.equal = FALSE)$p.value else NA
}

results <- df %>%
  group_by(m, n, nodes) %>%   # <-- nodes added here
  summarise(
    mean_alg0_ms = mean(runtime_ms[algorithm == 0]),
    mean_alg1_ms = mean(runtime_ms[algorithm == 1]),
    ratio = mean(runtime_ms[algorithm == 0]) / mean(runtime_ms[algorithm == 1]),
    shapiro_p_alg0 = safe_shapiro_p(runtime_ms[algorithm == 0]),
    shapiro_p_alg1 = safe_shapiro_p(runtime_ms[algorithm == 1]),
    p_value = safe_t_test_p(runtime_ms[algorithm == 0], runtime_ms[algorithm == 1]),
    .groups = "drop"
  ) %>%
  mutate(Significant = ifelse(p_value < 0.05, "YES", "NO"))
View(results)
results$p_value <- format(results$p_value, scientific = TRUE, digits = 3)

colnames(results) <- c(
  "M", "N", "Nodes", "Mean Alg0 (ms)", "Mean Alg1 (ms)", 
  "Ratio", "Shapiro p Alg0", "Shapiro p Alg1", "p-value", "Significant"
)

results %>%
  kable(format = "latex", booktabs = TRUE, 
        caption = "Performance comparison between Algorithm 0 and Algorithm 1 for each node count") %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))

# --------------------------
# 4. Test if ratio differs from 2
# --------------------------
lower_bound <- 1.9
upper_bound <- 2.1

sig_range_2_results <- df %>%
  group_by(m, n, nodes) %>%  # <-- nodes added here
  group_split() %>%
  map_dfr(function(group_key) {
    m_val <- unique(group_key$m)
    n_val <- unique(group_key$n)
    nodes_val <- unique(group_key$nodes)
    
    runtimes_alg0 <- filter(group_key, algorithm == 0) %>% pull(runtime_ms)
    runtimes_alg1 <- filter(group_key, algorithm == 1) %>% pull(runtime_ms)
    ratios <- runtimes_alg0 / runtimes_alg1
    
    mean_ratio <- mean(ratios, na.rm = TRUE)
    sd_ratio <- sd(ratios, na.rm = TRUE)
    n <- length(ratios)
    
    error_margin <- qt(0.975, df = n - 1) * sd_ratio / sqrt(n)
    ci_lower <- mean_ratio - error_margin
    ci_upper <- mean_ratio + error_margin
    
    p_lower <- t.test(ratios, mu = lower_bound, alternative = "greater")$p.value
    p_upper <- t.test(ratios, mu = upper_bound, alternative = "less")$p.value
    
    data.frame(
      m = m_val,
      n = n_val,
      nodes = nodes_val,
      mean_ratio = mean_ratio,
      sd_ratio = sd_ratio,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value_lower_vs_1.9 = format(p_lower, scientific = TRUE, digits = 3),
      p_value_upper_vs_2.1 = format(p_upper, scientific = TRUE, digits = 3),
      Significant = ifelse(p_lower < 0.05 && p_upper < 0.05, "YES", "NO")
    )
  })

colnames(sig_range_2_results) <- c(
  "Matrix m", "Matrix n", "Nodes", "Mean Ratio", "SD Ratio", 
  "CI Lower", "CI Upper", "P-Value Lower vs 1.9", "P-Value Upper vs 2.1", "Significant"
)

sig_range_2_results %>%
  kable(format = "latex", booktabs = TRUE, 
        caption = "Test if ratio differs from 2 for each node count") %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))


# --------------------------
# 5. Plot: Highlight "thin", "square", and "fat" shapes separately
# --------------------------

# Define shape categories (adjust thresholds as needed)
sig_range_2_results <- sig_range_2_results %>%
  mutate(
    shape_class = case_when(
      `Matrix m` < `Matrix n` ~ "thin",
      `Matrix m` == `Matrix n` ~ "square",
      `Matrix m` > `Matrix n` ~ "fat"
    )
  )

# Ensure proper numeric typing
sig_range_2_results$`Mean Ratio` <- as.numeric(sig_range_2_results$`Mean Ratio`)
sig_range_2_results$`CI Lower` <- as.numeric(sig_range_2_results$`CI Lower`)
sig_range_2_results$`CI Upper` <- as.numeric(sig_range_2_results$`CI Upper`)
sig_range_2_results$Nodes <- as.factor(sig_range_2_results$Nodes)
sig_range_2_results$`Matrix m` <- as.factor(sig_range_2_results$`Matrix m`)
sig_range_2_results$`Matrix n` <- as.factor(sig_range_2_results$`Matrix n`)

# Helper function for one plot
plot_shape_class <- function(data, highlight_class) {
  data <- data %>%
    mutate(highlight = ifelse(shape_class == highlight_class, "highlight", "other"))
  
  ggplot(data, aes(x = Nodes, y = `Mean Ratio`, group = interaction(`Matrix n`, shape_class))) +
    # background (gray) lines
    geom_point(data = subset(data, highlight == "other"),
               aes(group = `Matrix n`), color = "gray80", size = 2, alpha = 0.5) +
    
    # highlight lines
    geom_errorbar(
      data = subset(data, highlight == "highlight"),
      aes(ymin = `CI Lower`, ymax = `CI Upper`, color = `Matrix n`),
      width = 0.25, size = 0.8, alpha = 0.8
    ) +
    geom_point(
      data = subset(data, highlight == "highlight"),
      aes(color = `Matrix n`),
      size = 3, alpha = 0.9
    ) +
    
    # reference lines
    geom_hline(yintercept = 2, linetype = "dashed", color = "red", linewidth = 1.0) +
    geom_hline(yintercept = 1.9, linetype = "dotted", color = "gray60", linewidth = 0.8) +
    geom_hline(yintercept = 2.1, linetype = "dotted", color = "gray60", linewidth = 0.8) +
    
    # axes and styling
    #scale_y_continuous(limits = c(1.5, 2.5), breaks = seq(1.5, 2.5, 0.1)) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      title = paste("Runtime Ratio (A0/A1) — ", toupper(highlight_class)),
      subtitle = "Mean ratios with 95% CI; other shapes in gray",
      x = "Node Count",
      y = "Mean Runtime Ratio (A0/A1)",
      color = "Matrix n"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray30"),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 0, vjust = 0.5)
    )
}

# Create the three plots
plot_thin <- plot_shape_class(sig_range_2_results, "thin")
plot_square <- plot_shape_class(sig_range_2_results, "square")
plot_fat <- plot_shape_class(sig_range_2_results, "fat")

# Display them (optionally use patchwork to combine)
#library(patchwork)
plot_thin
plot_square
plot_fat
