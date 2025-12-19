library(kableExtra)
library(dplyr)
library(purrr)
library(kableExtra)


# --------------------------
# 1. Read all CSV files
# --------------------------
base_dir <- "../../experiments/raw"
folder <- file.path(base_dir, "2025-08-12")
files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)

df <- files %>% 
  lapply(read.csv, sep = ";", header = TRUE, stringsAsFactors = FALSE) %>% 
  bind_rows()

# --------------------------
# 2. Compute stats
# --------------------------

safe_shapiro_p <- function(x) {
  if(length(x) >= 3) shapiro.test(x)$p.value else NA
}

safe_t_test_p <- function(x, y) {
  if(length(x) >= 2 && length(y) >= 2) t.test(x, y, var.equal = FALSE)$p.value else NA
}


results <- df %>%
  group_by(m, n) %>%
  summarise(
    mean_alg0_ms = mean(runtime_ms[algorithm == 0]),
    mean_alg1_ms = mean(runtime_ms[algorithm == 1]),
    ratio = mean(runtime_ms[algorithm == 0]) / mean(runtime_ms[algorithm == 1]),
    # Shapiro-Wilk normality test p-values
    shapiro_p_alg0 = safe_shapiro_p(runtime_ms[algorithm == 0]),
    shapiro_p_alg1 = safe_shapiro_p(runtime_ms[algorithm == 1]),
    
    # t-test p-value (unequal variance)
    p_value = safe_t_test_p(runtime_ms[algorithm == 0], runtime_ms[algorithm == 1]),
  ) %>%
  ungroup() %>%
  mutate(Significant = ifelse(p_value < 0.05, "YES", "NO"))
View(results)



# --------------------------
# 3. Output LaTeX table
# --------------------------
results$p_value <- format(results$p_value, scientific = TRUE, digits = 3)
colnames(results) <- c("M", "N", "Mean Alg0 (ms)", "Mean Alg1 (ms)", 
                       "Ratio", "Shapiro p Alg0", "Shapiro p Alg1", "p-value", "Significant")
results %>%
  kable(format = "latex", booktabs = TRUE, 
        caption = "Performance comparison between Algorithm 0 and Algorithm 1") %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))

# ----
# 4. Test auf signifikant von 2 verschieden
# ----

lower_bound <- 1.9
upper_bound <- 2.1

sig_range_2_results <- df %>%
  group_by(m, n) %>%
  group_split() %>%
  map_dfr(function(group_key) {
    m_val <- unique(group_key$m)
    n_val <- unique(group_key$n)
    
    runtimes_alg0 <- filter(group_key, algorithm == 0) %>% pull(runtime_ms)
    runtimes_alg1 <- filter(group_key, algorithm == 1) %>% pull(runtime_ms)
    # Elementwise ratio
    ratios <- runtimes_alg0 / runtimes_alg1
    
    # Mean and SD
    mean_ratio <- mean(ratios, na.rm = TRUE)
    sd_ratio <- sd(ratios, na.rm = TRUE)
    n <- length(ratios)
    
    # 95% confidence interval for mean ratio using t-distribution
    error_margin <- qt(0.975, df = n - 1) * sd_ratio / sqrt(n)
    ci_lower <- mean_ratio - error_margin
    ci_upper <- mean_ratio + error_margin
    
    p_lower <- t.test(ratios, mu = lower_bound, alternative = "greater")$p.value
    p_upper <- t.test(ratios, mu = upper_bound, alternative = "less")$p.value
    
    data.frame(
      m = m_val,
      n = n_val,
      mean_ratio = mean_ratio,
      sd_ratio = sd_ratio,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value_lower_vs_1.9 = format(p_lower, scientific = TRUE, digits = 3),
      p_value_upper_vs_2.1 = format(p_upper, scientific = TRUE, digits = 3),
      Significant = ifelse(p_lower < 0.05 && p_upper < 0.05, "YES", "NO")
    )
  })
View(sig_range_2_results)

colnames(sig_range_2_results) <- c(
  "Matrix m", "Matrix n", "Mean Ratio", "SD Ratio", 
  "CI Lower", "CI Upper", "P-Value Lower vs 1.9", "P-Value Upper vs 2.1", "Significant"
)
sig_range_2_results %>%
  kable(format = "latex", booktabs = TRUE, 
        caption = "Test if ratio differs from 2") %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))
