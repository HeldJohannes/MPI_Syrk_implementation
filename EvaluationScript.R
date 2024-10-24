################################################################################
# Setup:
################################################################################

# add necessary libraries:
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

# set the directory to read form and the filename
dir <- "/Users/johannes/CLionProjects/MPI_Syrk_implentation/"
filename <- "results.csv"

# read the result data CSV
measurment_data <- read.csv(paste(dir, filename, sep = ""), sep = ";")

################################################################################
# Evaluation:
################################################################################

# group the date and summarise the runtime over all runs and 
# unite the matrix size to a string for better readability
avg_data <- measurment_data %>%
  arrange(algorithm, number_of_processors, matrix_size_m, matrix_size_n) %>%
  group_by(algorithm, number_of_processors, matrix_size_m, matrix_size_n) %>%
  summarise(avg_runtime = mean(run_time)) %>%
  unite(matrix_shape, c(matrix_size_m, matrix_size_n), sep = "x")

# plot the data using ggplot:
avg_data %>%
  ggplot(aes(x = matrix_shape, y = avg_runtime, colour = factor(number_of_processors))) +
  geom_point() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


# filter for matrix 10x... 
avg_data_2 <- measurment_data %>%
  filter(algorithm == 1) %>%
  arrange(algorithm, number_of_processors, matrix_size_m, matrix_size_n) %>%
  group_by(algorithm, number_of_processors, matrix_size_m, matrix_size_n) %>%
  summarise(avg_runtime = mean(run_time)) %>%
  unite(matrix_shape, c(matrix_size_m, matrix_size_n), sep = "x")

# plot filtered data:
plot1 <- avg_data_2 %>%
  ggplot(aes(x = matrix_shape, y = avg_runtime, colour = factor(number_of_processors))) +
  geom_point() + 
  ylim(0, 10) +
  scale_color_discrete(name = "number of processors", labels = c("1 Processor", "2 Processors", "3 Processors", "4 Processors", "5 Processors", "6 Processors")) +
  scale_y_continuous(trans='log10') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(title = "Rutime plot", subtitle =  "Triple for loop algorith", x = "Matrix", y = "Average Runtime [s]")
plot1

# filter for matrix 10x... 
avg_data_3 <- measurment_data %>%
  arrange(algorithm, number_of_processors, matrix_size_m, matrix_size_n) %>%
  group_by(algorithm, number_of_processors, matrix_size_m, matrix_size_n) %>%
  summarise(avg_runtime = mean(run_time)) %>%
  unite(matrix_shape, c(matrix_size_m, matrix_size_n), sep = "x")

# plot filtered data:
plot2 <- avg_data_3 %>%
  ggplot(aes(x = matrix_shape, y = avg_runtime, colour = factor(algorithm))) +
  geom_point() + 
  scale_y_continuous(trans = 'log10')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  labs(title = "Rutime plot", subtitle =  "Improved Triple for loop algorith", x = "Matrix", y = "Average Runtime [s]")
plot2

plot_grid(plot1, plot2, labels = "AUTO")

tmp <- avg_data_3 %>% group_by(number_of_processors, matrix_shape) %>% summarise(fact = first(avg_runtime, order_by = algorithm) / last(avg_runtime, order_by = algorithm))
tmp %>% ggplot(aes(x = matrix_shape, y = fact, colour = factor(number_of_processors))) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
