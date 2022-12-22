###########################
# Script for plotting the detection of all methods
# 
# created 03.11.2022
# monica.golumbeanu@unibas.ch
##########################

library(readxl)
library(tidyr)
library(ggplot2)
library(tidytext)

# Choose the file you want to plot:
# Read the data for all methods and long fragment in minority clone
# plot_data = read_excel("~/genotyping/test_data/All_Methods_V3_short.xlsx")
# plot_data_all = read_excel("~/genotyping/test_data/All_Methods_V3.xlsx", skip = 1)

# Read the data for all methods and short fragment in minority clone
plot_data = read_excel("~/genotyping/test_data/All_Methods_V2_appendix_short.xlsx")
plot_data_all = read_excel("~/genotyping/test_data/All_Methods_appendix_V2.xlsx", skip = 1)
#############

# Process the data table: remove the unnecessary rows
plot_data_all = as.data.frame(plot_data_all[-c(33:37), ])

# Transform the data frame to longer format to be able to plot
plot_data_long = plot_data %>% 
                  pivot_longer(cols = c("Detected", "Missing", "Not detected"), 
                               names_to = "bar_type", values_to = "value")

# We change the order of the methods so when it's plotted it doesn't order the methods alphabetically
plot_data_long$Method = factor(plot_data_long$Method, levels = c("Fast CE", "HR CE", "AmpSeq"))
# Same for the bar types
plot_data_long$bar_type = factor(plot_data_long$bar_type, levels = c("Not detected", "Missing", "Detected"))

# We define manually the y axes ticks
ratios_labels = rev(plot_data_all[, 2])

# Select only the detected bar
plot_data_long_detected = plot_data_long %>% filter(bar_type == "Detected")
plot_data_long_detected = plot_data_long_detected[order(plot_data_long_detected$value), ]
positions = plot_data_long_detected[order(plot_data_long_detected$value), "Marker"]

# Plotting
ggplot(plot_data_long_detected, aes(x = reorder_within(Marker, value, Method), y = value, fill = bar_type)) +
  geom_bar(stat="identity", show.legend = FALSE) + facet_wrap(~ Method, scales = "free_x") + 
  theme_minimal(base_size = 14) +
  scale_x_reordered() +
  # scale_x_discrete(drop = TRUE) + 
  scale_fill_manual(values = c("#74c476")) + 
  xlab("Marker") + ylab("Ratios") + labs(fill="") +
  scale_y_continuous(breaks = c(1:32),
                   labels = ratios_labels)

# Choose depending on which plot you generate
# ggsave("~/genotyping/plots/all_methods_V3_detected.pdf", width = 16, height = 6)
ggsave("~/genotyping/plots/CE_methods_V2_detected.pdf", width = 12, height = 6)

