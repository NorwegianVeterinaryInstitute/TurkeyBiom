# ABSTRACT
# This script summarise the sample
# data that we have on the 110 samples

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(funtools)
library(ggplot2)
library(patchwork)
library(here)

# Definitions ----
output_path <- paste(
  "/\\vetinst.no",
  "dfs-felles",
  "StasjonK",
  "FAG",
  "Tverrfaglig",
  "AMR",
  "FoU-aktiviteter & prosjekter",
  "12303_TurkeyBiom",
  "R",
  "Results",
  sep = "/"
)

input_path <- paste(
  "/\\vetinst.no",
  "dfs-felles",
  "StasjonK",
  "FAG",
  "Tverrfaglig",
  "AMR",
  "FoU-aktiviteter & prosjekter",
  "12303_TurkeyBiom",
  "R",
  "Data",
  sep = "/"
)


# Import data ----
data <- read_delim(
  paste(
    input_path,
    "00_metadata_samples.txt",
    sep = "/"),
  delim = "\t") 


metadata <- read_delim(
  here(
    input_path,
    "00_final_treatment_metadata.txt"
  ),
  delim = "\t"
)

# 01. Number of farms and samples per farm

all_farms <- data %>%
  count(eier_lokalitetnummer)

included_farms <- data %>%
  filter(saksnr %in% metadata$saksnr) %>%
  count(eier_lokalitetnummer)

# 01. Sampling date distribution ----
palette <- c("Male" = "#80b1d3",
             "Female" = "#fb8072")

p_samples <- ggplot(data, aes(eier_lokalitetnummer, fill = sex)) +
  geom_bar(color = "black") +
  labs(y = "N") +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

p_dates <- ggplot(data, aes(eier_lokalitetnummer,
                            sampling_date,
                            fill = sex,
                            group = eier_lokalitetnummer)) +
  geom_line() +
  geom_point(pch = 21,
             size = 3) +
  labs(y = "Sampling date") +
  scale_fill_manual(values = palette) +
  scale_y_date(date_breaks = "months",
               date_labels = "%b") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank())

p_density <- data %>%
  filter(!is.na(sex)) %>%
  ggplot(aes(sampling_date, fill = sex)) +
  geom_density(alpha = 0.8) +
  scale_fill_manual(values = palette) +
  scale_x_date(date_breaks = "months",
               date_labels = "%b") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  coord_flip()


p_all <- p_samples + guide_area() + p_dates + p_density +
  plot_layout(guides = "collect",
              heights = c(0.3, 1),
              widths = c(1, 0.3))

## Save plot
ggsave(
  paste(
    output_path,
    "01_sampling_date_distribution.png",
    sep = "/"
  ),
  p_all,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 25,
  width = 25
)
