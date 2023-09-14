# ABSTRACT
# THis script is used to plot the nonpareil
# curves for the mock and negative controls

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(tibble)
library(Nonpareil)

# Import data ----
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


mock_data <- paste0(
  input_path,
  "/mock_testing/nonpareil/MOCK.npo"
)

neg_data <- paste0(
  input_path,
  "/mock_testing/nonpareil/NEG.npo"
)

both_data <- c(mock_data,
               neg_data)

options(scipen = 999)


# Plot data ----
png(
  filename = paste0(
  output_path,
  "/mock_testing/nonpareil/04_curveplot.png"),
  width = 20, 
  height = 15, 
  units = "cm", 
  res = 600, 
  bg = "white"
  )

Nonpareil.curve.batch(both_data)

dev.off()

p_mock <- Nonpareil.curve(mock_data, plot.dispersion="sd", col="black")
p_neg <- Nonpareil.curve(neg_data, plot.dispersion="sd", col="black")

# Create summaries ----
summary_mock <- as.data.frame(summary.Nonpareil.Curve(p_mock)) %>%
  rownames_to_column("metric")

summary_neg <- as.data.frame(summary.Nonpareil.Curve(p_neg)) %>%
  rownames_to_column("metric")

summary_all <- left_join(summary_mock, summary_neg)

write_delim(
  summary_all,
  paste0(
    output_path,
    "/mock_testing/nonpareil/04_curve_summaries.txt"
  ),
  delim = "\t"
)
