# ABSTRACT
# This script presents the number of reads
# and basic QC of the shotgun metagenomic
# dataset

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(funtools)
library(ggplot2)
library(here)
library(readxl)
library(patchwork)

# 01. Import data ----
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
  sep = "/"
)

## Import data on number of reads
total_reads1 <- read_delim(
  here(
    input_path,
    "R",
    "Data",
    "number_of_reads",
    "runA.txt"
  ),
  delim = "\t",
  col_names = TRUE
)

total_reads2 <- read_delim(
  here(
    input_path,
    "R",
    "Data",
    "number_of_reads",
    "runB.txt"
  ),
  delim = "\t",
  col_names = TRUE
)

all_reads <- rbind(
  total_reads1,
  total_reads2
) %>%
  mutate(sample = sub("(.+)_S.+", "\\1", `Sample Name`),
         read = ifelse(grepl("_R1_",`Sample Name`), "R1", "R2")) %>%
  group_by(sample, read) %>%
  summarise(
    n_million_reads = sum(`M Seqs`)
  ) %>%
  ungroup()


## Import data on number of unclassified reads
unclassified_reads <- read_delim(
  here(
    input_path,
    "R",
    "Data",
    "misc",
    "kraken2_unclassified.txt"
  ),
  delim = "\t",
  col_names = FALSE
) %>%
  mutate(X2 = gsub(" ", "", X2),
         X2 = as.numeric(X2))

colnames(unclassified_reads) <- c(
  "id",
  "perc_fragments",
  "num_fragments_clade",
  "num_fragments_exact",
  "rank_code",
  "tax_id",
  "scientific_name"
  )

# 01. Calculate average number of total reads ----
stats_reads <- all_reads %>%
  filter(!sample %in% c("106-MOCK","107-NEG"),
         read == "R1") %>%
  summarise(
    mean = mean(n_million_reads),
    median = median(n_million_reads),
    sd = sd(n_million_reads),
    range = paste0(min(n_million_reads)," - ", max(n_million_reads))
  )

# 02. Calculate average unclassified reads----
unclassified_stats <- unclassified_reads %>%
  filter(! id %in% c("106-MOCK_pe_db1.bracken.kraken2.report.txt",
                     "107-NEG_pe_db1.bracken.kraken2.report.txt")) %>%
  summarise(
    avg_perc = mean(perc_fragments),
    median_perc = median(perc_fragments),
    sd_perc = sd(perc_fragments),
    avg_reads = mean(num_fragments_exact),
    median_reads = median(num_fragments_exact),
    sd_reads = sd(num_fragments_exact)
  )


