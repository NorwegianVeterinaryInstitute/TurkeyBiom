# ABSTRACT
# This script generates the input file for taxprofiler
# for the full shotgun dataset

# Libraries
library(tidyverse)
library(here)

# Import filename data ----
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

filenames <- read_delim(
  here(
    input_path,
    "read_filenames.txt"),
  delim = "\t",
  col_names = FALSE
)

# Generate samplesheet ----
samplesheet <- filenames %>%
  mutate(
    sample = sub(".+/(.+_S.+_L00.)_R.+", "\\1", X1),
    run_accession = paste0("run", sub(".+_L00(.)_.+", "\\1", X1))
  ) %>%
  mutate(
    sample = sub("(.+)_S.+", "\\1", sample),
    pair = sub(".+L00._(R.)_.+", "\\1", X1)
  ) %>%
  pivot_wider(names_from = "pair",
              values_from = "X1") %>%
  rename("fastq_1" = R1,
         "fastq_2" = R2) %>%
  mutate(instrument_platform = "ILLUMINA",
         fasta = "") %>%
  select(sample,
         run_accession,
         instrument_platform,
         fastq_1,
         fastq_2,
         fasta) %>%
  arrange(sample, run_accession)

# Write to file ----
write_delim(
  samplesheet,
  here(
    output_path,
    "/06_taxprofiler_samplesheet.txt"
  ),
  delim = ","
)
