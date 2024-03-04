# ABSTRACT
# This script is used to generate the
# clean phyloseq object for further
# analyses

# Libraries ----
library(dplyr)
library(readr)
library(phyloseq)
library(microViz)
library(here)
library(tibble)
library(lubridate)

# 01. Import data ----
## Define paths ----
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

## Import taxonomic assignments ----
bracken_data <- import_biom(
  here(
    input_path,
    "kraken_biom",
    "bracken_samples.biom"
  )
) %>%
  ### Remove remaining human annotations
  subset_taxa(Rank7 != "sapiens")

# 02. Add metadata to the phyloseq object ----
sample_data <- as.data.frame(bracken_data@otu_table) %>%
  t() %>%
  as.data.frame %>%
  rownames_to_column("id") %>%
  select(id)

metadata <- read_delim(
  here(
    input_path,
    "00_final_treatment_metadata.txt"
  ),
  delim = "\t"
) %>%
  mutate(sample_month = factor(month(sampling_date)))


seq_ids <- sample_data$id
metadata_ids <- metadata$seq_id

matches_seq_ids <- c()
matches_metadata_ids <- c()

for (i in seq_ids) {
  for (j in metadata_ids) {
    if (grepl(j, i)) {
      matches_seq_ids <- c(matches_seq_ids, i)
      matches_metadata_ids <- c(matches_metadata_ids, j)
    }
  }
}

no_of_reads <- as.data.frame(sample_sums(bracken_data)) %>%
  rename("reads" = `sample_sums(bracken_data)`) %>%
  rownames_to_column("id")

sample_df <- data.frame(id = matches_seq_ids,
                        seq_id = matches_metadata_ids) %>%
  left_join(metadata[, c(
    "seq_id",
    "sex",
    "monensin",
    "antibiotic",
    "antibiotic_treatments",
    "age_days",
    "sample_month"
  )]) %>%
  select(-seq_id) %>%
  mutate(sex_binary = ifelse(sex == "Female", 1, 0),
         monensin_antibiotic = ifelse(monensin == TRUE & 
                                        antibiotic == TRUE, 
                                      TRUE, FALSE),
         monensin_plot = monensin,
         antibiotic_plot = antibiotic) %>%
  left_join(no_of_reads) %>%
  column_to_rownames("id")

bracken_data@sam_data <- sample_data(sample_df)

# 03. Clean the taxonomic assignments ----
bracken_fixed <- bracken_data %>%
  tax_mutate(Rank1 = sub("k__", "", Rank1),
             Rank2 = sub("p__", "", Rank2),
             Rank3 = sub("c__", "", Rank3),
             Rank4 = sub("o__", "", Rank4),
             Rank5 = sub("f__", "", Rank5),
             Rank6 = sub("g__", "", Rank6),
             Rank7 = sub("s__", "", Rank7)) %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "current"
  )


# 04. Rarefy counts ----
set.seed(12345)
bracken_rarefied <- rarefy_even_depth(
  bracken_fixed
)

# 05. Identify lost taxa after rarefaction ----
old <- psmelt(bracken_fixed)
new <- psmelt(bracken_rarefied)

check_otus <- old %>%
  filter(!OTU %in% new$OTU) %>%
  select(OTU, contains("Rank")) %>%
  group_by(OTU) %>%
  distinct()

check_abundance <- old %>%
  filter(!OTU %in% new$OTU) %>%
  select(OTU, Abundance) %>%
  summarise(
    abundance = mean(Abundance),
    min = min(Abundance),
    max = max(Abundance)
  )

# 06. Save RDS ----
write_rds(
  bracken_fixed,
  here(
    output_path,
    "bracken_fixed.rds"
  )
)

write_rds(
  bracken_rarefied,
  here(
    output_path,
    "bracken_rarefied.rds"
  )
)

