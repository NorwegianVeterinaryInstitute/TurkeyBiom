# ABSTRACT
# This script looks into the taxa profiles for
# the shotgun mock and negative control

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(funtools)
library(ggplot2)
library(phyloseq)
library(vegan)
library(microbial)
library(here)
library(tibble)
library(patchwork)
library(microViz)

# 01. Import data ----
options(scipen = 999)

## Define paths
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

## Import taxonomic data
tax_data <- import_biom(
  here(
    input_path,
    "kraken_biom",
    "bracken_controls.biom"
  )
) %>%
  subset_taxa(Rank7 != "s__sapiens")

## Fix sample names
colnames(tax_data@otu_table) <- c("MOCK","NEG")


## Generate metadata
sample_metadata <- data.frame(
  sample = c("MOCK",
             "NEG"),
  group = c("ZymoBiomics Mock","Negative control")
)

SAM <- sample_data(
  sample_metadata %>%
    column_to_rownames("sample")
)

tax_data@sam_data <- SAM

# 02. Fix taxonomy ----
tax_table_fixed <- tax_data %>%
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

# Number of reads
reads_data <- as.data.frame(sample_sums(tax_table_fixed)) %>%
  rownames_to_column("sample") %>%
  rename("reads" = 2)

p_n_reads <- ggplot(reads_data, aes(sample, reads)) +
  geom_col(color = "black") +
  labs(y = "Number of reads") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank())

# Abundance
### Phylum level
p_phyla <- comp_barplot(
  tax_table_fixed,
  tax_level = "Rank2",
  sample_order = "bray",
  n_taxa = 8
) +
  labs(y = "Relative abundance (%)",
       fill = "Phyla") +
  theme(axis.ticks.x = element_blank())

### Genus level
genera_palette <- c(
  RColorBrewer::brewer.pal(
    8,
    "Set3"
    ),
  "grey80" 
  )

names(genera_palette) <-
  c(
    "Listeria",
    "Pseudomonas",
    "Saccharomyces",
    "Bacillus",
    "Faecalibacterium",
    "Bifidobacterium",
    "Alistipes",
    "Solanum",
    "Other"
  )

p_genera <- comp_barplot(
  tax_table_fixed,
  tax_level = "Rank6",
  sample_order = "bray",
  n_taxa = 8
) +
  labs(y = "Relative abundance (%)",
       fill = "Genera") +
  scale_fill_manual(values = genera_palette) +
  theme(axis.ticks.x = element_blank())


p_relabund_mock_neg <- p_n_reads + p_phyla + p_genera +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A")


ggsave(
  here(
    output_path,
    "05_relative_abundance_mock_neg.png"
  ),
  p_relabund_mock_neg,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 12,
  width = 25
)


# 03. Mock data ----

genus_list <- c(
  "Listeria",
  "Pseudomonas",
  "Bacillus",
  "Saccharomyces",
  "Escherichia",
  "Salmonella",
  "Lactobacillus",
  "Enterococcus",
  "Cryptococcus",
  "Staphylococcus"
)

original_mock_composition <- data.frame(
  group = c("Listeria",
            "Pseudomonas",
            "Bacillus",
            "Saccharomyces",
            "Escherichia",
            "Salmonella",
            "Lactobacillus",
            "Enterococcus",
            "Cryptococcus",
            "Staphylococcus"),
  Abundance_original = c(89.1,
                         8.9,
                         0.89,
                         0.89,
                         0.089,
                         0.089,
                         0.0089,
                         0.00089,
                         0.00089,
                         0.000089))

mock_data <- plot_data$`106-MOCK_pe_db1.bracken.kraken2.report_bracken_species` %>%
  mutate(group = ifelse(
    Rank6 %in% genus_list,
    Rank6, 
    "Other"
  )) %>%
  group_by(group) %>%
  summarise(Abundance = round(sum(Abundance), 5)) %>%
  ungroup() %>%
  left_join(original_mock_composition) %>%
  arrange(-Abundance_original)

write_delim(
  mock_data,
  here(
    output_path,
    "05_shotgun_mock_composition.txt"
  ),
  delim = "\t"
)


