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
  "C:",
  "Users",
  "vi1511",
  "temp",
  "data",
  "turkeybiom",
  sep = "/"
)

## Import taxonomic data
tax_data <- import_biom(
  here(
    input_path,
    "kraken_biom",
    "bracken.biom"
  )
)

## Generate metadata
sample_metadata <- data.frame(
  sample = c("106-MOCK_pe_run1_db1.bracken.kraken2.report_bracken_species",
             "107-NEG_pe_run1_db1.bracken.kraken2.report_bracken_species",
             "3-A5_pe_run1_db1.bracken.kraken2.report_bracken_species",
             "75-F3_pe_run1_db1.bracken.kraken2.report_bracken_species"),
  group = c("mock","neg","sample","sample")
)

SAM <- sample_data(
  sample_metadata %>%
    column_to_rownames("sample")
)

tax_data@sam_data <- SAM

# 02. Fix taxonomy ----
## We would like to merge at genus level, so
## we run tax_agg at that rank (Rank6) and 
## check if there are any problems
tax_agg(tax_data, rank = "Rank6")


## Problems were detected, these are fixed with
## tax_fix below
tax_table_fixed <- tax_data %>%
  tax_fix(
  min_length = 4,
  unknowns = c(
    "s__uncultured prokaryote",
    "s__uncultured organism",
    "s__uncultured microorganism",
    "s__synthetic construct",
    "s__unidentified plasmid",
    "s__uncultured marine organism",
    "s__unidentified",
    "s__Cloning vector pMSW107",
    "s__synthetic Caulobacter sp. 'ethensis'",
    "s__uncultured marine microorganism HF4000_APKG2H5",
    "s__unidentified microorganism",
    "s__eukaryotic synthetic construct",
    "s__Cloning vector pMSW105",
    "g__Alsophila",
    "g__Venturia"
    ),
  sep = " ", anon_unique = TRUE,
  suffix_rank = "current"
)

## Now we can aggregate correctly
taxdata_agg <- tax_agg(
  tax_table_fixed,
  rank = "Rank6"
)

## Calculate relative abundances
relabund_data <- transform_sample_counts(
  taxdata_agg, function(x) x*100 / sum(x)
  )

## Split data
plot_data <- psmelt(
  relabund_data
) %>%
  filter(Sample %in% c(
    "106-MOCK_pe_run1_db1.bracken.kraken2.report_bracken_species",
    "107-NEG_pe_run1_db1.bracken.kraken2.report_bracken_species")) %>%
  mutate(Rank6 = sub(".__", "", Rank6)) %>%
  select(Sample, Abundance, contains("Rank")) %>%
  split(f = .$Sample)



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

mock_data <- plot_data$`106-MOCK_pe_run1_db1.bracken.kraken2.report_bracken_species` %>%
  mutate(group = ifelse(
    Rank6 %in% genus_list,
    Rank6, 
    "Other"
  )) %>%
  group_by(group) %>%
  summarise(Abundance = round(sum(Abundance), 5)) %>%
  ungroup()

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
                0.000089)
) %>%
  left_join(mock_data)

write_delim(
  original_mock_composition,
  here(
    output_path,
    "05_shotgun_mock_composition.txt"
  ),
  delim = "\t"
)

## Identify which genera constitute the "Other" 
## category, likely contaminants which should
## be represented in the negative control.
## Aggregate this at a higher level
mock_other_phyla <- plot_data$`106-MOCK_pe_run1_db1.bracken.kraken2.report_bracken_species` %>%
  mutate(group = ifelse(
    Rank6 %in% genus_list,
    Rank6, 
    "Other"
  )) %>%
  filter(group == "Other") %>%
  group_by(Rank2) %>%
  summarise(Abundance_mock = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Rank2 = sub(".__", "", Rank2)) %>%
  mutate(total_perc = sum(Abundance_mock),
         new_abundance = Abundance_mock/total_perc * 100) %>%
  select(-c(Abundance_mock, total_perc)) %>%
  rename("Abundance_mock" = new_abundance) %>%
  arrange(-Abundance_mock)


# 04. Negative control ----
neg_data <- plot_data$`107-NEG_pe_run1_db1.bracken.kraken2.report_bracken_species` %>%
  group_by(Rank2) %>%
  summarise(Abundance_neg = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Rank2 = sub(".__", "", Rank2)) %>%
  arrange(-Abundance_neg)

## Compare the "other" category in the mock to the actual
## phyla identified in the negative control
mock_neg_comp <- mock_other_phyla %>%
  left_join(neg_data) %>%
  filter(Abundance_mock > 0 | Abundance_neg > 0) %>%
  mutate(group = ifelse(Abundance_mock > 0.1 | Abundance_neg > 0.1, Rank2, "Other")) %>%
  group_by(group) %>%
  summarise(Abundance_mock = sum(Abundance_mock),
            Abundance_neg = sum(Abundance_neg)) %>%
  ungroup() %>%
  pivot_longer(
    cols = -group,
    names_to = "type",
    values_to = "abundance"
  ) %>%
  mutate(type = sub("Abundance_", "", type))

p1 <- ggplot(mock_neg_comp, 
             aes(reorder(group, -abundance), abundance, color = type)) +
  geom_point(size = 3) +
  labs(x = "Phyla",
       y = "Relative abundance",
       color = "Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))

p2 <- ggplot(mock_neg_comp, aes(type, abundance, fill = group)) +
  geom_col(color = "white") +
  labs(y = "Relative abundance",
       fill = "Phyla") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank())


p3 <- ( p2 + guide_area() ) / p1

ggsave(
  here(
    output_path,
    "05_relative_abundance_contaminants.png"
  ),
  p3,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 20,
  width = 25
)

## Investigate major groups within the negative control ----
neg_data_major <- plot_data$`107-NEG_pe_run1_db1.bracken.kraken2.report_bracken_species` %>%
  filter(Rank2 %in% c("p__Firmicutes",
                      "p__Arthropoda",
                      "p__Bacteroidota")) %>%
  group_by(Rank1,Rank3) %>%
  summarise(Abundance_neg = sum(Abundance)) %>%
  mutate_at(vars(contains("Rank")),
            ~sub(".__", "", .))

mock_data_major <- plot_data$`106-MOCK_pe_run1_db1.bracken.kraken2.report_bracken_species` %>%
  filter(Rank2 %in% c("p__Firmicutes",
                      "p__Arthropoda",
                      "p__Bacteroidota")) %>%
  group_by(Rank1,Rank3) %>%
  summarise(Abundance_mock = sum(Abundance)) %>%
  mutate_at(vars(contains("Rank")),
            ~sub(".__", "", .))

major_data_all <- neg_data_major %>%
  left_join(mock_data_major) %>%
  pivot_longer(
    cols = c("Abundance_neg",
             "Abundance_mock"),
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  mutate(sample = sub("Abundance_", "", sample))


p4 <- ggplot(major_data_all, aes(Rank3, abundance, color = sample)) +
  geom_point(size = 3) +
  facet_grid(~Rank1, scales = "free_x") +
  labs(x = NULL,
       y = "Relative abundance",
       color = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.3))


ggsave(
  here(
    output_path,
    "05_relative_abundance_contaminants_major.png"
  ),
  p4,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 15,
  width = 20
)
