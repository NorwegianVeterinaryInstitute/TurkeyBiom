# ABSTRACT
# This script calculates the relative abundance and
# alpha/beta diversity of the 16S mock and test samples
# to compare the two DNA extraction kits

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

# 01. Import data ----
options(scipen = 999)

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
  "16S",
  sep = "/"
)

## Import metadata
metadata <- read_delim(
  here(
    input_path,
    "16S_metadata.txt"
  ),
  delim = ","
)


# 02.  Mock analysis ----
## Import ASV data
asv_names <- read_delim(
  here(
    input_path,
    "mock_analysis",
    "qiime2",
    "rel_abundance_tables",
    "rel-table-ASV_with-DADA2-tax.tsv"
  ),
  delim = "\t",
  col_names = TRUE
) %>%
  select(1:8) %>%
  rename("OTU" = ID) %>%
  column_to_rownames("OTU")

read_counts <- read_delim(
  here(
    input_path,
    "mock_analysis",
    "qiime2",
    "abundance_tables",
    "feature-table.tsv"
  ),
  delim = "\t",
  col_names = TRUE,
  skip = 1
) %>%
  rename("OTU" = `#OTU ID`) %>%
  column_to_rownames("OTU")

## Generate phyloseq object ----
OTU <- otu_table(read_counts, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(asv_names))
SAM <- sample_data(
  metadata %>% 
    filter(group == "mock") %>% 
    column_to_rownames("sample")
  )

abundance_data <- phyloseq(OTU,TAX,SAM)


## Normalize data ----
set.seed(12345)

norm_data <- rarefy_even_depth(
  abundance_data,
  rngseed = FALSE
)

## Calculate alpha diversity ----
alpha_div <- estimate_richness(norm_data,
                               measures = "Shannon") %>%
  rownames_to_column("sample") %>%
  left_join(metadata)

mock_otu_table <- otu_table(
  data.frame(
    row.names = c(
      "A",
      "B",
      "C",
      "D",
      "E",
      "F",
      "G",
      "H"
    ),
    original = c(
      89.1,
      8.9,
      0.89,
      0.089,
      0.089,
      0.0089,
      0.00089,
      0.000089
    )
  ),
  taxa_are_rows = TRUE
) 

mock_tax_table <- tax_table(as.matrix(data.frame(
  row.names = c("A",
                "B",
                "C",
                "D",
                "E",
                "F",
                "G",
                "H"),
  Genus = c(
    "Listeria",
    "Pseudomonas",
    "Bacillus",
    "Escherichia",
    "Salmonella",
    "Lactobacillus",
    "Enterococcus",
    "Staphylococcus"
  ),
  Species = c(
    "monocytogenes",
    "aeruginosa",
    "subtilis",
    "coli",
    "enterica",
    "fermentum",
    "faecalis",
    "aureus"
  )
)))

original_abundance <- phyloseq(mock_otu_table, mock_tax_table)
alpha_div_original <- estimate_richness(original_abundance, 
                                        measures = "Shannon") %>%
  rownames_to_column("sample") %>%
  left_join(metadata)

alpha_div_all <- rbind(alpha_div,
                       alpha_div_original) %>%
  mutate(kit = ifelse(is.na(kit), "original", kit))

p_alpha_mock <- ggplot(alpha_div_all, aes(kit, Shannon)) +
  geom_point(size = 4) +
  scale_y_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14))

## Calculate relative abundances ----
### Aggregate mock data to genus level
genus_data <- tax_glom(norm_data,
                       taxrank = "Genus")

original_mock_composition <- data.frame(
  Sample = rep("Original", 10),
  Genus = c("Listeria",
            "Pseudomonas",
            "Bacillus",
            "Saccharomyces",
            "Escherichia",
            "Salmonella",
            "Lactobacillus",
            "Enterococcus",
            "Cryptococcus",
            "Staphylococcus"),
  Abundance = c(89.1,
                8.9,
                0.89,
                0.89,
                0.089,
                0.089,
                0.0089,
                0.00089,
                0.00089,
                0.000089)
)

## Transform to relative abundance
relabund_data <- transform_sample_counts(
  genus_data, 
  function(x) x*100 / sum(x)
  )

tax_data <- psmelt(relabund_data) %>%
  select(Sample, Genus, Abundance) %>%
  mutate(Genus = sub("-Shigella", "", Genus))

plot_data <- rbind(original_mock_composition,
                   tax_data) %>%
  left_join(metadata, by = c("Sample" = "sample")) %>%
  mutate(kit = ifelse(is.na(kit), "original",kit)) %>%
  mutate(plot_group = ifelse(
    Genus %in% original_mock_composition$Genus, Genus, "Other"
  )) %>%
  mutate(plot_group = factor(plot_group,
                             levels = c("Listeria",
                                        "Pseudomonas",
                                        "Bacillus",
                                        "Saccharomyces",
                                        "Escherichia",
                                        "Salmonella",
                                        "Lactobacillus",
                                        "Enterococcus",
                                        "Cryptococcus",
                                        "Staphylococcus",
                                        "Other"))) %>%
  group_by(kit, plot_group) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup()

plot_data_2 <- plot_data %>%
  filter(plot_group != "Listeria")

data_summary <- plot_data %>%
  pivot_wider(names_from = "kit",
              values_from = "Abundance") %>%
  mutate_at(vars(c("original","purelink","qiagen")),
            ~round(., 4)) %>%
  mutate_at(vars(c("original","purelink","qiagen")),
            ~ifelse(is.na(.), "", .)) %>%
  rename("Genus" = plot_group,
         "Original" = original,
         "Purelink" = purelink,
         "Qiagen" = qiagen) %>%
  filter(! Genus %in% c("Cryptococcus","Saccharomyces"))
  

palette <- c(
  "Listeria" = "#80b1d3",
  "Pseudomonas" = "#fb8072",
  "Bacillus" = "#bebada",
  "Saccharomyces" = "#ffffb3",
  "Escherichia" = "#8dd3c7",
  "Salmonella" = "#fdb462",
  "Lactobacillus" = "#b3de69",
  "Enterococcus" = "#fccde5",
  "Cryptococcus" = "#bc80bd",
  "Staphylococcus" = "#ccebc5",
  "Other" = "#d9d9d9"
)


p_relabund_mock <- ggplot(plot_data, aes(kit, Abundance, fill = plot_group)) +
  geom_col() +
  scale_fill_manual(values = palette) +
  labs(y = "Relative abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank())

p_relabund_mock_detailed <- ggplot(plot_data_2, aes(kit, Abundance, fill = plot_group)) +
  geom_col() +
  scale_fill_manual(values = palette) +
  labs(y = "Relative abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 14))

## Merge plots ----

p_relabund_complete <- p_relabund_mock + 
  gridExtra::tableGrob(data_summary,rows = NULL) + 
  p_relabund_mock_detailed + 
  p_alpha_mock +
  plot_layout(heights = c(1,0.6)) +
  plot_annotation(tag_levels = "A")

ggsave(
  here(
    output_path,
    "02_16S_relabund_mock.png"
  ),
  p_relabund_complete,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 22,
  width = 25
)


# 03.  Sample analysis ----
## Generate phyloseq object ----
asv_names_samples <- read_delim(
  here(
    input_path,
    "sample_analysis",
    "qiime2",
    "rel_abundance_tables",
    "rel-table-ASV_with-DADA2-tax.tsv"
  ),
  delim = "\t",
  col_names = TRUE
) %>%
  select(1:8) %>%
  rename("OTU" = ID) %>%
  column_to_rownames("OTU")

read_counts_samples <- read_delim(
  here(
    input_path,
    "sample_analysis",
    "qiime2",
    "abundance_tables",
    "feature-table.tsv"
  ),
  delim = "\t",
  col_names = TRUE,
  skip = 1
) %>%
  rename("OTU" = `#OTU ID`) %>%
  column_to_rownames("OTU")

OTU_samples <- otu_table(read_counts_samples, taxa_are_rows = TRUE)
TAX_samples <- tax_table(as.matrix(asv_names_samples))
SAM_samples <- sample_data(
  metadata %>% 
    filter(group != "mock") %>% 
    column_to_rownames("sample")
)

abundance_data_samples <- phyloseq(OTU_samples,
                                   TAX_samples,
                                   SAM_samples)


## Normalize data ----
norm_data_samples <- rarefy_even_depth(
  abundance_data_samples,
  rngseed = FALSE
)


## Calculate alpha diversity ----
alpha_div_samples <- estimate_richness(norm_data_samples,
                               measures = "Shannon") %>%
  rownames_to_column("sample") %>%
  left_join(metadata) %>%
  mutate(sample = sub("Q|PL", "", sample))

p_alpha_samples <- ggplot(alpha_div_samples, 
                          aes(kit, 
                              Shannon, 
                              group = sample, 
                              color = kit,
                              shape = sample)) +
  geom_line(color = "black") +
  geom_point(size = 4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.position = "none")

## Calculate beta diversity ----
beta_bray <- ordinate(
  physeq = norm_data_samples,
  method = "MDS",
  distance = "bray",
  autotransform = FALSE
)

plot_data <- plot_ordination(
  norm_data_samples,
  beta_bray,
  color = "kit",
  justDF = TRUE
) %>%
  rownames_to_column("sample") %>%
  mutate(sample = sub("Q|PL", "", sample))

p_beta <- ggplot(plot_data, aes(Axis.1, Axis.2, color = kit, shape = sample)) +
  geom_point(size = 4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

## Calculate relative abundance ----
aggregated_sample_data <- tax_glom(norm_data_samples, 
                                   taxrank = "Phylum")

relabund_data_samples <- transform_sample_counts(
  aggregated_sample_data, 
  function(x) x*100 / sum(x)
)

plot_data_samples <- psmelt(relabund_data_samples)

summary_data_samples <- plot_data_samples %>%
  select(-kit) %>%
  mutate(Abundance = round(Abundance, 2)) %>%
  pivot_wider(names_from = "Sample",
              values_from = "Abundance") %>%
  select(Phylum,`370PL`,`370Q`,`371PL`,`371Q`,`386PL`,`386Q`)

p_relabund_samples <- ggplot(plot_data_samples, 
       aes(Sample, Abundance, fill = Phylum)) +
  geom_col() +
  labs(y = "Relative abundance") +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  facet_wrap(~kit, scales = "free_x")


## Merge plots and write table ----
p_samples <- (p_alpha_samples + p_beta) / p_relabund_samples +
  plot_layout(heights = c(0.6, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(
  here(
    output_path,
    "02_16S_samples.png"
  ),
  p_samples,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 20,
  width = 20
)

write_delim(
  summary_data_samples,
  here(
    output_path,
    "02_16S_relabundance_samples.txt"
  ),
  delim = "\t"
)
