# ABSTRACT
# This script calculates the relative abundances
# and generates plots of interest

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(phyloseq)
library(ggplot2)
library(here)
library(microViz)
library(purrr)
library(broom)

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

## Import phyloseq object ----
bracken_data <- read_rds(
  here(
    output_path,
    "bracken_rarefied.rds"
  )
)

options(scipen = 999)

# 02. Plot phylum abundance ----
## Calculate abundances ----
kingdom_data <- bracken_data %>%
  tax_agg("Rank1") %>%
  tax_transform(trans = "compositional") %>%
  psmelt() %>%
  group_by(OTU) %>%
  summarise(
    mean = mean(Abundance),
    median = median(Abundance),
    min = min(Abundance),
    max = max(Abundance),
    sd = sd(Abundance)
  ) %>%
  arrange(-mean) %>%
  mutate_at(vars(-OTU),
            ~round(.*100, 4))

phylum_data <- bracken_data %>%
  tax_agg("Rank2") %>%
  tax_transform(trans = "compositional") %>%
  psmelt() %>%
  group_by(OTU) %>%
  summarise(
    mean = mean(Abundance),
    median = median(Abundance),
    min = min(Abundance),
    max = max(Abundance),
    sd = sd(Abundance)
  ) %>%
  arrange(-mean) %>%
  mutate_at(vars(-OTU),
            ~round(.*100, 4))

phylum_data_sex <- bracken_data %>%
  tax_agg("Rank2") %>%
  tax_transform(trans = "compositional") %>%
  psmelt() %>%
  group_by(OTU,sex) %>%
  summarise(
    mean = mean(Abundance),
    median = median(Abundance),
    min = min(Abundance),
    max = max(Abundance),
    sd = sd(Abundance)
  ) %>%
  arrange(-mean) %>%
  mutate_at(vars(-c(OTU,sex)),
            ~round(.*100, 4))

## Wilcoxon rank test major phyla ----
phylum_data <- bracken_data %>%
  tax_agg("Rank2") %>%
  tax_transform(trans = "compositional") %>%
  psmelt()

stats_data_sex <- phylum_data %>%
  group_by(sex, OTU) %>%
  summarise(
    mean = mean(Abundance)
  ) %>%
  arrange(-mean) %>%
  mutate(group = ifelse(mean < 0.01, "Other", OTU)) %>%
  group_by(sex, group) %>%
  summarise(
    mean = sum(mean)
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = "sex",
              values_from = "mean") %>%
  arrange(-Female)

stats_data_monensin <- phylum_data %>%
  group_by(monensin, OTU) %>%
  summarise(
    mean = mean(Abundance)
  ) %>%
  arrange(-mean) %>%
  mutate(group = ifelse(mean < 0.01, "Other", OTU)) %>%
  group_by(monensin, group) %>%
  summarise(
    mean = sum(mean)
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = "monensin",
              values_from = "mean") %>%
  arrange(-`TRUE`)

phyla_list <- stats_data_sex$group[-6]

tests_sex <- lapply(
  phyla_list,
  function(x) {
    data <- phylum_data %>%
      filter(Rank2 == x)
    
    wilcox.test(
      formula = Abundance ~ sex,
      data = data,
      alternative = "two.sided"
    )
  }
) %>%
  map_dfr(tidy)


tests_monensin <- lapply(
  phyla_list,
  function(x) {
    data <- phylum_data %>%
      filter(Rank2 == x)
    
    wilcox.test(
      formula = Abundance ~ monensin,
      data = data,
      alternative = "two.sided"
    )
  }
) %>%
  map_dfr(tidy)

tests_sex$phylum <- phyla_list
tests_monensin$phylum <- phyla_list

tests_complete_sex <- tests_sex %>%
  select(phylum, statistic, p.value) %>%
  left_join(stats_data_sex, by = c("phylum" = "group")) %>%
  mutate(p_corr = p.adjust(p.value, method = "bonf")) %>%
  select(phylum, statistic, p.value, p_corr, Male, Female) %>%
  mutate_at(vars(-c(phylum, statistic)),
            ~round(., 3))

tests_complete_monensin <- tests_monensin %>%
  select(phylum, statistic, p.value) %>%
  left_join(stats_data_monensin, by = c("phylum" = "group")) %>%
  mutate(p_corr = p.adjust(p.value, method = "bonf")) %>%
  select(phylum, statistic, p.value, p_corr, `TRUE`, `FALSE`) %>%
  mutate_at(vars(-c(phylum, statistic)),
            ~round(., 3))

## Phylum plot ----
p_phylum_comp <- comp_barplot(
  bracken_data,
  tax_level = "Rank2",
  sample_order = "bray",
  n_taxa = 10
  ) +
  labs(fill = "Phylum",
       y = "Relative abundance") +
  facet_grid(~sex, scales = "free_x", space = "free") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(
  here(
    output_path,
    "08_phylum_comp.png"
  ),
  p_phylum_comp,
  device = "png",
  units = "cm", 
  dpi = 600,
  height = 15,
  width = 25
)


# 03. Plot genera abundance ----
p_genus_comp <- comp_barplot(
  bracken_data,
  tax_level = "Rank6",
  sample_order = "bray",
  n_taxa = 10
) +
  scale_fill_manual(values = c(
    RColorBrewer::brewer.pal(10, "Set3"), 
    "grey70"
    )) +
  labs(fill = "Genus",
       y = "Relative abundance") +
  facet_grid(~sex, scales = "free_x", space = "free") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(
  here(
    output_path,
    "08_genus_comp.png"
  ),
  p_genus_comp,
  device = "png",
  units = "cm", 
  dpi = 600,
  height = 15,
  width = 25
)


