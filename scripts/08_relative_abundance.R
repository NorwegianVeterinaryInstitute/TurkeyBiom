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
library(lsr)

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

group_names <- c(
  "monensin_antibiotic" = "M + P",
  "monensin" = "M",
  "antibiotic" = "P",
  "none" = "None"
)

## Import phyloseq object ----
bracken_data <- read_rds(
  here(
    output_path,
    "bracken_rarefied.rds"
  )
) %>%
  ps_mutate(group = case_when(
    monensin == TRUE & antibiotic == TRUE ~ "monensin_antibiotic",
    monensin == TRUE & antibiotic == FALSE ~ "monensin",
    monensin == FALSE & antibiotic == TRUE ~ "antibiotic",
    monensin == FALSE & antibiotic == FALSE ~ "none"
  )) %>%
  ps_mutate(group = factor(group,
                           levels = c("monensin_antibiotic",
                                      "antibiotic",
                                      "monensin",
                                      "none"),
                           ordered = TRUE),
            gender = ifelse(sex == "Female", "F", "M"))

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

phylum_data_group <- bracken_data %>%
  tax_agg("Rank2") %>%
  tax_transform(trans = "compositional") %>%
  psmelt() %>%
  group_by(OTU,group) %>%
  summarise(
    mean = mean(Abundance),
    median = median(Abundance),
    min = min(Abundance),
    max = max(Abundance),
    sd = sd(Abundance)
  ) %>%
  arrange(-mean) %>%
  mutate_at(vars(-c(OTU,group)),
            ~round(.*100, 4)) %>%
  filter(mean > 1) %>%
  select(OTU, group, mean) %>%
  pivot_wider(names_from = "group",
              values_from = "mean")

## ANOVA on top five phyla ----
phylum_data <- bracken_data %>%
  tax_agg("Rank2") %>%
  psmelt() %>%
  filter(OTU %in% phylum_data_group$OTU)

aov_phylum <- aov(
  Abundance ~ group*OTU + sex*OTU,
  data = phylum_data
)

eta_squared <- etaSquared(aov_phylum)

tukey_data <- TukeyHSD(
  aov_phylum,
  conf.level = 0.95
)

tukey_sex <- as.data.frame(tukey_data$`OTU:sex`) %>%
  tibble::rownames_to_column("id") %>%
  separate(id, sep = "-", into = c("comp1","comp2")) %>%
  separate(comp1, sep = ":", into = c("OTU1","sex1")) %>%
  separate(comp2, sep = ":", into = c("OTU2","sex2")) %>%
  filter(sex1 != sex2) %>%
  filter(OTU1 == OTU2) %>%
  filter(`p adj` <= 0.05)

tukey_group <- as.data.frame(tukey_data$`group:OTU`) %>%
  tibble::rownames_to_column("id") %>%
  separate(id, sep = "-", into = c("comp1","comp2")) %>%
  separate(comp1, sep = ":", into = c("group1","OTU1")) %>%
  separate(comp2, sep = ":", into = c("group2","OTU2")) %>%
  filter(group1 != group2) %>%
  filter(OTU1 == OTU2) %>%
  filter(`p adj` <= 0.05)
  

### No significant differences detected


## Phylum plot ----
p_phylum_comp <- comp_barplot(
  bracken_data,
  tax_level = "Rank2",
  sample_order = "bray",
  n_taxa = 5
  ) +
  geom_point(aes(shape = sex),
            y = -0.01,
            size = 1.5) +
  labs(fill = "Phylum",
       shape = "Host sex",
       y = "Relative abundance") +
  facet_grid(
    ~group, 
    scales = "free_x", 
    space = "free", 
    labeller = as_labeller(group_names)
    ) +
  guides(
    fill = guide_legend(override.aes = list(shape = NA))
    ) +
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
  width = 30
)

# 03. Top genera table ----

genera_data <- bracken_data %>%
  tax_agg("Rank6") %>%
  tax_transform(trans = "compositional") %>%
  psmelt()

top_genera <- genera_data %>%
  mutate(Abundance = round(Abundance, 3)) %>%
  group_by(group, Rank6) %>%
  summarise(
    mean_abundance = round(mean(Abundance), 3)
  ) %>%
  ungroup() %>%
  slice_max(order_by = mean_abundance,
            n = 5,
            by = group) %>%
  pivot_wider(
    names_from = "group",
    values_from = "mean_abundance",
    values_fill = 0
  )

write_delim(
  top_genera,
  here(
    output_path,
    "08_relabund_stats_genera.txt"
  ),
  delim = "\t"
)
