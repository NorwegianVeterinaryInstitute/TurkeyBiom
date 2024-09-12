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
library(tibble)

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

## Import metadata -----
metadata <- read_delim(
  here(
    input_path,
    "00_final_treatment_metadata.txt"
  ),
  delim = "\t"
)

## Import phyloseq object ----
bracken_data <- read_rds(
  here(
    output_path,
    "bracken_rarefied.rds"
  )
)

metadata_ps <- bracken_data@sam_data
metadata_ps$id <- rownames(metadata_ps)

metadata_ps_clean <- as_tibble(metadata_ps) %>%
  mutate(group = case_when(
    monensin == TRUE & antibiotic == TRUE ~ "monensin_antibiotic",
    monensin == TRUE & antibiotic == FALSE ~ "monensin",
    monensin == FALSE & antibiotic == TRUE ~ "antibiotic",
    monensin == FALSE & antibiotic == FALSE ~ "none"
  )) %>%
  mutate(group = factor(group,
                        levels = c("monensin_antibiotic",
                                   "antibiotic",
                                   "monensin",
                                   "none"))) %>%
  mutate(seq_id = sub(".+-([A-Z].+-.+)_pe_db1..+", "\\1", id)) %>%
  left_join(metadata[,c("seq_id","farm_id","flock_id")], by = "seq_id") %>%
  select(id,
         sex,
         age_days,
         sample_month,
         group,
         farm_id,
         flock_id) %>%
  mutate(sex = factor(sex),
         farm_id = factor(farm_id),
         flock_id = factor(flock_id)) %>%
  column_to_rownames("id")

bracken_data@sam_data <- sample_data(metadata_ps_clean)


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
  Abundance ~ group*OTU  + sex*OTU + farm_id*OTU,
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
  separate(comp1, sep = ":", into = c("OTU1","group1")) %>%
  separate(comp2, sep = ":", into = c("OTU2","group2")) %>%
  filter(group1 != group2) %>%
  filter(OTU1 == OTU2) %>%
  filter(`p adj` <= 0.05)

tukey_group <- as.data.frame(tukey_data$`group:OTU`) %>%
  tibble::rownames_to_column("id") %>%
  separate(id, sep = "-", into = c("comp1","comp2")) %>%
  separate(comp1, sep = ":", into = c("group1","OTU1")) %>%
  separate(comp2, sep = ":", into = c("group2","OTU2")) %>%
  filter(group1 != group2) %>%
  filter(OTU1 == OTU2) %>%
  filter(`p adj` <= 0.05) %>%
  select(OTU1, group1, OTU2, group2, diff, lwr, upr, `p adj`)

tukey_farm <- as.data.frame(tukey_data$`OTU:farm_id`) %>%
  tibble::rownames_to_column("id") %>%
  separate(id, sep = "-", into = c("comp1","comp2")) %>%
  separate(comp1, sep = ":", into = c("OTU1","group1")) %>%
  separate(comp2, sep = ":", into = c("OTU2","group2")) %>%
  filter(group1 != group2) %>%
  filter(OTU1 == OTU2) %>%
  filter(`p adj` <= 0.05) %>%
  select(OTU1, group1, OTU2, group2, diff, lwr, upr, `p adj`)

tukey_all <- rbind(
  tukey_sex,
  tukey_group,
  tukey_farm
) %>%
  select(OTU1, group1, group2, `p adj`) %>%
  mutate(`p adj` = round(`p adj`, 4))


write_delim(
  tukey_all,
  here(
    output_path,
    "08_phylum_tukey.txt"
  ),
  delim = "\t"
)


## Phylum plot ----
p_phylum_comp <- comp_barplot(
  bracken_data,
  tax_level = "Rank2",
  sample_order = "bray",
  n_taxa = 5
  ) +
  geom_point(aes(color = sex),
            y = -0.01,
            size = 2) +
  scale_color_manual(
    values = c(
      "Female" = "grey80",
      "Male" = "black"
    )
  ) +
  labs(fill = "Phylum",
       color = "Host sex",
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
  mutate(mean_abundance = mean_abundance * 100) %>%
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
