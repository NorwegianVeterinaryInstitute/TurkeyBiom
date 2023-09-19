# ABSTRACT
# This script calculates the alpha and beta
# diversity for the 105 turkey caecal samples
# based on the bracken reports

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)
library(vegan)
library(phyloseq)
library(patchwork)
library(microViz)
library(here)

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
  "R",
  "Data",
  sep = "/"
)

bracken_data <- import_biom(
  here(
    input_path,
    "kraken_biom",
    "bracken_pluspfp.biom"
  )
)

## Add metadata to the phyloseq object
sample_data <- as.data.frame(bracken_data@otu_table) %>%
  t() %>%
  as.data.frame %>%
  rownames_to_column("id") %>%
  select(id)

metadata <- read_delim(
  here(
    input_path,
    "00_metadata_samples.txt"
  ),
  delim = "\t"
)

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

sample_df <- data.frame(id = matches_seq_ids,
                        seq_id = matches_metadata_ids) %>%
  left_join(metadata[,c("seq_id","sex","farm_type","age_days")]) %>%
  select(-seq_id) %>%
  column_to_rownames("id")
  
bracken_data@sam_data <- sample_data(sample_df)

bracken_fixed <- bracken_data %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )

# 02. Alpha diversity ----
alpha_div <- estimate_richness(
  bracken_fixed,
  measures = c("Shannon","Simpson","Observed")
  ) %>%
  rownames_to_column("id") %>%
  mutate(id = sub("_pe_pluspfp.bracken.kraken2.report_bracken_species", "", id),
         id = sub("X.+([A-Z].+\\..+)", "\\1", id),
         id = sub("\\.", "-", id)) %>%
  left_join(metadata[,c("seq_id","sex","age_days","farm_type")],
            by = c("id" = "seq_id"))

p_alpha_sex <- ggplot(alpha_div, aes(sex, Shannon, fill = sex)) +
  stat_boxplot(geom = "errorbar", width = 0.4) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "none")

p_alpha_sex_point <- ggplot(
  alpha_div, aes(age_days, Shannon, color = sex, shape = farm_type)
  ) +
  geom_point(size = 4) +
  labs(x = "Age (days)",
       color = NULL,
       shape = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

p_alpha <- p_alpha_sex + p_alpha_sex_point +
  plot_layout(widths = c(0.5, 1))


ggsave(
  here(
    output_path,
    "07_alpha_div_sex.png"
  ),
  p_alpha,
  device = "png",
  units = "cm", 
  dpi = 600,
  height = 20,
  width = 25
)

# 03. Beta diversity ----
p_beta <- bracken_fixed %>%
  tax_filter(min_prevalence = 2.5 / 100, verbose = FALSE) %>%
  tax_transform(rank = "Rank6", trans = "compositional") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(
    alpha = 0.6, size = 2, color = "sex", shape = "farm_type",
    plot_taxa = 1:8, tax_vec_length = 0.5,
    tax_lab_style = tax_lab_style(
      type = "text", max_angle = 90, aspect_ratio = 0.7,
      size = 3, fontface = "bold"
    ),
  ) +
  theme_classic(12) +
  coord_fixed(0.7, clip = "off") +
  stat_ellipse(aes(color = sex)) +
  xlim(-0.3, 0.6)


p_relabund <- bracken_fixed %>%
  comp_barplot(
    tax_level = "Rank6", n_taxa = 12, bar_width = 0.7, sample_order = "bray",
    bar_outline_colour = "white", bar_outline_width = 0.001
  ) +
  facet_grid(~sex, scales = "free_x", space = "free") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank())

p_all <- p_beta / p_relabund + plot_layout(guides = "collect", ncol = 1,
                                           heights = c(1, 0.7))


ggsave(
  here(
    output_path,
    "07_beta_div_sex.png"
  ),
  p_all,
  device = "png",
  units = "cm", 
  dpi = 600,
  height = 25,
  width = 22
)

ord_explore(bracken_fixed)
