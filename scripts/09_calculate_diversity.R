# ABSTRACT
# This script calculates the diversity indices, runs
# statistics, and generates plots of interest

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggsignif)
library(here)
library(microViz)
library(tibble)
library(patchwork)
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

## Import phyloseq object ----
bracken_data <- read_rds(
  here(
    output_path,
    "bracken_rarefied.rds"
  )
)

# 02. Define palettes ----
palette_sex <- c(
  "Male" = "#80b1d3",
  "Female" = "#fb8072"
)

palette_group <- c(
  "monensin_antibiotic" = "#fdb462",
  "antibiotic" = "#8dd3c7",
  "monensin" = "#fb8072",
  "none" = "#bebada"
)

## Set metadata
metadata <- bracken_data@sam_data
metadata$id <- rownames(metadata)

metadata <- as_tibble(metadata) %>%
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
                                   "none"),
                        ordered = TRUE))

group_names <- c(
  "monensin_antibiotic" = "M + P",
  "monensin" = "M",
  "antibiotic" = "P",
  "none" = "None"
)

options(scipen = 999)

# 02. Alpha diversity ----
## Calculate alpha diversity ----
alpha_div <- bracken_data %>%
  estimate_richness(measures = c("Shannon","Chao1","Observed")) %>%
  cbind(metadata)

## Generate summary of Alpha diversity values ----
alpha_summary <- alpha_div %>%
  mutate(Shannon = round(Shannon, 3),
         Chao1 = round(Chao1, 3)) %>%
  group_by(group) %>%
  summarise(
    n = n(),
    shannon_mean = mean(Shannon),
    shannon_range = paste0(min(Shannon), " - ", max(Shannon)),
    chao1_mean = mean(Chao1),
    chao1_range = paste0(min(Chao1), " - ", max(Chao1))
  )

## Statistical analyses ----
### ANOVA on linear model controlling for sex
alpha_model <- aov(
  Shannon ~ group + sex, 
  data = alpha_div
)

eta_squared <- etaSquared(alpha_model)

TukeyHSD(
  alpha_model,
  conf.level = 0.95
  )

## Plot results ----
p_alpha_div <- ggplot(alpha_div, aes(group, Shannon, fill = group)) +
  stat_boxplot(geom = "errorbar", 
               width = 0.4, 
               linewidth = 1) +
  geom_boxplot(outlier.shape = NA,
               linewidth = 1) +
  geom_jitter(width = 0.05,
              alpha = 0.5,
              size = 3) +
  geom_signif(y_position = 6.1,
              xmin = 1,
              xmax = 4,
              annotation = "p = 0.0347",
              tip_length = 0.02,
              size = 0.8,
              textsize = 4) +
  scale_x_discrete(labels = group_names) +
  scale_y_continuous(limits = c(4.4,6.2)) +
  scale_fill_manual(values = palette_group) +
  labs(y = "Shannon diversity") +
  guides(fill = guide_legend(override.aes=list(shape = 15))) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

ggsave(
  here(
    output_path,
    "09_alpha_div.png"
  ),
  p_alpha_div,
  device = "png",
  units = "cm", 
  dpi = 600,
  height = 12,
  width = 20
)

# 03. Beta diversity ----
### Calculate beta diversity ----
beta_ord <- ordinate(
  physeq = bracken_data,
  method = "NMDS", 
  distance = "bray",
  autotransform = FALSE,
  k = 2,
  try = 30
  )

### Statistical analyses ----
sample_beta_data <- as(sample_data(bracken_data), "data.frame") %>%
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
                                   "none"),
                        ordered = TRUE))

adonis_test <- adonis2(
  phyloseq::distance(
    bracken_data,
    method = "bray"
  ) ~ group + sex,
  data = sample_beta_data,
  method = "bray",
  permutations = 9999)

disp_test_groups <- betadisper(
    phyloseq::distance(
      bracken_data,
      method = "bray"
    ),
    type = "median",
    group = sample_beta_data$group
  )

disp_test_sex <- betadisper(
  phyloseq::distance(
    bracken_data,
    method = "bray"
  ),
  type = "median",
  group = sample_beta_data$sex
)

TukeyHSD(disp_test_groups,
         conf.level = 0.95)

TukeyHSD(disp_test_sex,
         conf.level = 0.95)


### Plot NMDS ----
#### NMDS dotplot
beta_plot_data <- as.data.frame(beta_ord$points) %>%
  rownames_to_column("id") %>%
  left_join(metadata)

p_nmds_beta <- ggplot(beta_plot_data, aes(MDS1, MDS2, color = group)) +
  geom_point(size = 3) +
  labs(title = "A") +
  scale_color_manual(values = palette_group,
                     labels = group_names) +
  theme_bw() +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 7)) +
  xlim(-0.5, 0.53)


#### Boxplot of distances
boxplot_data <- as.data.frame(
  as.matrix(
    distance(
      bracken_data,
      method = "bray"
    )
  )
) %>%
  rownames_to_column("id1") %>%
  pivot_longer(cols = -id1,
               names_to = "id2",
               values_to = "dist") %>%
  filter(id1 != id2) %>%
  left_join(metadata[,c("id", "sex", "group")],
            by = c("id1" = "id")) %>%
  rename("sex1" = sex,
         "group1" = group) %>%
  left_join(metadata[,c("id", "sex", "group")],
            by = c("id2" = "id")) %>%
  rename("sex2" = sex,
         "group2" = group)

p_bray_boxplot <- boxplot_data %>%
  filter(group1 == group2) %>%
  ggplot(aes(group1, dist, fill = group1)) +
  stat_boxplot(geom = "errorbar", width = 0.4, 
               position = position_dodge(0.8)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(y = "Bray-Curtis distance",
       title = "B") +
  scale_fill_manual(values = palette_group) +
  scale_x_discrete(labels = group_names) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

## Plot PCA ----
pca_data <- bracken_data %>%
  tax_select("Bacteria", 
             ranks_searched = "Rank1", 
             strict_matches = TRUE) %>%
  tax_agg("Rank6") %>%
  ord_calc(method = "PCA") %>%
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
                        ordered = TRUE))

p_pca <- pca_data %>%
  ord_plot( 
    size = 2.5, 
    color = "group",
    plot_taxa = 1:10,
    plot_samples = TRUE,
    tax_vec_length = 0.6,
    auto_caption = NA,
    tax_lab_style = tax_lab_style(
      type = "text", 
      max_angle = 90, 
      aspect_ratio = 1,
      size = 2, 
      fontface = "italic"
    )
  ) +
  theme_classic() +
  scale_color_manual(values = palette_group,
                     labels = group_names) +
  labs(color = NULL,
       title = "C") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  xlim(-600, 2500)

p_scree <- bracken_data %>%
  tax_select("Bacteria", 
             ranks_searched = "Rank1", 
             strict_matches = TRUE) %>%
  tax_agg("Rank6") %>%
  ord_calc(method = "PCA") %>%
  ord_get() %>%
  plot_scree() +
  scale_x_discrete(limits = paste0("PC", 1:6)) +
  labs(y = "Eigenvalues",
       x = "Principal components",
       title = "D") +
  theme_bw() +
  theme(panel.grid = element_blank())


p_beta <- p_nmds_beta + p_bray_boxplot + p_pca + p_scree +
  plot_layout(ncol = 2,
              widths = c(0.7, 0.3),
              heights = c(0.4, 0.6))

ggsave(
  here(
    output_path,
    "09_beta_div.png"
  ),
  p_beta,
  device = "png",
  units = "cm", 
  dpi = 600,
  height = 16,
  width = 20
)


# 04. Canonical correspondence analysis ----
### TODO: Need to run this correctly! ----
cca_data <- bracken_data %>%
  tax_select("Bacteria", 
             ranks_searched = "Rank1", 
             strict_matches = TRUE) %>%
  tax_agg(rank = "Rank6") %>%
  ord_calc(
    constraints = c("monensin"),
    method = "CCA"
  )

anova(cca_data@ord)
anova(cca_data@ord, by = "axis")

p_cca <- cca_data %>% ord_plot(
  size = 3,
  color = "monensin_plot",
  plot_taxa = 1:20,
  tax_vec_length = 2.5,
  auto_caption = NA,
  tax_lab_style = tax_lab_style(
    type = "text",
    max_angle = 90,
    aspect_ratio = 0.6,
    size = 2,
    fontface = "italic"
  )
) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = palette_monensin) +
  labs(color = "Monensin") +
  xlim(-3, 5) +
  ylim(-5, 6) +
  coord_fixed(ratio = 0.6)

ggsave(
  here(
    output_path,
    "09_cca.png"
  ),
  p_cca,
  device = "png",
  units = "cm", 
  dpi = 600,
  height = 15,
  width = 20
)
