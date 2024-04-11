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
library(pairwiseAdonis)

# Set options ----
options(scipen = 999)
set.seed(12345)

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
### Calculate bray-curtis distance ----
bray_dists <- phyloseq::distance(
  bracken_data,
  method = "bray"
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

#### PERMANOVA ----
adonis_test <- adonis2(
  bray_dists ~ group + sex,
  data = sample_beta_data,
  method = "bray",
  permutations = 9999)

#### Dispersion tests ----
disp_test_groups <- betadisper(
    bray_dists,
    type = "median",
    group = sample_beta_data$group
  )

disp_test_sex <- betadisper(
  bray_dists,
  type = "median",
  group = sample_beta_data$sex
)

#### Pairwise PERMANOVA ----
pairwise.adonis(
  bray_dists,
  sample_beta_data$group
)

pairwise.adonis(
  bray_dists,
  sample_beta_data$sex
)

### Boxplot of distances ----
boxplot_data <- as.data.frame(
  as.matrix(
    bray_dists
  )
) %>%
  rownames_to_column("id1") %>%
  pivot_longer(cols = -id1,
               names_to = "id2",
               values_to = "dist") %>%
  filter(id1 != id2) %>%
  left_join(metadata[,c("id", "sex", "group")],
            by = c("id1" = "id")) %>%
  dplyr::rename("sex1" = sex,
         "group1" = group) %>%
  left_join(metadata[,c("id", "sex", "group")],
            by = c("id2" = "id")) %>%
  dplyr::rename("sex2" = sex,
         "group2" = group)

p_bray_boxplot <- boxplot_data %>%
  filter(group1 == group2) %>%
  ggplot(aes(group1, dist, fill = group1)) +
  stat_boxplot(geom = "errorbar", width = 0.4, 
               position = position_dodge(0.8)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_signif(y_position = c(0.62, 0.64),
              xmin = c(3,2),
              xmax = c(4,3),
              annotation = c("*","*"),
              tip_length = 0.02,
              size = 0.6,
              textsize = 4,
              vjust = 0.5) +
  labs(y = "Bray-Curtis distance",
       title = "B") +
  scale_fill_manual(values = palette_group) +
  scale_x_discrete(labels = group_names) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")



# 04. Constrained analysis ----
## Prepare data ----
ord_data <- bracken_data %>%
  # select the bacterial fraction of the data
  tax_select("Bacteria", 
             ranks_searched = "Rank1", 
             strict_matches = TRUE) %>%
  # aggregate to genus-level
  tax_agg(rank = "Rank6") %>%
  # transform to compositional
  tax_transform(trans = "compositional", keep_counts = TRUE) %>%
  # filter out low-abundant taxa
  tax_filter(
    min_sample_abundance = 0.01,
    tax_level = "Rank6",
    use_counts = FALSE,
    prev_detection_threshold = 0
  ) %>%
  # revert to counts
  ps_get(counts = TRUE)

sampledf <- data.frame(sample_data(bracken_data)) %>%
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
                                      "none")),
         sex = factor(sex)) %>%
  select(sex, group)

vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

otutable <- vegan_otu(ord_data)

## Run distance-based constrained analysis ----
ord_all <- capscale(
  otutable ~ ., 
  distance = "bray", 
  data = sampledf,
  autotransform = FALSE,
  sqrt.dist = TRUE
  )

ord_group <- capscale(
  otutable ~ group, 
  distance = "bray", 
  data = sampledf,
  autotransform = FALSE,
  sqrt.dist = TRUE
)

ord_sex <- capscale(
  otutable ~ sex, 
  distance = "bray", 
  data = sampledf,
  autotransform = FALSE,
  sqrt.dist = TRUE
)

ord_interaction <- capscale(
  otutable ~ sex * group, 
  distance = "bray", 
  data = sampledf,
  autotransform = FALSE,
  sqrt.dist = TRUE
)

## Check significance of models ----
anova.cca(ord_all)
anova.cca(ord_group)
anova.cca(ord_sex)
anova.cca(ord_interaction)

anova.cca(ord_all, by = "axis")
anova.cca(ord_all, by = "terms")
anova.cca(ord_all, by = "margin")

## Check if there are unimportant terms in the model
spec_rda <- step(ord_all, scope = formula(ord_all), test = "perm")
summary(spec_rda)
## No change to model

## Inspect collinearity, exclude variables with a VIF > 20
vif.cca(ord_all)
## No variables > 20

## Get R-squared for the whole model
RsquareAdj(ord_all)

## Check significance of factors ----
anova.cca(ord_all, by = "onedf", permutations = 999)
anova.cca(ord_group, by = "onedf", permutations = 999)
anova.cca(ord_sex, by = "onedf", permutations = 999)
anova.cca(ord_interaction, by = "onedf", permutations = 999)
### Group "Antibiotic" is not significant
### Interaction model not significant for the interaction terms

## Fit variables to model ----
### Extract the scores from the RDA for plotting
site.scrs <- as.data.frame(scores(ord_all, display = "sites"))
site.scrs <- cbind(site.scrs, sex = sampledf$sex)
site.scrs <- cbind(site.scrs, group = sampledf$group)

### Fit environmental variables to the RDA
otu_env_fit <- envfit(ord_all, sampledf, permutations = 999)

env.scrs <- as.data.frame(scores(otu_env_fit, display = "factors")) %>%
  rownames_to_column("env.variables") %>%
  filter(env.variables != "groupantibiotic") %>%
  mutate(env.variables = case_when(
    env.variables == "sexMale" ~ "Male",
    env.variables == "sexFemale" ~ "Female",
    env.variables == "groupmonensin_antibiotic" ~ "M+P",
    env.variables == "groupmonensin" ~ "M",
    env.variables == "groupnone" ~ "None",
    TRUE ~ env.variables
  ))

### Fit taxa to the RDA
otu_taxa_fit <- envfit(ord_all, otutable, permutations = 999)

spp.scrs <- as.data.frame(scores(otu_taxa_fit, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs <- cbind(spp.scrs, pval = otu_taxa_fit$vectors$pvals)
spp.scrs <- cbind(spp.scrs, r = otu_taxa_fit$vectors$r)
sig.spp.scrs <- subset(spp.scrs, pval <= 0.05) %>%
  mutate(r = round(r, 3))

write_delim(
  sig.spp.scrs,
  here(
    output_path,
    "09_envfit_taxa_sign.txt"
  ),
  delim = "\t"
)

sig.spp.scrs <- slice_max(sig.spp.scrs, n = 10, order_by = r)

## Plot the results ----
p_rda <- ggplot(site.scrs, aes(x = CAP1, y = CAP2))+ 
  geom_hline(yintercept = 0,
             color = "grey80") +
  geom_vline(xintercept = 0,
             color = "grey80") +
  geom_point(aes(colour = group,
                 shape = sex),
             size = 3) + 
  stat_ellipse(aes(color = group),
               type = "t",
               level = 0.6) +
  geom_segment(data = sig.spp.scrs, 
               aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               colour = "grey40",
               alpha = 0.5,
               lwd = 0.4,
               inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = sig.spp.scrs, 
            aes(x = CAP1, y = CAP2, label = Species), 
            size = 4,
            fontface = "italic",
            min.segment.length = Inf,
            max.overlaps = 25) +
  geom_segment(data = env.scrs,
               aes(x = 0, xend = CAP1, y = 0, yend = CAP2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "#377eb8",
               lwd = 1) +
  geom_text(data = env.scrs,
            aes(x = CAP1, y = CAP2, label = env.variables,
                angle = (180/pi) * atan(CAP2/CAP1),
                hjust = (1 - 2 * sign(CAP1)) / 3),
            cex = 4,
            color = "#377eb8") +
  labs(title = "A",
       x = "CAP1 (53.14%)",
       y = "CAP2 (33.85%)") +
  scale_color_manual(values = palette_group, labels = group_names) +
  scale_fill_manual(values = palette_group) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

# 05. Abundance of specific taxa ----
genera_data <- bracken_data %>%
  tax_transform(trans = "compositional") %>%
  tax_agg(rank = "Rank6") %>%
  tax_select(sig.spp.scrs$Species, 
             ranks_searched = "Rank6", 
             strict_matches = TRUE) %>%
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
                           ordered = TRUE)) %>%
  psmelt()


tests <- lapply(
  sig.spp.scrs$Species, function(x) {
    genera_data %>%
      filter(Rank6 == x) %>%
      aov(Abundance ~ group, data = .) %>%
      broom::tidy() %>%
      mutate(ref = x)
  }
) %>%
  bind_rows() %>%
  filter(!is.na(p.value)) %>%
  filter(p.value < 0.05)

stat_bifido <- genera_data %>%
  filter(Rank6 == "Bifidobacterium") %>%
  aov(Abundance ~ group, data = .)

stat_megaspaera <- genera_data %>%
  filter(Rank6 == "Megasphaera") %>%
  aov(Abundance ~ group, data = .)

stat_megamonas <- genera_data %>%
  filter(Rank6 == "Megamonas") %>%
  aov(Abundance ~ group, data = .)

TukeyHSD(stat_bifido)
TukeyHSD(stat_megaspaera)
TukeyHSD(stat_megamonas)

## Create plots
p_bifido <- genera_data %>%
  filter(Rank6 == "Bifidobacterium") %>%
  ggplot(aes(group, Abundance, fill = group)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot() +
  geom_signif(y_position = c(0.32, 0.35),
              xmin = c(1,1),
              xmax = c(2,4),
              annotation = c("*","*"),
              tip_length = 0.02,
              size = 0.6,
              textsize = 5) +
  labs(title = "C",
       subtitle = "*Bifidobacterium*",
       y = "Relative abundance") +
  scale_fill_manual(values = palette_group) +
  scale_y_continuous(limits = c(0, 0.38)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.subtitle = ggtext::element_markdown(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

p_megasphaera <- genera_data %>%
  filter(Rank6 == "Megasphaera") %>%
  ggplot(aes(group, Abundance, fill = group)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot() +
  geom_signif(y_position = c(0.07, 0.09),
              xmin = c(2,3),
              xmax = c(3,4),
              annotation = c("*","*"),
              tip_length = 0.04,
              size = 0.6,
              textsize = 5) +
  labs(subtitle = "*Megasphaera*",
       y = "Relative abundance") +
  scale_fill_manual(values = palette_group) +
  scale_y_continuous(limits = c(0, 0.38)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.subtitle = ggtext::element_markdown(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

p_megamonas <- genera_data %>%
  filter(Rank6 == "Megamonas") %>%
  ggplot(aes(group, Abundance, fill = group)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot() +
  geom_signif(y_position = c(0.12, 0.14),
              xmin = c(2,3),
              xmax = c(3,4),
              annotation = c("*","*"),
              tip_length = 0.04,
              size = 0.6,
              textsize = 5) +
  labs(subtitle = "*Megamonas*",
       y = "Relative abundance") +
  scale_fill_manual(values = palette_group) +
  scale_y_continuous(limits = c(0, 0.38)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.subtitle = ggtext::element_markdown(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")


rest_plots <- lapply(
  sig.spp.scrs$Species[-c(1,9,10)],
  function(x) {
    genera_data %>%
      filter(Rank6 == x) %>%
      ggplot(aes(group, Abundance, fill = group)) +
      stat_boxplot(geom = "errorbar", width = 0.5) +
      geom_boxplot() +
      labs(subtitle = paste0("*", x, "*"),
           y = "Relative abundance") +
      scale_fill_manual(values = palette_group) +
      scale_y_continuous(limits = c(0, 0.38)) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            plot.subtitle = ggtext::element_markdown(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none")
  }
)

# 06. Merge plots ----

all_plots <- append(rest_plots, list(p_megasphaera, p_megamonas), after = 7)
all_plots <- append(all_plots, list(p_bifido), after = 0)


plot_design <- c(
  "AAAB
   AAAB
   CDEF
  "
)

p_beta <- p_rda + 
  p_bray_boxplot + 
  p_bifido +
  p_megasphaera +
  p_megamonas +
  guide_area() +
  plot_layout(design = plot_design,
              guides = "collect")

ggsave(
  here(
    output_path,
    "09_rda_ord.png"
  ),
  p_beta,
  device = "png",
  units = "cm", 
  dpi = 600,
  height = 25,
  width = 27
)

