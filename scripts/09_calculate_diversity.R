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

## Import metadata
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

# 02. Define palettes ----
palette_sex <- c(
  "Male" = "#80b1d3",
  "Female" = "#b2df8a"
)

palette_group <- c(
  "monensin_antibiotic" = "#fdb462",
  "antibiotic" = "#8dd3c7",
  "monensin" = "#fb8072",
  "none" = "#bebada"
)

## Set metadata
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
         flock_id = factor(flock_id),
         sample_month = factor(sample_month)) %>%
  column_to_rownames("id")

bracken_data@sam_data <- sample_data(metadata_ps_clean)

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
  cbind(metadata_ps_clean)

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
### ANOVA on linear model controlling for sex, farm_id, and flock_id
alpha_model <- aov(
  Shannon ~ group + sex + farm_id + flock_id, 
  data = alpha_div
)

summary(alpha_model)

## None significant, remove flock_id
alpha_model2 <- aov(
  Shannon ~ group + sex + farm_id, 
  data = alpha_div
)

summary(alpha_model2)
## All three significant

alpha_model3 <- aov(
  Shannon ~ group * sex * farm_id, 
  data = alpha_div
)

summary(alpha_model3)
### Interaction not significant, farm_id not significant
### Interaction model not relevant for alpha diversity

## Get R2 value from model
summary(
  lm(
    Shannon ~ group + sex + farm_id,
    data = alpha_div
    )
)

eta_squared <- etaSquared(alpha_model2)

alpha_tukey <- TukeyHSD(
  alpha_model2,
  conf.level = 0.95
  )

alpha_tukey_farm <- as.data.frame(alpha_tukey$farm_id)


## Plot results ----
p_alpha_div_sex <- ggplot(alpha_div, aes(sex, Shannon)) +
  stat_boxplot(geom = "errorbar", 
               width = 0.4, 
               linewidth = 1) +
  geom_boxplot(outlier.shape = NA,
               linewidth = 1,
               fill = "grey40") +
  geom_jitter(width = 0.05,
              alpha = 0.5,
              size = 3) +
  geom_signif(y_position = 6,
              xmin = 1,
              xmax = 2,
              annotation = "*",tip_length = 0.015,
              size = 0.8,
              textsize = 7) +
  scale_y_continuous(limits = c(4.4,6.3)) +
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

p_alpha_div_group <- ggplot(alpha_div, aes(group, Shannon, fill = group)) +
  stat_boxplot(geom = "errorbar", 
               width = 0.4, 
               linewidth = 1) +
  geom_boxplot(outlier.shape = NA,
               linewidth = 1) +
  geom_jitter(width = 0.05,
              alpha = 0.5,
              size = 3) +
  geom_signif(y_position = c(6,6.1,6.2),
              xmin = c(1,1,3),
              xmax = c(2,4,4),
              annotation = c("*","*","*"),
              tip_length = 0.015,
              size = 0.8,
              textsize = 7) +
  scale_x_discrete(labels = group_names) +
  scale_y_continuous(limits = c(4.4,6.3)) +
  scale_fill_manual(values = palette_group) +
  guides(fill = guide_legend(override.aes=list(shape = 15))) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))


p_alpha_div <- p_alpha_div_sex + p_alpha_div_group

ggsave(
  here(
    output_path,
    "09_alpha_div.png"
  ),
  p_alpha_div,
  device = "png",
  units = "cm", 
  dpi = 600,
  height = 16,
  width = 22
)

# 03. Beta diversity ----
### Calculate bray-curtis distance ----
bray_dists <- phyloseq::distance(
  bracken_data,
  method = "bray"
)

### Statistical analyses ----
sample_beta_data <- as(sample_data(bracken_data), "data.frame")

#### PERMANOVA ----
## Full model
adonis_test <- adonis2(
  bray_dists ~ group + sex + farm_id + flock_id,
  data = sample_beta_data,
  method = "bray",
  permutations = 9999)

## Model without flock_id
adonis_test2 <- adonis2(
  bray_dists ~ group + sex + farm_id,
  data = sample_beta_data,
  method = "bray",
  permutations = 9999)

## Model with interaction
adonis_test3 <- adonis2(
  bray_dists ~ group * sex * farm_id,
  data = sample_beta_data,
  method = "bray",
  permutations = 9999)

### Interaction terms not significant

#### Dispersion tests ----
disp_test_groups <- anova(
  betadisper(
    bray_dists,
    type = "median",
    group = sample_beta_data$group
    )
)

disp_test_sex <- anova(
  betadisper(
    bray_dists,
    type = "median",
    group = sample_beta_data$sex
  )
)

disp_test_farm <- anova(
  betadisper(
    bray_dists,
    type = "median",
    group = sample_beta_data$farm_id
  )
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

pairwise.adonis(
  bray_dists,
  sample_beta_data$farm_id
)

### Boxplot of distances ----
metadata_beta <- metadata_ps_clean %>%
  rownames_to_column("id")

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
  left_join(metadata_beta[,c("id", "sex", "group")],
            by = c("id1" = "id")) %>%
  dplyr::rename("sex1" = sex,
         "group1" = group) %>%
  left_join(metadata_beta[,c("id", "sex", "group")],
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
  scale_y_continuous(limits = c(0, 0.65)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

p_bray_boxplot_sex <- boxplot_data %>%
  filter(sex1 == sex2) %>%
  ggplot(aes(sex1, dist)) +
  stat_boxplot(geom = "errorbar", width = 0.4, 
               position = position_dodge(0.8)) +
  geom_boxplot(position = position_dodge(0.8),
               fill = "grey40") +
  geom_signif(y_position = 0.61,
              xmin = 1,
              xmax = 2,
              annotation = "**",
              tip_length = 0.02,
              size = 0.6,
              textsize = 4,
              vjust = 0.5) +
  labs(y = "Bray-Curtis distance",
       title = "C") +
  scale_y_continuous(limits = c(0, 0.65)) +
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
  select(sex, group, farm_id, flock_id, age_days) %>%
  mutate(sex = factor(sex),
         group = factor(group, ordered = FALSE),
         farm_id = factor(farm_id),
         flock_id = factor(flock_id))

vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

otutable <- vegan_otu(ord_data)

## Run distance-based constrained analysis ----
ord <- capscale(
  otutable ~ sex + group + farm_id, 
  distance = "bray", 
  data = sampledf,
  autotransform = FALSE,
  sqrt.dist = TRUE
  )

## Check significance of models ----
ord_summary <- summary(ord)

### Check significance
anova.cca(ord)
anova.cca(ord, by = "axis")
anova.cca(ord, by = "terms")
anova.cca(ord, by = "onedf")

## Inspect collinearity, exclude variables with a VIF > 20
vif.cca(ord)
## No variables > 20

## Get R-squared for the whole model
RsquareAdj(ord)

## Fit variables to model ----
### Extract the scores from the RDA for plotting
site.scrs <- as.data.frame(scores(ord, display = "sites"))
site.scrs <- cbind(site.scrs, sex = sampledf$sex)
site.scrs <- cbind(site.scrs, group = sampledf$group)
site.scrs <- cbind(site.scrs, farm_id = sampledf$farm_id)

### Fit environmental variables to the RDA
otu_env_fit <- envfit(ord, sampledf, permutations = 999)

env.scrs <- as.data.frame(scores(otu_env_fit, display = "factors")) %>%
  rownames_to_column("env.variables") %>%
  mutate(env.variables = case_when(
    env.variables == "sexMale" ~ "Male",
    env.variables == "sexFemale" ~ "Female",
    env.variables == "groupmonensin_antibiotic" ~ "M+P",
    env.variables == "groupmonensin" ~ "M",
    env.variables == "groupantibiotic" ~ "P",
    env.variables == "groupnone" ~ "None"
  )) %>%
  filter(!is.na(env.variables))

### Fit taxa to the RDA
otu_taxa_fit <- envfit(ord, otutable, permutations = 999)

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

sig.spp.scrs <- slice_max(sig.spp.scrs, n = 15, order_by = r)

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
       x = paste0("CAP1: ",
                  round(ord_summary$cont$importance[2,"CAP1"]*100, 2),
                  "%"),
       y = paste0("CAP2: ",
                  round(ord_summary$cont$importance[2,"CAP2"]*100, 2),
                  "%")) +
  scale_color_manual(values = palette_group, labels = group_names) +
  scale_fill_manual(values = palette_group) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

# 05. Merge plots ----

plot_design <- "
  AAB
  AAC
"

p_beta <- p_rda + 
  p_bray_boxplot +
  p_bray_boxplot_sex +
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

