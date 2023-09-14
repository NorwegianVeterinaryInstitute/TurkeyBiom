# ABSTRACT
# This script presents the number of reads
# and basic QC of the shotgun metagenomic
# dataset

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(funtools)
library(ggplot2)
library(here)
library(readxl)
library(patchwork)

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
  sep = "/"
)

reads_run1 <- read_delim(
  here(
    input_path,
    "R",
    "Data",
    "number_of_reads",
    "runA.txt"
  ),
  delim = "\t"
)

reads_run2 <- read_delim(
  here(
    input_path,
    "R",
    "Data",
    "number_of_reads",
    "runB.txt"
  ),
  delim = "\t"
)

metadata <- read_delim(
  here(
    input_path,
    "R",
    "Data",
    "00_metadata_samples.txt"
  ),
  delim = "\t"
)

dna_ids <- read_xlsx(
  here(
    input_path,
    "Lab",
    "TurkeyBiom_DNA.xlsx"
  )
) %>%
  select(
    Journalnummer,
    `BOKS nr`,
    `BoksKoord.`,
    `Vekt (mg)`
  ) %>%
  rename("saksnr" = Journalnummer,
         "weight" = `Vekt (mg)`) %>%
  mutate(sample = paste0(`BoksKoord.`, "-", `BOKS nr`),
         saksnr = sub("-1$", "", saksnr))

seq_metadata <- read_xlsx(
  here(
    input_path,
    "Lab",
    "sample_table_96_template_TurkeyBiom.xlsx"
  )
) %>%
  slice(-1) %>%
  select(4, 6:8) %>%
  rename(
    "sample" = `Sample Name (only letters, numbers and hyphen allowed, MAX 16 characters)`,
    "concentration" = `Conc. (ng/µl)`,
    "A260_280" = `A260/280`,
    "A260_230" = `A260/230`
    ) %>%
  left_join(dna_ids) %>%
  select(saksnr, sample, weight, concentration, A260_280, A260_230) %>%
  left_join(metadata[,c("saksnr","age_days","sex")])

write_delim(
  seq_metadata,
  here(
    output_path,
    "03_seq_metadata.tsv"
  ),
  delim = "\t"
)

# 02. Merge read counts ----
reads_data <- rbind(
  reads_run1,
  reads_run2
) %>%
  mutate(sample = sub("_S.+_L00.+", "", `Sample Name`),
         lane = sub(".+_(L00.)_R.+", "\\1", `Sample Name`)) %>%
  rename("n_seqs" = `M Seqs`) %>%
  select(sample, lane, n_seqs) %>%
  pivot_wider(names_from = "lane",
              values_from = "n_seqs",
              values_fn = func_paste,
              values_fill = "0") %>%
  mutate_at(vars(contains("L00")),
            ~as.numeric(.)) %>%
  mutate(total_reads = L002 + L003 + L004) %>%
  select(sample, total_reads) %>%
  mutate(sample = sub(".+-([A,B,C,D,E,F,G,H,M,N].+)", "\\1", sample)) %>%
  left_join(seq_metadata)

plot_data <- filter(
  reads_data,
  !is.na(sex)
)
  

# 03. Run stats ----
test1 <- cor.test(plot_data$total_reads,
                  plot_data$concentration,
                  alternative = "two.sided",
                  method = "spearman")

test2 <- cor.test(plot_data$total_reads,
                  plot_data$A260_280,
                  alternative = "two.sided",
                  method = "spearman")

test3 <- cor.test(plot_data$total_reads,
                  plot_data$weight,
                  alternative = "two.sided",
                  method = "spearman")

test4 <- cor.test(plot_data$total_reads,
                  plot_data$A260_230,
                  alternative = "two.sided",
                  method = "spearman")

# 04. Plot data ----
palette <- c("Male" = "#80b1d3",
             "Female" = "#fb8072")

p1 <- ggplot(plot_data, aes(total_reads, concentration)) +
  geom_point(aes(color = sex)) +
  geom_smooth(method = "lm",
              fullrange = TRUE,
              color = "black",
              se = FALSE,
              alpha = 0.2,
              linewidth = 0.2) +
  ggtitle(paste0("Rho: ",
                 round(test1$estimate, 3),
                 ", p = ",
                 round(test1$p.value, 3))) +
  scale_color_manual(values = palette) +
  labs(x = "Millon reads",
       y = "DNA concentration (ng/µl)",
       color = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank())

p2 <- ggplot(plot_data, aes(total_reads, A260_280)) +
  geom_point(aes(color = sex)) +
  geom_smooth(method = "lm",
              fullrange = TRUE,
              color = "black",
              se = FALSE,
              alpha = 0.2,
              linewidth = 0.2) +
  ggtitle(paste0("Rho: ",
                 round(test2$estimate, 3),
                 ", p = ",
                 round(test2$p.value, 3))) +
  scale_color_manual(values = palette) +
  labs(x = "Millon reads",
       y = "DNA concentration (ng/µl)",
       color = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank())

p3 <- ggplot(plot_data, aes(total_reads, weight)) +
  geom_point(aes(color = sex)) +
  geom_smooth(method = "lm",
              fullrange = TRUE,
              color = "black",
              se = FALSE,
              alpha = 0.2,
              linewidth = 0.2) +
  ggtitle(paste0("Rho: ",
                 round(test3$estimate, 3),
                 ", p = ",
                 round(test3$p.value, 3))) +
  scale_color_manual(values = palette) +
  labs(x = "Millon reads",
       y = "DNA concentration (ng/µl)",
       color = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank())

p4 <- ggplot(plot_data, aes(total_reads, A260_230)) +
  geom_point(aes(color = sex)) +
  geom_smooth(method = "lm",
              fullrange = TRUE,
              color = "black",
              se = FALSE,
              alpha = 0.2,
              linewidth = 0.2) +
  ggtitle(paste0("Rho: ",
                 round(test4$estimate, 3),
                 ", p = ",
                 round(test4$p.value, 3))) +
  scale_color_manual(values = palette) +
  labs(x = "Millon reads",
       y = "DNA concentration (ng/µl)",
       color = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank())

p_all <- p1 + p2 + p3 + p4 +
  plot_layout(guides = "collect")

ggsave(
  here(
    output_path,
    "03_reads_and_metadata_corr.png"
  ),
  p_all,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 22,
  width = 25
)

median_values <- plot_data %>%
  select(sex, total_reads) %>%
  group_by(sex) %>%
  summarise(
    median = median(total_reads),
    mean = mean(total_reads)
  ) %>%
  ungroup()

p_nreads <- ggplot(plot_data, aes(
  reorder(sample, -total_reads), total_reads, fill = sex)
  ) +
  geom_col(color = "black") +
  geom_hline(aes(yintercept = median),
             data = median_values) +
  labs(x = "Samples",
       y = "Million reads") +
  scale_fill_manual(values = palette) +
  facet_grid(~sex, scales = "free_x", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

p_reads_sex <- ggplot(plot_data, aes(sex, total_reads, fill = sex)) +
  stat_boxplot(geom = "errorbar", width = 0.4) +
  geom_boxplot() +
  labs(x = "Sex",
       y = "Million reads") +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank())

p_reads_all <- p_nreads / (p_reads_sex | plot_spacer()) +
  plot_layout()

ggsave(
  here(
    output_path,
    "03_number_of_reads_nsc.png"
  ),
  p_reads_all,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 18,
  width = 22
)
