# ABSTRACT
# This script will run DESeq2 differential abundance
# analysis on the bacterial fraction of the data and 
# generate plots

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(phyloseq)
library(ggplot2)
library(here)
library(patchwork)
library(microViz)
library(DESeq2)
library(RColorBrewer)
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

## Import phyloseq object ----
### Import the non-rarefied data
bracken_data <- read_rds(
  here(
    output_path,
    "bracken_fixed.rds"
  )
) %>%
  ## Filter out low-abundance taxa
  tax_transform("compositional") %>%
  tax_filter(
    min_sample_abundance = 0.01, tax_level = "unique", use_counts = FALSE,
    prev_detection_threshold = 0
  ) %>%
  ps_get(counts = TRUE) %>%
  ## Filter to only include bacteria
  tax_select("Bacteria", ranks_searched = "Rank1", strict_matches = TRUE) %>%
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
                                      "none")),
            sex = factor(sex)) %>%
  ps_select(sex, group) %>%
  ps_mutate(group = relevel(group, ref = "none"))

# 02. Run DESeq2 ----
deseq_data <- phyloseq_to_deseq2(bracken_data, design =  ~ group + sex)

deseq_res <- DESeq(
  deseq_data,
  test = "Wald",
  fitType = "mean"
)

## Function to extract the results from the DESeq object
extract_deseq_results <- function(x, tax_data, name = NULL, pval = NULL) {
  dds_res <- results(
    x,
    name = name,
    cooksCutoff = FALSE
  )
  
  if (!is.null(pval)) {
    dds_res_sig <- dds_res[which(dds_res$padj <= pval), ]
  } else {
    dds_res_sig <- dds_res
  }
  
  dds_res_clean <- cbind(
    as(dds_res_sig, "data.frame"),
    as(tax_data[rownames(dds_res_sig),], 
       "matrix")
  ) %>%
    mutate(species = paste0(Rank6, " ", Rank7)) %>%
    select(-c(Rank3, Rank4, Rank5, Rank7)) %>%
    mutate_at(vars(pvalue, padj, baseMean, log2FoldChange, lfcSE, stat),
              ~ round(., 4)) %>%
    select(Rank1, Rank2, Rank6, species, everything())
  
  return(dds_res_clean)
}

## Apply the function
results_deseq <- lapply(resultsNames(deseq_res), function(x) {
  extract_deseq_results(
    deseq_res, 
    tax_table(bracken_data), 
    name = x
  )
}
)

names(results_deseq) <- resultsNames(deseq_res)

final_results <- bind_rows(results_deseq, .id = "ref") %>%
  filter(ref != "Intercept")

## Filter out non-significant hits
deseq_sign_results <- final_results %>%
  filter(padj < 0.05) %>%
  mutate(abs_lfc = abs(log2FoldChange)) %>%
  arrange(ref, -abs_lfc)


## Create list of results on genera of interest
## from the diversity analyses (RDA)
genera_of_interest <- c(
  "Megasphaera",
  "Megamonas",
  "Bifidobacterium",
  "Clostridium",
  "Dysosmobacter",
  "Flintibacter",
  "Eubacterium",
  "Enterocloster",
  "Flavoniflactor",
  "Faecalibacterium"
)

deseq_genera <- final_results %>%
  filter(Rank6 %in% genera_of_interest)

## Write results to file
write_delim(
  final_results,
  here(
    output_path,
    "10_DESeq_groups.txt"
  ),
  delim = "\t"
)

write_delim(
  deseq_sign_results,
  here(
    output_path,
    "10_DESeq_groups_sign.txt"
  ),
  delim = "\t"
)


write_delim(
  deseq_genera,
  here(
    output_path,
    "10_DESeq_genera.txt"
  ),
  delim = "\t"
)


# 03. Plot results ----

label_vector <- labeller(
  .cols = c(
    "group_monensin_antibiotic_vs_none" = "M+P vs. None",
    "group_monensin_vs_none" = "M vs. None",
    "group_antibiotic_vs_none" = "P vs. None",
    "sex_Male_vs_Female" = "Male vs. Female"
  )
)

col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

p_deseq <- ggplot(deseq_sign_results, aes(log2FoldChange, species, fill = log2FoldChange)) +
  geom_errorbar(aes(xmin = log2FoldChange - lfcSE,
                    xmax = log2FoldChange + lfcSE),
                width = 0.5) +
  geom_point(pch = 21,
             size = 4) +
  scale_fill_gradientn(colors = brewer.pal(9, "Spectral"),
                       guide = guide_colourbar(
                         ticks = TRUE,
                         nbin = 50,
                         barheight = 10,
                         label = TRUE,
                         barwidth = 0.5
                       )) +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic", size = 8),
        panel.grid = element_line(color = col_grid),
        axis.title.y = element_blank()) +
  facet_wrap(~ref, nrow = 1, labeller = label_vector)


ggsave(
  here(
    output_path,
    "10_deseq2.png"
  ),
  p_deseq,
  device = "png",
  units = "cm", 
  dpi = 600,
  height = 25,
  width = 25
)
