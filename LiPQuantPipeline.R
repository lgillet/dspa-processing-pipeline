
# LiPQuant pipeline
# not finished 


library(yaml)
library(protti)
library(tidyr)
library(dplyr)
library(missForest)
library(doParallel)
library(ggplot2)
library(gridExtra)
library(ggrepel)

registerDoParallel(cores = 10) 

# Read the YAML file
args <- commandArgs(trailingOnly = TRUE)
# Expect the first argument to be the path to the YAML file
yaml_file <- args[1]
params <- yaml::read_yaml(yaml_file)

group_id <- params$group_id
input_file <- params$input_file
experiment_ids <- params$experiment_id
treatment <- params$treatment
ref_condition <- params$ref_condition
comparisons <- params$comparison
output_dir <- params$output_dir

roup_folder_path <- file.path(output_dir, group_id)

# bug: adapted "./preprocessed" to output_dir
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(group_folder_path)) {
  dir.create(group_folder_path, recursive = TRUE)
  cat("Created group folder:", group_folder_path, "\n")
}

# Log outputs to check for errors 
# Redirect output to a log file
logfile_dir <- file.path(output_dir, group_id, "processing_log2.txt")
sink(logfile_dir, append = TRUE)

# Print or use the parameters
print(paste("Input File:", input_file))
print(paste("Treatment:", treatment))
print(paste("Ref Condtion:", ref_condition))
print(paste("Output Directory:", output_dir))

# Log R session information, including loaded packages
cat("Logging R session info:\n")
sessionInfo() 

# Create a list to store the plots
plot_list <- list()
plot_list2 <- list()
plot_list3 <- list()
plot_list4 <- list()

# load file
DIA_raw <- read_protti(input_file)

# ------------------------------------------------------------------------------
# Preprocessing 
# ------------------------------------------------------------------------------

DIA_raw$intensity_log2 <- log2(DIA_raw$fg_ms2raw_quantity)
DIA_raw$condrep <- paste(DIA_raw$r_condition, DIA_raw$r_replicate, sep = "_")

DIA_raw_norm <- protti::normalise(
  DIA_raw,
  sample = r_file_name,
  intensity_log2 = intensity_log2,
  method = "median"
)

# bug: included parsing of the protein accession when more than 1 to return the first one alphabetically (to help fetch_uniprot) 
DIA_clean <- DIA_raw_norm %>%
  dplyr::filter(intensity_log2 > 10) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(pg_protein_accessions2 = ifelse(base::grepl(";", pg_protein_accessions, fixed = FALSE), base::sort(base::strsplit(pg_protein_accessions, ";", fixed = TRUE)[[1]])[1], pg_protein_accessions)) %>% 
  dplyr::ungroup()

DIA_clean$fg_id <- paste0(DIA_clean$fg_labeled_sequence, DIA_clean$fg_charge)
unis <- unique(DIA_clean$pg_protein_accessions2) # make vector for fetch_uniprot

## Load data from uniprot and join with DIA dataframe
uniprot <-
  protti::fetch_uniprot(
    unis,
    columns =  c(
      "protein_name",
      "gene_names",
      "length",
      "sequence",
      "xref_pdb",
      "go_f",
      "go_p",
      "go_c"
    )
  ) 

DIA_clean_uniprot <- DIA_clean %>%
  left_join(uniprot, by = c("pg_protein_accessions2" = "accession")) %>% 
  find_peptide(sequence, pep_stripped_sequence) %>%
  assign_peptide_type(aa_before, last_aa, aa_after) %>%
  distinct() %>% 
  calculate_sequence_coverage(protein_sequence = sequence, peptides = pep_stripped_sequence) %>% 
  dplyr::mutate(normalised_intensity = 2^normalised_intensity_log2)

DIA_clean_uniprot_complete <- DIA_clean_uniprot %>% 
  distinct(r_file_name, fg_id, normalised_intensity_log2, eg_modified_peptide, pep_stripped_sequence, pg_protein_accessions, gene_names, go_f, r_condition, start, end) %>% 
  tidyr::complete(nesting(r_file_name, r_condition), nesting(pg_protein_accessions, gene_names, go_f, fg_id, eg_modified_peptide, pep_stripped_sequence, start, end))

imputed <- impute_randomforest(
  DIA_clean_uniprot_complete,
  sample = r_file_name,
  grouping = fg_id,
  intensity_log2 = normalised_intensity_log2,
  retain_columns = c("eg_modified_peptide", "pep_stripped_sequence", "pg_protein_accessions", 
                     "gene_names", "go_f", "r_condition", "start", "end"),
  parallelize = "variables"
)
# ------------------------------------------------------------------------------
# QC plots
# ------------------------------------------------------------------------------

# bug: for QCs, consider plotting first the plots in the order of processing: DIA_raw first, then DIA_raw_norm and then DIA_clean_uniprot!
# bug: relabeled the plots accordingly:

plot_list[[1]] <- protti::qc_ids(
  data = DIA_raw, 
  sample = condrep, 
  grouping = pep_grouping_key, 
  condition = r_condition, 
  intensity = fg_quantity
)+ ggtitle('Precursor ID count per sample')

# bug: added a plot for the number of protein groups
plot_list[[2]] <- protti::qc_ids(
  data = DIA_raw, 
  sample = condrep, 
  grouping = pg_protein_accessions, 
  condition = r_condition, 
  intensity = pg_quantity
)+ ggtitle('Protein ID count per sample')


plot_list[[3]] <- protti::qc_intensity_distribution(
  data = DIA_raw_norm,
  sample = condrep,
  grouping = pep_grouping_key,
  intensity_log2 = intensity_log2,
  plot_style = "histogram"
) + ggtitle('Overall log2 Intensity distribution before normalisation')

plot_list[[4]] <- protti::qc_intensity_distribution(
  data = DIA_raw_norm,
  sample = condrep,
  grouping = pep_grouping_key,
  intensity_log2 = normalised_intensity_log2,
  plot_style = "boxplot"
) + ggtitle('Run intensities after normalisation')

plot_list[[5]] <- protti::qc_cvs(
  data = DIA_clean_uniprot,
  grouping = fg_id,
  condition = r_condition,
  intensity = normalised_intensity,
  plot_style = "violin",
  plot = TRUE
)


## Missed cleavages
# bug: Swaped the count and intensity to have the plot in the same order as missed cleavages
plot_list[[6]] <- protti::qc_missed_cleavages(
  data = DIA_clean_uniprot,
  sample = condrep,
  grouping = fg_id,
  missed_cleavages = pep_nr_of_missed_cleavages,
  intensity = normalised_intensity,
  method = "count",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[7]] <- protti::qc_missed_cleavages(
  data = DIA_clean_uniprot,
  sample = condrep,
  grouping = fg_id,
  missed_cleavages = pep_nr_of_missed_cleavages,
  intensity = normalised_intensity,
  method = "intensity",
  plot = TRUE,
  interactive = FALSE
)


plot_list[[8]] <- protti::qc_peptide_type(
  DIA_clean_uniprot,
  condrep,
  fg_id,
  pep_type,
  intensity = normalised_intensity,
  method = "count",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[9]] <- protti::qc_peptide_type(
  DIA_clean_uniprot,
  condrep,
  fg_id,
  pep_type,
  intensity = normalised_intensity,
  method = "intensity",
  plot = TRUE,
  interactive = FALSE
)

## Principal component analysis (PCA)
plot_list[[10]] <-DIA_clean_uniprot %>%
  protti::qc_pca(
    sample = condrep, 
    grouping = pep_grouping_key, 
    intensity = intensity_log2, 
    condition = r_condition
  )

## corelation_map
# bug: replaced r_filename with r_condition for consistency reasons with other plots
# bug: added [[4]] to get directly the grob from the pheatmap
plot_list[[11]] <- protti::qc_sample_correlation(
  data = DIA_clean_uniprot,
  sample = condrep,
  grouping = fg_id,
  intensity_log2 = intensity_log2,
  condition = r_condition
)[[4]]

# bug: added the QC plot to check the intensities before and after imputation
plot_list[[12]] <- imputed %>%
  dplyr::rename(imputed_intensity_log2 = normalised_intensity_log2) %>% 
  left_join(distinct(DIA_clean_uniprot_complete, r_file_name, fg_id, normalised_intensity_log2), by = c("r_file_name", "fg_id")) %>% 
  dplyr::select(imputed_intensity_log2, normalised_intensity_log2) %>%
  pivot_longer(cols = everything(),
               names_to = "imputed",
               values_to = "intensity") %>%
  dplyr::filter(!is.na(intensity)) %>%
  dplyr::mutate(imputed = factor(imputed, levels = c("imputed_intensity_log2", "normalised_intensity_log2"))) %>% 
  ggplot(aes(intensity, fill = imputed)) +
  labs(title = "Histogram of intensities before and after imputation (Log2)",
       x = "Log2 Intensity",
       y = "Frequency",
       fill = "Type") +
  geom_histogram(
    binwidth = 0.5,
    color = "black",
    position = "identity"
  ) +
  scale_fill_manual(values = protti_colours[c(2,1)]) +
  theme_bw() + 
  coord_cartesian(xlim = c(5, 30))

# Save QC plots
# bug: changed to ggsave for cleaner output
output_qc_pdf <- file.path(group_folder_path, "qc_plots.pdf")
ggsave(
  filename = output_qc_pdf, 
  plot = marrangeGrob(plot_list, nrow=1, ncol=1), 
  width = 8, height = 6
)

# ------------------------------------------------------------------------------
# LiP-Quant Analysis
# ------------------------------------------------------------------------------
# Adding dose-response analysis using `parallel_fit_drc_4p`

# Ensure a numeric concentration column exists for dose-response fitting
DIA_clean_uniprot_summed_protti$conc_frag <- as.numeric(gsub("_mM", "", DIA_clean_uniprot_summed_protti$r_condition))

# Perform parallel dose-response fitting
lipquant_results <- protti::parallel_fit_drc_4p(
  data = DIA_clean_uniprot_summed_protti,
  sample = r_file_name,
  grouping = pep_stripped_sequence,
  intensity_log2 = normalised_intensity_log2,
  conc_frag = conc_frag,
  filter = "post",
  replicate_completeness = 0.7,
  condition_completeness = 0.5,
  correlation_cutoff = 0.8,
  log_logarithmic = TRUE,
  retain_columns = c("pg_protein_accessions", "start", "end", "pep_type", "gene_names", "go_f")
)

# Save LiP-Quant results
lipquant_output_file <- file.path(group_folder_path, "lipquant_results.csv")
write.csv(lipquant_results, lipquant_output_file)


lipquant_plots <- protti::plot_drc_results(
  lipquant_results,
  grouping = pep_stripped_sequence,
  conc_frag = conc_frag,
  intensity_log2 = normalised_intensity_log2
)

# Save diagnostic plots for dose-response fitting
lipquant_plot_pdf <- file.path(group_folder_path, "lipquant_dose_response_plots.pdf")
pdf(lipquant_plot_pdf, width = 7, height = 5)
lapply(lipquant_plots, print)
dev.off()


# copy yaml file into the output as well
yaml_file_path <- file.path(group_folder_path, "params.yaml")
file.copy(yaml_file, yaml_file_path)

# Stop redirecting output to the log file
sink()
