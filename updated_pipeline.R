# load arguments
# use bash: Rscript pipeline.R params.yaml

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

#####
# debugging: load the files manually and fix the OSX path to Windows:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
yaml_file <- "params_GRP000002.yaml"
params <- yaml::read_yaml("./params_GRP000002.yaml")
params <- lapply(params, function(x) gsub("/Volumes/biol_bc_picotti_1/", "Y:/", x, fixed = FALSE))
params$output_dir <- paste(params$output_dir, "_ludo", sep = "")
#####

# bug: insert an extra return at the end of the yaml file otherwise get a warning: 
# In readLines(file, warn = readLines.warn) :
# incomplete final line found on 'params_template.yaml'

# Access the parameters
group_id <- params$group_id
input_file <- params$input_file
input_file_tryptic_control <- params$input_file_tryptic_control
experiment_ids <- params$experiment_id
treatment <- params$treatment
ref_condition <- params$ref_condition
comparisons <- params$comparison
output_dir <- params$output_dir


if (!is.null(params$input_file_tryptic_control)) {
  input_file_tryptic_control <- params$input_file_tryptic_control
  tryptic_control_data <- protti::read_protti(input_file_tryptic_control)
} else {
  input_file_tryptic_control <- NULL
  tryptic_control_data <- NULL
}

group_folder_path <- file.path(output_dir, group_id)

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

# bug: adapted LiP Quant to LiP => make a conditional (if) section specific if LiP Quant (aka dose response) dataset
# ------------------------------------------------------------------------------
# LiP 
# ------------------------------------------------------------------------------

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

# bug: removed the following 2 lines: conrep already in those dataframes
# DIA_clean_uniprot$condrep <- paste(DIA_clean_uniprot$r_condition, DIA_clean_uniprot$r_replicate, sep = "_")
# DIA_raw$condrep <- paste(DIA_raw$r_condition, DIA_raw$r_replicate, sep = "_")

# bug: removed the following 2 lines: not used afterwards
# proteins_identified <- uniprot %>%
#   distinct(accession)

# Imputation on precursor level 
# BUG: imputation incomplete since started from the dataframe with incomplete/missing rows => complete first!!
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

write.table(imputed, "./preprocessed_ludo/GRP000002/imputed.tsv", sep = "\t", row.names= FALSE, quote = FALSE)
save.image("./all_data.RData")
# load("./all_data.RData")

# sum up precursors to peptide level and keep only one entry per pep_stripped_sequence
# bug: do not lump to the stripped_peptide_sequence but rather to the modified_peptide_sequence; remove accordingly "pep_stripped_sequemce" from the retained columns
DIA_clean_uniprot_summed_protti <- protti::calculate_protein_abundance(
  imputed,
  sample = r_file_name,
  protein_id = eg_modified_peptide,
  precursor = fg_id,
  peptide = eg_modified_peptide,
  intensity_log2 = normalised_intensity_log2,
  min_n_peptides = 1,
  method = "sum",
  for_plot = FALSE,
  retain_columns = c("pg_protein_accessions",
                     "gene_names", "go_f", "r_condition", "start", "end")
)

dia_clean_file <- file.path(group_folder_path, paste0("dia_clean_uniprot.tsv"))
# bug: changed to write.table to remove quotes and row names
write.csv(DIA_clean_uniprot_summed_protti, dia_clean_file)
write.table(DIA_clean_uniprot_summed_protti, dia_clean_file, sep = "\t", row.names= FALSE, quote = FALSE)

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

# bug: kind of useless: after normalization => flat line?! Consider removing?
# plot_list[[1]] <- protti::qc_median_intensities(
#   data = DIA_raw_norm,
#   sample = r_file_name,
#   grouping = pep_grouping_key,
#   intensity = normalised_intensity_log2,
#   plot = TRUE,
#   interactive = FALSE
# )

# bug: this is not a plot but a table => Not super usefull?! Consider removing?
# plot_list[[2]] <- protti::qc_cvs(
#   data = DIA_clean_uniprot,
#   grouping = fg_id,
#   condition = r_condition,
#   intensity = fg_quantity,
#   plot = FALSE
# )

# bug: Kind of duplicated with the violin plot?! Consider removing?
# plot_list[[3]] <- protti::qc_cvs(
#   data = DIA_clean_uniprot,
#   grouping = fg_id,
#   condition = r_condition,
#   intensity = fg_quantity,
#   plot_style = "density",
#   plot = TRUE
# )

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

# output_qc_pdf <- file.path(group_folder_path, "qc_plots.pdf")
# pdf(output_qc_pdf, width = 9, height = 5)  # Open the PDF device
# # Loop through the list and print each ggplot object
# lapply(plot_list, print)
# dev.off()  

# ------------------------------------------------------------------------------
# Analysis ...
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# differential abundance and go term
# ------------------------------------------------------------------------------

# Data analysis
## Volcano plots peptide level
# Iterate over each comparison and associated experiment ID

# bug: replace pep_stripped_sequence by eg_modified_peptide

for (i in seq_along(comparisons)) {
  comparison_filter <- comparisons[[i]]
  experiment_id <- experiment_ids[[i]]
  
  comparison_parts <- strsplit(comparison_filter, "_vs_")[[1]] 
  print(comparison_parts)
  
  filtered_data <- DIA_clean_uniprot_summed_protti %>%
    dplyr::filter(r_condition %in% comparison_parts)
  
  Volcano_input <- filtered_data %>%
    unique() %>%
    protti::assign_missingness(
      sample = r_file_name,
      condition = r_condition,
      grouping = eg_modified_peptide,
      intensity = normalised_intensity_log2,
      ref_condition = ref_condition,
      retain_columns = all_of(c("pg_protein_accessions", "r_file_name", "r_condition", 
                                "normalised_intensity_log2", 
                                "go_f", "start", "end"))
    )
  
  df_diff_abundance <- protti::calculate_diff_abundance(
    data = Volcano_input,
    sample = r_file_name,
    condition = r_condition,
    grouping = eg_modified_peptide,
    intensity_log2 = normalised_intensity_log2,
    missingness = missingness,
    comparison = comparison,
    method = "moderated_t-test",
    retain_columns = all_of(c("pg_protein_accessions", "eg_modified_peptide", 
                              "comparison", "go_f", "start", "end"))
  )

  # bug: output some important plots to pdf and results to the log file: 
  candidates <- df_diff_abundance %>% 
    dplyr::filter(adj_pval < 0.05 & abs(diff) > 1) %>% 
    left_join(dplyr::distinct(DIA_clean_uniprot, pg_protein_accessions, gene_names), by = "pg_protein_accessions") %>% 
    rowwise() %>% 
    dplyr::mutate(gene = sort(strsplit(gene_names, " ", fixed = TRUE)[[1]])[1]) %>% 
    ungroup() %>% 
    dplyr::mutate(label = paste(gene, "@", start, "-", end, sep = "")) %>% 
    dplyr::mutate(significant = TRUE)
    
  print(paste("Number of significantly changing candidates: ", length(unique(candidates$eg_modified_peptide)), " peptide(s) belonging to ", length(unique(candidates$pg_protein_accessions)), " protein(s)."), sep = "")
  print(distinct(candidates, pg_protein_accessions, eg_modified_peptide, label))
  
  df_diff_abundance <- df_diff_abundance %>%
    left_join(distinct(candidates, eg_modified_peptide, label, significant), by = "eg_modified_peptide") %>% 
    left_join(dplyr::distinct(DIA_clean_uniprot, pg_protein_accessions, gene_names, sequence, length, coverage), by = "pg_protein_accessions") %>% 
    rowwise() %>% 
    dplyr::mutate(gene = sort(strsplit(gene_names, " ", fixed = TRUE)[[1]])[1]) 
  
  # Save Differential Abundance results
  diff_abundance_file <- file.path(
    group_folder_path, 
    paste0("differential_abundance_", experiment_id, "_", comparison_filter, ".csv")
  )
  # bug: changed to write.table to remove quotes and row names
  # write.csv(df_diff_abundance, diff_abundance_file)
  write.table(df_diff_abundance, diff_abundance_file, sep = "\t", row.names= FALSE, quote = FALSE)
  
  # plot the corresponding volcano plots
  plot_list2[[1]] <- volcano_plot(
    data = df_diff_abundance,
    grouping = eg_modified_peptide,
    log2FC = diff,
    significance = pval,
    method = "significant",
    # target_column = pg_protein_accessions,
    # target = "P62942",
    x_axis_label = "log2(fold change) Rapamycin treated vs. untreated",
    significance_cutoff = c(0.05, "adj_pval") 
  ) +
    geom_text_repel(data = df_diff_abundance, size = 3, aes(x = diff, y = -log10(pval), label = label), min.segment.length = unit(0, 'lines'), nudge_y = 0.1)
  
  
  # plot the distribution of the p-values 
  plot_list2[[2]] <- pval_distribution_plot(data = df_diff_abundance,
                         grouping = eg_modified_peptide,
                         pval = pval
  )  
  
  output_stats_pdf <- file.path(group_folder_path, "statistical_analysis_plots.pdf")
  
  ggsave(
    filename = output_stats_pdf, 
    plot = marrangeGrob(plot_list2, nrow=1, ncol=1), 
    width = 8, height = 6
  )
  
  # plot the profile plots for the candidates: 
  candidate_summed <- DIA_clean_uniprot_summed_protti %>% 
    dplyr::filter(pg_protein_accessions %in% sort(unique(candidates$pg_protein_accessions))) %>% 
    left_join(distinct(DIA_clean, r_file_name, condrep), by = "r_file_name") %>% 
    left_join(distinct(candidates, eg_modified_peptide, significant), by = "eg_modified_peptide") %>% 
    dplyr::mutate(significant = ifelse(is.na(significant), FALSE, significant)) %>% 
    dplyr::mutate(line_type = ifelse(significant == TRUE, 'solid', 'dotted'))
      
  plot_list3 <- peptide_profile_plot(
    data = candidate_summed,
    sample = condrep,
    peptide = eg_modified_peptide,
    intensity_log2 = normalised_intensity_log2,
    grouping = pg_protein_accessions,
    targets = sort(unique(candidates$pg_protein_accessions)),
    protein_abundance_plot = FALSE
  ) #+
    # scale_linetype_manual(values=candidate_summed$line_type)
  
  output_profile_pdf <- file.path(group_folder_path, "candidates_profile_plots.pdf")
  
  ggsave(
    filename = output_profile_pdf, 
    plot = marrangeGrob(plot_list3, nrow=1, ncol=1), 
    width = 8, height = 6
  )
  
  plot_list4 <- woods_plot(
    data = df_diff_abundance,
    fold_change = diff,
    start_position = start,
    end_position = end,
    protein_length = length,
    coverage = coverage,
    colouring = adj_pval,
    protein_id = pg_protein_accessions,
    targets = sort(unique(candidates$pg_protein_accessions)), 
    facet = FALSE,
    fold_change_cutoff = 1,
    highlight = significant
  )
  
  output_woods_pdf <- file.path(group_folder_path, "candidates_woods_plots.pdf")
  
  ggsave(
    filename = output_woods_pdf, 
    plot = marrangeGrob(plot_list4, nrow=1, ncol=1), 
    width = 8, height = 6
  )
  
    # Calculate GO term enrichment with error handling
  tryCatch({
    df_diff_abundance_significant <- df_diff_abundance %>%
      dplyr::mutate(significant = ifelse(!is.na(adj_pval) & adj_pval < 0.05, TRUE, FALSE)) 
    
    df_go_term <- protti::calculate_go_enrichment(
      data = df_diff_abundance_significant,
      protein_id = pg_protein_accessions,
      go_annotations_uniprot = go_f,
      is_significant = significant,
      min_n_detected_proteins_in_process = 3,
      plot = FALSE
    )
    
    # Save GO Term enrichment results
    go_term_file <- file.path(group_folder_path, paste0("go_term_", experiment_id, "_", comparison_filter, ".csv"))
    # bug: changed to write.table to remove quotes and row names
    # write.csv(df_go_term, go_term_file)
    write.table(df_go_term, go_term_file, sep = "\t", row.names= FALSE, quote = FALSE)
  }, error = function(e) {
    message(paste("Error in GO term enrichment for comparison", comparison_filter, ":", e))
  })
}

# copy yaml file into the output as well
yaml_file_path <- file.path(group_folder_path, "params.yaml")
file.copy(yaml_file, yaml_file_path)

# Stop redirecting output to the log file
sink()

# bug: stop here in case of LiP analysis only: 

# # ------------------------------------------------------------------------------
# # LiP-Quant Analysis
# # ------------------------------------------------------------------------------
# # Adding dose-response analysis using `parallel_fit_drc_4p`
# 
# # Ensure a numeric concentration column exists for dose-response fitting
# DIA_clean_uniprot_summed_protti$conc_frag <- as.numeric(gsub("_mM", "", DIA_clean_uniprot_summed_protti$r_condition))
# 
# # Perform parallel dose-response fitting
# lipquant_results <- protti::parallel_fit_drc_4p(
#   data = DIA_clean_uniprot_summed_protti,
#   sample = r_file_name,
#   grouping = pep_stripped_sequence,
#   intensity_log2 = normalised_intensity_log2,
#   conc_frag = conc_frag,
#   filter = "post",
#   replicate_completeness = 0.7,
#   condition_completeness = 0.5,
#   correlation_cutoff = 0.8,
#   log_logarithmic = TRUE,
#   retain_columns = c("pg_protein_accessions", "start", "end", "pep_type", "gene_names", "go_f")
# )
# 
# # Save LiP-Quant results
# lipquant_output_file <- file.path(group_folder_path, "lipquant_results.csv")
# write.csv(lipquant_results, lipquant_output_file)
# 
# 
# lipquant_plots <- protti::plot_drc_results(
#   lipquant_results,
#   grouping = pep_stripped_sequence,
#   conc_frag = conc_frag,
#   intensity_log2 = normalised_intensity_log2
# )
# 
# # Save diagnostic plots for dose-response fitting
# lipquant_plot_pdf <- file.path(group_folder_path, "lipquant_dose_response_plots.pdf")
# pdf(lipquant_plot_pdf, width = 7, height = 5)
# lapply(lipquant_plots, print)
# dev.off()
# 
# 
# # ------------------------------------------------------------------------------
# # TRYPTIC CONTROL
# # ------------------------------------------------------------------------------
# 
# # if tryptic control present
# # Tryptic Control (if provided)
# if (!is.null(input_file_tryptic_control)) {
#   tryptic_control_data <- protti::read_protti(input_file_tryptic_control)
#   tryptic_control_data$intensity_log2 <- log2(tryptic_control_data$pg_quantity)
#   
#   # Preprocess tryptic control data
#   tryptic_clean <- tryptic_control_data %>%
#    dplyr::filter(intensity_log2 > 10) %>%
#     protti::normalise(
#       sample = r_file_name,
#       intensity_log2 = intensity_log2,
#       method = "median"
#     )
#   
#   # Fetch UniProt annotations for tryptic control
#   unis_tryptic <- unique(tryptic_clean$pg_protein_accessions)
#   uniprot_tryptic <- protti::fetch_uniprot(unis_tryptic, columns = c("protein_name", "gene_names", "length"))
#   
#   tryptic_clean_uniprot <- tryptic_clean %>%
#     left_join(uniprot_tryptic, by = c("pg_protein_accessions" = "accession"))
#   
#   # Save preprocessed tryptic control data
#   write.csv(tryptic_clean_uniprot, file.path(group_folder_path, "tryptic_control_clean.csv"))
# }
# 
# 
# # copy yaml file into the output as well
# yaml_file_path <- file.path(group_folder_path, "params.yaml")
# file.copy(yaml_file, yaml_file_path)
# 
# # Stop redirecting output to the log file
# sink()
