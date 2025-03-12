# load arguments
# use bash: Rscript pipeline.R params.yaml

# LiP MS pipeline
# with optional Tryptic control

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
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# yaml_file <- "PXD022297_Valentina_Osmotic.yaml"
# params <- yaml::read_yaml("Y:/LiP-Atlas/preprocessed/20250124_142925_PXD022297_Valentina_Osmotic_LiP/PXD022297_Valentina_Osmotic.yaml")
# params <- lapply(params, function(x) gsub("G:/Biognosys/Spectronaut/results/", "Y:/LiP-Atlas/preprocessed/", x, fixed = FALSE))
# params$output_dir <- paste(params$output_dir, "_ValeOmso_Ludo", sep = "")
#####

# Access the parameters
group_id <- params$group_id
input_file <- params$input_file
input_file_tryptic_control <- params$input_file_tryptic_control
experiment_ids <- params$experiment_id
treatment <- params$treatment
ref_condition <- params$ref_condition
comparisons <- params$comparison
output_dir <- params$output_dir

# BUG: is group_id and experiment_id in the yaml correct?? GRP00001 and LiP00001 resp. 

group_folder_path <- file.path(output_dir, group_id)

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
plot_list5 <- list()
plot_list6 <- list()
plot_list7 <- list()

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

# imputed <- read.csv("C:/Users/lgillet/Documents/IMSB2/github_lgillet/dspa-processing-pipeline/preprocessed_ValeOmso_Ludo/GRP000001/imputed.tsv", sep = "\t", header = TRUE)

imputed <- impute_randomforest(
  DIA_clean_uniprot_complete,
  sample = r_file_name,
  grouping = fg_id,
  intensity_log2 = normalised_intensity_log2,
  retain_columns = c("eg_modified_peptide", "pep_stripped_sequence", "pg_protein_accessions",
                     "gene_names", "go_f", "r_condition", "start", "end"),
  parallelize = "variables"
)

imputed_file <- file.path(group_folder_path, paste0("imputed.tsv"))
write.table(imputed, imputed_file, sep = "\t", row.names= FALSE, quote = FALSE)
save.image("./all_data.RData")

# sum up precursors to peptide level and keep only one entry per pep_stripped_sequence
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
write.table(DIA_clean_uniprot_summed_protti, dia_clean_file, sep = "\t", row.names= FALSE, quote = FALSE)

# ------------------------------------------------------------------------------
# QC plots
# ------------------------------------------------------------------------------


plot_list[[1]] <- protti::qc_ids(
  data = DIA_raw, 
  sample = condrep, 
  grouping = pep_grouping_key, 
  condition = r_condition, 
  intensity = fg_quantity
)+ ggtitle('Precursor ID count per sample')

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

## correlation_map
plot_list[[11]] <- protti::qc_sample_correlation(
  data = DIA_clean_uniprot,
  sample = condrep,
  grouping = fg_id,
  intensity_log2 = intensity_log2,
  condition = r_condition
)[[4]]

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
output_qc_pdf <- file.path(group_folder_path, "qc_plots.pdf")
ggsave(
  filename = output_qc_pdf, 
  plot = marrangeGrob(plot_list, nrow=1, ncol=1), 
  width = 8, height = 6
)

# bug: if TrP also present then make another set of QC plots for the TrP

# ------------------------------------------------------------------------------
# TRYPTIC CONTROL
# ------------------------------------------------------------------------------

# Tryptic Control (if provided)
if (!is.null(input_file_tryptic_control)) {
  
  input_file_tryptic_control <- params$input_file_tryptic_control
  tryptic_control_data <- protti::read_protti(input_file_tryptic_control)
  # BUG: do the calculations (normalisation etc. at the fg_id/fg_quantity level, NOT at the pg_quantity!!)
  tryptic_control_data$intensity_log2 <- log2(tryptic_control_data$fg_ms2raw_quantity)
  tryptic_control_data$condrep <- paste(tryptic_control_data$r_condition, tryptic_control_data$r_replicate, sep = "_")
  
  # Preprocess tryptic control data  
  tryptic_norm <- protti::normalise(
    tryptic_control_data,
    sample = r_file_name,
    intensity_log2 = intensity_log2,
    method = "median"
  )
  
  tryptic_clean <- tryptic_norm %>%
    dplyr::filter(intensity_log2 > 10) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(pg_protein_accessions2 = ifelse(base::grepl(";", pg_protein_accessions, fixed = FALSE), base::sort(base::strsplit(pg_protein_accessions, ";", fixed = TRUE)[[1]])[1], pg_protein_accessions)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(normalised_intensity = 2^normalised_intensity_log2)

  tryptic_clean$fg_id <- paste0(tryptic_clean$fg_labeled_sequence, tryptic_clean$fg_charge)
  
  # Fetch UniProt annotations for tryptic control
  unis_tryptic <- unique(tryptic_clean$pg_protein_accessions2)

  ## Load data from uniprot and join with DIA dataframe
  uniprot_tryptic <-
    protti::fetch_uniprot(
      unis_tryptic,
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
  
  tryptic_clean_uniprot <- tryptic_clean %>%
    left_join(uniprot_tryptic, by = c("pg_protein_accessions2" = "accession")) %>% 
    find_peptide(sequence, pep_stripped_sequence) %>%
    assign_peptide_type(aa_before, last_aa, aa_after) %>%
    distinct() %>% 
    calculate_sequence_coverage(protein_sequence = sequence, peptides = pep_stripped_sequence)
  
  tryptic_clean_uniprot_complete <- tryptic_clean_uniprot %>% 
    distinct(r_file_name, fg_id, normalised_intensity_log2, eg_modified_peptide, pep_stripped_sequence, pg_protein_accessions, pg_protein_accessions2, gene_names, go_f, r_condition, start, end) %>% 
    tidyr::complete(nesting(r_file_name, r_condition), nesting(pg_protein_accessions, pg_protein_accessions2, gene_names, go_f, fg_id, eg_modified_peptide, pep_stripped_sequence, start, end))
  
  # Save preprocessed tryptic control data
  write.table(tryptic_clean_uniprot_complete, file.path(group_folder_path, "tryptic_control_clean.tsv"), sep = "\t", row.names= FALSE, quote = FALSE)

  
  # calculate protein abundance
  tryptic_protein_abundance <- calculate_protein_abundance(
    tryptic_clean_uniprot,
    sample = r_file_name,
    protein_id = pg_protein_accessions,
    precursor = fg_id,
    peptide = eg_modified_peptide,
    intensity_log2 = normalised_intensity_log2,
    method = "iq",
    for_plot = FALSE,
    retain_columns = c("pg_protein_accessions",
                       "gene_names", "go_f", "r_condition", "condrep", "pg_protein_accessions2")
  )
  
  
  tryptic_clean_file <- file.path(group_folder_path, paste0("tryptic_uniprot.tsv"))

  write.table(tryptic_protein_abundance, tryptic_clean_file, sep = "\t", row.names= FALSE, quote = FALSE)
  
  
  
  # make the QC plots as for LiP
  plot_list2[[1]] <- protti::qc_ids(
    data = tryptic_control_data, 
    sample = condrep, 
    grouping = pep_grouping_key, 
    condition = r_condition, 
    intensity = fg_quantity
  )+ ggtitle('Tryptic control: Precursor ID count per sample')
  
  plot_list2[[2]] <- protti::qc_ids(
    data = tryptic_control_data, 
    sample = condrep, 
    grouping = pg_protein_accessions, 
    condition = r_condition, 
    intensity = pg_quantity
  )+ ggtitle('Tryptic control: Protein ID count per sample')
  
  plot_list2[[3]] <- protti::qc_intensity_distribution(
    data = tryptic_norm,
    sample = condrep,
    grouping = fg_id,
    intensity_log2 = intensity_log2,
    plot_style = "histogram"
  ) + ggtitle('Tryptic control: Overall log2 Intensity distribution before normalisation')
  
  plot_list2[[4]] <- protti::qc_intensity_distribution(
    data = tryptic_norm,
    sample = condrep,
    grouping = pep_grouping_key,
    intensity_log2 = normalised_intensity_log2,
    plot_style = "boxplot"
  ) + ggtitle('Tryptic control: Run intensities after normalisation')
  
  plot_list2[[5]] <- protti::qc_cvs(
    data = tryptic_clean_uniprot,
    grouping = fg_id,
    condition = r_condition,
    intensity = normalised_intensity,
    plot_style = "violin",
    plot = TRUE
  ) + ggtitle('Tryptic control: Violon plot')
  
  ## Missed cleavages
  plot_list2[[6]] <- protti::qc_missed_cleavages(
    data = tryptic_clean_uniprot,
    sample = condrep,
    grouping = fg_id,
    missed_cleavages = pep_nr_of_missed_cleavages,
    intensity = normalised_intensity,
    method = "count",
    plot = TRUE,
    interactive = FALSE
  ) + ggtitle('Tryptic control: number of missed cleavages')
  
  plot_list2[[7]] <- protti::qc_missed_cleavages(
    data = tryptic_clean_uniprot,
    sample = condrep,
    grouping = fg_id,
    missed_cleavages = pep_nr_of_missed_cleavages,
    intensity = normalised_intensity,
    method = "intensity",
    plot = TRUE,
    interactive = FALSE
  ) + ggtitle('Tryptic control: intensity of missed cleavages')
  
  plot_list2[[8]] <- protti::qc_peptide_type(
    tryptic_clean_uniprot,
    condrep,
    fg_id,
    pep_type,
    intensity = normalised_intensity,
    method = "count",
    plot = TRUE,
    interactive = FALSE
  ) + ggtitle('Tryptic control: number of half tryptics')
  
  plot_list2[[9]] <- protti::qc_peptide_type(
    tryptic_clean_uniprot,
    condrep,
    fg_id,
    pep_type,
    intensity = normalised_intensity,
    method = "intensity",
    plot = TRUE,
    interactive = FALSE
  ) + ggtitle('Tryptic control: intensity of half tryptics')
  
  ## Principal component analysis (PCA)
  plot_list2[[10]] <- tryptic_protein_abundance %>%
    protti::qc_pca(
      sample = condrep, 
      grouping = pg_protein_accessions, 
      intensity = normalised_intensity_log2, 
      condition = r_condition
    ) + ggtitle('Tryptic control: PCA at protein level')
  
  ## correlation_map
  plot_list2[[11]] <- protti::qc_sample_correlation(
    data = tryptic_protein_abundance,
    sample = condrep,
    grouping = pg_protein_accessions,
    intensity_log2 = normalised_intensity_log2,
    condition = r_condition
  )[[4]] 
  
  # Save QC plots
  output_tryptic_qc_pdf <- file.path(group_folder_path, "tryptic_controls_qc_plots.pdf")
  ggsave(
    filename = output_tryptic_qc_pdf, 
    plot = marrangeGrob(plot_list2, nrow=1, ncol=1), 
    width = 8, height = 6
  )
  
}

  
i <- 1

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
    method = "t-test",
    retain_columns = all_of(c("pg_protein_accessions", "eg_modified_peptide", 
                              "comparison", "go_f", "start", "end"))
  )

  
  if (!is.null(input_file_tryptic_control)) {
    diff_trp <- tryptic_protein_abundance %>%
      assign_missingness(
        sample = r_file_name,
        condition = r_condition,
        intensity = normalised_intensity_log2,
        grouping = pg_protein_accessions,
        ref_condition = "CTR_Trp",
        retain_columns = all_of(c("pg_protein_accessions2"))
      ) %>%
      calculate_diff_abundance(
        sample = r_file_name,
        condition = r_condition,
        grouping = pg_protein_accessions,
        intensity_log2 = normalised_intensity_log2,
        comparison = comparison,
        method = "t-test",
        retain_columns = all_of(c("pg_protein_accessions2"))
      )
    
    Trp_candidates <- diff_trp %>% 
      dplyr::filter(adj_pval < 0.05 & abs(diff) > 1) %>% 
      left_join(dplyr::distinct(tryptic_protein_abundance, pg_protein_accessions2, gene_names), by = "pg_protein_accessions2") %>% 
      rowwise() %>% 
      dplyr::mutate(gene = sort(strsplit(gene_names, " ", fixed = TRUE)[[1]])[1]) %>% 
      ungroup() %>% 
      dplyr::mutate(significant = TRUE)
    
    print(paste("Tryptic control: Number of significantly changing candidates: ", length(unique(Trp_candidates$pg_protein_accessions)), " protein(s)."), sep = "")
    print(distinct(Trp_candidates, gene))
    
    diff_trp <- diff_trp %>%
      left_join(distinct(Trp_candidates, pg_protein_accessions2, significant, gene), by = "pg_protein_accessions2") %>% 
      left_join(dplyr::distinct(tryptic_protein_abundance, pg_protein_accessions2, gene_names), by = "pg_protein_accessions2")

    
    # BUG: did the yaml parameter changed? 
    # before: experiment_id => experiment_ids
    # before: comparison_filter => comparisons
    
    # Save Differential Abundance results
    diff_trp_file <- file.path(
      group_folder_path, 
      paste0("trp_differential_abundance_", experiment_ids, "_", comparisons, ".tsv")
    )
    
    write.table(diff_trp, diff_trp_file, sep = "\t", row.names= FALSE, quote = FALSE)

    # plot the corresponding volcano plots
    plot_list3[[1]] <- volcano_plot(
      data = diff_trp,
      grouping = pg_protein_accessions2,
      log2FC = diff,
      significance = pval,
      method = "significant",
      x_axis_label = "log2(fold change) treated vs. untreated",
      significance_cutoff = c(0.05, "adj_pval") 
    ) +
      geom_text_repel(data = diff_trp, 
                      size = 3, 
                      aes(x = diff, y = -log10(pval), label = gene), 
                      min.segment.length = unit(0, 'lines'), nudge_y = 0.1)
    
    
    # plot the distribution of the p-values 
    plot_list3[[2]] <- pval_distribution_plot(data = diff_trp,
                                              grouping = pg_protein_accessions2,
                                              pval = pval
    )  
    
    trp_output_stats_pdf <- file.path(group_folder_path, "Trp_statistical_analysis_plots.pdf")
    
    ggsave(
      filename = trp_output_stats_pdf, 
      plot = marrangeGrob(plot_list3, nrow=1, ncol=1), 
      width = 8, height = 6
    )
    
    # plot the profile plots for the candidates: 
    trp_candidate_summed <- tryptic_clean_uniprot %>% 
      dplyr::filter(pg_protein_accessions2 %in% sort(unique(Trp_candidates$pg_protein_accessions2))) %>% 
      left_join(distinct(Trp_candidates, pg_protein_accessions2, gene, significant), by = "pg_protein_accessions2") %>% 
      dplyr::mutate(significant = ifelse(is.na(significant), FALSE, significant)) %>% 
      dplyr::mutate(line_type = ifelse(significant == TRUE, 'solid', 'dotted'))
    
    plot_list4 <- peptide_profile_plot(
      data = trp_candidate_summed,
      sample = condrep,
      peptide = fg_id,
      intensity_log2 = normalised_intensity_log2,
      grouping = gene,
      targets = sort(unique(Trp_candidates$gene)),
      protein_abundance_plot = FALSE
    ) 
    
    output_profile_pdf <- file.path(group_folder_path, "Trp_candidates_profile_plots.pdf")
    
    ggsave(
      filename = output_profile_pdf, 
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
      go_term_file <- file.path(group_folder_path, paste0("go_term_", experiment_id, "_", comparison_filter, ".tsv"))
      write.table(df_go_term, go_term_file, sep = "\t", row.names= FALSE, quote = FALSE)
    }, error = function(e) {
      message(paste("Error in GO term enrichment for comparison", comparison_filter, ":", e))
    })
    
    # perform TrP protein correction on LiP:
    diff_trp2 <- diff_trp %>% 
      rowwise() %>% 
      dplyr::mutate(comparison = gsub("_Trp", "_LiP", x = comparison, fixed = TRUE)) %>% 
      ungroup()
    
    corrected <- correct_lip_for_abundance(
      lip_data = df_diff_abundance,
      trp_data = diff_trp2,
      protein_id = pg_protein_accessions,
      grouping = eg_modified_peptide,
      comparison = comparison, 
      diff = diff,
      n_obs = n_obs,
      std_error = std_error,
      p_adj_method = "BH",
      retain_columns = c("missingness"),
      method = "satterthwaite"
    )

    # bug: output some important plots to pdf and results to the log file: 
    candidates <- corrected %>% 
      dplyr::filter(adj_pval < 0.05 & abs(adj_diff) > 1) %>% 
      left_join(dplyr::distinct(DIA_clean_uniprot, eg_modified_peptide, pg_protein_accessions2, gene_names, coverage, length, start, end), by = "eg_modified_peptide") %>% 
      rowwise() %>% 
      dplyr::mutate(gene = sort(strsplit(gene_names, " ", fixed = TRUE)[[1]])[1]) %>% 
      ungroup() %>% 
      dplyr::mutate(label = paste(gene, "@", start, "-", end, sep = "")) %>% 
      dplyr::mutate(significant = TRUE)
      
    print(paste("Number of significantly changing candidates: ", length(unique(candidates$eg_modified_peptide)), " peptide(s) belonging to ", length(unique(candidates$pg_protein_accessions)), " protein(s)."), sep = "")
    print(distinct(candidates, pg_protein_accessions, eg_modified_peptide, label))
    
    corrected <- corrected %>%
      left_join(distinct(candidates, eg_modified_peptide, start, end, label, coverage, length, significant), by = "eg_modified_peptide") 
      
    # Save Differential Abundance results
    corrected_file <- file.path(
      group_folder_path, 
      paste0("differential_abundance_", experiment_id, "_", comparison_filter, ".tsv")
    )

    write.table(corrected, corrected_file, sep = "\t", row.names= FALSE, quote = FALSE)
  
    # plot the corresponding volcano plots
    plot_list5[[1]] <- volcano_plot(
      data = corrected,
      grouping = eg_modified_peptide,
      log2FC = adj_diff,
      significance = pval,
      method = "significant",
      x_axis_label = "log2(fold change) treated vs. untreated",
      significance_cutoff = c(0.05, "adj_pval") 
    ) +
      geom_text_repel(data = corrected, 
                      size = 3, 
                      aes(x = adj_diff, y = -log10(pval), label = label), 
                      min.segment.length = unit(0, 'lines'), nudge_y = 0.1)
  
  
    # plot the distribution of the p-values 
    plot_list5[[2]] <- pval_distribution_plot(data = corrected,
                           grouping = eg_modified_peptide,
                           pval = pval
    )  
  
    output_stats_pdf <- file.path(group_folder_path, "LiP_corrected_statistical_analysis_plots.pdf")
  
    ggsave(
      filename = output_stats_pdf, 
      plot = marrangeGrob(plot_list5, nrow=1, ncol=1), 
      width = 8, height = 6
    )
  
    # plot the profile plots for the candidates: 
    candidate_summed <- DIA_clean_uniprot_summed_protti %>% 
      dplyr::filter(pg_protein_accessions %in% sort(unique(candidates$pg_protein_accessions))) %>% 
      left_join(distinct(DIA_clean, r_file_name, condrep), by = "r_file_name") %>% 
      left_join(distinct(candidates, eg_modified_peptide, significant), by = "eg_modified_peptide") %>% 
      dplyr::mutate(significant = ifelse(is.na(significant), FALSE, significant)) %>% 
      dplyr::mutate(eg_modified_peptide2 = ifelse(significant == TRUE, paste("*", eg_modified_peptide, sep = ""), eg_modified_peptide))
      
    plot_list6 <- peptide_profile_plot(
      data = candidate_summed,
      sample = condrep,
      peptide = eg_modified_peptide2,
      intensity_log2 = normalised_intensity_log2,
      grouping = pg_protein_accessions,
      targets = sort(unique(candidates$pg_protein_accessions)),
      protein_abundance_plot = FALSE
    ) 
  
    output_profile_pdf <- file.path(group_folder_path, "LiP_corrected_candidates_profile_plots.pdf")
    
    ggsave(
      filename = output_profile_pdf, 
      plot = marrangeGrob(plot_list6, nrow=1, ncol=1), 
      width = 8, height = 6
    )
  
    plot_list7 <- woods_plot(
      data = corrected,
      fold_change = adj_diff,
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
  
    output_woods_pdf <- file.path(group_folder_path, "LiP_corrected_candidates_woods_plots.pdf")
  
    ggsave(
      filename = output_woods_pdf, 
      plot = marrangeGrob(plot_list7, nrow=1, ncol=1), 
      width = 8, height = 6
    )
  
    # Calculate GO term enrichment with error handling
    tryCatch({
      corrected_significant <- corrected %>%
        dplyr::mutate(significant = ifelse(!is.na(adj_pval) & adj_pval < 0.05, TRUE, FALSE)) 
      
      df_go_term <- protti::calculate_go_enrichment(
        data = corrected_significant,
        protein_id = pg_protein_accessions,
        go_annotations_uniprot = go_f,
        is_significant = significant,
        min_n_detected_proteins_in_process = 3,
        plot = FALSE
      )
    
      # Save GO Term enrichment results
      go_term_file <- file.path(group_folder_path, paste0("go_term_", experiment_id, "_", comparison_filter, ".tsv"))
      write.table(df_go_term, go_term_file, sep = "\t", row.names= FALSE, quote = FALSE)
    }, error = function(e) {
      message(paste("Error in GO term enrichment for comparison", comparison_filter, ":", e))
    })
    
  } else {
    
    ## LiP analysis only when no TrP control
    
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
      paste0("differential_abundance_", experiment_id, "_", comparison_filter, ".tsv")
    )
    
    write.table(df_diff_abundance, diff_abundance_file, sep = "\t", row.names= FALSE, quote = FALSE)
    
    # plot the corresponding volcano plots
    plot_list5[[1]] <- volcano_plot(
      data = df_diff_abundance,
      grouping = eg_modified_peptide,
      log2FC = diff,
      significance = pval,
      method = "significant",
      x_axis_label = "log2(fold change) treated vs. untreated",
      significance_cutoff = c(0.05, "adj_pval") 
    ) +
      geom_text_repel(data = df_diff_abundance, 
                      size = 3, 
                      aes(x = diff, y = -log10(pval), label = label), 
                      min.segment.length = unit(0, 'lines'), nudge_y = 0.1)
    
    
    # plot the distribution of the p-values 
    plot_list5[[2]] <- pval_distribution_plot(data = df_diff_abundance,
                                              grouping = eg_modified_peptide,
                                              pval = pval
    )  
    
    output_stats_pdf <- file.path(group_folder_path, "statistical_analysis_plots.pdf")
    
    ggsave(
      filename = output_stats_pdf, 
      plot = marrangeGrob(plot_list5, nrow=1, ncol=1), 
      width = 8, height = 6
    )
    
    # plot the profile plots for the candidates: 
    candidate_summed <- DIA_clean_uniprot_summed_protti %>% 
      dplyr::filter(pg_protein_accessions %in% sort(unique(candidates$pg_protein_accessions))) %>% 
      left_join(distinct(DIA_clean, r_file_name, condrep), by = "r_file_name") %>% 
      left_join(distinct(candidates, eg_modified_peptide, significant), by = "eg_modified_peptide") %>% 
      dplyr::mutate(significant = ifelse(is.na(significant), FALSE, significant)) %>% 
      dplyr::mutate(line_type = ifelse(significant == TRUE, 'solid', 'dotted'))
    
    plot_list6 <- peptide_profile_plot(
      data = candidate_summed,
      sample = condrep,
      peptide = eg_modified_peptide,
      intensity_log2 = normalised_intensity_log2,
      grouping = pg_protein_accessions,
      targets = sort(unique(candidates$pg_protein_accessions)),
      protein_abundance_plot = FALSE
    ) 
    
    output_profile_pdf <- file.path(group_folder_path, "candidates_profile_plots.pdf")
    
    ggsave(
      filename = output_profile_pdf, 
      plot = marrangeGrob(plot_list6, nrow=1, ncol=1), 
      width = 8, height = 6
    )
    
    plot_list7 <- woods_plot(
      data = corrected,
      fold_change = adj_diff,
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
      plot = marrangeGrob(plot_list7, nrow=1, ncol=1), 
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
      go_term_file <- file.path(group_folder_path, paste0("go_term_", experiment_id, "_", comparison_filter, ".tsv"))
      write.table(df_go_term, go_term_file, sep = "\t", row.names= FALSE, quote = FALSE)
    }, error = function(e) {
      message(paste("Error in GO term enrichment for comparison", comparison_filter, ":", e))
    })
  }
    
}


# copy yaml file into the output as well
yaml_file_path <- file.path(group_folder_path, "params.yaml")
file.copy(yaml_file, yaml_file_path)

# Stop redirecting output to the log file
sink()
