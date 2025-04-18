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
library(magrittr)


registerDoParallel(cores = 10) 

# Read the YAML file
args <- commandArgs(trailingOnly = TRUE)
# Expect the first argument to be the path to the YAML file
yaml_file <- args[1]
params <- yaml::read_yaml(yaml_file)

# Access the parameters
group_id <- params$group_id
input_file <- params$input_file
input_file_tryptic_control <- params$input_file_tryptic_control
experiment_ids <- params$experiment_id
treatment <- params$treatment
ref_condition <- params$ref_condition
comparisons <- params$comparison
output_dir <- params$output_dir

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
logfile_dir <- file.path(output_dir, group_id, "processing_log.txt")
sink(logfile_dir, append = TRUE)

# Log R session information, including loaded packages
cat("Logging R session info:\n")
sessionInfo() 

# Create a list to store the plots
plot_list <- list()
plot_list2 <- list()

# load file
df <- read_protti(input_file)

# ------------------------------------------------------------------------------
# LiP 
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Preprocessing 
# ------------------------------------------------------------------------------

df$intensity_log2 <- log2(df$fg_ms2raw_quantity)
df$condrep <- paste(df$r_condition, df$r_replicate, sep = "_")

plot_list[[1]] <- protti::qc_ids(
  data = df, 
  sample = condrep, 
  grouping = pep_grouping_key, 
  condition = r_condition, 
  intensity = fg_quantity
)+ ggtitle('Precursor ID count per sample')

plot_list[[2]] <- protti::qc_ids(
  data = df, 
  sample = condrep, 
  grouping = pg_protein_accessions, 
  condition = r_condition, 
  intensity = pg_quantity
)+ ggtitle('Protein ID count per sample')

# ------------------------------------------------------------------------------
# Normalise
# ------------------------------------------------------------------------------

df %<>% protti::normalise(
  sample = r_file_name,
  intensity_log2 = intensity_log2,
  method = "median"
) 

plot_list[[3]] <- protti::qc_intensity_distribution(
  data = df,
  sample = condrep,
  grouping = pep_grouping_key,
  intensity_log2 = intensity_log2,
  plot_style = "histogram"
) + ggtitle('Overall log2 Intensity distribution before normalisation')

plot_list[[4]] <- protti::qc_intensity_distribution(
  data =df,
  sample = condrep,
  grouping = pep_grouping_key,
  intensity_log2 = normalised_intensity_log2,
  plot_style = "boxplot"
) + ggtitle('Run intensities after normalisation')

df %<>% 
  dplyr::filter(intensity_log2 > 10) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(pg_protein_accessions2 = ifelse(base::grepl(";", pg_protein_accessions, fixed = FALSE), 
    base::sort(base::strsplit(pg_protein_accessions, ";", fixed = TRUE)[[1]])[1], pg_protein_accessions)) %>% 
  dplyr::ungroup()

df$fg_id <- paste0(df$fg_labeled_sequence, df$fg_charge)
unis <- unique(df$pg_protein_accessions2) # make vector for fetch_uniprot

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

df %<>% 
  left_join(uniprot, by = c("pg_protein_accessions2" = "accession")) %>% 
  find_peptide(sequence, pep_stripped_sequence) %>%
  assign_peptide_type(aa_before, last_aa, aa_after) %>%
  distinct() %>% 
  calculate_sequence_coverage(protein_sequence = sequence, peptides = pep_stripped_sequence) %>% 
  dplyr::mutate(normalised_intensity = 2^normalised_intensity_log2)

plot_list[[5]] <- protti::qc_cvs(
  data = df,
  grouping = fg_id,
  condition = r_condition,
  intensity = normalised_intensity,
  plot_style = "violin",
  plot = TRUE
)

## Missed cleavages
plot_list[[6]] <- protti::qc_missed_cleavages(
  data = df,
  sample = condrep,
  grouping = fg_id,
  missed_cleavages = pep_nr_of_missed_cleavages,
  intensity = normalised_intensity,
  method = "count",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[7]] <- protti::qc_missed_cleavages(
  data = df,
  sample = condrep,
  grouping = fg_id,
  missed_cleavages = pep_nr_of_missed_cleavages,
  intensity = normalised_intensity,
  method = "intensity",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[8]] <- protti::qc_peptide_type(
  df,
  condrep,
  fg_id,
  pep_type,
  intensity = normalised_intensity,
  method = "count",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[9]] <- protti::qc_peptide_type(
  df,
  condrep,
  fg_id,
  pep_type,
  intensity = normalised_intensity,
  method = "intensity",
  plot = TRUE,
  interactive = FALSE
)

## Principal component analysis (PCA)
plot_list[[10]] <- df %>%
  protti::qc_pca(
    sample = condrep, 
    grouping = pep_grouping_key, 
    intensity = intensity_log2, 
    condition = r_condition
  )

## correlation_map
plot_list[[11]] <- protti::qc_sample_correlation(
  data = df,
  sample = condrep,
  grouping = fg_id,
  intensity_log2 = intensity_log2,
  condition = r_condition
)[[4]]


df %<>% 
  distinct(r_file_name, fg_id, normalised_intensity_log2, eg_modified_peptide, pep_stripped_sequence, pg_protein_accessions, gene_names, go_f, r_condition, start, end) %>% 
  tidyr::complete(nesting(r_file_name, r_condition), nesting(pg_protein_accessions, gene_names, go_f, fg_id, eg_modified_peptide, pep_stripped_sequence, start, end))

df %<>% impute_randomforest(
  sample = r_file_name,
  grouping = fg_id,
  intensity_log2 = normalised_intensity_log2,
  retain_columns = c("eg_modified_peptide", "pep_stripped_sequence", "pg_protein_accessions",
                     "gene_names", "go_f", "r_condition", "start", "end"),
  parallelize = "variables"
)

imputed_file <- file.path(group_folder_path, paste0("imputed.tsv"))
write.table(df, imputed_file, sep = "\t", row.names= FALSE, quote = FALSE)
save.image("./all_data.RData")

plot_list[[12]] <- df %>%
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

# sum up precursors to peptide level and keep only one entry per pep_stripped_sequence
df %<>% protti::calculate_protein_abundance(
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
write.table(df, dia_clean_file, sep = "\t", row.names= FALSE, quote = FALSE)

# ------------------------------------------------------------------------------
# Save QC plots
# ------------------------------------------------------------------------------

# Save QC plots
output_qc_pdf <- file.path(group_folder_path, "qc_plots.pdf")
ggsave(
  filename = output_qc_pdf, 
  plot = marrangeGrob(plot_list, nrow=1, ncol=1), 
  width = 8, height = 6
)

# ------------------------------------------------------------------------------
# TRYPTIC CONTROL
# ------------------------------------------------------------------------------

# Tryptic Control (if provided)
if (!is.null(input_file_tryptic_control)) {
  
  input_file_tryptic_control <- params$input_file_tryptic_control
  df_tryptic <- protti::read_protti(input_file_tryptic_control)
  df_tryptic$intensity_log2 <- log2(df_tryptica$fg_ms2raw_quantity)
  df_tryptic$condrep <- paste(df_tryptic$r_condition, df_tryptic$r_replicate, sep = "_")
  
  # make the QC plots as for LiP
  plot_list2[[1]] <- protti::qc_ids(
    data = df_tryptic, 
    sample = condrep, 
    grouping = pep_grouping_key, 
    condition = r_condition, 
    intensity = fg_quantity
  )+ ggtitle('Tryptic control: Precursor ID count per sample')
  
  plot_list2[[2]] <- protti::qc_ids(
    data = df_tryptic, 
    sample = condrep, 
    grouping = pg_protein_accessions, 
    condition = r_condition, 
    intensity = pg_quantity
  )+ ggtitle('Tryptic control: Protein ID count per sample')
  
  # Preprocess tryptic control data  
  df_tryptic %<>% protti::normalise(
    sample = r_file_name,
    intensity_log2 = intensity_log2,
    method = "median"
  )
  
  plot_list2[[3]] <- protti::qc_intensity_distribution(
    data = df_tryptic,
    sample = condrep,
    grouping = fg_id,
    intensity_log2 = intensity_log2,
    plot_style = "histogram"
  ) + ggtitle('Tryptic control: Overall log2 Intensity distribution before normalisation')
  
  plot_list2[[4]] <- protti::qc_intensity_distribution(
    data = df_tryptic,
    sample = condrep,
    grouping = pep_grouping_key,
    intensity_log2 = normalised_intensity_log2,
    plot_style = "boxplot"
  ) + ggtitle('Tryptic control: Run intensities after normalisation')
  
  df_tryptic %<>%
    dplyr::filter(intensity_log2 > 10) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(pg_protein_accessions2 = ifelse(base::grepl(";", pg_protein_accessions, fixed = FALSE), 
      base::sort(base::strsplit(pg_protein_accessions, ";", fixed = TRUE)[[1]])[1], pg_protein_accessions)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(normalised_intensity = 2^normalised_intensity_log2)
  
  df_tryptic$fg_id <- paste0(df_tryptic$fg_labeled_sequence, df_tryptic$fg_charge)
  
  # Fetch UniProt annotations for tryptic control
  unis_tryptic <- unique(df_tryptic$pg_protein_accessions2)

  ## Load data from uniprot and join with DIA dataframe
  uniprot_tryptic <-
    protti::fetch_uniprot(
      df_tryptic,
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
  
  df_tryptic %<>%
    left_join(uniprot_tryptic, by = c("pg_protein_accessions2" = "accession")) %>% 
    find_peptide(sequence, pep_stripped_sequence) %>%
    assign_peptide_type(aa_before, last_aa, aa_after) %>%
    distinct() %>% 
    calculate_sequence_coverage(protein_sequence = sequence, peptides = pep_stripped_sequence)
  
  plot_list2[[5]] <- protti::qc_cvs(
    data = df_tryptic,
    grouping = fg_id,
    condition = r_condition,
    intensity = normalised_intensity,
    plot_style = "violin",
    plot = TRUE
  ) + ggtitle('Tryptic control: Violon plot')
  
  ## Missed cleavages
  plot_list2[[6]] <- protti::qc_missed_cleavages(
    data = df_tryptic,
    sample = condrep,
    grouping = fg_id,
    missed_cleavages = pep_nr_of_missed_cleavages,
    intensity = normalised_intensity,
    method = "count",
    plot = TRUE,
    interactive = FALSE
  ) + ggtitle('Tryptic control: number of missed cleavages')
  
  plot_list2[[7]] <- protti::qc_missed_cleavages(
    data = df_tryptic,
    sample = condrep,
    grouping = fg_id,
    missed_cleavages = pep_nr_of_missed_cleavages,
    intensity = normalised_intensity,
    method = "intensity",
    plot = TRUE,
    interactive = FALSE
  ) + ggtitle('Tryptic control: intensity of missed cleavages')
  
  plot_list2[[8]] <- protti::qc_peptide_type(
    df_tryptic,
    condrep,
    fg_id,
    pep_type,
    intensity = normalised_intensity,
    method = "count",
    plot = TRUE,
    interactive = FALSE
  ) + ggtitle('Tryptic control: number of half tryptics')
  
  plot_list2[[9]] <- protti::qc_peptide_type(
    df_tryptic,
    condrep,
    fg_id,
    pep_type,
    intensity = normalised_intensity,
    method = "intensity",
    plot = TRUE,
    interactive = FALSE
  ) + ggtitle('Tryptic control: intensity of half tryptics')
  
  df_tryptic %<>%
    distinct(r_file_name, fg_id, normalised_intensity_log2, eg_modified_peptide, pep_stripped_sequence, pg_protein_accessions, pg_protein_accessions2, gene_names, go_f, r_condition, start, end) %>% 
    tidyr::complete(nesting(r_file_name, r_condition), nesting(pg_protein_accessions, pg_protein_accessions2, gene_names, go_f, fg_id, eg_modified_peptide, pep_stripped_sequence, start, end))
  
  # Save preprocessed tryptic control data
  write.table(df_tryptic, file.path(group_folder_path, "tryptic_control_clean.tsv"), sep = "\t", row.names= FALSE, quote = FALSE)

  # calculate protein abundance
  df_tryptic %<>% calculate_protein_abundance(
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
  write.table(df_tryptic, tryptic_clean_file, sep = "\t", row.names= FALSE, quote = FALSE)
  
  ## Principal component analysis (PCA)
  plot_list2[[10]] <- df_tryptic %>%
    protti::qc_pca(
      sample = condrep, 
      grouping = pg_protein_accessions, 
      intensity = normalised_intensity_log2, 
      condition = r_condition
    ) + ggtitle('Tryptic control: PCA at protein level')
  
  ## correlation_map
  plot_list2[[11]] <- protti::qc_sample_correlation(
    data = df_tryptic,
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


for (i in seq_along(comparisons)) {
  comparison_filter <- comparisons[[i]]
  experiment_id <- experiment_ids[[i]]
  
  comparison_parts <- strsplit(comparison_filter, "_vs_")[[1]] 
  
  df_filtered <- df %>%
    dplyr::filter(r_condition %in% comparison_parts)
  
  df_diff <- df_filtered  %>%
    unique() %>%
    protti::assign_missingness(
      sample = r_file_name,
      condition = r_condition,
      grouping = eg_modified_peptide,
      intensity = normalised_intensity_log2,
      ref_condition = ref_condition,
      retain_columns = all_of(c("pg_protein_accessions", "r_file_name", "r_condition", 
                                "normalised_intensity_log2", 
                                "go_f", "start", "end")))%>%
    protti::calculate_diff_abundance(
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
    df_trp_filtered <- df_trp %>%
      dplyr::filter(r_condition %in% comparison_parts)
    
    df_trp_filtered_diff <- df_trp_filtered  %>%
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
    
    Trp_candidates <- df_trp_filtered_diff %>% 
      dplyr::filter(adj_pval < 0.05 & abs(diff) > 1) %>% 
      left_join(dplyr::distinct(tryptic_protein_abundance, pg_protein_accessions2, gene_names), by = "pg_protein_accessions2") %>% 
      rowwise() %>% 
      dplyr::mutate(gene = sort(strsplit(gene_names, " ", fixed = TRUE)[[1]])[1]) %>% 
      ungroup() %>% 
      dplyr::mutate(significant = TRUE)

    df_trp_filtered_diff %<>%
      left_join(distinct(Trp_candidates, pg_protein_accessions2, significant, gene), by = "pg_protein_accessions2") %>% 
      left_join(dplyr::distinct(tryptic_protein_abundance, pg_protein_accessions2, gene_names), by = "pg_protein_accessions2")
    
    # Save Differential Abundance results
    diff_trp_file_path <- file.path(
      group_folder_path, 
      paste0("trp_differential_abundance_", experiment_ids, "_", comparisons, ".tsv")
    )
    write.table( df_trp_filtered_diff, diff_trp_file_path, sep = "\t", row.names= FALSE, quote = FALSE)
    
    # perform TrP protein correction on LiP:
    df_trp_filtered_diff %<>%
      rowwise() %>% 
      dplyr::mutate(comparison = gsub("_Trp", "_LiP", x = comparison, fixed = TRUE)) %>% 
      ungroup()
    
    df_diff<- correct_lip_for_abundance(
      lip_data = df_diff,
      trp_data =  df_trp_filtered_diff,
      protein_id = pg_protein_accessions,
      grouping = eg_modified_peptide,
      comparison = comparison, 
      diff = diff,
      n_obs = n_obs,
      std_error = std_error,
      p_adj_method = "BH",
      retain_columns = all_of(c("missingness", "go_f")),
      method = "satterthwaite"
    )

    # bug: output some important plots to pdf and results to the log file: 
    candidates <- df_diff %>% 
      dplyr::filter(adj_pval < 0.05 & abs(adj_diff) > 1) %>% 
      left_join(dplyr::distinct(DIA_clean_uniprot, eg_modified_peptide, pg_protein_accessions2, gene_names, coverage, length, start, end), by = "eg_modified_peptide") %>% 
      rowwise() %>% 
      dplyr::mutate(gene = sort(strsplit(gene_names, " ", fixed = TRUE)[[1]])[1]) %>% 
      ungroup() %>% 
      dplyr::mutate(label = paste(gene, "@", start, "-", end, sep = "")) %>% 
      dplyr::mutate(significant = TRUE)
    
    df_diff %<>%
      left_join(distinct(candidates, eg_modified_peptide, start, end, label, coverage, length, significant), by = "eg_modified_peptide") 
      
    # Save Differential Abundance results
    corrected_file <- file.path(
      group_folder_path, 
      paste0("differential_abundance_", experiment_id, "_", comparison_filter, ".tsv")
    )
    write.table(df_diff, corrected_file, sep = "\t", row.names= FALSE, quote = FALSE)
  
  } else {
    
    ## LiP analysis only when no TrP control
    candidates <- df_diff %>% 
      dplyr::filter(adj_pval < 0.05 & abs(diff) > 1) %>% 
      left_join(dplyr::distinct(DIA_clean_uniprot, pg_protein_accessions, gene_names), by = "pg_protein_accessions") %>% 
      rowwise() %>% 
      dplyr::mutate(gene = sort(strsplit(gene_names, " ", fixed = TRUE)[[1]])[1]) %>% 
      ungroup() %>% 
      dplyr::mutate(label = paste(gene, "@", start, "-", end, sep = "")) %>% 
      dplyr::mutate(significant = TRUE)
    
    df_diff %<>%
      left_join(distinct(candidates, eg_modified_peptide, label, significant), by = "eg_modified_peptide") %>% 
      left_join(dplyr::distinct(pg_protein_accessions, gene_names, sequence, length, coverage), by = "pg_protein_accessions") %>% 
      rowwise() %>% 
      dplyr::mutate(gene = sort(strsplit(gene_names, " ", fixed = TRUE)[[1]])[1]) 
    
    # Save Differential Abundance results
    diff_abundance_file <- file.path(
      group_folder_path, 
      paste0("differential_abundance_", experiment_id, "_", comparison_filter, ".tsv")
    )
    write.table(df_diff, diff_abundance_file, sep = "\t", row.names= FALSE, quote = FALSE)
    
    
  }
  # plot the profile plots for the candidates: 
  candidate_summed <- df_diff %>% 
    dplyr::filter(pg_protein_accessions %in% sort(unique(candidates$pg_protein_accessions))) %>% 
    left_join(distinct(DIA_clean, r_file_name, condrep), by = "r_file_name") %>% 
    left_join(distinct(candidates, eg_modified_peptide, significant), by = "eg_modified_peptide") %>% 
    dplyr::mutate(significant = ifelse(is.na(significant), FALSE, significant)) %>% 
    dplyr::mutate(line_type = ifelse(significant == TRUE, 'solid', 'dotted'))
  
  peptide_profile_plot <- protti::peptide_profile_plot(
    data = candidate_summed,
    sample = condrep,
    peptide = eg_modified_peptide,
    intensity_log2 = normalised_intensity_log2,
    grouping = pg_protein_accessions,
    targets = sort(unique(candidates$pg_protein_accessions)),
    protein_abundance_plot = FALSE
  ) 
  
  output_profile_pdf <- file.path(group_folder_path, 
                                  paste0("candidates_profile_plots_", experiment_id, "_", comparison_filter, ".pdf"))
  ggsave(
    filename = output_profile_pdf, 
    plot = marrangeGrob(peptide_profile_plot, nrow=1, ncol=1), 
    width = 8, height = 6
  )
  
  woods_plot <- protti::woods_plot(
    data = df_diff,
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
  
  output_woods_pdf <- file.path(
    group_folder_path, 
    paste0("LiP_corrected_candidates_woods_plots_", experiment_id, "_", comparison_filter, ".pdf"))
  ggsave(
    filename = output_woods_pdf, 
    plot = marrangeGrob(woods_plot, nrow=1, ncol=1), 
    width = 8, height = 6
  )
  
  tryCatch({
    df_go_term <- df_diff %>%
      dplyr::mutate(significant = ifelse(!is.na(adj_pval) & adj_pval < 0.05, TRUE, FALSE)) 
      protti::calculate_go_enrichment(
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
  
}


# copy yaml file into the output as well
yaml_file_path <- file.path(group_folder_path, "params.yaml")
file.copy(yaml_file, yaml_file_path)

# Stop redirecting output to the log file
sink()
