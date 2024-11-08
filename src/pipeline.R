
# Comments
# Add package to every function
# for instance protti::normalise() instead of normalise()
# add parameter to every function 
# for instance assign_peptide_type(aa_before = aa_before, last_aa= last_aa, 
# aa_after = aa_after) 
# the previous version is not wrong but keeping this format ensures that the 
# code is more
# readable and avoids certain errors when some of the packages get updated

# Do we need the volcano plot or should we just stick with the differential 
# abundace results? 
# do you want to save it to make it downloadable on the website? 
# the target parameter and xlabel parameter in your volcano plot function are right
# now specific do your dataset so you would need to give it to the pipeline as
# parameter instead of hardcoding it in the script 

# the ref condition needs to be a parameter as well

# we need some code for the go term analysis 

# i figured that it would be better to have a YAML- file to pass the parameters
# also i would suggest that we save further information about the experiment in 
# the yaml so i can save in the database in the lipatlas

# load arguments
# use bash: Rscript pipeline.R

library(yaml)
library(protti)
library(dplyr)
library(missForest)
library(doParallel)

registerDoParallel(cores = 4) 

# Read the YAML file
args <- commandArgs(trailingOnly = TRUE)
# Expect the first argument to be the path to the YAML file
yaml_file <- args[1]
params <- yaml::read_yaml(yaml_file)

# Access the parameters
group_id <- params$group_id
input_file <- params$input_file
experiment_ids <- params$experiment_id
treatment <- params$treatment
ref_condition <- params$ref_condition
comparisons <- params$comparison
#taxonomy_id <- params$taxonomy_id
output_dir <- params$output_dir


group_folder_path <- file.path(output_dir, group_id)

if (!dir.exists("./preprocessed")) {
  dir.create("./preprocessed", recursive = TRUE)
}

# Now create the experiment folder
if (!dir.exists(group_folder_path)) {
  dir.create(group_folder_path, recursive = TRUE)
  cat("Created group folder:", group_folder_path, "\n")
}


# Log outputs to check for errors 
# Redirect output to a log file
logfile_dir <- file.path(output_dir, group_id, "processing_log.txt")
sink(logfile_dir, append = TRUE)

# Print or use the parameters
print(paste("Input File:", input_file))
print(paste("Experiment IDs:", experiment_ids))
print(paste("Treatment:", treatment))
print(paste("Ref Condtion:", ref_condition))
print(paste("Output Directory:", output_dir))


# Log R session information, including loaded packages
cat("Logging R session info:\n")
sessionInfo() 

# Check if the folder exists, and create it if it doesn't

# Create a list to store the plots
plot_list <- list()

# load file
DIA_raw <- read_protti(input_file)

# ------------------------------------------------------------------------------
# Preprocessing 
# ------------------------------------------------------------------------------

DIA_raw$intensity_log2 <- log2(DIA_raw$fg_quantity)

DIA_raw_norm <- protti::normalise(
  DIA_raw,
  sample = r_file_name,
  intensity_log2 = intensity_log2,
  method = "median"
)

DIA_clean <- DIA_raw_norm %>%
  dplyr::filter(fg_quantity > 1000) %>%
  dplyr::filter(pep_is_proteotypic == T)

DIA_clean$fg_id <- paste0(DIA_clean$fg_labeled_sequence, DIA_clean$fg_charge)

unis <- unique(DIA_clean$pg_protein_accessions) # make vector for fetch_uniprot

## Load data from uniprot and join with DIA dataframe
uniprot <- # download protein information from UniProt
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
  left_join(uniprot, by = c("pg_protein_accessions" = "accession")) %>% 
  find_peptide(sequence, pep_stripped_sequence) %>%
  assign_peptide_type(aa_before, last_aa, aa_after) %>%
  distinct()

DIA_clean_uniprot$condrep <- paste(DIA_clean_uniprot$r_condition, DIA_clean_uniprot$r_replicate, sep = "_")

proteins_identified <- uniprot %>%
  distinct(accession)


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# This might be interesting as an output for us?

# fwrite(proteins_identified, file = "identified_proteins_CaM_LiP.csv", sep = ",")

##Imputation on precursor lvl 
imputed <- impute_randomforest(
  DIA_clean_uniprot,
  sample = r_file_name,
  grouping = fg_id,
  intensity_log2 = normalised_intensity_log2,
  retain_columns = c("pep_stripped_sequence", "pg_protein_accessions", "gene_names", "go_f", "r_condition", "start", "end"),
  parallelize = "variables"
)

##sum up precursors to peptide lvl and keep only one entry per pep_stripped_sequence
###code here needs to be adjusted for protti function

DIA_clean_uniprot_summed_protti <- protti::calculate_protein_abundance(
  imputed,
  sample = r_file_name,
  protein_id = pep_stripped_sequence,
  precursor = fg_id,
  peptide = pep_stripped_sequence,
  intensity_log2 = normalised_intensity_log2,
  min_n_peptides = 1,
  method = "sum",
  for_plot = FALSE,
  retain_columns = c("pep_stripped_sequence", "pg_protein_accessions", "gene_names", "go_f", "r_condition", "start", "end")
)

# ------------------------------------------------------------------------------
# QC
# ------------------------------------------------------------------------------

plot_list[[1]] <- protti::qc_median_intensities(
  data = DIA_raw_norm,
  sample = r_file_name,
  grouping = pep_grouping_key,
  intensity = normalised_intensity_log2,
  plot = TRUE,
  interactive = FALSE
)

#QC

##CV
plot_list[[2]] <- protti::qc_cvs(
  data = DIA_clean_uniprot,
  grouping = fg_id,
  condition = r_condition,
  intensity = fg_quantity,
  plot = FALSE
)

plot_list[[3]] <- protti::qc_cvs(
  data = DIA_clean_uniprot,
  grouping = fg_id,
  condition = r_condition,
  intensity = fg_quantity,
  plot_style = "density",
  plot = TRUE
)

plot_list[[4]] <- protti::qc_cvs(
  data = DIA_clean_uniprot,
  grouping = fg_id,
  condition = r_condition,
  intensity = fg_quantity,
  plot_style = "violin",
  plot = TRUE
)

## Intensity distribution
#Intensity distributions are plotted for the whole dataset.

plot_list[[5]] <- protti::qc_intensity_distribution(
  data = DIA_clean_uniprot,
  sample = condrep,
  grouping = pep_grouping_key,
  intensity_log2 = intensity_log2,
  plot_style = "histogram"
)

## Missed cleavages
plot_list[[6]] <- protti::qc_missed_cleavages(
  data = DIA_clean_uniprot,
  sample = condrep,
  grouping = fg_id,
  missed_cleavages = pep_nr_of_missed_cleavages,
  intensity = fg_quantity,
  method = "intensity",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[7]] <- protti::qc_missed_cleavages(
  data = DIA_clean_uniprot,
  sample = condrep,
  grouping = fg_id,
  missed_cleavages = pep_nr_of_missed_cleavages,
  intensity = fg_quantity,
  method = "count",
  plot = TRUE,
  interactive = FALSE
)


## Peptide types

#Peptide type (tryptic, semi-tryptic, non-tryptic) distributions are even throughout the different samples.

plot_list[[8]] <- qc_peptide_type(
  DIA_clean_uniprot,
  condrep,
  fg_id,
  pep_type,
  intensity = fg_quantity,
  method = "count",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[9]] <- qc_peptide_type(
  DIA_clean_uniprot,
  condrep,
  fg_id,
  pep_type,
  intensity = fg_quantity,
  method = "intensity",
  plot = TRUE,
  interactive = FALSE
)


## Number of peptide IDs per sample
#The numbers of identified peptides are consistent throughout the different samples.
DIA_raw$condrep <- paste(DIA_raw$r_condition, DIA_raw$r_replicate, sep = "_")

plot_list[[10]] <- qc_ids(
  data = DIA_raw, 
  sample = condrep, 
  grouping = pep_grouping_key, 
  condition = r_condition, 
  intensity = fg_quantity
)

## Principal component analysis (PCA)
plot_list[[11]] <-DIA_clean_uniprot %>%
  qc_pca(
    sample = condrep, 
    grouping = pep_grouping_key, 
    intensity = intensity_log2, 
    condition = r_condition
  )

## corelation_map
plot_list[[12]] <- qc_sample_correlation(
  data = DIA_clean_uniprot,
  sample = r_file_name,
  grouping = fg_id,
  intensity_log2 = intensity_log2,
  condition = r_condition
)

# ------------------------------------------------------------------------------
# Analysis ...
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# differential abundance
# ------------------------------------------------------------------------------

# Data analysis
## Volcano plots peptide lvl
# Iterate over each comparison and associated experiment ID
for (i in seq_along(comparisons)) {
  comparison_filter <- comparisons[[i]]
  experiment_id <- experiment_ids[[i]]
  
  comparison_parts <- strsplit(comparison_filter, "_vs_")[[1]] 
  
  filtered_data <- DIA_clean_uniprot_summed_protti %>%
    dplyr::filter(r_condition %in% comparison_parts)
  
  Volcano_input <- filtered_data %>%
    unique() %>%
    protti::assign_missingness(
      sample = r_file_name,
      condition = r_condition,
      grouping = pep_stripped_sequence,
      intensity = normalised_intensity_log2,
      ref_condition = ref_condition,
      retain_columns = all_of(c("pg_protein_accessions", "r_file_name", "r_condition", 
                                "normalised_intensity_log2", 
                                  "pep_stripped_sequence", "go_f", "start", "end"))
    )
  
  df_diff_abundance <- protti::calculate_diff_abundance(
    data = Volcano_input,
    sample = r_file_name,
    condition = r_condition,
    grouping = pep_stripped_sequence,
    intensity_log2 = normalised_intensity_log2,
    missingness = missingness,
    comparison = comparison,
    method = "moderated_t-test",
    retain_columns = all_of(c("pg_protein_accessions", "pep_stripped_sequence", 
                              "comparison", "go_f", "start", "end"))
  )

  # Subset data for the specific comparison
  df_diff_abundance_subset <- df_diff_abundance %>%
    dplyr::filter(comparison == comparison_filter)

  # Save Differential Abundance results
  diff_abundance_file <- file.path(group_folder_path, paste0("differential_abundance_", experiment_id, "_", comparison, ".csv"))
  write.csv(df_diff_abundance_subset, diff_abundance_file)
  
  # Filter significant proteins for GO term analysis
  df_diff_abundance_significant <- df_diff_abundance_subset %>%
    dplyr::mutate(significant = ifelse(adj_pval < 0.05, TRUE, FALSE)) %>%
    dplyr::filter(!is.na(go_f), significant == TRUE)
  
  # Calculate GO term enrichment for the current comparison
  df_go_term <- protti::calculate_go_enrichment(
    data = df_diff_abundance_significant,
    protein_id = pg_protein_accessions,
    go_annotations_uniprot = go_f,
    is_significant = significant,
    plot = FALSE,
    plot_cutoff = "pval 0.01"
  )
  
  # Save GO Term enrichment results
  go_term_file <- file.path(group_folder_path, paste0("go_term_", experiment_id, "_", comparison, ".csv"))
  write.csv(df_go_term, go_term_file)
}


# ------------------------------------------------------------------------------
# Save everything 
# ------------------------------------------------------------------------------

# Save QC plots
output_qc_pdf <- file.path(group_folder_path, "qc_plots.pdf")
pdf(output_qc_pdf, width = 7, height = 5)  # Open the PDF device
# Loop through the list and print each ggplot object
lapply(plot_list, print)
dev.off()  

# saving a file with the metadata for the experiment would be great as well 

# copy yaml file into the output as well
yaml_file_path <- file.path(group_folder_path, "params.yaml")
file.copy(yaml_file, yaml_file_path)

# Stop redirecting output to the log file
sink()