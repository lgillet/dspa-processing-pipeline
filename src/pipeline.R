
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

registerDoParallel(cores = 6)

# Read the YAML file
params <- yaml::read_yaml("param_files/params_LIP00010.yaml")

# Access the parameters
input_file <- params$input_file
experiment_id <- params$experiment_id
treatment <- params$treatment
dose <- params$dose
ref_condition <- params$ref_condition
comparison <- params$comparison
#taxonomy_id <- params$taxonomy_id
output_dir <- params$output_dir


experiment_folder_path <- file.path(output_dir, experiment_id)

if (!dir.exists("/preprocessed")) {
  dir.create("/preprocessed", recursive = TRUE)
}

# Now create the experiment folder
if (!dir.exists(experiment_folder_path)) {
  dir.create(experiment_folder_path, recursive = TRUE)
  cat("Created experiment folder:", experiment_folder_path, "\n")
}


# Log outputs to check for errors 
# Redirect output to a log file
logfile_dir <- file.path(output_dir, experiment_id, "processing_log.txt")
sink(logfile_dir, append = TRUE)

# Print or use the parameters
print(paste("Input File:", input_file))
print(paste("Experiment ID:", experiment_id))
print(paste("Treatment:", treatment))
print(paste("Dose:", dose))
print(paste("Ref Condtion:", ref_condition))
print(paste("Output Directory:", output_dir))


# Log R session information, including loaded packages
cat("Logging R session info:\n")
sessionInfo() 

# Check if the folder exists, and create it if it doesn't

# Create a list to store the plots
plot_list <- list()

# ------------------------------------------------------------------------------
# Preprocessing 
# ------------------------------------------------------------------------------

DIA_raw$intensity_log2 <- log2(DIA_raw$fg_ms2raw_quantity)

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
  fetch_uniprot(
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
  left_join(uniprot, by = c("pg_protein_accessions" = "accession")) %>% # rejoin with annotations
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
  r_file_name,
  fg_id,
  normalised_intensity_log2,
  retain_columns = names(DIA_clean_uniprot),
  parallelize = "variables"
)

##sum up precursors to peptide lvl and keep only one entry per pep_stripped_sequence
###code here needs to be adjusted for protti function

DIA_clean_uniprot_summed_protti <- calculate_protein_abundance(
  imputed,
  sample = r_file_name,
  protein_id = pep_stripped_sequence,
  precursor = fg_id,
  peptide = pep_stripped_sequence,
  intensity_log2 = normalised_intensity_log2,
  min_n_peptides = 1,
  method = "sum",
  for_plot = FALSE,
  retain_columns = names(DIA_clean_uniprot)
)

# ------------------------------------------------------------------------------
# QC
# ------------------------------------------------------------------------------

plot_list[[1]] <- qc_median_intensities(
  DIA_raw_norm,
  r_file_name,
  pep_grouping_key,
  normalised_intensity_log2,
  plot = TRUE,
  interactive = FALSE
)

#QC

##CV
plot_list[[2]] <- qc_cvs(
  data = DIA_clean_uniprot,
  grouping = fg_id,
  condition = r_condition,
  intensity = fg_quantity,
  plot = FALSE
)

plot_list[[3]] <- qc_cvs(
  data = DIA_clean_uniprot,
  grouping = fg_id,
  condition = r_condition,
  intensity = fg_quantity,
  plot_style = "density",
  plot = TRUE
)

plot_list[[4]] <- qc_cvs(
  data = DIA_clean_uniprot,
  grouping = fg_id,
  condition = r_condition,
  intensity = fg_quantity,
  plot_style = "violin",
  plot = TRUE
)

## Intensity distribution
#Intensity distributions are plotted for the whole dataset.

plot_list[[5]] <- qc_intensity_distribution(
  DIA_clean_uniprot,
  condrep,
  pep_grouping_key,
  intensity_log2,
  plot_style = "histogram"
)

## Missed cleavages
plot_list[[6]] <- qc_missed_cleavages(
  DIA_clean_uniprot,
  condrep,
  fg_id,
  pep_nr_of_missed_cleavages,
  fg_quantity,
  method = "intensity",
  plot = TRUE,
  interactive = FALSE
)

plot_list[[7]] <- qc_missed_cleavages(
  DIA_clean_uniprot,
  condrep,
  fg_id,
  pep_nr_of_missed_cleavages,
  fg_quantity,
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
  DIA_raw, 
  condrep, 
  pep_grouping_key, 
  condition = r_condition, 
  intensity = fg_quantity
)

## Principal component analysis (PCA)
plot_list[[11]] <-DIA_clean_uniprot %>%
  qc_pca(
    condrep, 
    pep_grouping_key, 
    intensity_log2, 
    r_condition
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
Volcano_input <- unique_pep_filt %>%
  unique() %>%
  assign_missingness(
    r_file_name,
    r_condition,
    pep_stripped_sequence,
    normalised_intensity_log2,
    ref_condition = ref_condition, # cannot be hardcoded
    retain_columns = c(pg_protein_accessions, pep_stripped_sequence, start, end, pep_type, sequence, length))

Volcano_input = Volcano_input %>%
  dplyr::filter(is.na(sequence)== F)

#differential analysis
df_diff_abundance <- calculate_diff_abundance(
  data = Volcano_input,
  r_file_name,
  r_condition,
  pep_stripped_sequence,
  normalised_intensity_log2,
  missingness,
  comparison = comparison,
  ref_condition = ref_condition, # cannot be hardcoded 
  method = "moderated_t-test",
  retain_columns = c(pg_protein_accessions, pep_stripped_sequence, start, end, pep_type, sequence, length, comparison))

# ------------------------------------------------------------------------------
# GO terms
# ------------------------------------------------------------------------------

go <- fetch_go(taxonomy_id)

df_go_term <- calculate_go_enrichment(
  data,
  protein_id = accession,
  go_annotations_uniprot = go_f,
  is_significant = significant,
  plot = FALSE,
  plot_cutoff = "pval 0.01"
)

# ------------------------------------------------------------------------------
# Save everything 
# ------------------------------------------------------------------------------

# Save QC plots
output_qc_pdf <- file.path(output_dir, experiment_folder, "qc_plots.pdf")
pdf(output_qc_pdf, width = 7, height = 5)  # Open the PDF device
# Loop through the list and print each ggplot object
lapply(plot_list, print)
dev.off()  

# Save Differential Abundance 
differential_abundance_file_path <- file.path(output_dir, experiment_folder, "differential_abundance.csv")
write.csv(df_diff_abundance, differential_abundance_file_path)

# Save GO Term
go_term_file_path <- file.path(output_dir, experiment_folder, "go_term.csv")
write.csv(df_go_term, go_term_file_path)

# saving a file with the metadata for the experiment would be great as well 

# copy yaml file into the output as well
yaml_file_path <- file.path(output_dir, experiment_folder, "params.yaml")
file.copy("param_files/params_LIP00010.yaml", yaml_file_path)

# Stop redirecting output to the log file
sink()