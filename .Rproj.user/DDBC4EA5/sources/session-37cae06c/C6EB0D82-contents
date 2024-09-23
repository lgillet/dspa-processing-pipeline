ludovic_assign_missingness <- function(data, file_name, r_condition, r_condrep, precursor, intensity, peptide = pep_stripped_sequence) {
  # Reshape data: Pivot to wide format with file_name as columns and precursor as row names
  data_wide <- data %>%
    dplyr::select({{file_name}}, {{precursor}}, {{intensity}}, {{peptide}}) %>%  # Include peptide here
    tidyr::pivot_wider(names_from = {{file_name}}, values_from = {{intensity}}) 
  
  DIA_clean_unique <- DIA_clean %>%
    dplyr::distinct(r_file_name, Condition, Replicate)
  
  # Reshape back to long format to keep NA values, add r_condition and r_condrep
  data_long <- data_wide %>%
    tidyr::pivot_longer(cols = -c(fg_id, {{peptide}}), names_to = "r_file_name", values_to = "intensity_log2") %>%  # Include peptide
    left_join(DIA_clean_unique, by = "r_file_name")
  
  # Group by precursor and condition to determine missingness
  data_missing_type <- data_long %>%
    dplyr::group_by({{precursor}}, {{r_condition}}, {{peptide}}) %>%  # Group by peptide as well
    # Calculate missing count and replicate count per group
    dplyr::mutate(missing = sum(is.na({{intensity}}))) %>%
    dplyr::mutate(repl = dplyr::n()) %>%
    rowwise() %>%
    # Determine NA type: MNAR, MAR, or complete
    dplyr::mutate(missingness = ifelse(is.na({{intensity}}),
                                       ifelse(repl - missing == 0, "MNAR", "MAR"), "complete")) %>%
    # Flag missing values
    dplyr::mutate(is_missing = is.na({{intensity}})) %>%
    # Select relevant columns and include peptide
    dplyr::select({{r_condrep}}, {{precursor}}, {{peptide}}, {{intensity}}, {{file_name}}, missingness, is_missing) %>%
    ungroup()
  
  return(data_missing_type)
}
