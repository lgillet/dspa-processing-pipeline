#' Label-free peptide quantification
#'
#' Determines relative peptide abundances from ion quantification. 
#'
#' @param data a data frame that contains at least the input variables.
#' @param sample a character column in the \code{data} data frame that contains the sample name.
#' @param peptide a character column in the \code{data} data frame that contains peptide sequences.
#' This column is needed to filter for proteins with at least 3 unique peptides. This can equate
#' to more than three precursors. The quantification is done on the precursor level.
#' @param precursor a character column in the \code{data} data frame that contains precursors.
#' @param intensity_log2 a numeric column in the \code{data} data frame that contains log2
#' transformed precursor intensities.
#' @param method a character value specifying with which method peptide quantities should be
#' calculated. Possible options include \code{"sum"}, which takes the sum of all precursor
#' intensities as the protein abundance. Another option is \code{"iq"}, which performs protein
#' quantification based on a maximal peptide ratio extraction algorithm that is adapted from the
#' MaxLFQ algorithm of the MaxQuant software. Functions from the
#' \href{https://academic.oup.com/bioinformatics/article/36/8/2611/5697917}{\code{iq}} package are
#' used. Default is \code{"iq"}.
#' @param for_plot a logical value indicating whether the result should be only protein intensities
#' or protein intensities together with precursor intensities that can be used for plotting using
#' \code{peptide_profile_plot()}. Default is \code{FALSE}.
#' @param retain_columns a vector indicating if certain columns should be retained from the input
#' data frame. Default is not retaining additional columns \code{retain_columns = NULL}. Specific
#' columns can be retained by providing their names (not in quotations marks, just like other
#' column names, but in a vector).
#'
#' @return If \code{for_plot = FALSE}, protein abundances are returned, if \code{for_plot = TRUE}
#' also precursor intensities are returned in a data frame. The later output is ideal for plotting
#' with \code{peptide_profile_plot()} and can be filtered to only include protein abundances.
#'
#' @import dplyr
#' @import progress
#' @importFrom tidyr complete pivot_wider drop_na
#' @importFrom rlang .data := !! ensym as_name enquo
#' @importFrom tibble column_to_rownames as_tibble rownames_to_column
#' @importFrom magrittr %>%
#' @importFrom purrr map map2_df discard pluck
#' @export
#'
#' @examples
#' \donttest{
#' # Create example data
#' data <- data.frame(
#'   sample = c(
#'     rep("S1", 6),
#'     rep("S2", 6),
#'     rep("S1", 2),
#'     rep("S2", 2)
#'   ),,
#'   precursor = c(
#'     rep(c("A1", "A2", "B1", "B2", "C1", "D1"), 2),
#'     rep(c("E1", "F1"), 2)
#'   ),
#'   peptide = c(
#'     rep(c("A", "A", "B", "B", "C", "D"), 2),
#'     rep(c("E", "F"), 2)
#'   ),
#'   intensity = c(
#'     rnorm(n = 6, mean = 15, sd = 2),
#'     rnorm(n = 6, mean = 21, sd = 1),
#'     rnorm(n = 2, mean = 15, sd = 1),
#'     rnorm(n = 2, mean = 15, sd = 2)
#'   )
#' )
#'
#' data
#'
#' # Calculate peptide abundances
#' peptide_abundance <- calculate_peptide_abundance(
#'   data,
#'   sample = sample,
#'   precursor = precursor,
#'   peptide = peptide,
#'   intensity_log2 = intensity,
#'   method = "sum",
#'   for_plot = FALSE
#' )
#'
#' peptide_abundance
#'
#' # abundances that can be used in a peptide profile plot
#' complete_abundances <- calculate_peptide_abundance(
#'   data,
#'   sample = sample,
#'   precursor = precursor,
#'   peptide = peptide,
#'   intensity_log2 = intensity,
#'   method = "sum",
#'   for_plot = TRUE
#' )
#'
#' complete_abundances
#' }
calculate_peptide_abundance <- function(data,
                                        sample,
                                        peptide,
                                        precursor,
                                        intensity_log2,
                                        min_n_peptides = 3,
                                        method = "sum",
                                        for_plot = FALSE,
                                        retain_columns = NULL) {
  . <- NULL
  
  # Filter out any proteins with less than 3 peptides
  input <- data %>%
    dplyr::ungroup() %>%
    dplyr::distinct(
      {{ sample }},
      {{ precursor }},
      {{ peptide }},
      {{ intensity_log2 }}
    ) %>%
    tidyr::drop_na() %>%
    dplyr::group_by({{ peptide }}, {{ sample }}) %>%
    dplyr::ungroup()
  
  if (method == "sum") {
    result <- input %>%
      dplyr::group_by({{ sample }}, {{ peptide }}) %>%
      dplyr::summarise({{ intensity_log2 }} := log2(sum(2^{{ intensity_log2 }})), .groups = "drop")
    
    if (missing(retain_columns) & for_plot == FALSE) {
      return(result)
    }
    
    combined <- result %>%
      dplyr::mutate({{ precursor }} := "preptide_intensity") %>%
      dplyr::bind_rows(input)
    
    if (missing(retain_columns) & for_plot == TRUE) {
      return(combined)
    }
  }
  
  # needs to be checked
  if (method == "iq") {
    if (!requireNamespace("iq", quietly = TRUE)) {
      message("Package \"iq\" is needed for this function to work. Please install it.", call. = FALSE)
      return(invisible(NULL))
    }
    pb <- progress::progress_bar$new(
      total = length(unique(dplyr::pull(input, {{ protein_id }}))),
      format = "Preparing data [:bar] :current/:total (:percent) :eta"
    )
    
    input <- input %>%
      dplyr::distinct({{ sample }}, {{ peptide }}, {{ precursor }}, {{ intensity_log2 }}) %>%
      tidyr::complete(!!rlang::ensym(sample), nesting(!!rlang::ensym(precursor), !!rlang::ensym(peptide))) %>%
      split(dplyr::pull(., {{ peptide }})) %>%
      purrr::map(.f = ~ {
        pb$tick()
        .x %>%
          dplyr::select(-{{ peptide }}) %>%
          tidyr::pivot_wider(names_from = {{ sample }}, values_from = {{ intensity_log2 }}) %>%
          tibble::column_to_rownames(rlang::as_name(rlang::enquo(precursor))) %>%
          as.matrix()
      })
    
    pb <- progress::progress_bar$new(
      total = length(input),
      format = "Applying maximal peptide ratio extraction algorithm [:bar] :current/:total (:percent) :eta"
    )
    
    combined <- input %>%
      purrr::map2_df(
        .y = names(.),
        .f = ~ {
          pb$tick()
          iq::maxLFQ(.x) %>%
            purrr::pluck("estimate") %>%
            matrix(
              ncol = ncol(.x),
              nrow = 1,
              dimnames = list("peptide_intensity", colnames(.x))
            ) %>%
            rbind(.x) %>%
            tibble::as_tibble(rownames = NA) %>%
            tibble::rownames_to_column(var = rlang::as_name(rlang::enquo(precursor))) %>%
            tidyr::pivot_longer(-{{ precursor }},
                                names_to = rlang::as_name(rlang::enquo(sample)),
                                values_to = rlang::as_name(rlang::enquo(intensity_log2))
            ) %>%
            dplyr::mutate({{ protein_id }} := .y)
        }
      ) %>%
      tidyr::drop_na()
    
    if (missing(retain_columns) & for_plot == TRUE) {
      return(combined)
    }
    
    result <- combined %>%
      dplyr::filter({{ precursor }} == "peptide_intensity") %>%
      dplyr::select(-{{ precursor }})
  }
  
  if (!missing(retain_columns)) {
    peptide_intensity_retain <- data %>%
      dplyr::select(
        !!enquo(retain_columns),
        colnames(combined)[!colnames(combined) %in%
                             c(
                               rlang::as_name(rlang::enquo(intensity_log2)),
                               rlang::as_name(rlang::enquo(precursor))
                             )]
      ) %>%
      dplyr::distinct() %>%
      dplyr::mutate({{ precursor }} := "peptide_intensity")
  }
  
  if (!missing(retain_columns) & for_plot == FALSE) {
    result <- data %>%
      dplyr::select(
        !!enquo(retain_columns),
        colnames(result)[!colnames(result) %in%
                           c(rlang::as_name(rlang::enquo(intensity_log2)))]
      ) %>%
      dplyr::distinct() %>%
      dplyr::right_join(result, by = colnames(result)[!colnames(result) %in%
                                                        c(rlang::as_name(rlang::enquo(intensity_log2)))])
    
    return(result)
  }
  if (!missing(retain_columns) & for_plot == TRUE) {
    combined <- data %>%
      dplyr::select(
        !!enquo(retain_columns),
        colnames(combined)[!colnames(combined) %in%
                             c(
                               rlang::as_name(rlang::enquo(intensity_log2))
                             )]
      ) %>%
      dplyr::distinct() %>%
      dplyr::bind_rows(peptide_intensity_retain) %>%
      dplyr::right_join(combined, by = colnames(combined)[!colnames(combined) %in%
                                                            c(
                                                              rlang::as_name(rlang::enquo(intensity_log2))
                                                            )])
    
    return(combined)
  }
  
  return(result)
}