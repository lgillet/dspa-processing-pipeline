impute_without_comparison <- function(data,
                   sample,
                   grouping,
                   intensity_log2,
                   condition,
                   missingness = missingness,
                   noise = NULL,
                   method = "ludovic",
                   skip_log2_transform_error = FALSE,
                   retain_columns = NULL) {
  noise_missing <- missing(noise) # check if argument noise was provided or not
  result <- data %>%
    dplyr::distinct(
      {{ sample }},
      {{ grouping }},
      {{ intensity_log2 }},
      {{ condition }},
      {{ missingness }},
      {{ noise }}
    ) %>%
    dplyr::group_by({{ grouping }}, {{ condition }}) %>%
    dplyr::mutate(mean = mean({{ intensity_log2 }}, na.rm = TRUE)) %>%
    dplyr::mutate(sd = stats::sd({{ intensity_log2 }}, na.rm = TRUE)) %>%
    dplyr::group_by({{ grouping }}) %>%
    dplyr::mutate(min = min({{ intensity_log2 }}, na.rm = TRUE)) %>%
    dplyr::mutate(noise_mean = ifelse(noise_missing, NA, mean({{ noise }}, na.rm = TRUE))) %>%
    dplyr::mutate(sd = mean(unique(.data$sd), na.rm = TRUE)) %>%
    dplyr::group_by({{ grouping }}, {{ sample }}) %>%
    dplyr::mutate(
      impute =
        calculate_imputation(
          min = .data$min,
          noise = .data$noise_mean,
          mean = mean,
          sd = .data$sd,
          missingness = {{ missingness }},
          method = method,
          skip_log2_transform_error = skip_log2_transform_error
        )
    ) %>%
    dplyr::group_by({{ grouping }}, {{ sample }}, {{ missingness }}) %>%
    dplyr::mutate(impute = .data$impute[1]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(imputed_intensity = ifelse(is.na({{ intensity_log2 }}) == TRUE,
                                             .data$impute,
                                             {{ intensity_log2 }}
    )) %>%
    dplyr::mutate(imputed = is.na({{ intensity_log2 }}) & !is.na(.data$imputed_intensity)) %>%
    dplyr::select(-c("impute", "mean", "sd", "min", "noise_mean")) %>%
    dplyr::arrange(factor({{ grouping }}, levels = unique(stringr::str_sort({{ grouping }}, numeric = TRUE))))
  
  if (missing(retain_columns)) {
    return(result)
  } else {
    result <- data %>%
      dplyr::select(
        !!enquo(retain_columns),
        colnames(result)[!colnames(result) %in% c(
          "imputed_intensity",
          "imputed"
        )]
      ) %>%
      dplyr::distinct() %>%
      dplyr::right_join(result, by = colnames(result)[!colnames(result) %in% c(
        "imputed_intensity",
        "imputed"
      )]) %>%
      dplyr::arrange(factor({{ grouping }}, levels = unique(stringr::str_sort({{ grouping }}, numeric = TRUE))))
    return(result)
  }
}

calculate_imputation <-
  function(min = NULL,
           noise = NULL,
           mean = NULL,
           sd,
           missingness = c("MNAR", "MAR"),
           method = c("ludovic", "noise"),
           skip_log2_transform_error = FALSE) {
    if ((ifelse(is.na(ifelse(is.null(min), 0, min) > 40),
                FALSE,
                ifelse(is.null(min), 0, min) > 40
    ) |
    ifelse(is.na(ifelse(is.null(mean), 0, mean) > 40),
           FALSE,
           ifelse(is.null(mean), 0, mean) > 40
    ) |
    ifelse(is.na(ifelse(is.null(noise), 0, noise) > 40),
           FALSE,
           ifelse(is.null(noise), 0, noise) > 40
    )) &
    skip_log2_transform_error == FALSE) {
      stop(strwrap("Input intensities seem not to be log2 transformed. If they are and you want
                   to proceed set the skip_log2_transform_error argument to TRUE. Notice that
                   this function does not give correct results for non-log2 transformed data.",
                   prefix = "\n", initial = ""
      ))
    }
    if (!(missingness %in% c("MNAR", "MAR"))) {
      return(NA)
    }
    if (method == "ludovic") {
      if (missingness == "MNAR") {
        result <- suppressWarnings(stats::rnorm(1, mean = min - 3, sd = sd))
      }
      if (missingness == "MAR") {
        result <- suppressWarnings(stats::rnorm(1, mean = mean, sd = sd))
      }
    }
    if (method == "noise") {
      if (missingness == "MNAR") {
        result <- suppressWarnings(stats::rnorm(1, mean = noise, sd = sd))
      }
      if (missingness == "MAR") {
        result <- suppressWarnings(stats::rnorm(1, mean = mean, sd = sd))
      }
    }
    result
  }