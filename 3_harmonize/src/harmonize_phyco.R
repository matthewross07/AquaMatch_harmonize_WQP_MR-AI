#' @title Harmonize phycocyanin data from the Water Quality Portal
#'
#' @description
#' Harmonizes raw phycocyanin records retrieved from the WQP after
#' pre-harmonization (`clean_wqp_data()`), separating incompatible measurement
#' scales (`Phycocyanin` extractive laboratory `ug/L`, `Phycocyanin (probe)`
#' optical sensor `ug/L`, `Phycocyanin relative fluorescence` `RFU`/`IVFU`, and
#' `Phycocyanin cell density` `#/mL`), filtering blanks, periphyton, mislabeled
#' phycoerythrin, unit-misfiled batches, and non-mean daily summary statistics,
#' assigning QA/QC flags (`mdl_flag`, `approx_flag`, `greater_flag`,
#' `depth_flag`, `field_flag`, `misc_flag`), tiering analytical methods (`0` =
#' Restrictive, `1` = Narrowed, `2` = Inclusive), aggregating simultaneous
#' observations, and exporting summary tables and diagnostic figures.
#'
#' @param raw_phyco Data frame containing pre-harmonized phycocyanin records.
#' @param p_codes Data frame of USGS parameter codes (`p3_p_codes`).
#'
#' @returns
#' A named list of file paths to outputs generated during harmonization:
#' `phyco_param_drop_record_path`, `phyco_tiering_record_path`,
#' `phyco_misc_flag_table_path`, `phyco_grouped_preagg_path`,
#' `phyco_harmonized_path`, and `compiled_drops_path`.
#'
harmonize_phyco <- function(raw_phyco, p_codes){
  
  # Starting values for dataset
  starting_data <- tibble(
    step = "phyco harmonization",
    reason = "Starting dataset",
    short_reason = "Start",
    number_dropped = 0,
    n_rows = nrow(raw_phyco),
    order = 0
  )
  
  # Minor data prep ---------------------------------------------------------
  
  # Produce a bar chart of the current CharacteristicName counts
  stack_chart <- raw_phyco %>%
    count(CharacteristicName) %>%
    arrange(desc(CharacteristicName)) %>%
    mutate(cume_pos = cumsum(n) - n / 2) %>%
    ggplot(aes(x = 1, y = n)) +
    geom_bar(aes(fill = CharacteristicName),
             color = "black", position = "stack", stat = "identity",
             alpha = 0.85) +
    geom_text_repel(
      aes(
        label = paste0(CharacteristicName, ": n = ", comma(n)),
        y = cume_pos
      ),
      hjust = 0,
      force_pull = 0,
      nudge_x = 1,
      min.segment.length = 0.75,
      direction = "y",
      show.legend = FALSE
    ) +
    xlim(c(0, 10)) +
    scale_fill_viridis_d() +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    xlab(NULL) +
    ylab("Cumulative record count") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  ggsave(filename = "3_harmonize/out/phyco_stacked_char_names.png",
         plot = stack_chart, width = 6.5, height = 7.5, units = "in",
         device = "png")
  
  # Grab the column names of the incoming dataset
  raw_names <- names(raw_phyco)
  
  # Join USGS p-codes, recover missing measure units from detection limit
  # metadata, and filter media, blanks, periphyton, non-target pigments,
  # known unit-misfiled batches, and non-mean daily summary statistics
  lab_method_pattern <- paste0(
    c("Otsuki", "Kasinak", "Zimba", "Horvath", "Sarada", "Bennett",
      "Phycobiliproteins,wf", "Phycocyanin, wf", "MCGLYNN METHOD-PC",
      "Phycocyanin in Water by Fluorometry", "2338\\.22A",
      "PHYCOCYANIN-MDEQ-WQ", "^EPA 445$", "freeze[- ]thaw",
      "GF/F filter"),
    collapse = "|"
  )

  probe_indicator_pattern <- paste0(
    c("probe", "sonde", "sensor", "in[- ]situ", "\\bEXO", "\\bYSI\\b",
      "\\b6131\\b", "\\b6600\\b", "ProDSS", "Cyclops", "Trilogy", "SCUFA",
      "AquaFluor", "Hydrolab", "Eureka", "Manta", "Aqua ?TROLL", "FluoroProbe",
      "Data Logger", "Continuous Time Series", "Field Msr/Obs"),
    collapse = "|"
  )

  phyco_narrowed <- raw_phyco %>%
    left_join(x = ., y = p_codes, by = c("USGSPCode" = "parm_cd")) %>%
    mutate(
      ResultMeasure.MeasureUnitCode = coalesce(
        na_if(trimws(ResultMeasure.MeasureUnitCode), ""),
        na_if(trimws(DetectionQuantitationLimitMeasure.MeasureUnitCode), "")
      ),
      is_water_media = (ActivityMediaSubdivisionName %in%
                          c("Surface Water", "Water", "Estuary") |
                          is.na(ActivityMediaSubdivisionName)) &
        !grepl("Blank", coalesce(ActivityTypeCode, ""), ignore.case = TRUE) &
        !grepl("Periphyton", coalesce(SampleCollectionEquipmentName, ""),
               ignore.case = TRUE),
      is_non_phyco = grepl(
        pattern = "BGA-PE|phycoerythrin|allophycocyanin|chlorophyll fluorescence",
        x = paste(ResultCommentText,
                  ResultAnalyticalMethod.MethodName,
                  parameter_name_description),
        ignore.case = TRUE
      ),
      # CA_BVR 2021-2022 submitted integer cell-density counts (#/mL) under
      # ug/L prior to updating WQX unit codes to #/mL in 2023
      is_mislabeled_unit_batch = OrganizationIdentifier == "CA_BVR" &
        year %in% c(2021, 2022) &
        ResultMeasure.MeasureUnitCode %in% c("ug/L", "ug/l"),
      is_valid_stat = is.na(StatisticalBaseCode) |
        StatisticalBaseCode %in% c("Mean"),
      is_mass_unit = ResultMeasure.MeasureUnitCode %in% c(
        "ug/l", "ug/L", "mg/m3", "ppb", "mg/l", "mg/L", "ppm",
        "ug/ml", "ug/mL", "mg/cm3"
      ),
      is_lab_assay = is_mass_unit &
        !(CharacteristicName %in% "Phycocyanin (probe)") &
        !(USGSPCode %in% c("32319", "95203")) &
        (
          (USGSPCode %in% "52905") |
            grepl(
              pattern = lab_method_pattern,
              x = paste(ResultAnalyticalMethod.MethodIdentifier,
                        ResultAnalyticalMethod.MethodName,
                        ResultCommentText),
              ignore.case = TRUE
            ) |
            (!grepl(
              pattern = probe_indicator_pattern,
              x = paste(CharacteristicName, ActivityTypeCode,
                        SampleCollectionEquipmentName,
                        SampleCollectionMethod.MethodName,
                        ResultAnalyticalMethod.MethodName),
              ignore.case = TRUE
            )) |
            (grepl("Water Bottle",
                   coalesce(SampleCollectionEquipmentName, ""),
                   ignore.case = TRUE) &
               !grepl("probe|sonde|sensor",
                      paste(ResultAnalyticalMethod.MethodName,
                            SampleCollectionMethod.MethodName),
                      ignore.case = TRUE))
        ),
      parameter = case_when(
        is_lab_assay ~ "Phycocyanin",
        is_mass_unit & !is_lab_assay ~ "Phycocyanin (probe)",
        ResultMeasure.MeasureUnitCode %in% c("RFU", "IVFU") ~
          "Phycocyanin relative fluorescence",
        ResultMeasure.MeasureUnitCode %in% c("#/mL", "#/ml", "cells/mL") ~
          "Phycocyanin cell density",
        .default = "Unmapped unit"
      )
    ) %>%
    rowid_to_column(., "index")
  
  if (any(is.na(phyco_narrowed$parameter))) {
    stop("Unexpected NA values generated when classifying phycocyanin parameters.")
  }
  
  phyco <- phyco_narrowed %>%
    filter(is_water_media, !is_non_phyco, !is_mislabeled_unit_batch, is_valid_stat) %>%
    select(-is_water_media, -is_non_phyco, -is_mislabeled_unit_batch,
           -is_valid_stat, -is_mass_unit, -is_lab_assay)
  
  param_drop_record_out_path <- "3_harmonize/out/phyco_param_drop_record.csv"
  
  write_csv(
    x = phyco_narrowed %>%
      filter(!(is_water_media & !is_non_phyco &
                 !is_mislabeled_unit_batch & is_valid_stat)) %>%
      count(CharacteristicName, ResultAnalyticalMethod.MethodName,
            USGSPCode, ResultAnalyticalMethod.MethodIdentifier,
            ResultAnalyticalMethod.MethodIdentifierContext,
            ResultMeasure.MeasureUnitCode, ResultSampleFractionText) %>%
      arrange(desc(n)),
    file = param_drop_record_out_path
  )
  
  # Plot counts across parameter unit families (prior to dropping unmapped units)
  stack_param_chart <- phyco %>%
    count(parameter) %>%
    arrange(desc(parameter)) %>%
    mutate(cume_pos = cumsum(n) - n / 2) %>%
    ggplot(aes(x = 1, y = n)) +
    geom_bar(aes(fill = parameter),
             color = "black", position = "stack", stat = "identity",
             alpha = 0.85) +
    geom_text_repel(
      aes(
        label = paste0(parameter, ": n = ", comma(n)),
        y = cume_pos
      ),
      hjust = 0,
      force_pull = 0,
      nudge_x = 1,
      min.segment.length = 0.75,
      direction = "y",
      show.legend = FALSE
    ) +
    xlim(c(0, 10)) +
    scale_fill_viridis_d() +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    xlab(NULL) +
    ylab("Cumulative record count") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  ggsave(filename = "3_harmonize/out/phyco_stacked_param_names.png",
         plot = stack_param_chart, width = 6.5, height = 7.5, units = "in",
         device = "png")
  
  dropped_media <- tibble(
    step = "phyco harmonization",
    reason = "Filtered for water media, phycocyanin pigment targets, and valid stats",
    short_reason = "Target water media & PC filters",
    number_dropped = nrow(raw_phyco) - nrow(phyco),
    n_rows = nrow(phyco),
    order = 1
  )
  
  rm(raw_phyco, phyco_narrowed)
  gc()
  
  
  # Document and remove fail language & disqualifying QC codes ---------------
  
  # Clear stale fail pie charts before generating current run diagnostics
  unlink(list.files("3_harmonize/out", pattern = "^phyco_.*_fail_pie\\.png$",
                    full.names = TRUE))
  
  # Fail phrases checked across comment, condition, and original value fields.
  # Bare "unstable" is excluded so valid records noting unrelated DO probe drift
  # remain intact, while biofouling pre-clean audits and non-deployed checks
  # are removed per USGS TM 1-D3.
  fail_text <- c(
    "beyond accept", "cancelled", "contaminat", "error", "fail", "improper",
    "instrument down", "interference", "invalid", "no result", "no sample",
    "no test", "not accept", "outside of accept", "problem", "questionable",
    "suspect", "unable", "violation", "reject", "no data", "value extrapolated",
    "criteria exceed", "biased", "parameter not required", "not visited",
    "warm sample", "broken bottle", "broken cooler", "received broken",
    "broken sensor", "time not documented", "wiper stopped",
    "wiper not parked", "sonde out of water", "sonde stopped recording",
    "not logging", "not recording", "not deployed", "audit at streambank",
    "preclean", "pre-clean", "readings will not settle", "reading high on bga",
    "stopped working", "not working", "off scale", "censored"
  )
  
  disqual_codes <- c("FQC", "R", "SUS", "QCI")
  
  fail_counts <- list(
    "ActivityCommentText", "ResultLaboratoryCommentText",
    "ResultCommentText", "ResultMeasureValue_original",
    "ResultDetectionConditionText"
  ) %>%
    set_names() %>%
    map(
      .x = .,
      .f = ~ {
        col_name <- .x
        map_df(
          .x = fail_text,
          .f = ~ {
            hit_count <- phyco %>%
              filter(grepl(pattern = .x,
                           x = !!sym(col_name),
                           ignore.case = TRUE)) %>%
              nrow()
            tibble(
              word = .x,
              record_count = hit_count
            )
          }
        ) %>%
          filter(record_count > 0)
      }
    ) %>%
    keep(~ nrow(.) > 0)
  
  qual_fail_counts <- phyco %>%
    filter(coalesce(trimws(MeasureQualifierCode), "") %in% disqual_codes) %>%
    count(word = trimws(MeasureQualifierCode), name = "record_count")
  
  if (nrow(qual_fail_counts) > 0) {
    fail_counts$MeasureQualifierCode <- qual_fail_counts
  }
  
  walk2(
    .x = fail_counts,
    .y = names(fail_counts),
    .f = ~ ggsave(
      filename = paste0("3_harmonize/out/phyco_", .y, "_fail_pie.png"),
      plot = plot_fail_pie(dataset = .x, col_name = .y, text_size = 2.6),
      width = 7, height = 7, units = "in", device = "png"
    )
  )
  
  phyco_fails_removed <- phyco %>%
    filter(
      !(coalesce(trimws(MeasureQualifierCode), "") %in% disqual_codes),
      if_all(
        .cols = c(ActivityCommentText, ResultLaboratoryCommentText,
                  ResultCommentText, ResultMeasureValue_original,
                  ResultDetectionConditionText),
        .fns = ~ !grepl(
          pattern = paste0(fail_text, collapse = "|"),
          x = .x,
          ignore.case = TRUE
        )
      )
    )
  
  print(paste0(
    "Rows removed due to fail-related language or QC qualifiers: ",
    nrow(phyco) - nrow(phyco_fails_removed)
  ))
  
  dropped_fails <- tibble(
    step = "phyco harmonization",
    reason = "Dropped rows containing fail-related language or QC qualifier codes",
    short_reason = "Fails, etc.",
    number_dropped = nrow(phyco) - nrow(phyco_fails_removed),
    n_rows = nrow(phyco_fails_removed),
    order = 2
  )
  
  
  # Clean up MDLs -----------------------------------------------------------
  
  non_detect_text <- "non-detect|not detect|non detect|undetect|below|Present <QL"
  
  mdl_updates <- phyco_fails_removed %>%
    filter(
      is.na(ResultMeasureValue) |
        (ResultMeasureValue == 0 &
           (coalesce(trimws(MeasureQualifierCode), "") == "U" |
              grepl(non_detect_text,
                    coalesce(ResultDetectionConditionText, ""),
                    ignore.case = TRUE)))
    ) %>%
    mutate(
      mdl_vals = ifelse(
        test = ((is.na(ResultMeasureValue_original) | ResultMeasureValue == 0) &
                  (coalesce(trimws(MeasureQualifierCode), "") == "U" |
                     grepl(non_detect_text, ResultLaboratoryCommentText, ignore.case = TRUE) |
                     grepl(non_detect_text, ResultCommentText, ignore.case = TRUE) |
                     grepl(non_detect_text, ResultDetectionConditionText, ignore.case = TRUE))) |
          (grepl(non_detect_text, ResultMeasureValue_original, ignore.case = TRUE) |
             ResultMeasureValue_original %in% c("ND", "nd")),
        yes = DetectionQuantitationLimitMeasure.MeasureValue,
        no = ifelse(
          test = grepl("[0-9]", ResultMeasureValue_original) &
            grepl("<", ResultMeasureValue_original),
          yes = str_replace_all(
            string = ResultMeasureValue_original,
            pattern = c("\\<" = "", "\\*" = "", "\\=" = "")
          ),
          no = NA
        )
      ),
      mdl_units = ifelse(
        !is.na(mdl_vals),
        coalesce(
          na_if(trimws(DetectionQuantitationLimitMeasure.MeasureUnitCode), ""),
          na_if(trimws(ResultMeasure.MeasureUnitCode), "")
        ),
        ResultMeasure.MeasureUnitCode
      ),
      mdl_num = suppressWarnings(as.numeric(mdl_vals)),
      half = ifelse(!is.na(mdl_num) & mdl_num > 0, mdl_num / 2, NA_real_)
    )
  
  mdl_updates$std_value <- with(mdl_updates, suppressWarnings(runif(nrow(mdl_updates), 0, half)))
  mdl_updates$std_value[is.nan(mdl_updates$std_value)] <- NA
  
  mdl_updates <- mdl_updates %>%
    select(index, std_value, mdl_vals, mdl_units) %>%
    filter(!is.na(std_value))
  
  print(paste(
    round((nrow(mdl_updates)) / nrow(phyco_fails_removed) * 100, 1),
    "% of samples had values listed as being below a detection limit"
  ))
  
  phyco_mdls_added <- phyco_fails_removed %>%
    left_join(x = ., y = mdl_updates, by = "index") %>%
    mutate(
      harmonized_value = ifelse(index %in% mdl_updates$index, std_value, ResultMeasureValue),
      harmonized_unit_raw = ifelse(index %in% mdl_updates$index, mdl_units, ResultMeasure.MeasureUnitCode),
      mdl_flag = case_when(
        index %in% mdl_updates$index ~ 1,
        (!(index %in% mdl_updates$index) &
           DetectionQuantitationLimitMeasure.MeasureValue <= ResultMeasureValue) |
          (!(index %in% mdl_updates$index) &
             is.na(DetectionQuantitationLimitMeasure.MeasureValue)) ~ 0,
        DetectionQuantitationLimitMeasure.MeasureValue > ResultMeasureValue ~ 2,
        .default = 0
      )
    )
  
  dropped_mdls <- tibble(
    step = "phyco harmonization",
    reason = "Dropped rows while cleaning MDLs",
    short_reason = "Clean MDLs",
    number_dropped = nrow(phyco_fails_removed) - nrow(phyco_mdls_added),
    n_rows = nrow(phyco_mdls_added),
    order = 3
  )
  
  
  # Clean up approximated values --------------------------------------------
  
  approx_text <- "result approx|RESULT IS APPROX|value approx"
  
  phyco_approx <- phyco_mdls_added %>%
    filter(
      !index %in% mdl_updates$index,
      is.na(ResultMeasureValue) &
        grepl("[0-9]", ResultMeasureValue_original) &
        (grepl(approx_text, ResultLaboratoryCommentText, ignore.case = TRUE) |
           grepl(approx_text, ResultCommentText, ignore.case = TRUE) |
           grepl(approx_text, ResultDetectionConditionText, ignore.case = TRUE))
    )
  
  phyco_approx$approx_value <- suppressWarnings(
    as.numeric(str_replace_all(phyco_approx$ResultMeasureValue_original, c("\\*" = "")))
  )
  phyco_approx$approx_value[is.nan(phyco_approx$approx_value)] <- NA
  
  phyco_approx <- phyco_approx %>%
    select(approx_value, index)
  
  print(paste(
    round((nrow(phyco_approx)) / nrow(phyco_mdls_added) * 100, 3),
    "% of samples had values listed as approximated"
  ))
  
  phyco_approx_added <- phyco_mdls_added %>%
    left_join(x = ., y = phyco_approx, by = "index") %>%
    mutate(
      harmonized_value = ifelse(index %in% phyco_approx$index, approx_value, harmonized_value),
      approx_flag = ifelse(index %in% phyco_approx$index, 1, 0)
    )
  
  dropped_approximates <- tibble(
    step = "phyco harmonization",
    reason = "Dropped rows while cleaning approximate values",
    short_reason = "Clean approximates",
    number_dropped = nrow(phyco_mdls_added) - nrow(phyco_approx_added),
    n_rows = nrow(phyco_approx_added),
    order = 4
  )
  
  
  # Clean up "greater than" values ------------------------------------------
  
  greater_vals <- phyco_approx_added %>%
    filter((!index %in% mdl_updates$index) & (!index %in% phyco_approx$index)) %>%
    filter(is.na(ResultMeasureValue)) %>%
    mutate(
      greater_value = case_when(
        grepl("[0-9]", coalesce(ResultMeasureValue_original, "")) &
          grepl(">", coalesce(ResultMeasureValue_original, "")) ~
          suppressWarnings(as.numeric(
            str_replace_all(
              ResultMeasureValue_original,
              c("\\>" = "", "\\*" = "", "\\=" = "")
            )
          )),
        grepl("above|exceed|greater",
              coalesce(ResultDetectionConditionText, ""),
              ignore.case = TRUE) &
          !is.na(DetectionQuantitationLimitMeasure.MeasureValue) ~
          suppressWarnings(as.numeric(DetectionQuantitationLimitMeasure.MeasureValue)),
        .default = NA_real_
      ),
      greater_units = coalesce(
        na_if(trimws(DetectionQuantitationLimitMeasure.MeasureUnitCode), ""),
        harmonized_unit_raw
      )
    ) %>%
    filter(!is.na(greater_value) & !is.nan(greater_value)) %>%
    select(greater_value, greater_units, index)
  
  print(paste(
    round((nrow(greater_vals)) / nrow(phyco_approx_added) * 100, 9),
    "% of samples had values listed as being above a detection limit//greater than"
  ))
  
  phyco_harmonized_values <- phyco_approx_added %>%
    left_join(x = ., y = greater_vals, by = "index") %>%
    mutate(
      harmonized_value = ifelse(index %in% greater_vals$index, greater_value, harmonized_value),
      harmonized_unit_raw = ifelse(index %in% greater_vals$index, greater_units, harmonized_unit_raw),
      greater_flag = ifelse(index %in% greater_vals$index, 1, 0)
    )
  
  dropped_greater_than <- tibble(
    step = "phyco harmonization",
    reason = "Dropped rows while cleaning 'greater than' values",
    short_reason = "Greater thans",
    number_dropped = nrow(phyco_approx_added) - nrow(phyco_harmonized_values),
    n_rows = nrow(phyco_harmonized_values),
    order = 5
  )
  
  rm(phyco)
  gc()
  
  
  # Remove remaining NAs and negative values --------------------------------
  
  phyco_no_na <- phyco_harmonized_values %>%
    filter(
      !is.na(harmonized_value),
      !is.nan(harmonized_value),
      harmonized_value >= 0
    )
  
  dropped_na <- tibble(
    step = "phyco harmonization",
    reason = "Dropped unresolved NAs and negative values",
    short_reason = "Unresolved NAs",
    number_dropped = nrow(phyco_harmonized_values) - nrow(phyco_no_na),
    n_rows = nrow(phyco_no_na),
    order = 6
  )
  
  rm(phyco_harmonized_values, phyco_approx_added, phyco_mdls_added,
     phyco_fails_removed)
  gc()
  
  
  # Harmonize value units ---------------------------------------------------
  
  unit_conversion_table <- tribble(
    ~ResultMeasure.MeasureUnitCode, ~harmonized_units,               ~conversion,
    "mg/l",                         "ug/L",                          1000,
    "mg/L",                         "ug/L",                          1000,
    "ppm",                          "ug/L",                          1000,
    "ug/l",                         "ug/L",                          1,
    "ug/L",                         "ug/L",                          1,
    "mg/m3",                        "ug/L",                          1,
    "ppb",                          "ug/L",                          1,
    "mg/cm3",                       "ug/L",                          1000000,
    "ug/ml",                        "ug/L",                          1000,
    "ug/mL",                        "ug/L",                          1000,
    "RFU",                          "RFU",                           1,
    "IVFU",                         "IVFU",                          1,
    "#/mL",                         "#/mL",                          1,
    "#/ml",                         "#/mL",                          1,
    "cells/mL",                     "#/mL",                          1
  )
  
  unit_table_out_path <- "3_harmonize/out/phyco_unit_table.csv"
  write_csv(x = unit_conversion_table, file = unit_table_out_path)
  
  # Plot dropped unit codes before inner_join removes unconvertible rows
  tryCatch({
    dropped_unit_counts <- phyco_no_na %>%
      anti_join(
        x = .,
        y = unit_conversion_table,
        by = c("harmonized_unit_raw" = "ResultMeasure.MeasureUnitCode")
      ) %>%
      transmute(ResultMeasure.MeasureUnitCode = harmonized_unit_raw) %>%
      count(ResultMeasure.MeasureUnitCode, name = "record_count")
    if (nrow(dropped_unit_counts) > 0) {
      ggsave(
        filename = "3_harmonize/out/phyco_unit_drop_pie.png",
        plot = plot_unit_pie(dropped_unit_counts),
        width = 6, height = 6, units = "in", device = "png"
      )
    }
  }, error = function(e) NULL)
  
  converted_units_phyco <- phyco_no_na %>%
    inner_join(
      x = .,
      y = unit_conversion_table,
      by = c("harmonized_unit_raw" = "ResultMeasure.MeasureUnitCode")
    ) %>%
    mutate(
      harmonized_value = harmonized_value * conversion,
      parameter = case_when(
        parameter == "Unmapped unit" & harmonized_units == "ug/L" &
          ((USGSPCode %in% "52905") |
             grepl(lab_method_pattern,
                   paste(ResultAnalyticalMethod.MethodIdentifier,
                         ResultAnalyticalMethod.MethodName,
                         ResultCommentText),
                   ignore.case = TRUE)) ~ "Phycocyanin",
        parameter == "Unmapped unit" & harmonized_units == "ug/L" ~ "Phycocyanin (probe)",
        harmonized_units %in% c("RFU", "IVFU") ~ "Phycocyanin relative fluorescence",
        harmonized_units == "#/mL" ~ "Phycocyanin cell density",
        .default = parameter
      )
    ) %>%
    select(-harmonized_unit_raw)
  
  if (any(converted_units_phyco$parameter == "Unmapped unit")) {
    stop("Unresolved parameter labels remain after unit conversion.")
  }
  
  print(paste0(
    "Rows removed while harmonizing units: ",
    nrow(phyco_no_na) - nrow(converted_units_phyco)
  ))
  
  dropped_harmonization <- tibble(
    step = "phyco harmonization",
    reason = "Dropped rows while harmonizing units",
    short_reason = "Harmonize units",
    number_dropped = nrow(phyco_no_na) - nrow(converted_units_phyco),
    n_rows = nrow(converted_units_phyco),
    order = 7
  )
  
  
  # Clean and flag depth data -----------------------------------------------
  
  recode_depth_na_phyco <- converted_units_phyco %>%
    mutate(
      across(
        .cols = c(ActivityDepthHeightMeasure.MeasureValue,
                  ResultDepthHeightMeasure.MeasureValue,
                  ActivityTopDepthHeightMeasure.MeasureValue,
                  ActivityBottomDepthHeightMeasure.MeasureValue),
        .fns = ~ if_else(
          condition = .x %in% c("NA", "999", "-999", "9999", "-9999", "-99", "99", "NaN"),
          true = NA_character_,
          false = .x
        )
      )
    )
  
  depth_unit_conversion_table <- tibble(
    depth_units = c("in", "ft", "feet", "cm", "m", "meters"),
    depth_conversion = c(0.0254, 0.3048, 0.3048, 0.01, 1, 1)
  )
  
  converted_depth_units_phyco <- recode_depth_na_phyco %>%
    left_join(
      x = .,
      y = depth_unit_conversion_table,
      by = c("ActivityDepthHeightMeasure.MeasureUnitCode" = "depth_units")
    ) %>%
    mutate(
      harmonized_activity_depth_value = suppressWarnings(
        as.numeric(ActivityDepthHeightMeasure.MeasureValue) * depth_conversion
      )
    ) %>%
    select(-depth_conversion) %>%
    left_join(
      x = .,
      y = depth_unit_conversion_table,
      by = c("ResultDepthHeightMeasure.MeasureUnitCode" = "depth_units")
    ) %>%
    mutate(
      harmonized_result_depth_value = suppressWarnings(
        as.numeric(ResultDepthHeightMeasure.MeasureValue) * depth_conversion
      )
    ) %>%
    select(-depth_conversion) %>%
    left_join(
      x = .,
      y = depth_unit_conversion_table,
      by = c("ActivityTopDepthHeightMeasure.MeasureUnitCode" = "depth_units")
    ) %>%
    mutate(
      harmonized_top_depth_value = suppressWarnings(
        as.numeric(ActivityTopDepthHeightMeasure.MeasureValue) * depth_conversion
      ),
      harmonized_top_depth_unit = "m"
    ) %>%
    select(-depth_conversion) %>%
    left_join(
      x = .,
      y = depth_unit_conversion_table,
      by = c("ActivityBottomDepthHeightMeasure.MeasureUnitCode" = "depth_units")
    ) %>%
    mutate(
      harmonized_bottom_depth_value = suppressWarnings(
        as.numeric(ActivityBottomDepthHeightMeasure.MeasureValue) * depth_conversion
      ),
      harmonized_bottom_depth_unit = "m"
    )
  
  harmonized_depth_phyco <- converted_depth_units_phyco %>%
    rowwise() %>%
    mutate(
      harmonized_discrete_depth_value = case_when(
        !is.na(harmonized_activity_depth_value) &
          is.na(harmonized_result_depth_value) ~ harmonized_activity_depth_value,
        is.na(harmonized_activity_depth_value) &
          !is.na(harmonized_result_depth_value) ~ harmonized_result_depth_value,
        (!is.na(harmonized_activity_depth_value) &
           !is.na(harmonized_result_depth_value)) &
          harmonized_activity_depth_value != harmonized_result_depth_value ~ mean(
            c(harmonized_activity_depth_value, harmonized_result_depth_value)
          ),
        harmonized_activity_depth_value == harmonized_result_depth_value ~ harmonized_activity_depth_value,
        .default = NA_real_
      ),
      harmonized_discrete_depth_unit = "m"
    ) %>%
    ungroup()
  
  flagged_depth_phyco <- harmonized_depth_phyco %>%
    mutate(
      depth_flag = case_when(
        is.na(harmonized_discrete_depth_value) &
          is.na(harmonized_top_depth_value) &
          is.na(harmonized_bottom_depth_value) ~ 0,
        !is.na(harmonized_discrete_depth_value) &
          !is.na(harmonized_top_depth_value) &
          !is.na(harmonized_bottom_depth_value) ~ 3,
        (!is.na(harmonized_top_depth_value) |
           !is.na(harmonized_bottom_depth_value)) &
          is.na(harmonized_discrete_depth_value) ~ 2,
        !is.na(harmonized_discrete_depth_value) &
          is.na(harmonized_top_depth_value) &
          is.na(harmonized_bottom_depth_value) ~ 1,
        ((!is.na(harmonized_activity_depth_value) | !is.na(harmonized_result_depth_value)) &
           !is.na(harmonized_top_depth_value)) |
          ((!is.na(harmonized_activity_depth_value) | !is.na(harmonized_result_depth_value)) &
             !is.na(harmonized_bottom_depth_value)) ~ 3,
        .default = NA_real_
      )
    ) %>%
    select(-c(harmonized_activity_depth_value, harmonized_result_depth_value,
              depth_conversion))
  
  depth_check_table <- flagged_depth_phyco %>%
    mutate(
      three_cols_present = if_else(
        !is.na(harmonized_discrete_depth_value) &
          !is.na(harmonized_top_depth_value) &
          !is.na(harmonized_bottom_depth_value),
        true = 1, false = 0
      ),
      only_discrete = if_else(
        !is.na(harmonized_discrete_depth_value) &
          is.na(harmonized_top_depth_value) &
          is.na(harmonized_bottom_depth_value),
        true = 1, false = 0
      ),
      only_top = if_else(
        is.na(harmonized_discrete_depth_value) &
          !is.na(harmonized_top_depth_value) &
          is.na(harmonized_bottom_depth_value),
        true = 1, false = 0
      ),
      only_bottom = if_else(
        is.na(harmonized_discrete_depth_value) &
          is.na(harmonized_top_depth_value) &
          !is.na(harmonized_bottom_depth_value),
        true = 1, false = 0
      ),
      fully_integrated = if_else(
        is.na(harmonized_discrete_depth_value) &
          !is.na(harmonized_top_depth_value) &
          !is.na(harmonized_bottom_depth_value),
        true = 1, false = 0
      ),
      no_depths = if_else(
        is.na(harmonized_discrete_depth_value) &
          is.na(harmonized_top_depth_value) &
          is.na(harmonized_bottom_depth_value),
        true = 1, false = 0
      ),
      discrete_partial_integ = if_else(
        (!is.na(harmonized_discrete_depth_value) &
           !is.na(harmonized_top_depth_value) &
           is.na(harmonized_bottom_depth_value)) |
          (!is.na(harmonized_discrete_depth_value) &
             is.na(harmonized_top_depth_value) &
             !is.na(harmonized_bottom_depth_value)),
        true = 1, false = 0
      )
    ) %>%
    count(three_cols_present, only_discrete, discrete_partial_integ,
          only_top, only_bottom, fully_integrated, no_depths, depth_flag,
          harmonized_units) %>%
    arrange(depth_flag, harmonized_units)
  
  write_csv(x = depth_check_table, file = "3_harmonize/out/phyco_depth_check_table.csv")
  
  depth_counts <- flagged_depth_phyco %>%
    mutate(depth_agg_flag = case_when(
      depth_flag == 1 & harmonized_discrete_depth_value <= 1 ~ "<=1m",
      depth_flag == 1 & harmonized_discrete_depth_value <= 5 ~ "<=5m",
      depth_flag == 1 & harmonized_discrete_depth_value > 5 ~ ">5m",
      depth_flag == 2 & harmonized_bottom_depth_value <= 1 ~ "<=1m",
      depth_flag == 2 & harmonized_bottom_depth_value <= 5 ~ "<=5m",
      depth_flag == 2 & harmonized_bottom_depth_value > 5 ~ ">5m",
      .default = "No or inconsistent depth"
    )) %>%
    count(depth_agg_flag, harmonized_units)
  
  write_csv(x = depth_counts, file = "3_harmonize/out/phyco_depth_counts.csv")
  
  print(paste0(
    "Rows removed due to unmatched depth units: ",
    nrow(converted_units_phyco) - nrow(flagged_depth_phyco)
  ))
  
  dropped_depths <- tibble(
    step = "phyco harmonization",
    reason = "Dropped rows while cleaning depths",
    short_reason = "Clean depths",
    number_dropped = nrow(converted_units_phyco) - nrow(flagged_depth_phyco),
    n_rows = nrow(flagged_depth_phyco),
    order = 8
  )
  
  rm(converted_units_phyco, recode_depth_na_phyco, converted_depth_units_phyco,
     harmonized_depth_phyco)
  gc()
  
  
  # Filter and tier analytical methods --------------------------------------
  
  # Tier 0 (Restrictive): Verified extractive laboratory phycocyanin assays in
  #   ug/L with documented aqueous extraction and cell-lysis/spectral protocols
  #   (Otsuki et al. 1994, Zimba 2012, Horvath et al. 2013, Kasinak et al. 2015)
  # Tier 1 (Narrowed): Documented laboratory optical assays lacking verified
  #   extraction comparability, OR documented in situ optical sensor platforms/
  #   USGS sensor p-codes across probe ug/L, RFU/IVFU, and #/mL
  # Tier 2 (Inclusive): Missing/generic methods, acetone Chl-a methods (EPA 445)
  #   that denature phycobiliproteins, or unverified 0.5 mg/L LRL HPLC entries
  
  documented_sensor_pattern <- paste0(
    c("\\bYSI\\b", "\\bEXO[0-9]*\\b", "\\b6131\\b", "\\b6600\\b", "\\bProDSS\\b",
      "\\bTurner\\b", "\\bCyclops\\b", "\\bSCUFA\\b", "\\bTrilogy\\b",
      "\\bAquaFluor\\b", "\\bEureka\\b", "\\bManta\\b", "\\bHydrolab\\b",
      "\\bAqua ?TROLL\\b", "\\bFluoroProbe\\b", "Multiparameter.*Sonde",
      "^FL035$", "^ALG02$", "YSI - BGA PC"),
    collapse = "|"
  )
  
  # Remove clearly misassigned non-phycocyanin analytical methods (e.g., 4500-P B phosphorus)
  phyco_relevant_methods <- flagged_depth_phyco %>%
    filter(
      is.na(ResultAnalyticalMethod.MethodName) |
        !grepl("4500[- ]?P|Phosphorus",
               ResultAnalyticalMethod.MethodName,
               ignore.case = TRUE)
    )
  
  tiered_methods_phyco <- phyco_relevant_methods %>%
    mutate(
      is_chl_method = grepl(
        "445",
        paste(ResultAnalyticalMethod.MethodIdentifier,
              ResultAnalyticalMethod.MethodName),
        ignore.case = TRUE
      ),
      is_unverified_hplc = coalesce(ResultAnalyticalMethod.MethodIdentifier, "") ==
        "PHYCOCYANIN-MDEQ-WQ",
      is_extractive_lab = parameter == "Phycocyanin" &
        !is_chl_method &
        !is_unverified_hplc &
        (
          (USGSPCode %in% "52905") |
            (ResultAnalyticalMethod.MethodIdentifier %in%
               c("SP012", "FL083", "USGS PHYCOCYANIN", "2338.22A")) |
            grepl(
              pattern = "Otsuki|Kasinak|Zimba|Horvath|Sarada|fbu096|GF/F filter",
              x = paste(ResultAnalyticalMethod.MethodName,
                        MethodDescriptionText,
                        ResultCommentText),
              ignore.case = TRUE
            )
        ),
      is_doc_lab_optical = parameter == "Phycocyanin" &
        !is_extractive_lab &
        !is_chl_method &
        !is_unverified_hplc &
        (
          (ResultAnalyticalMethod.MethodIdentifier %in%
             c("MCGLYNN METHOD-PC", "FLUOROMETRY")) |
            grepl(
              pattern = "Spectromet|Spectrophoto|Fluoromet",
              x = coalesce(ResultAnalyticalMethod.MethodName, ""),
              ignore.case = TRUE
            )
        ),
      is_documented_sensor = parameter != "Phycocyanin" &
        !is_chl_method &
        (
          (USGSPCode %in% c("32315", "32319", "32321", "32329", "32335",
                            "95201", "95202", "95203", "95204")) |
            grepl(
              pattern = documented_sensor_pattern,
              x = paste(ResultAnalyticalMethod.MethodIdentifier,
                        ResultAnalyticalMethod.MethodName,
                        SampleCollectionMethod.MethodName,
                        if_else(
                          grepl("YSI - BGA PC", coalesce(ResultCommentText, ""),
                                ignore.case = TRUE),
                          ResultCommentText,
                          ""
                        )),
              ignore.case = TRUE,
              perl = TRUE
            )
        ),
      tier = case_when(
        is_extractive_lab ~ 0,
        is_doc_lab_optical | is_documented_sensor ~ 1,
        .default = 2
      )
    ) %>%
    select(-is_chl_method, -is_unverified_hplc, -is_extractive_lab,
           -is_doc_lab_optical, -is_documented_sensor)
  
  tiering_record <- tiered_methods_phyco %>%
    group_by(
      parameter, CharacteristicName, ResultAnalyticalMethod.MethodName,
      USGSPCode, ResultAnalyticalMethod.MethodIdentifier,
      ResultAnalyticalMethod.MethodIdentifierContext,
      ResultMeasure.MeasureUnitCode, ResultSampleFractionText, tier
    ) %>%
    summarize(
      min_value = min(harmonized_value, na.rm = TRUE),
      max_value = max(harmonized_value, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
  
  tiering_record_out_path <- "3_harmonize/out/phyco_tiering_record.csv"
  write_csv(x = tiering_record, file = tiering_record_out_path)
  
  if (nrow(phyco_relevant_methods) != nrow(tiered_methods_phyco)) {
    stop("Rows were lost during analytical method tiering. This is not expected.")
  }
  
  dropped_methods <- tibble(
    step = "phyco harmonization",
    reason = "Dropped rows while filtering and tiering analytical methods",
    short_reason = "Analytical methods",
    number_dropped = nrow(flagged_depth_phyco) - nrow(tiered_methods_phyco),
    n_rows = nrow(tiered_methods_phyco),
    order = 9
  )
  
  
  # Flag field methods ------------------------------------------------------
  
  # Check sampling modality consistency against assigned parameter & tier:
  # 0 = discrete grab/bottle/Kemmerer/Van Dorn/filter for lab Phycocyanin
  #     (Tiers 0-1), or in situ sonde / discrete non-dissolved subsample for
  #     documented optical sensors (Tier 1);
  # 1 = depth/width-integrated composite or pump sample (homogenizes surface
  #     scums and can shear filaments), or modality/fraction mismatch
  #     (Probe/Sensor on lab Phycocyanin, Dissolved fraction on optical sensor);
  # 2 = inconclusive (Tier 2 or Unknown/Other/NA collection metadata).
  composite_pump_pattern <- paste0(
    c("integrated", "composite", "multiple verticals", "equal width",
      "\\bewi\\b", "vinyl tube", "bag sampler", "pump"),
    collapse = "|"
  )
  discrete_lab_pattern <- paste0(
    c("grab", "bottle", "kemmerer", "van dorn", "bucket", "bail",
      "jar", "filter", "GF/F", "wading", "point", "water sampler"),
    collapse = "|"
  )
  discrete_probe_pattern <- paste0(
    c("probe", "sonde", "sensor", "logger", "time series",
      "field msr/obs", "point", "grab", "kemmerer", "van dorn",
      "weighted[- ]bottle", "bail"),
    collapse = "|"
  )
  
  field_flagged_phyco <- tiered_methods_phyco %>%
    mutate(
      coll_text_all = paste(
        ActivityTypeCode,
        SampleCollectionMethod.MethodName,
        SampleCollectionEquipmentName,
        ResultCommentText
      ),
      coll_text_meta = paste(
        ActivityTypeCode,
        SampleCollectionMethod.MethodName,
        SampleCollectionEquipmentName
      ),
      is_comp_pump = grepl(
        composite_pump_pattern, coll_text_meta,
        ignore.case = TRUE, perl = TRUE
      ),
      field_flag = case_when(
        tier == 2 ~ 2,
        parameter == "Phycocyanin" &
          grepl("Probe/Sensor", coalesce(SampleCollectionEquipmentName, ""),
                ignore.case = TRUE) ~ 1,
        parameter == "Phycocyanin" & is_comp_pump ~ 1,
        parameter == "Phycocyanin" &
          grepl(discrete_lab_pattern, coll_text_all,
                ignore.case = TRUE, perl = TRUE) ~ 0,
        parameter == "Phycocyanin" ~ 2,
        parameter != "Phycocyanin" &
          (ResultSampleFractionText %in% "Dissolved") ~ 1,
        parameter != "Phycocyanin" & is_comp_pump ~ 1,
        parameter != "Phycocyanin" &
          grepl(discrete_probe_pattern, coll_text_meta,
                ignore.case = TRUE, perl = TRUE) &
          !(coalesce(SampleCollectionEquipmentName, "") == "Unknown" &
              coalesce(ActivityTypeCode, "") %in% c("Unknown", "Sample-Routine")) ~ 0,
        .default = 2
      )
    ) %>%
    select(-coll_text_all, -coll_text_meta, -is_comp_pump)
  
  print(paste0(
    "Rows removed while assigning field flags: ",
    nrow(tiered_methods_phyco) - nrow(field_flagged_phyco)
  ))
  
  dropped_field <- tibble(
    step = "phyco harmonization",
    reason = "Dropped rows while assigning field flags",
    short_reason = "Field flagging",
    number_dropped = nrow(tiered_methods_phyco) - nrow(field_flagged_phyco),
    n_rows = nrow(field_flagged_phyco),
    order = 10
  )
  
  
  # Miscellaneous flag ------------------------------------------------------
  
  # misc_flag definitions for phycocyanin:
  # 0 = discrete/instantaneous observation within linear calibration bounds
  # 1 = daily "Mean" statistical summary (held out from replicate averaging)
  # 2 = extreme hypertrophic / sensor saturation-ceiling exceedance
  #     (> 1,000 ug/L for extractive lab Phycocyanin, > 100 ug/L for in situ
  #     Phycocyanin (probe), > 100 RFU, > 200 IVFU, or > 280,000 #/mL); held out
  #     from simultaneous replicate averaging matching AquaMatch convention
  misc_flagged_phyco <- field_flagged_phyco %>%
    mutate(
      misc_flag = case_when(
        (parameter == "Phycocyanin" & harmonized_value > 1000) |
          (parameter == "Phycocyanin (probe)" & harmonized_value > 100) |
          (harmonized_units == "RFU" & harmonized_value > 100) |
          (harmonized_units == "IVFU" & harmonized_value > 200) |
          (harmonized_units == "#/mL" & harmonized_value > 280000) ~ 2,
        StatisticalBaseCode %in% "Mean" ~ 1,
        .default = 0
      )
    )
  
  misc_flag_table_out_path <- "3_harmonize/out/phyco_misc_flag_table.csv"
  write_csv(
    x = misc_flagged_phyco %>%
      count(parameter, harmonized_units, misc_flag),
    file = misc_flag_table_out_path
  )
  
  dropped_misc <- tibble(
    step = "phyco harmonization",
    reason = "Dropped rows while assigning misc flags",
    short_reason = "Misc flagging",
    number_dropped = nrow(field_flagged_phyco) - nrow(misc_flagged_phyco),
    n_rows = nrow(misc_flagged_phyco),
    order = 11
  )
  
  
  # Unrealistic values ------------------------------------------------------
  
  # Apply modality-specific physical upper bounds:
  # - Extractive lab Phycocyanin > 50,000 ug/L (preserves diluted bloom/scum
  #   extracts up to 41,000 ug/L under misc_flag = 2)
  # - In situ Phycocyanin (probe) > 5,000 ug/L (above commercial fluorometer
  #   saturation ceilings; removes unit-entry errors such as 19,900 ug/L)
  # - Relative fluorescence > 2,000 RFU/IVFU; cell density > 10,000,000 #/mL
  # - Depths > 592 m (maximum US lake depth, Crater Lake OR)
  realistic_phyco <- misc_flagged_phyco %>%
    filter(
      !(parameter == "Phycocyanin" & harmonized_value > 50000),
      !(parameter == "Phycocyanin (probe)" & harmonized_value > 5000),
      !(parameter == "Phycocyanin relative fluorescence" & harmonized_value > 2000),
      !(parameter == "Phycocyanin cell density" & harmonized_value > 10000000),
      harmonized_top_depth_value <= 592 | is.na(harmonized_top_depth_value),
      harmonized_bottom_depth_value <= 592 | is.na(harmonized_bottom_depth_value),
      harmonized_discrete_depth_value <= 592 | is.na(harmonized_discrete_depth_value)
    )
  
  dropped_unreal <- tibble(
    step = "phyco harmonization",
    reason = "Dropped rows with unrealistic values",
    short_reason = "Unrealistic values",
    number_dropped = nrow(misc_flagged_phyco) - nrow(realistic_phyco),
    n_rows = nrow(realistic_phyco),
    order = 12
  )
  
  
  # Aggregate simultaneous records ------------------------------------------
  
  grouped_phyco <- realistic_phyco %>%
    filter(misc_flag == 0) %>%
    group_by(
      parameter, OrganizationIdentifier, MonitoringLocationIdentifier,
      MonitoringLocationTypeName, ResolvedMonitoringLocationTypeName,
      ActivityStartDate, ActivityStartTime.Time,
      ActivityStartTime.TimeZoneCode, harmonized_tz,
      harmonized_local_time, harmonized_utc, ActivityStartDateTime,
      harmonized_top_depth_value, harmonized_top_depth_unit,
      harmonized_bottom_depth_value, harmonized_bottom_depth_unit,
      harmonized_discrete_depth_value, harmonized_discrete_depth_unit,
      depth_flag, mdl_flag, approx_flag, greater_flag, tier, field_flag,
      misc_flag, harmonized_units
    ) %>%
    mutate(subgroup_id = cur_group_id())
  
  grouped_phyco_out_path <- "3_harmonize/out/phyco_harmonized_grouped.feather"
  
  grouped_phyco %>%
    select(
      all_of(c(raw_names,
               "parameter_code", "group_name", "parameter_name_description",
               "subgroup_id")),
      group_cols(),
      harmonized_value
    ) %>%
    write_feather(path = grouped_phyco_out_path)
  
  no_simul_phyco_partial <- grouped_phyco %>%
    group_by(subgroup_id, .add = TRUE) %>%
    summarize(
      harmonized_row_count = n(),
      harmonized_value_sd = sd(harmonized_value),
      harmonized_value = mean(harmonized_value),
      lon = unique(lon)[1],
      lat = unique(lat)[1],
      datum = unique(datum)[1]
    ) %>%
    mutate(
      harmonized_value_cv = if_else(
        harmonized_row_count > 1 & harmonized_value == 0,
        0,
        harmonized_value_sd / harmonized_value
      )
    ) %>%
    ungroup() %>%
    select(-harmonized_value_sd) %>%
    relocate(
      c(subgroup_id, harmonized_row_count, harmonized_units,
        harmonized_value, harmonized_value_cv, lat, lon, datum),
      .after = misc_flag
    )
  
  no_simul_phyco <- no_simul_phyco_partial %>%
    bind_rows(
      .,
      realistic_phyco %>%
        filter(misc_flag != 0) %>%
        select(any_of(names(no_simul_phyco_partial)))
    )
  
  rm(grouped_phyco)
  gc()
  
  # Post-aggregation diagnostic plots
  no_simul_phyco_tier_label <- no_simul_phyco %>%
    mutate(
      tier_label = case_when(
        tier == 0 ~ "Restrictive (Tier 0)",
        tier == 1 ~ "Narrowed (Tier 1)",
        tier == 2 ~ "Inclusive (Tier 2)"
      ),
      tier_label = factor(
        x = tier_label,
        levels = c("Restrictive (Tier 0)", "Narrowed (Tier 1)", "Inclusive (Tier 2)"),
        ordered = TRUE
      ),
      param_label = case_when(
        parameter == "Phycocyanin" ~ "Phycocyanin (lab)",
        parameter == "Phycocyanin (probe)" ~ "Phycocyanin (probe)",
        parameter == "Phycocyanin relative fluorescence" ~ "PC Rel. Fluor.",
        parameter == "Phycocyanin cell density" ~ "PC Cell Density",
        .default = parameter
      ),
      facet_label = paste0(param_label, " (", harmonized_units, ")")
    )
  
  # 1. Harmonized value distributions
  tier_dists <- no_simul_phyco_tier_label %>%
    select(parameter, facet_label, tier_label, harmonized_value) %>%
    filter(harmonized_value > 0) %>%
    ggplot() +
    geom_histogram(aes(harmonized_value, fill = parameter),
                   color = "black", alpha = 0.85, bins = 35) +
    facet_wrap(facet_label ~ tier_label, ncol = 3, scales = "free") +
    xlab("Harmonized values") +
    ylab("Record count") +
    ggtitle(
      label = "Distribution of positive harmonized phycocyanin values by parameter, unit, and tier"
    ) +
    scale_x_log10(label = label_scientific()) +
    scale_y_continuous(label = label_number(scale_cut = cut_short_scale())) +
    scale_fill_viridis_d() +
    theme_bw() +
    theme(
      strip.text = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      plot.margin = margin(10, 22, 10, 10),
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 2))
  
  ggsave(filename = "3_harmonize/out/phyco_tier_dists_postagg.png",
         plot = tier_dists,
         width = 9.5, height = 9.5, units = "in", device = "png")
  
  # 2. Harmonized CVs
  cv_plot_data <- no_simul_phyco_tier_label %>%
    select(param_label, tier_label, harmonized_value_cv) %>%
    filter(!is.na(harmonized_value_cv)) %>%
    mutate(plot_value = harmonized_value_cv + 0.001)
  
  tier_cv_dist <- cv_plot_data %>%
    ggplot() +
    geom_histogram(aes(plot_value), color = "black", fill = "white", bins = 30) +
    facet_wrap(param_label ~ tier_label, scales = "free_y", ncol = 3) +
    xlab(expression("Harmonized coefficient of variation (+ 0.001, " ~ log[10] ~ " scale)")) +
    ylab("Record count") +
    ggtitle(
      label = paste0(
        "Distribution of harmonized CVs (n = ",
        comma(nrow(cv_plot_data)),
        ")"
      ),
      subtitle = "0.001 added to each CV prior to log10 transformation"
    ) +
    scale_x_log10(label = label_scientific()) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 8),
      plot.margin = margin(10, 16, 10, 10)
    )
  
  ggsave(filename = "3_harmonize/out/phyco_tier_cv_dists_postagg.png",
         plot = tier_cv_dist,
         width = 9.5, height = 7.5, units = "in", device = "png")
  
  # 3. Maps
  plot_tier_maps(dataset = no_simul_phyco, custom_width = 9.5, custom_height = 10.5,
                 n_bins = 25, param_name = "phyco", flip_facets = TRUE,
                 legend_position = "bottom")
  
  # 4. Time charts
  plot_time_charts(dataset = no_simul_phyco, custom_width = 9, custom_height = 7,
                   year_seq = 5, param_name = "phyco", legend_position = "bottom",
                   scale_type = "free_y", n_cols = 2)
  
  # 5. Depth distributions
  plot_depth_hist <- function(df, depth_col, out_file) {
    sub_df <- df %>%
      filter(!is.na(.data[[depth_col]]))
    if (nrow(sub_df) == 0) {
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste("No non-missing values in", depth_col)) +
        theme_void()
      ggsave(filename = out_file, plot = p,
             width = 7, height = 4, units = "in", device = "png")
      return(invisible(NULL))
    }
    p <- sub_df %>%
      ggplot() +
      geom_histogram(aes(.data[[depth_col]], fill = tier_label),
                     color = "black", bins = 30) +
      facet_grid(
        cols = vars(ResolvedMonitoringLocationTypeName),
        rows = vars(param_label),
        scales = "free_y"
      ) +
      scale_fill_viridis_d("Tier", direction = -1) +
      xlab(paste0(depth_col, ", m")) +
      ylab("Record count") +
      ggtitle(paste0(depth_col, " distribution by parameter and location type")) +
      theme_bw() +
      theme(legend.position = "bottom")
    ggsave(filename = out_file, plot = p,
           width = 8.5, height = 7.5, units = "in", device = "png")
  }
  
  plot_depth_hist(no_simul_phyco_tier_label, "harmonized_top_depth_value",
                  "3_harmonize/out/phyco_tier_top_depth_dist_postagg.png")
  plot_depth_hist(no_simul_phyco_tier_label, "harmonized_bottom_depth_value",
                  "3_harmonize/out/phyco_tier_bottom_depth_dist_postagg.png")
  plot_depth_hist(no_simul_phyco_tier_label, "harmonized_discrete_depth_value",
                  "3_harmonize/out/phyco_tier_discrete_depth_dist_postagg.png")
  
  rm(no_simul_phyco_tier_label, cv_plot_data)
  gc()
  
  print(paste0(
    "Rows removed while aggregating simultaneous records: ",
    nrow(realistic_phyco) - nrow(no_simul_phyco)
  ))
  
  dropped_simul <- tibble(
    step = "phyco harmonization",
    reason = "Dropped rows while aggregating simultaneous records",
    short_reason = "Simultaneous records",
    number_dropped = nrow(realistic_phyco) - nrow(no_simul_phyco),
    n_rows = nrow(no_simul_phyco),
    order = 13
  )
  
  
  # Export ------------------------------------------------------------------
  
  compiled_dropped <- bind_rows(
    starting_data, dropped_media, dropped_fails, dropped_mdls,
    dropped_approximates, dropped_greater_than, dropped_na,
    dropped_harmonization, dropped_depths, dropped_methods,
    dropped_field, dropped_misc, dropped_unreal, dropped_simul
  )
  
  documented_drops_out_path <- "3_harmonize/out/phyco_harmonize_dropped_metadata.csv"
  write_csv(x = compiled_dropped, file = documented_drops_out_path)
  
  data_out_path <- "3_harmonize/out/phyco_harmonized_final.csv"
  write_csv(no_simul_phyco, data_out_path)
  
  print(paste0("Final number of records: ", nrow(no_simul_phyco)))
  
  return(list(
    phyco_param_drop_record_path = param_drop_record_out_path,
    phyco_tiering_record_path = tiering_record_out_path,
    phyco_misc_flag_table_path = misc_flag_table_out_path,
    phyco_grouped_preagg_path = grouped_phyco_out_path,
    phyco_harmonized_path = data_out_path,
    compiled_drops_path = documented_drops_out_path
  ))
}
