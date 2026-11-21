# Phycocyanin specific targets list for the harmonization step

p3_phyco_targets_list <- list(
  
  # Pre-harmonization data prep ---------------------------------------------
  
  # All columns in p2_wqp_data_aoi_phyco are of class character. Coerce select
  # columns back to numeric, retaining original entries in "_original" columns.
  tar_target(
    name = p3_wqp_data_aoi_formatted_phyco,
    command = format_columns(p2_wqp_data_aoi_phyco),
    format = "feather"
  ),
  
  # Join in column containing site type info
  tar_target(
    name = p3_wqp_data_aoi_sitetype_phyco,
    command = left_join(
      x = p3_wqp_data_aoi_formatted_phyco,
      y = p1_wqp_inventory_aoi_phyco %>%
        select(OrganizationIdentifier, MonitoringLocationIdentifier,
               ResolvedMonitoringLocationTypeName, CharacteristicName,
               OrganizationFormalName, ProviderName, MonitoringLocationTypeName) %>%
        distinct()
    )
  ),
  
  
  # Pre-harmonization -------------------------------------------------------
  
  # Time and time zone fills
  tar_target(
    name = p3_wqp_data_aoi_date_time_phyco,
    command = fill_date_time(dataset = p3_wqp_data_aoi_sitetype_phyco,
                             site_data = p2_site_counts_phyco),
    packages = c("tidyverse", "lutz", "sf", "sfheaders")
  ),
  
  tar_target(
    name = p3_wqp_data_aoi_ready_phyco,
    command = clean_wqp_data(wqp_data = p3_wqp_data_aoi_date_time_phyco,
                             char_names_crosswalk = p1_char_names_crosswalk_phyco,
                             site_data = p2_site_counts_phyco %>%
                               distinct(MonitoringLocationIdentifier, CharacteristicName,
                                        .keep_all = TRUE),
                             wqp_metadata = p1_wqp_inventory_aoi_phyco),
    packages = c("tidyverse", "feather")
  ),
  
  # Connect cleaned data output to the pipeline
  tar_file_read(
    name = p3_cleaned_wqp_data_phyco,
    command = p3_wqp_data_aoi_ready_phyco$wqp_data_clean_path,
    read = read_feather(path = !!.x),
    packages = "feather"
  ),
  
  
  # Harmonization -----------------------------------------------------------
  
  tar_target(
    name = p3_phyco_harmonized,
    command = harmonize_phyco(raw_phyco = p3_cleaned_wqp_data_phyco,
                              p_codes = p3_p_codes),
    packages = c("tidyverse", "feather", "ggrepel", "scales", "snakecase",
                 "sf", "sfheaders")
  ),
  
  # Record of rows excluded in early parameter/metadata filtering
  tar_file_read(
    name = p3_phyco_param_drop_record,
    command = p3_phyco_harmonized$phyco_param_drop_record_path,
    read = read_csv(file = !!.x, guess_max = Inf, show_col_types = FALSE)
  ),
  
  tar_file_read(
    name = p3_phyco_tiering_record,
    command = p3_phyco_harmonized$phyco_tiering_record_path,
    read = read_csv(file = !!.x, guess_max = Inf, show_col_types = FALSE)
  ),
  
  tar_file_read(
    name = p3_phyco_misc_flag_table,
    command = p3_phyco_harmonized$phyco_misc_flag_table_path,
    read = read_csv(file = !!.x, guess_max = Inf, show_col_types = FALSE)
  ),
  
  tar_file_read(
    name = p3_phyco_preharm_drops,
    command = p3_wqp_data_aoi_ready_phyco$compiled_drops_path,
    read = read_csv(file = !!.x, guess_max = Inf, show_col_types = FALSE)
  ),
  
  tar_file_read(
    name = p3_phyco_harmonize_drops,
    command = p3_phyco_harmonized$compiled_drops_path,
    read = read_csv(file = !!.x, guess_max = Inf, show_col_types = FALSE)
  ),
  
  # Compiled drop log across pre-harmonization and harmonization
  tar_target(
    name = p3_documented_drops_phyco,
    command = bind_rows(
      p3_phyco_preharm_drops,
      p3_phyco_harmonize_drops
    )
  ),
  
  # Export compiled drop log for downstream documentation
  tar_target(
    name = p3_documented_drops_phyco_drive_file,
    command = {
      p0_check_phyco_drive
      export_single_file(target = p3_documented_drops_phyco,
                         drive_path = p0_phyco_output_path,
                         stable = p0_harmonization_config$phyco_use_stable,
                         google_email = p0_harmonization_config$google_email,
                         date_stamp = p0_harmonization_config$phyco_stable_date,
                         local_drive = isTRUE(p0_harmonization_config$local_drive))
    },
    packages = c("tidyverse", "googledrive"),
    error = "stop"
  ),
  
  # Harmonized phycocyanin data containing grouping IDs for simultaneous
  # records, prior to aggregation (excludes misc_flag != 0)
  tar_file_read(
    name = p3_phyco_preagg_grouped,
    command = p3_phyco_harmonized$phyco_grouped_preagg_path,
    read = read_feather(path = !!.x),
    packages = "feather"
  ),
  
  # Harmonized phycocyanin data after simultaneous record aggregation
  tar_file_read(
    name = p3_phyco_agg_harmonized,
    command = p3_phyco_harmonized$phyco_harmonized_path,
    read = read_csv(file = !!.x, guess_max = Inf, show_col_types = FALSE)
  ),
  
  # Create a copy of the csv in feather format
  tar_file_read(
    name = p3_phyco_agg_harmonized_feather,
    command = {
      out_path <- gsub(x = p3_phyco_harmonized$phyco_harmonized_path,
                       pattern = ".csv",
                       replacement = ".feather")
      
      write_feather(x = p3_phyco_agg_harmonized,
                    path = out_path)
      
      out_path
    },
    packages = c("targets", "feather"),
    read = read_feather(path = !!.x)
  ),
  
  # Export
  tar_target(
    name = p3_phyco_agg_harmonized_feather_drive_file,
    command = {
      p0_check_phyco_drive
      export_single_file(target = p3_phyco_agg_harmonized_feather,
                         drive_path = p0_phyco_output_path,
                         stable = p0_harmonization_config$phyco_use_stable,
                         google_email = p0_harmonization_config$google_email,
                         date_stamp = p0_harmonization_config$phyco_stable_date,
                         local_drive = isTRUE(p0_harmonization_config$local_drive))
    },
    packages = c("tidyverse", "googledrive"),
    error = "stop"
  ),
  
  
  # Site info ---------------------------------------------------------------
  
  # Generate site metadata after harmonization is complete
  tar_file_read(
    name = p3_phyco_harmonized_site_info,
    command = {
      phyco_sites <- get_site_info(dataset = p3_phyco_preagg_grouped)
      
      out_path <- "3_harmonize/out/phyco_harmonized_site_info.feather"
      
      phyco_sites %>%
        write_feather(path = out_path)
      
      out_path
    },
    read = read_feather(path = !!.x),
    packages = c("tidyverse", "dataRetrieval", "feather")
  ),
  
  # Export
  tar_target(
    name = p3_phyco_site_info_drive_file,
    command = {
      p0_check_phyco_drive
      export_single_file(target = p3_phyco_harmonized_site_info,
                         drive_path = p0_phyco_output_path,
                         stable = p0_harmonization_config$phyco_use_stable,
                         google_email = p0_harmonization_config$google_email,
                         date_stamp = p0_harmonization_config$phyco_stable_date,
                         local_drive = isTRUE(p0_harmonization_config$local_drive))
    },
    packages = c("tidyverse", "googledrive"),
    error = "stop"
  ),
  
  # Get file IDs ------------------------------------------------------------
  
  tar_file_read(
    name = p3_phyco_drive_ids,
    command = get_file_ids(google_email = p0_harmonization_config$google_email,
                           drive_folder = p0_phyco_output_path,
                           file_path = "3_harmonize/out/phyco_drive_ids.csv",
                           depend = c(p3_phyco_agg_harmonized_feather_drive_file,
                                      p3_phyco_site_info_drive_file,
                                      p3_documented_drops_phyco_drive_file),
                           local_drive = isTRUE(p0_harmonization_config$local_drive)),
    read = read_csv(file = !!.x, guess_max = Inf, show_col_types = FALSE),
    packages = c("tidyverse", "googledrive")
  )
  
)
