#' @title Retrieve a dataset from Google Drive
#' 
#' @description
#' A function to retrieve a dataset from Google Drive after it has been uploaded
#' in a previous step.
#' 
#' @param target A string containing the name of the target to be retrieved.
#' 
#' @param id_df A dataframe containing the cols `name` and `id` to use for
#' retrieving the uploaded dataset from Google Drive.
#' 
#' @param local_folder A string specifying the folder where the file should be
#' downloaded to.
#' 
#' @param stable Logical value. If TRUE, look for file in the "stable" subfolder
#' in Google Drive. If FALSE, use the path as provided by the user.
#' 
#' @param google_email A string containing the gmail address to use for
#' Google Drive authentication.
#' 
#' @param file_type A string giving the file extension to be used. (".rds" or 
#' ".feather")
#' 
#' @param stable_date A string containing an eight-digit date (i.e., in
#' ISO 8601 "basic" format: YYYYMMDD) that should be used to identify the
#' correct file version on Google Drive.
#' 
#' @param local_drive Logical; if TRUE (or if the resolved ID is a local file
#' path), copy directly from the local filesystem mirror instead of calling
#' Google Drive.
#' 
#' @returns 
#' The dataset after being downloaded and read into the pipeline from Google Drive
#' or the local filesystem mirror.
#' 
retrieve_data <- function(target, id_df, local_folder, stable, 
                          google_email, file_type = ".rds", stable_date,
                          local_drive = FALSE){
  
  # Local file download location
  local_path <- file.path(local_folder, paste0(target, file_type))
  
  # Make the directory if it does not exist:
  if(!dir.exists(local_folder)) {dir.create(local_folder, recursive = TRUE)}
  
  # Get file contents of the Google Drive folder specified. If stable == TRUE
  # then append stable_date to match the dated snapshot filename.
  if(stable){
    
    file_name <- paste0(target, "_", stable_date, file_type)
    
  } else {
    
    file_name <- paste0(target, file_type)
    
  }
  
  # Filter the contents to the file requested and obtain its ID or local path
  raw_id <- id_df %>%
    filter(name == file_name) %>%
    pull(id)
  
  if (length(raw_id) == 0) {
    stop(sprintf("Could not find file '%s' in id_df.", file_name))
  }
  raw_id <- raw_id[1]
  
  if (isTRUE(local_drive) || file.exists(raw_id)) {
    if (!file.exists(raw_id)) {
      stop(sprintf("Local mirror file not found: %s", raw_id))
    }
    file.copy(from = raw_id, to = local_path, overwrite = TRUE)
  } else {
    if (is.null(google_email) || !nzchar(google_email)) {
      drive_deauth()
    } else {
      drive_auth(google_email)
    }
    
    # Run the download
    drive_download(file = as_id(raw_id),
                   path = local_path,
                   overwrite = TRUE)
  }
  
  # Read dataset into pipeline
  if(file_type == ".rds"){
    
    read_rds(local_path)
    
  } else if(file_type == ".feather"){
    
    read_feather(local_path)   
    
  } else {
    
    stop("file_type does not appear to be either .rds or .feather.")
    
  }
  
}
