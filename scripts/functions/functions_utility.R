#### Title: Utility functions
#### Author: María Bueno Álvez
#### Description: script collecting general functions for data manipulation
#### Last edited : 12/08/2024

# Utility packages
library(tidyverse)
library(impute)

# Function to import data frames
import_df <- function(file_path) {
  
  # Determine file extension from file path
  file_extension <- tools::file_ext(file_path)
  
  df <- switch(tolower(file_extension),
               csv = readr::read_delim(file_path),
               tsv = utils::read.delim(file_path, stringsAsFactors = FALSE),
               txt = utils::read.table(file_path, header = TRUE, stringsAsFactors = FALSE),
               rda = { load(file_path); get(ls()[1]) },
               rds = readRDS(file_path),
               xlsx = readxl::read_excel(file_path, guess_max=10000000),
               parquet = arrow::read_parquet(file_path),
               stop("Unsupported file type: ", file_extension))
  
  df <- tibble::as_tibble(df)
  return(df)
}

# Functions to save results 
savepath <- 
  function(savename) { 
    result_folder <- paste0("results/", Sys.Date())
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    
    return(savename)
    
  }

savepath_data <- 
  function(folder, savename) { 
    result_folder <- paste0("data/processed/", folder)
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    
    return(savename)
    
  }
