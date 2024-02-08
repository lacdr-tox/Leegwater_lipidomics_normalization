# Functions to calculate intraclass correlation
# The input data frame has an original_name, a sample, a 
# Sometimes all headers are in row 1, but sometimes they are split in row 1 and 3
# Input:
# - The lipidomics data frame to use as starting material. Ideally, this is a
# data frame openend with the read_lipidomics.R script
# - The lipid name for which the ICC will be calculated
# Output:
# - The ICC value

library(tidyverse)
library(irr)
# 1. subset lipidomics df to keep only lipid
.subset_lipidomics <- function(df, lipid, metadata = NA, log10Transform = FALSE, 
                               replicate_column_name = "biological_replicate", group_column_name = "group"){
  # Metadata should have unique sample names in Original_name
  if(typeof(metadata) == "list"){
    metadata <- dplyr::select(metadata, Original_name, all_of(replicate_column_name), all_of(group_column_name))
    df <- df %>%
      inner_join(metadata, by = "Original_name")
  }
  
  df_subset <- df %>% 
    rename(replicate = replicate_column_name,
           group = group_column_name) %>%
    dplyr::select(group, replicate, all_of(lipid)) %>%
    mutate(replicate = str_glue("replicate_{replicate}")) %>%
    pivot_wider(id_cols = group, names_from = replicate, values_from = all_of(lipid)) %>%
    column_to_rownames("group")
  
  if(log10Transform){
    df_subset <- mutate(df_subset, across(everything(), ~log10(.x)))
  }
  return(df_subset)
}

calculate_one_ICC <- function(df,lipid, metadata, log10Transform){
  result <- icc(df, model = "twoway", type = "agreement", unit = "single")
  if(is.na(result$value)){
      result$value <- 0
      warning(str_glue("ICC cannot be calculated for {lipid}. Set to 0.\n"))
    }
  return(result$value)
}

.calculate_ICC <- function(df,lipid = "not specified"){
  return(calculate_one_ICC(df, lipid))
}
  
.subset_and_calculate_ICC <- function(lipid, df, metadata = NA, log10Transform = FALSE,
                                      replicate_column_name = "biological_replicate", group_column_name = "group"){
  df <- .subset_lipidomics(df, lipid, metadata, log10Transform,replicate_column_name, group_column_name)
  result <- calculate_one_ICC(df, lipid)
  return(result)
}

calculate_all_ICC <- function(df, metadata = NA, log10Transform = FALSE, 
                              replicate_column_name = "biological_replicate", group_column_name = "group"){
  # Lipids need to be in columns, and a sample column can be there, but it should
  # be called Original_name. 
  all_lipids <- data.frame("lipid" = setdiff(colnames(df), "Original_name"),
                           "ICC" = NaN) %>%
    column_to_rownames("lipid")
  for(lipid in rownames(all_lipids)){
    nr <- .subset_and_calculate_ICC(lipid,
      df = df, 
      metadata = metadata, 
      log10Transform = log10Transform,
      replicate_column_name = replicate_column_name,
      group_column_name = group_column_name)
    all_lipids[lipid, "ICC"] <- nr
  }
  
  return(all_lipids)
}

calculate_ICC_for_each_dataset <- function(df_list, metadata, 
                                           replicate_column_name = "biological_replicate",
                                           group_column_name = "group"){
  # Returns a data frame with lipids as rows and ICC values as columns. 
  # Input is a list of data frames for which ICCs need to be calculated. 
  # This code assumes (does not check) that all data frames have the same lipids as columns and samples as rows, in the same order
  all_ICC_vals <- lapply(df_list, 
                         function(x) calculate_all_ICC(x, 
                                                       metadata,
                                                       replicate_column_name = replicate_column_name,
                                                       group_column_name = group_column_name
                         ))
  names(all_ICC_vals)[1:3]
  
  all_ICC_vals_df <- bind_cols(all_ICC_vals)
  colnames(all_ICC_vals_df) <- names(all_ICC_vals)
  return(all_ICC_vals_df)
}

################# Unit tests-like tools in this script #########################
.calculate_ICC_omics_unittests <- function(){
  # Normalization
  small_df <- data.frame("Original_name" = c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6"),
                         "lipid1" = c(1.0, 1.1, 1.3, 4.0, 4.2, 3.9),
                         "lipid2" = c(1.0, 1.1, 1.3, 1.0, 1.1, 1.3),
                         "lipid3" = c(1.0, 1.1, 1.3, 1.2, 1.3, 1.5),
                         "lipid4" = c(0.3, 1.8, 2.5, 9.2, 1.1, 2.0))
  small_metadata <- data.frame("Original_name" = small_df$Original_name,
                      "sample" = c("group1", "group1", "group1", "group2", "group2", "group2"),
                      "group" = c("group1", "group1", "group1", "group2", "group2", "group2"),
                      "replicate" = c("1", "2", "3", "1", "2", "3"),
                      "biological_replicate" = c("1", "2", "3", "1", "2", "3"))
  
  # Check if one lipid can be extracted and put in a column per group
  .subset_lipidomics(small_df, "lipid1", small_metadata, 
                                 replicate_column_name = "biological_replicate", group_column_name = "group")
  
  # Test four options of calculating a ICC for one lipid
  .subset_lipidomics(small_df, "lipid1", small_metadata) %>%
    calculate_one_ICC() # should be very high
  .subset_lipidomics(small_df, "lipid2", small_metadata) %>%
    calculate_one_ICC() # should be very low
  .subset_lipidomics(small_df, "lipid3", small_metadata) %>%
    calculate_one_ICC() # should be medium
  .subset_lipidomics(small_df, "lipid4", small_metadata) %>%
    calculate_one_ICC() # should be random-ish
  
  # Calculate all ICC at once
  calculate_all_ICC(small_df, small_metadata)
}