# Normalize a lipidomics data frame by a metric that is given
# Input should be a data frame with lipids as columns and samples as rows.
# Output is the data frame divided by the specified normalization df
# This should contain a column with some value for each sample.
library(tidyverse)
library(matrixStats) # for rowMedians

source(here::here("code/functions/quotNorm.R"))

.quality_check_normalizationdata <- function(df, norm_values, column_to_use){
  # Many quality checks just in case I break something by accident
  
  if(!"Original_name" %in% colnames(df)){
    stop("Lipidomics data frame needs samples in a column called Original_name")
  } 
  if(!"Original_name" %in% colnames(norm_values)){
    stop("The data frame with normalization values needs samples in a column called Original_name")
  } 
  
  if(nrow(df) != nrow(norm_values)){
    warning("Amount of normalization values or lipids differs. Only overlapping Original_name will be used")
    print("Only in lipidomics df:")
    print(setdiff(df$Original_name, norm_values$Original_name))
    print("Only in normalization df:")
    print(setdiff(norm_values$Original_name, df$Original_name))
    
  }
  
  if(ncol(norm_values) != 2){
    if(ncol(norm_values) == 1){
      stop("Either normalization value or Original_name is missing")
    }
    if(is.na(column_to_use) | !column_to_use %in% colnames(norm_values)){
      stop("Not clear how to normalize. Please specify column in column_to_use")
    } 
  } else{
    column_to_use <- setdiff(colnames(norm_values), "Original_name")
  }
  return(column_to_use)
}

.prepare_df <- function(df){
  if("batch" %in% colnames(df) | "NMC_name" %in% colnames(df)){
    df <- df %>% dplyr::select(-any_of(c("batch", "NMC_name")))
  }
  return(df)
}

.prepare_norm_values <- function(norm_values, column_to_use){
  norm_values <- dplyr::select(norm_values, Original_name, all_of(column_to_use))
  norm_values <- column_to_rownames(norm_values, "Original_name")
  colnames(norm_values) <- c("value_for_normalization")
  norm_values <- rownames_to_column(norm_values, "Original_name")
  return(norm_values)
}

.divide_by_norm_value <- function(df, norm_values){
  # Divide each number by the value_for_normalization in norm_values
  # samples are first matched by Original_name
  norm_df <- inner_join(norm_values, df, by = "Original_name") %>%
    rowwise() %>%
    mutate(across(where(is.double) & !c(Original_name, value_for_normalization),
                  ~.x/value_for_normalization)) %>%
    dplyr::select(-value_for_normalization)
  return(norm_df)
}

.get_normalization_values <- function(df, norm_value){
  norm_values <- data.frame(Original_name = df$Original_name)
  if(norm_value == "median"){
    norm_values$value_for_normalization <- df %>%
      dplyr::select(-Original_name) %>%
      as.matrix() %>%
      rowMedians()
  } else if(norm_value == "sum"){
    norm_values$value_for_normalization <- df %>%
      dplyr::select(-Original_name) %>%
      rowSums()
  } else{
    stop("Unknown normalization. Choose sum, median or provide a data frame")
  }
  return(norm_values)
}
normalize_lipidomics <- function(df, norm_values, column_to_use = NA){
  # If norm_values has more columns than Original_name and one other, 
  # use column_to_use to specify which one
  if(typeof(norm_values) != "character"){
    column_to_use <- .quality_check_normalizationdata(df, norm_values, column_to_use)
    norm_values <- .prepare_norm_values(norm_values, column_to_use)
    df <- .prepare_df(df)
    norm_df <- .divide_by_norm_value(df, norm_values)
  } else if( norm_values != "median" & norm_values != "sum"){
    stop("Norm values not recognized")
  } else{
    df <- .prepare_df(df)
    norm_values <- .get_normalization_values(df, norm_values)
    norm_df <- .divide_by_norm_value(df, norm_values)
  }
  
  print("Done with normalization")
  return(norm_df)
}

########################## Currently used functions ###########################
normalize_lipidomics_quotNorm <- function(df, to_return = "df"){
  # to_return can be the data frame or the dilution factors per sample
  if("Original_name" %in% colnames(df)){
    df <- column_to_rownames(df, "Original_name")
  }
  quotnorm_answer <- quotNorm(df)
  df <- quotnorm_answer$X %>%
    data.frame() %>%
    rownames_to_column("Original_name")
  dilutions <- quotnorm_answer$dilution
  
  if(to_return == "df") {return(df)} else
    if(to_return == "dilutions") {return(dilutions)} else
    {stop("Not sure what to return. Choose df or dilutions")}
}

normalize_lipidomics_quotNorm_structural_lipids <- function(df, structural_lipids){
  # Use quotient normalization based on a subset of lipids - vector called structural lipids - and apply calculated dilution factor to full data frame.
  rownames_in_col1 <- F
  if("Original_name" %in% colnames(df)){
    df <- column_to_rownames(df, "Original_name")
    rownames_in_col1 <- T
  }
  subset_structural_lipids <- dplyr::select(df, all_of(structural_lipids))
  print(str_glue("Normalization will be based on {ncol(subset_structural_lipids)} out of {ncol(df)} lipids."))
  
  dilutions <- normalize_lipidomics_quotNorm(subset_structural_lipids, to_return = "dilutions")
  
  # Divide each sample by its dilution factor
  df_normalized <- df / dilutions
  
  if(rownames_in_col1){
    df_normalized <- rownames_to_column(df_normalized, "Original_name")
  }
  return(df_normalized)
}

log2transform_lipidomics <- function(df){
  if("Original_name" %in% colnames(df)){
    df <- column_to_rownames(df, "Original_name")
  }
  
  df <- df %>% 
    mutate(across(everything(), ~log2(.x))) %>%
    rownames_to_column("Original_name")
  
  return(df)
}

cubeRootTransform_lipidomics <- function(df){
    if("Original_name" %in% colnames(df)){
      df <- column_to_rownames(df, "Original_name")
    }
    df <- df %>% 
      mutate(across(everything(), ~.x^(1/3))) %>%
      rownames_to_column("Original_name")
    
    return(df)
}

auto_scale_lipids <- function(df){
  if("Original_name" %in% colnames(df)){
    df <- column_to_rownames(df, "Original_name")
  }
  df <- scale(df, center = TRUE, scale = TRUE) %>%
    data.frame() %>% rownames_to_column("Original_name")
  return(df)
}

replace_missing_values_by_10percent <- function(df){
  if("Original_name" %in% colnames(df)){
    df <- column_to_rownames(df, "Original_name")
  }
  df <- df %>% 
    mutate(across(everything(),~replace(., . == 0, min(.[.>0], na.rm = TRUE)/10))) %>%
    rownames_to_column("Original_name")
  return(df)
}

save_normalized_df <- function(df, output_dir, LCMS_mode = "all", data_subset = "all", suffix = "unknown_normalization"){
  if(!"Original_name" %in% colnames(df)){
    df <- rownames_to_column(df, "Original_name")
  }
  file_name <- file.path(output_dir, data_subset, str_glue("{data_subset}_{LCMS_mode}_{suffix}.csv"))
  print(str_glue("Saving as {data_subset}_{LCMS_mode}_{suffix}.csv"))
  write_csv(df, file_name)
  return(df)
}

################# Unit tests-like tools in this script #########################
.normalize_lipidomics_unittests <- function(){
  # Normalization
  small_df <- data.frame("Original_name" = c("sample1", "sample2", "sample3", "sample4", "sample5"),
                                    "lipid1" = c(0.03584259, 0.001513433, 0.02783237, 0.048012362, 0.5105104),
                                    "lipid2" = c(0.01014100, 0.001287854, 0.01381689, 0.006427615, 0.1673311),
                                    "lipid3" = c(0.08781706, 0.008241205, 0.40488115, 0.493580496, 2.5723126),
                                    "lipid4" = c(0.10120294, 0.007295504, 0.34559514, 0.291994430, 2.4659894))
  small_df2 <- data.frame("Original_name" = c("sample1", "sample2", "sample3", "sample4", "sample5"),
                         "lipid1" = c(1, 0, 1, 0, 1),
                         "lipid2" = c(1, 2, 3, 4, 5),
                         "lipid3" = c(6, 7, 8, 9, 10),
                         "lipid4" = c(3, 2, 3, 2, 3))
  
  normalize_lipidomics_PQN(small_df)
  normalize_lipidomics_quotNorm(small_df)
  
  # Transformation
  log2transform_lipidomics(small_df)
  cubeRootTransform_lipidomics(small_df)
  
  # Scaling
  pareto_scale_lipids(small_df)
  pareto_scale_lipids(small_df2)
  auto_scale_lipids(small_df)
  auto_scale_lipids(small_df2)
}