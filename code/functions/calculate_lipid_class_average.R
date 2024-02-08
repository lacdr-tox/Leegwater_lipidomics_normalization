# Function to calculate the average of a lipid class per sample.

# Input is a data frame with lipids as columns, and the first column name as sample identifiers
# Also a lipid_info data frame with individual lipids and their class is needed
# For this, column names need to be "class" or "subclass" and "internal_identifiers".
# Class averages will be calculated based on the class colum, or subclass column if subclass=T

calculate_lipid_class_average <- function(df, lipid_info, subclass = F){
  
  if(colnames(df)[1] != "Original_name"){
    warning("The first column name is not Original_name. Are you sure this contains sample names?")
    print("Renaming first column to Original_name...")
    colnames(df)[1] <- "Original_name"
  }
  
  if(subclass == F){
    df_average <- df %>% 
      pivot_longer(!Original_name, names_to = "internal_identifiers", values_to = "value") %>%
      left_join(dplyr::select(lipid_info, internal_identifiers, class), by = "internal_identifiers") %>%
      group_by(Original_name, class) %>%
      summarise(class_mean = mean(value),
              class_sd = sd(value)) %>%
      ungroup()
  } else{
    df_average <- df %>% 
      pivot_longer(!Original_name, names_to = "internal_identifiers", values_to = "value") %>%
      left_join(dplyr::select(lipid_info, internal_identifiers, subclass), by = "internal_identifiers") %>%
      group_by(Original_name, subclass) %>%
      summarise(subclass_mean = mean(value),
                subclass_sd = sd(value)) %>%
      ungroup()
  }
 
  return(df_average)
}
