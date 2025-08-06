# To duplicate the T=0 control samples (and metadata if needed).
# This can be useful for example for plots where you want to plot two lines
# One for control and one for treatment over time.
library(here)

if(file.exists(here("code/functions/transpose_tibble.R"))) source(here("code/functions/transpose_tibble.R")) else
  stop("file transpose_tibble.R not found. Source manually.")

duplicate_control_T0_samples <- function(df, samples_as_rows = T){
  # If samples are rows, check for Original_name. If samples are column, check
  # for a lipid column or set lipids to row names yourself
  if(samples_as_rows & !"Original_name" %in% colnames(df)){
    df <- rownames_to_column(df, "Original_name")
  }
  if(!samples_as_rows){
    if("lipid" %in% colnames(df)){ 
      df <- column_to_rownames(df, "lipid")
    }
    df <- transpose_tibble(df) %>%
      rownames_to_column("Original_name")
  }
  
  extra_data <- df %>%
    filter(str_detect(Original_name, "ctrl-0")) %>%
    mutate(Original_name = str_replace(Original_name, "ctrl-0", "2DG-0"))
  
  # Merge data and return
  df <- bind_rows(df, extra_data)
  
  if(!samples_as_rows){
    df <- df %>%
      column_to_rownames("Original_name") %>%
      transpose_tibble()
  }
  return(df)
}

duplicate_control_T0_samples_long <- function(long_df, id_col = "Original_name"){
  extra_data <- long_df %>%
    dplyr::rename(Original_name = all_of(id_col)) %>%
    filter(str_detect(Original_name, "ctrl-0|_0h|Ctrl 0h")) %>%
    mutate(Original_name = str_replace(Original_name, "ctrl-0|_0h", "2DG-0")) %>%
    mutate(Original_name = str_replace(Original_name, "Ctrl 0h", "2DG 0h")) %>%
    dplyr::rename(!!id_col := "Original_name")
  if("Treatment" %in% colnames(extra_data)){
    extra_data <- extra_data %>%
      mutate(Treatment = str_replace(Treatment, "Ctrl|ctrl", "2DG"))
  }
  # Merge data and return
  long_df <- bind_rows(long_df, extra_data)
  return(long_df)
}

duplicate_control_T0_metadata <- function(metadata_df){
  extra_metadata <- metadata_df %>%
    filter(str_detect(Original_name, "ctrl-0|_0h")) %>%
    mutate(Original_name = str_replace(Original_name, "ctrl-0|_0h", "2DG-0"),
           if("group" %in% colnames(metadata_df)){group = str_replace(group, "ctrl-0|_0h", "2DG-0")},
           if("Name_in_original_data" %in% colnames(metadata_df)){Name_in_original_data = str_replace(Name_in_original_data, "ctrl-0|_0h", "2DG-0")},
           if("Treatment_group" %in% colnames(metadata_df)){Treatment_group = "2DG_0"},
           Treatment = "2DG") 
  metadata_df <- bind_rows(metadata_df, extra_metadata)
  return(metadata_df)
}

################# Unit tests-like tools in this script #########################
.duplicate_control_T0_samples_unittests <- function(){
  # Normalization
  small_df <- data.frame("Original_name" = c("168_MDA-MB-231_ctrl-0", "169_MDA-MB-231_ctrl-6", "170_MDA-MB-231_ctrl-24", "172_MDA-MB-231_2DG-6", "173_MDA-MB-231_2DG-24"),
    "lipid1" = c(0.03189661, 0.16462307, 0.278610327, 0.4578655, 0.15865607),
    "lipid2" = c(-0.06938994, -0.03057414, 0.001784968, 0.3086332, -0.07931443),
    "lipid3" = c(0.03391186, 0.14951118, 0.040897301, 0.6095089, -0.01200846),
    "lipid4" = c(-0.28277344, -0.32142634, -0.300435926, -0.2057648, -0.30387240))
  small_df2 <- small_df %>% column_to_rownames("Original_name") %>% transpose_tibble()
  
  # This should work
  duplicate_control_T0_samples(small_df)
  duplicate_control_T0_samples(small_df2, samples_as_rows = F)
  # and this should fail/be odd/break
  duplicate_control_T0_samples(small_df, samples_as_rows = F)
  duplicate_control_T0_samples(small_df2)
  
  # For the long format
  small_df_long <- small_df %>% 
    pivot_longer(!Original_name, names_to = "ID", values_to = "val")
  small_df_long2 <- dplyr::rename(small_df_long, other_name = Original_name)
  duplicate_control_T0_samples_long(small_df_long) %>% dim()
  duplicate_control_T0_samples_long(small_df_long2, id_col = "other_name") %>% dim()
  
  # Note that it does not generate a different sampleID or subject!
  small_metadata <- data.frame("Original_name" = c("168_MDA-MB-231_ctrl-0", "169_MDA-MB-231_ctrl-6", "170_MDA-MB-231_ctrl-24", "172_MDA-MB-231_2DG-6"),
    "sampleID" = c("168", "169", "170", "172"),
    "group" = c("MDA-MB-231_ctrl-0", "MDA-MB-231_ctrl-6", "MDA-MB-231_ctrl-24", "MDA-MB-231_2DG-6"),
    "Name_in_original_data" = c("168_MDA-MB-231_ctrl-0", "169_MDA-MB-232_ctrl-6", "170_MDA-MB-233_ctrl-24", "172_MDA-MB-235_2DG-6"),
    "sampleName" = c("MDA-MB-231_0_1", "MDA-MB-231_6_1", "MDA-MB-231_24_1", "MDA-MB-231_6_1"),
    "cell_line" = c("MDA-MB-231", "MDA-MB-231", "MDA-MB-231", "MDA-MB-231"),
    "Treatment_group" = c("ctrl_0", "ctrl_6", "ctrl_24", "2DG_6"),
    "Treatment" = c("ctrl", "ctrl", "ctrl", "2DG"),
    "Time" = c(0, 6, 24, 6),
    "Subject" = c("S2", "S2", "S2", "S3"),
    "biological_replicate" = c("1", "1", "1", "1"),
    "batch_MSMS" = c(1, 1, 1, 1),
    "Study" = c("treatment_2DG", "treatment_2DG", "treatment_2DG", "treatment_2DG"))
  
  duplicate_control_T0_metadata(small_metadata)
  
}