# These functions scale the time course of a sample series to start from T0 at 0.
# For each lipid, the T0 point is set to 0 this value is subtracted (or divided?) from other points in the time series. 
# Note that the T0=0 needs to be shared for treatment and control to have this code work. If there are two T0's, the code does not work.

.subtract_T0 <- function(lipid_value, control_T0value){
  return(lipid_value - control_T0value)
}


scale_lipid_to_T0is0 <- function(lipid_set, group_name = NA, lipid_column = NA, amount_of_controls =1){
  # Set the sample name column to Original_name and as row name
  if(!"Original_name" %in% colnames(lipid_set)){
    lipid_set <- rownames_to_column(lipid_set, "Original_name") 
  }
  
  # Find the T0 sample
  if("Treatment_group" %in% colnames(lipid_set)){
    control_T0sample <- lipid_set %>%
      filter(str_detect(Treatment_group, "0|zero")) %>% pull(Original_name)
    if(amount_of_controls == 2) control_T0sample <- control_T0sample[1]
    if(length(control_T0sample) != 1){
      stop("Cannot determine control sample in column Treatment_group")
    }
  } else if("Time" %in% colnames(lipid_set)){
    control_T0sample <- lipid_set %>%
      filter(str_detect(Time, "0") | Time == 0) %>% pull(Original_name)
    if(amount_of_controls == 2) control_T0sample <- control_T0sample[1]
    if(length(control_T0sample) != 1){
      stop("Cannot determine control sample in column Time")
    }
  } else{
    control_T0sample <- lipid_set %>%
      filter(str_detect(Original_name, "ctrl-0|ctrl_0|0h")) %>% pull(Original_name)
    if(amount_of_controls == 2) control_T0sample <- control_T0sample[1]
    if(length(control_T0sample) != 1){
      stop("Cannot determine control sample in column Original_name")
    }
  }
  
  # Find the lipid column
  if(is.na(lipid_column)){
    if("lipid" %in% colnames(lipid_set)){
      lipid_column <- "lipid"
    } else if("value" %in% colnames(lipid_set)){
      lipid_column <- "value"
    } else {
      lipid_column <- setdiff(colnames(lipid_set), c("Original_name", "Time", "Treatment_group"))
      if(length(lipid_column) != 1){
        stop(str_glue("Cannot find column with lipid values to correct. Options are: {lipid_column}"))
      }
    }
  }
  
  # Find the T0 value and subtract
  control_T0value <- filter(lipid_set, Original_name == control_T0sample) %>% pull(lipid_column)
  
  scaled_set <- lipid_set %>%
    mutate(across(all_of(lipid_column), ~.x - control_T0value))
  return(scaled_set)
}

scale_multiple_lipids_to_T0is0 <- function(long_df, group_name = "sample", 
            lipid_id_column = "identifiers", lipid_value_column = "value",
                                           amount_of_controls = 1){
  if(!lipid_id_column %in% colnames(long_df)){
    stop(str_glue("Could not find lipid_id_column {lipid_id_column}. Setting to NA"))
  }
  
  result <- long_df %>%
    group_by(across(all_of(c(group_name, lipid_id_column)))) %>%
    group_modify(scale_lipid_to_T0is0, 
                 lipid_column = lipid_value_column, 
                 amount_of_controls = amount_of_controls)
  return(result)
}

################# Unit tests-like tools in this script #########################
.scale_lipid_to_T0is0_unittests <- function(){
  # Normalization
  small_df <- data.frame("Original_name" = c("168_MDA-MB-231_ctrl-0", "169_MDA-MB-231_ctrl-6", "170_MDA-MB-231_ctrl-24", 
                                                                      "172_MDA-MB-231_2DG-6", "173_MDA-MB-231_2DG-24"),
                         "Treatment_group" = c("ctrl_0", "ctrl_6", "ctrl_24", "2DG_6", "2DG_24"),
                         "lipid1" = c(0.03189661, 0.16462307, 0.278610327, 0.4578655, 0.15865607))
  
  .subtract_T0(0.16462307, 0.03189661)
  
  scale_lipid_to_T0is0(small_df)
  
  bigger_df <- data.frame("Original_name" = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"),
                          "Treatment_group" = c("ctrl_0", "ctrl_6", "ctrl_24", "2DG_6", "2DG_24",
                                                "ctrl_0", "ctrl_6", "ctrl_24", "2DG_6", "2DG_24"),
                          "sample" = c("A", "A", "A", "A", "A", "B", "B", "B", "B", "B"),
                          "lipid" = c(1, 2, 3, 4, 5, 1, 2, 2, 3, 3))
  bigger_df %>%
    group_by(sample) %>%
    group_modify(scale_lipid_to_T0is0)
  
  even_bigger_df <- data.frame("Original_name" = rep(c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"), 2),
                          "Treatment_group" = rep(c("ctrl_0", "ctrl_6", "ctrl_24", "2DG_6", "2DG_24"), 4),
                          "sample" = rep(c(rep("A", 5), rep("B", 5)),2),
                          "identifiers" = c(rep("lipid1", 10), rep("lipid2", 10)),
                          "value" = c(1,2,3,4,5, 1,2,2,3,3, 8,6,6,9,9, 3,7,3,2,8))
  
  scale_multiple_lipids_to_T0is0(even_bigger_df, group_name = "sample", 
                                 lipid_id_column = "identifiers", lipid_value_column = "value", 
                                 amount_of_controls = 1)
  
}