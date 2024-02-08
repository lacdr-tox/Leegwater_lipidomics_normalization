# Function to calculate the median signal per sample
# Input should have samples of interest as rows and lipids of interest as columns
# Optional: calculate the log2 value of the median after obtaining the median per sample

calculate_median_signal <- function(df, log2transform = F){
  # Set column to rownames
  if("Original_name" %in% colnames(df)){
    data_frame_without_rownames <- T
    df <- column_to_rownames(df, "Original_name") 
    
  } else {
    data_frame_without_rownames <- F
  }
  
  # Calculate median
  median_signal_per_sample <- df %>%
    apply(1, median, na.rm=T) %>%
    as.data.frame() %>%
    rename("median" = 1)
  
  # Log2 transform if needed
  if(log2transform){
    median_signal_per_sample <- median_signal_per_sample %>%
    mutate(log2median = log2(median))
  }
  
  # Set rownames to columns if needed
  if(data_frame_without_rownames){
    median_signal_per_sample <- rownames_to_column(median_signal_per_sample, "Original_name")
  }
  return(median_signal_per_sample)
}
