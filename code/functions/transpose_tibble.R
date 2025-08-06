# This is a function to transpose a tibble. It works if the first column contains rownames or if rownames are in rownames.
# The output has the transposed rownames in rownames
# Written by Hanneke Leegwater (2024).

transpose_tibble <- function(df, rowname_col = NA){
  # rowname_col = the name of the column that contains rownames
  if(is.na(rowname_col)){
    df <- rownames_to_column(df, "rn")
  } else{
    df <- rename(df, rn = rowname_col)
  }
  df <- df %>%
    pivot_longer(!rn, names_to = "name") %>%
    pivot_wider(id_cols = name, names_from = rn) %>%
    column_to_rownames("name")
  return(df)
}
