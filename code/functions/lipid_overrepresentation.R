# Calculate the overrepresentation of a lipid property in a list of lipids
# Input is the number of features of interest and the total lipid list length
# Output is an overrepresentation number and a statistic

# Reference blogpost https://statsandr.com/blog/fisher-s-exact-test-in-r-independence-test-for-a-small-sample/ accessed on 05-04-2024
# And post on the difference between an odds ratio and the result of the fisher exact test https://stats.stackexchange.com/questions/54530/why-do-odds-ratios-from-formula-and-rs-fisher-test-differ-which-one-should-one
# Idea for the code and choice of statistics is Panther, the default for gene ontology enrichment in gene lists (https://pantherdb.org/webservices/go/overrep.jsp)
# Help for getting nice filters from https://stackoverflow.com/questions/70355426/referring-to-columns-and-variables-with-the-same-name-in-dplyr-filter

count_pattern_in_limma_result <- function(limma_result_table, pattern, direction, adj.P.Val = 0.1, id_col = "internal_identifiers", exact_match = F){
  # Remove empty/rows with NA
  df <- limma_result_table %>%
    filter(!is.na(.data[[id_col]]))
  # Find the pattern in the column of interest.
  if(exact_match){
    df <- df %>%
      mutate(contains_pattern = .data[[id_col]] == pattern)
  } else{
    df <- df %>%
      mutate(contains_pattern = str_detect(.data[[id_col]], pattern))
  }
  df <- df %>% mutate(
           in_filter = (adj.P.Val < .env$adj.P.Val & direction == .env$direction))
  
  df_filter <- df %>%
    filter(adj.P.Val < .env$adj.P.Val, direction == .env$direction)
  
  # Make a table
  nr_of_hits <- sum(df_filter$contains_pattern)
  features_in_total_list <- sum(df$contains_pattern)
  hit_list_length <- nrow(df_filter)
  total_list_length <- nrow(df)
  
  contingency_table <- data.frame(
    "yes" = c(nr_of_hits, features_in_total_list),
    "no" = c((hit_list_length - nr_of_hits), (total_list_length - features_in_total_list)),
    row.names = c("hits", "total")
  )
  
  return(contingency_table)
}

calc_overrepresentation <- function(contingency_table){
  # Perform fisher test and get odds ratio and p value
  result <- fisher.test(contingency_table)
  p_value_fisher <- result$p.value
  conditional_odds_ratio <- unname(result$estimate)
  return(list(p_value_fisher = p_value_fisher, 
              conditional_odds_ratio = conditional_odds_ratio))
}

lipid_overrepresentation <- function(limma_result_table, pattern, direction, adj.P.Val = 0.1, id_col = "internal_identifiers", exact_match = F){
  # Create an overrepresentation data frame with relevant information
  
  contingency_table <- count_pattern_in_limma_result(limma_result_table, pattern, direction, adj.P.Val, id_col, exact_match)
  df <- calc_overrepresentation(contingency_table) %>%
    data.frame()
  rownames(df) <- c(pattern)
  
  df$nr_of_hits <- contingency_table["hits", "yes"]
  df$hit_list_length <- sum(contingency_table["hits",])
  df$features_in_total_list <- contingency_table["total", "yes"]
  df$total_list_length <- sum(contingency_table["total",])
  
  df <- relocate(df, conditional_odds_ratio)
  
  return(df)
  
}

calculate_column_overrepresentations <- function(limma_result_table, id_col = "internal_identifiers", direction, adj.P.Val = 0.1, exact_match = F){
  print(str_glue("***Enriched {id_col}***"))
  patterns <- levels(limma_result_table[[id_col]])
  df <- map_dfr(patterns, 
                function(x) lipid_overrepresentation(limma_result_table, x, direction, adj.P.Val = adj.P.Val, id_col = id_col, exact_match = exact_match)) 
  
  df$padj <- p.adjust(df$p_value_fisher, method = "bonferroni") 
  df <- df %>% arrange(p_value_fisher) %>%
    arrange(padj) %>% relocate(padj, .after = conditional_odds_ratio)
  
  return(df)
}

lipid_overrepresentation_from_hits <- function(nr_of_hits, hit_list_length, 
                                     features_in_total_list, total_list_length,
                                     show_mosaic = F, return_pval = T){
  # Check if there are any hits
  if(nr_of_hits == 0){
    if(return_pval){
      return(list("conditional_odds_ratio" = 0, 
                  "pval_fisher" = NA))
    } else{
      return(0)
    }
  }
  # Make a table
  contingency_table <- data.frame(
    "yes_feature" = c(nr_of_hits, features_in_total_list),
    "no_feature" = c((hit_list_length - nr_of_hits), (total_list_length - features_in_total_list)),
    row.names = c("hits", "total")
  )
  res <- calc_overrepresentation(contingency_table)
  
  if(show_mosaic){
    p <- mosaicplot(contingency_table,
               main = "Mosaic plot",
               color = TRUE
    )
    print(p)
  }
  
  if(return_pval){
    return(res)
  } else{
    return(res$conditional_odds_ratio)
  }
}


################# Unit tests-like tools in this script #########################
.lipid_overrepresentation_unittests <- function(){
  count_pattern_in_limma_result(
    limma_result_table = tibble(
      "internal_identifiers" = c("lipid_A1", "lipid_A2", "lipid_AA2"),
      "class" = c("A", "A", "AA"),
      "adj.P.Val" = c(0.01, 0.2, 0.01)
    ),
    pattern = "A",
    direction = "up",
    id_col = "class",
    exact_match = T  
  )
  lipid_overrepresentation(
    limma_result_table = tibble(
        "internal_identifiers" = c("lipid_A1", "lipid_A2", "lipid_AA2"),
        "class" = c("A", "A", "AA"),
        "adj.P.Val" = c(0.01, 0.2, 0.01)
      ),
      pattern = "A",
      direction = "up",
      id_col = "class",
      exact_match = T  
  )
  
  # #library(testthat)
  # #test_that("calculation works", {
  #   expect_equal(lipid_overrepresentation(5, 10, 50, 100, return_pval = F), 1) # expect 1
  #   expect_lt(lipid_overrepresentation(5, 10, 80, 100, return_pval = F), 1) # expect less than 1
  #   expect_gt(lipid_overrepresentation(8, 10, 50, 100, return_pval = F), 1) # expect more than 1
  #   expect_equal(lipid_overrepresentation(20, 40, 100, 400, return_pval = F), 2) # expect two
  #   expect_equal(lipid_overrepresentation(3, 6, 5, 6, return_pval = F), 0.6) # expect 0.6
  #   expect_equal(lipid_overrepresentation(0, 6, 5, 6, return_pval = F), 0) # expect 0
  # })
  
  # More real examples
  lipid_overrepresentation(8, 40, 10, 200,
                           show_mosaic = T)
  # and one in my dataset
  nr_of_hits <- 19
  hit_list_length <- 19+17
  features_in_total_list <- 67
  total_list_length <- 892
  
  lipid_overrepresentation(nr_of_hits, hit_list_length, features_in_total_list, total_list_length,
                           show_mosaic = F)
  lipid_overrepresentation(nr_of_hits, hit_list_length, features_in_total_list, total_list_length,
                           show_mosaic = T)
  
  contingency_table2 <- data.frame(
    "yes_feature" = c(8, 40),
    "no_feature" = c(10, 100),
    row.names = c("hits", "total")
  )
  calc_overrepresentation(contingency_table2)
}