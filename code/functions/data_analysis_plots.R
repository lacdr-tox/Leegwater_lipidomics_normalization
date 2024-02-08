# Plot df info
# This set of functions can plot a variety of things for data frames.

# A density plot per sample, showing the distribution of lipids for each biological replicate.
plot_density_per_individual_sample_group <- function(df, data_subset, metadata){
  # metadata should contain the groups of interest
  # example: plot_density_per_individual_sample_group(data_cell, "cell_lines)
  # output: a ggplot object with facet_wrap per sample group
  if(data_subset == "cell_lines"){
    p <- .density_per_cell_line(df, metadata)
  } else if(data_subset == "treatment_2DG") {
    p <- .density_per_treatment(df, metadata)
  } else{ stop("data_subset type not recognized.")}
  return(p)
}

.density_per_cell_line <- function(df, metadata){
  p <- df %>%
    left_join(select(metadata, Original_name, cell_line, biological_replicate), by = "Original_name") %>%
    pivot_longer(!c(Original_name, cell_line, biological_replicate), names_to = "lipid") %>%
    ggplot(aes(x = value, group = Original_name)) +
    geom_density(aes(color = biological_replicate)) +
    theme_classic() +
    facet_wrap(~cell_line)
  return(p)
}

.density_per_treatment <- function(df, metadata){
  p <- df %>%
    left_join(select(metadata, Original_name, cell_line, Treatment, Time, biological_replicate), by = "Original_name") %>%
    unite(Timepoint, Treatment, Time, remove = FALSE) %>%
    pivot_longer(!c(Original_name, cell_line, Treatment, Time, Timepoint, biological_replicate), names_to = "lipid") %>%
    ggplot(aes(x = value, group = Original_name)) +
    geom_density(aes(color = biological_replicate)) +
    theme_classic() +
    facet_wrap(~cell_line + Timepoint)
  return(p)
}

#Define function to plot a boxplot based on a data frame. It shows samples on the X axis and the signal distribution per sample on the Y axis.
make_sample_boxplot <- function(df, title){
  p <- df %>%
    pivot_longer(!Original_name) %>%
    ggplot(aes(x = Original_name, y = value)) + 
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank() 
    ) +
    labs(x = "Sample", y = "signal", title = title)
  return(p)
}

# A function to have one overlay of a density distribution
make_sample_densityplot <- function(df, title){
  p <- df %>%
    pivot_longer(!Original_name, names_to = "lipid") %>%
    ggplot(aes(x = value, group = Original_name)) +
    geom_density(linetype="dashed", color = "#808080") +
    theme_classic() +
    ggtitle(title)
  return(p)
}

# A function to compare the old and new distribution overlay of all samples
compare_densities <- function(df, title, original_df = lipidomics_cell_lines_log2){
  p_new <- make_sample_densityplot(df, title = title)
  p_old <- make_sample_densityplot(original_df, title = "Log2 transformed only")
  return(p_old + p_new)
}

plot_PCA <- function(df, metadata, LCMS_mode, data_subset){
  # Plot a PCA using pca in mixOmics without scaling
  
  # Prepare the data frames
  if("Original_name" %in% colnames(df)){
    df <- column_to_rownames(df, "Original_name")
  }
  
  # Calculate PCA
  pca_obj <- mixOmics::pca(df, ncomp = ifelse(ncol(df) >=8, 8, ncol(df)))
  plot_list <- list()
  plot_list$scree_plot <- plot_scree_plot(pca_obj, LCMS_mode, data_subset, ncomp = ifelse(ncol(df) >=8, 8, ncol(df)))
  plot_list$group <- plot_PCA_df(pca_obj, metadata, "group", LCMS_mode, data_subset, add_line = T)
  plot_list$PCA_cell_line <- plot_PCA_df(pca_obj, metadata, "cell_line", LCMS_mode, data_subset)
  plot_list$PCA_biological_replicate <- plot_PCA_df(pca_obj, metadata, "biological_replicate", LCMS_mode, data_subset, legend = T)
  
  if(data_subset == "cell_lines"){
  plot_list$PCA_batch_MSMS <- plot_PCA_df(pca_obj, metadata, "batch_MSMS", LCMS_mode, data_subset)
  }
  if(data_subset == "treatment_2DG"){
  plot_list$PCA_treatment <- plot_PCA_df(pca_obj, metadata, "Treatment", LCMS_mode, data_subset, legend = T)
  plot_list$PCA_timepoint <- plot_PCA_df(pca_obj, metadata, "Time", LCMS_mode, data_subset, legend = T)
  }
  
  return(plot_list)
}
get_PCA_df <- function(pca_mixomics_object){
  # To extract the principal components from a PCA mixOmics object
  df <- data.frame(pca_mixomics_object$variates$X) %>%
    rownames_to_column("Original_name")
  return(df)
}

plot_scree_plot <- function(pca_obj, LCMS_mode, data_subset, ncomp = 8){
  # Extract explained variance from a mixOmics PCA object and plot it
  var_explained_df<- data.frame("var_explained" = pca_obj$prop_expl_var$X) %>%
    rownames_to_column("PC") %>%
    mutate(var_explained = round(var_explained, 2)*100)
  
  p <-var_explained_df %>%
    ggplot(aes(x=PC,y=var_explained, group=1))+
    geom_point(size=4)+
    geom_line()+
    theme_classic() +
    scale_y_continuous(limits = c(0, NA)) +
    labs(title=str_glue("Scree plot: PCA for {data_subset} {LCMS_mode}"), 
         y = "variance explained (%)", x = NULL)
  return(p)
}

plot_PCA_df <- function(df, metadata, color_metadata, LCMS_mode = "all", data_subset = "all", 
                        legend = FALSE, add_line = F, comp = c(1,2)){
  
  # Set x and y axis
  PCx <- str_glue("PC{comp[1]}")
  PCy <- str_glue("PC{comp[2]}")
  
  # Turn mixomics pca object into data frame for ggplot
  if(class(df) == "pca"){
    PCx_var <- round(df$prop_expl_var$X[[PCx]], 2) *100
    x_label <- str_glue("{PCx} ({PCx_var}% expl. var)")
    PCy_var <- round(df$prop_expl_var$X[[PCy]],2) *100
    y_label <- str_glue("{PCy} ({PCy_var}% expl. var)")
    df <- get_PCA_df(df)
  } else{
    x_label <- PCx
    y_label <- PCy
  }
  
  df <- left_join(df, metadata, by = "Original_name")
  p <- ggplot(df, aes(x = .data[[PCx]], y = .data[[PCy]], fill = .data[[color_metadata]], color = .data[[color_metadata]])) +
    geom_point() +
    theme_classic() +
    theme(legend.position = ifelse(legend, "bottom", "none")) +
    ggtitle(str_glue("PCA for {data_subset} {LCMS_mode} colored by {color_metadata}")) +
    labs(x = x_label, y = y_label) 
  if(add_line){
    p <- p + geom_line(aes(group = group))
  }
  return(p)
}

print_report_plots <- function(df, LCMS_mode, data_subset, metadata, save_plots = F){
  # Prepare metadata for plotting
  if(!"Original_name" %in% colnames(metadata)){
    metadata <- rownames_to_column(metadata, "Original_name")
  }
  if(!"Original_name" %in% colnames(df)){
    df <- rownames_to_column(df, "Original_name")
  }
  if(nrow(metadata) != nrow(df)){
    metadata <- filter(metadata, Original_name %in% df$Original_name)
  }
  
  # Plot a variety of plots to test normalization quality
  plot_list <- list()
  plot_list$sample_densityplot <- make_sample_densityplot(df, title = str_glue("Density plot for {data_subset} {LCMS_mode}"))
  plot_list$sample_boxplot <- make_sample_boxplot(df, title = str_glue("Sample boxplot for {data_subset} {LCMS_mode}"))
  plot_list$density_per_individual_sample_group <- plot_density_per_individual_sample_group(df, data_subset, metadata)
  plot_list <- append(plot_list, plot_PCA(df, metadata, LCMS_mode, data_subset))
  
  # Either return list of plot objects or the original df
  if(save_plots){
    return(plot_list)
  } else{
    for(p in plot_list){
      plot(p)
    }
    return(df)
  }
}

save_report_plots <- function(df, output_dir, LCMS_mode, data_subset, metadata, suffix, one_file = F){
  # Save a variety of plots
  plot_list <- print_report_plots(df, LCMS_mode, data_subset, metadata, save_plots = T)
  
  if(one_file){
  # Save in one file
    datestamp <- str_remove_all(Sys.Date(), "-")
    file_name <- file.path(output_dir, str_glue("{datestamp}_plots_{data_subset}_{LCMS_mode}_{suffix}.pdf"))
    pdf(file_name)
    for(plot_name in names(plot_list)){
      p <- plot_list[[plot_name]]
      plot(p)
    }
    dev.off()
  } else{
  # Save in multiple files
    output_dir <- file.path(output_dir, str_glue("{datestamp}_plots_{data_subset}_{LCMS_mode}"))
    for(plot_name in names(plot_list)){
      p <- plot_list[[plot_name]]
      file_name <- file.path(output_dir, str_glue("{datestamp}_plot_{data_subset}_{LCMS_mode}_{plot_name}.png"))
      ggsave(file_name, p)
    }
  }
  
  return(df)
}
