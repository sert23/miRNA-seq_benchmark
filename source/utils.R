# change column names
fix_names <- function(input_cols){
  return(gsub("\\_R1\\|grp1", "", input_cols))
}

define_groups <- function(input_df) {
  groups <- list()
  
  for (x in c(20, 40, 60, 80)) {
    current_group <- as.character(x)
    columns <- grep(paste0("^", current_group, "(?!.*_2$)"), colnames(input_df), perl = TRUE)
    members <- colnames(input_df)[columns]
    groups[[paste0(current_group, "S")]] <- members
  }
  
  columns <- grep("^[1-5]L(?!.*_.*_)", colnames(input_df), perl = TRUE)
  groups[["0S"]] <- colnames(input_df)[columns]
  
  columns <- grep("^[1-5]S(?!.*_2$)", colnames(input_df), perl = TRUE)
  groups[["100S"]] <- colnames(input_df)[columns]
  
  return(groups)
}

get_group <- function(entry, group_list) {
  entry <- trimws(entry)
  
  # Special case
  if (entry == "1L_2L_3L_4L_5L") return("0S")
  
  # First: replace '_2' with '_1'
  entry_alt1 <- sub("_2$", "_1", entry)
  for (name in names(group_list)) {
    if (entry_alt1 %in% group_list[[name]]) return(name)
  }
  
  # Second: remove '_2' altogether
  entry_alt2 <- sub("_2$", "", entry)
  if (entry_alt2 != entry_alt1) {
    for (name in names(group_list)) {
      if (entry_alt2 %in% group_list[[name]]) return(name)
    }
  }
  
  return(NA_character_)
}

# Function to check and remove directories with less than 15 comparisons
cleanup_directories <- function(method_dir) {
  # List all parameter combination directories
  combo_dirs <- list.dirs(method_dir, recursive = FALSE)
  
  for (combo_dir in combo_dirs) {
    # Count the number of subdirectories (pairwise comparisons)
    subdirs <- list.dirs(combo_dir, recursive = FALSE)
    
    # If there are fewer than 15 comparisons, remove the directory
    if (length(subdirs) < 15) {
      print(paste("Removing incomplete parameter combination:", combo_dir))
      unlink(combo_dir, recursive = TRUE)
    }
  }
}

check_monotonicity <- function(x) {
  if (all(diff(x) > 0)) {
    return("increasing")
  } else if (all(diff(x) < 0)) {
    return("decreasing")
  } else {
    return("neither")
  }
}

# Function to get method names and aggregate .tsv files for a specific comparison
aggregate_tsv_files_old <- function(base_path, comparison_name) {
  # Point 1: Get the names of all methods (first-level subdirectories)
  method_names <- list.dirs(path = base_path, full.names = FALSE, recursive = FALSE)
  # print(method_names)
  # Initialize a list to store data from each method
  aggregated_data <- list()
  
  # Iterate over methods
  for (method in method_names) {
    method_path <- file.path(base_path, method)
    
    # Get level 2 subdirectories (parameter combinations)
    parameter_combinations <- list.dirs(path = method_path, full.names = TRUE, recursive = FALSE)
    # print(parameter_combinations)
    for (param_path in parameter_combinations) {
      # Check if the specified comparison exists (level 3 subdirectory)
      comparison_path <- file.path(param_path, comparison_name)
      # print(comparison_path)
      if (dir.exists(comparison_path)) {
        # Get the .tsv files in the comparison directory
        tsv_files <- list.files(path = comparison_path, pattern = "\\.tsv$", full.names = TRUE)
        
        # Read each .tsv file and add it to the aggregated data
        for (tsv_file in tsv_files) {
          tsv_data <- read.delim(tsv_file, header = TRUE, stringsAsFactors = FALSE)
          tsv_data$Method <- method  # Add the method name as a column
          tsv_data$Parameter <- basename(param_path)  # Add the parameter combination as a column
          tsv_data$Comparison <- comparison_name  # Add the comparison name as a column
          
          # Append to the aggregated data list
          aggregated_data[[length(aggregated_data) + 1]] <- tsv_data
        }
      }
    }
  }
  
  # Combine all dataframes in the list into a single dataframe
  final_data <- do.call(rbind, aggregated_data)
  final_data$FC_direction <- ifelse(final_data$log2FC > 0, "up", "down")
  
  
  return(final_data)
}

aggregate_tsv_files_old <- function(base_path, comparison_name) {
  # Point 1: Get the names of all methods (first-level subdirectories)
  method_names <- list.dirs(path = base_path, full.names = FALSE, recursive = FALSE)
  method_names <- method_names[!grepl("^_", method_names)]  # ← ignore method dirs that start with "_"
  # Initialize a list to store data from each method
  aggregated_data <- list()
  
  # Iterate over methods
  for (method in method_names) {
    method_path <- file.path(base_path, method)
    
    # Get level 2 subdirectories (parameter combinations)
    parameter_combinations <- list.dirs(path = method_path, full.names = TRUE, recursive = FALSE)
    for (param_path in parameter_combinations) {
      # Check if the specified comparison exists (level 3 subdirectory)
      comparison_path <- file.path(param_path, comparison_name)
      if (dir.exists(comparison_path)) {
        # Get the .tsv files in the comparison directory
        tsv_files <- list.files(path = comparison_path, pattern = "\\.tsv$", full.names = TRUE)
        
        # Read each .tsv file and add it to the aggregated data
        for (tsv_file in tsv_files) {
          tsv_data <- read.delim(tsv_file, header = TRUE, stringsAsFactors = FALSE)
          
          # Rename columns if needed
          if ("method"  %in% colnames(tsv_data)) {
            names(tsv_data)[names(tsv_data) == "method"]  <- "DE_method"
          }
          if ("logFC"   %in% colnames(tsv_data)) {
            names(tsv_data)[names(tsv_data) == "logFC"]   <- "log2FC"
          }
          
          # Ensure an 'average' column exists
          if (!"average"    %in% colnames(tsv_data)) {
            tsv_data$average    <- 0
          }
          
          # Ensure a 'parameters' column exists
          if (!"parameters" %in% colnames(tsv_data)) {
            tsv_data$parameters <- "default"
          }
          
          tsv_data$Method     <- method                   # Add the method name
          tsv_data$Parameter  <- basename(param_path)     # Add the parameter combination
          tsv_data$Comparison <- comparison_name          # Add the comparison name
          
          # Append to the aggregated data list
          aggregated_data[[length(aggregated_data) + 1]] <- tsv_data
        }
      }
    }
  }
  
  # Combine all dataframes in the list into a single dataframe
  final_data <- do.call(rbind, aggregated_data)
  final_data$FC_direction <- ifelse(final_data$log2FC > 0, "up", "down")
  
  return(final_data)
}

aggregate_tsv_files <- function(base_path, comparison_name) {
  # Point 1: Get the names of all methods (first-level subdirectories)
  method_names <- list.dirs(path = base_path, full.names = FALSE, recursive = FALSE)
  method_names <- method_names[!grepl("^_", method_names)]  # ← ignore method dirs that start with "_"
  
  # Initialize a list to store data from each method
  aggregated_data <- list()
  
  # Iterate over methods
  for (method in method_names) {
    method_path <- file.path(base_path, method)
    
    # Get level-2 subdirectories (parameter combinations)
    parameter_combinations <- list.dirs(path = method_path, full.names = TRUE, recursive = FALSE)
    for (param_path in parameter_combinations) {
      # Check if the specified comparison exists (level-3 subdirectory)
      comparison_path <- file.path(param_path, comparison_name)
      if (!dir.exists(comparison_path)) next
      
      # Process each .tsv in that comparison folder
      tsv_files <- list.files(path = comparison_path, pattern = "\\.tsv$", full.names = TRUE)
      for (tsv_file in tsv_files) {
        tsv_data <- read.delim(tsv_file, header = TRUE, stringsAsFactors = FALSE)
        
        # Standardize column names
        if ("method"  %in% colnames(tsv_data)) names(tsv_data)[names(tsv_data) == "method"]  <- "DE_method"
        if ("logFC"   %in% colnames(tsv_data)) names(tsv_data)[names(tsv_data) == "logFC"]   <- "log2FC"
        
        # Ensure required columns
        if (!"average"    %in% colnames(tsv_data)) tsv_data$average    <- 0
        if (!"parameters" %in% colnames(tsv_data)) tsv_data$parameters <- "default"
        
        # Annotate with metadata
        tsv_data$Method     <- method
        tsv_data$Parameter  <- basename(param_path)
        tsv_data$Comparison <- comparison_name
        
        # Collect
        aggregated_data[[length(aggregated_data) + 1]] <- tsv_data
      }
    }
  }
  
  # Combine all and compute initial FC direction
  final_data <- do.call(rbind, aggregated_data)
  final_data$FC_direction <- ifelse(final_data$log2FC > 0, "up", "down")
  
  # Flip log2FC (and direction) for specified methods
  reverse_methods <- c("NBSR", "miRglmm")
  rev_idx <- final_data$Method %in% reverse_methods
  
  final_data$log2FC[rev_idx] <- -final_data$log2FC[rev_idx]
  final_data$FC_direction[rev_idx] <- ifelse(
    final_data$log2FC[rev_idx] > 0, 
    "up", 
    ifelse(final_data$log2FC[rev_idx] < 0, "down", final_data$FC_direction[rev_idx])
  )
  
  return(final_data)
}

# revese FC direction of potentially wrong methods
reverse_fc <- function(df, method_names) {
  # Identify rows to flip (changed: allow a vector of names)
  idx <- df$Method %in% method_names
  
  # Flip log2FC
  df$log2FC[idx] <- -df$log2FC[idx]
  
  # Reverse FC_direction
  df$FC_direction[idx] <- ifelse(
    df$log2FC[idx] > 0, 
    "up", 
    ifelse(df$log2FC[idx] < 0, "down", df$FC_direction[idx])
  )
  
  df
}


# Define the function to process and create ROC data
create_ROC_data <- function(ground_df, comparison_df, padj_threshold = 0.05) {
  # Check if input dataframes are non-empty
  if (nrow(ground_df) == 0 || nrow(comparison_df) == 0) {
    stop("Input dataframes are empty. Please provide valid data.")
  }
  
  # Check required columns
  required_columns <- c("Method", "miRNA", "padj", "pval")
  if (!all(required_columns %in% colnames(ground_df))) {
    stop("ground_df is missing required columns: Method, miRNA, padj, pval.")
  }
  if (!all(required_columns %in% colnames(comparison_df))) {
    stop("comparison_df is missing required columns: Method, miRNA, padj, pval.")
  }
  
  methods <- unique(ground_df$Method)
  roc_data <- list()
  
  for (method in methods) {
    # Subset the data for the current method
    ground_truth <- subset(ground_df, Method == method)
    test_data <- subset(comparison_df, Method == method)
    test_data_clean <- test_data[
      !is.na(test_data$log2FC) & !is.na(test_data$pval),
    ]
    print(method)
    
    # Skip methods with no data
    if (nrow(ground_truth) == 0 || nrow(test_data) == 0) {
      warning(paste("No data for method:", method))
      next
    }
    
    # Define true DE and non-DE based on padj threshold
    true_DE <- ground_truth$miRNA[ground_truth$padj < padj_threshold]
    true_no_DE <- ground_truth$miRNA[ground_truth$padj >= padj_threshold]
    print("here")
    # Check if thresholds exist
    thresholds <- sort(unique(test_data$pval))
    if (length(thresholds) == 0) {
      warning(paste("No thresholds found for method:", method))
      next
    }
    
    # Compute sensitivity and specificity
    for (threshold in thresholds) {
      predicted_DE <- test_data$miRNA[test_data$pval <= threshold]
      predicted_no_DE <- test_data$miRNA[test_data$pval > threshold]
      
      TP <- length(intersect(predicted_DE, true_DE))
      FP <- length(intersect(predicted_DE, true_no_DE))
      FN <- length(intersect(predicted_no_DE, true_DE))
      TN <- length(intersect(predicted_no_DE, true_no_DE))
      # print("hete 2")
      # print(TP)
      # print(FN)
      sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
      specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
      
      # Append results to roc_data list
      roc_data[[length(roc_data) + 1]] <- data.frame(
        Method = method,
        Threshold = threshold,
        Sensitivity = sensitivity,
        Specificity = specificity
      )
    }
  }
  
  # Combine all ROC data into a single dataframe
  if (length(roc_data) > 0) {
    roc_df <- do.call(rbind, roc_data)
    
    # Adjust specificity to avoid overlap at 1
    epsilon <- 1e-7
    specificity_1_indices <- which(roc_df$Specificity == 1)
    for (i in seq(length(specificity_1_indices) - 1, 1)) {
      roc_df$Specificity[specificity_1_indices[i]] <-
        roc_df$Specificity[specificity_1_indices[i]] - epsilon * (length(specificity_1_indices) - i)
    }
    return(roc_df)
  } else {
    warning("No ROC data generated. Check input data and thresholds.")
    return(NULL)
  }
}

# Plot ROC curve
plot_ROC <- function(roc_df, output_file, input_title = "") {
  gg <- ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity, color = Method)) +
    geom_line(linewidth = 1) +  # Use `linewidth` instead of `size`
    labs(x = "1 - Specificity (False Positive Rate)", 
         y = "Sensitivity (True Positive Rate)", 
         title = paste0("ROC Curve by Method ", input_title)) +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  # Save the plot to a file
  ggsave(output_file, plot = gg, width = 10, height = 7, dpi = 300)
  
  return(gg)
}

# Function to generate an UpSet plot using ComplexUpset
create_upset_plot_old <- function(input_df, padj_threshold = 0.05, output_file = "upset_plot.png") {
  
  # Identify all unique miRNAs and methods
  all_miRNAs <- unique(input_df$miRNA)
  all_methods <- unique(input_df$Method)
  
  # Filter for significant miRNAs based on padj threshold
  df_filtered <- input_df %>%
    filter(padj < padj_threshold) %>%
    select(miRNA, Method)
  
  # Convert to wide format (binary matrix)
  df_wide <- df_filtered %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = Method, values_from = present, values_fill = 0)
  
  # Ensure all miRNAs are included, even if absent from all methods
  df_complete <- expand.grid(miRNA = all_miRNAs) %>%
    left_join(df_wide, by = "miRNA") %>%
    replace(is.na(.), 0)  # Fill missing values with 0
  
  # Convert to data frame
  df_complete <- as.data.frame(df_complete)
  
  # Check if there is data to plot
  if (nrow(df_complete) == 0) {
    message("No significant miRNAs found for the given threshold.")
    return(NULL)
  }
  
  # Generate the UpSet plot
  p <- upset(
    df_complete,
    colnames(df_complete)[-1],  # Use all Method columns
    name = "Method Overlap",
    width_ratio = 0.2,
    min_size = 1
  ) + ggtitle(paste("UpSet Plot (padj <", padj_threshold, ")"))
  
  # Save plot to file
  ggsave(output_file, p, width = 8, height = 6, dpi = 300)
  
  # Return the plot object
  return(p)
}
create_upset_plot <- function(input_df, padj_threshold = 0.05, output_file = "upset_plot.png") {
  
  # Identify all unique miRNAs
  all_miRNAs <- unique(input_df$miRNA)
  
  # Filter for significant miRNAs and create a combined method_direction label
  df_filtered <- input_df %>%
    filter(padj < padj_threshold) %>%                   # Same as before: filter by padj
    mutate(method_direction = paste(Method, FC_direction, sep = "_")) %>%  # <-- NEW: Combine Method and FC_direction
    select(miRNA, method_direction) %>%                 # <-- UPDATED: Select the new label instead of just Method
    distinct()
  
  # Convert to wide format (binary matrix)
  df_wide <- df_filtered %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = method_direction, values_from = present, values_fill = 0)  # <-- UPDATED: Pivot using method_direction
  
  # Ensure all miRNAs are included, even if absent from all sets
  df_complete <- expand.grid(miRNA = all_miRNAs) %>%
    left_join(df_wide, by = "miRNA") %>%
    replace(is.na(.), 0)
  
  # Convert to data frame
  df_complete <- as.data.frame(df_complete)
  
  # Check if there is data to plot
  if (nrow(df_complete) == 0) {
    message("No significant miRNAs found for the given threshold.")
    return(NULL)
  }
  
  # Generate the UpSet plot using the method_direction columns
  p <- upset(
    df_complete,
    colnames(df_complete)[-1],
    name = "Method_Direction Overlap",  # <-- UPDATED: Changed title to reflect method_direction
    width_ratio = 0.2,
    min_size = 1
  ) + ggtitle(paste("UpSet Plot (padj <", padj_threshold, ")"))
  
  # Save plot to file
  ggsave(output_file, p, width = 8, height = 6, dpi = 300)
  
  return(p)
}
create_upset_plot <- function(input_df,
                              padj_threshold   = 0.05,
                              min_intersections = 2,          
                              output_file      = "upset_plot.png") {
  
  # Identify all unique miRNAs
  all_miRNAs <- unique(input_df$miRNA)
  
  # Filter for significant miRNAs and create a combined method_direction label
  df_filtered <- input_df %>%
    filter(padj < padj_threshold) %>%
    mutate(method_direction = paste(Method, FC_direction, sep = "_")) %>%
    select(miRNA, method_direction) %>%
    distinct()
  
  # Convert to wide format (binary matrix)
  df_wide <- df_filtered %>%
    mutate(present = 1) %>%
    pivot_wider(
      names_from   = method_direction,
      values_from  = present,
      values_fill  = 0
    )
  
  # Ensure all miRNAs are included
  df_complete <- expand.grid(miRNA = all_miRNAs) %>%
    left_join(df_wide, by = "miRNA") %>%
    replace(is.na(.), 0) %>%
    as.data.frame()
  
  # Check for data
  if (nrow(df_complete) == 0) {
    message("No significant miRNAs found for the given threshold.")
    return(NULL)
  }
  
  # Generate the UpSet plot with the user‐defined minimum intersection size
  p <- upset(
    df_complete,
    colnames(df_complete)[-1],
    name       = "Method_Direction Overlap",
    width_ratio = 0.2,
    min_size   = min_intersections         # Use the new parameter here
  ) +
    ggtitle(paste0("UpSet Plot (padj < ", padj_threshold, 
                   "; min intersection = ", min_intersections, ")"))
  
  # Save to file
  ggsave(output_file, p, width = 8, height = 7, dpi = 300)
  
  return(p)
}



# Function to generate an UpSet plot using ComplexUpset using parameters for sets

create_upset_plot_by_parameters_old <- function(input_df, padj_threshold = 0.05, output_file = "upset_plot.png") {
  
  # Identify all unique miRNAs and parameters
  all_miRNAs <- unique(input_df$miRNA)
  all_parameters <- unique(input_df$parameters)
  
  # Filter for significant miRNAs based on padj threshold
  df_filtered <- input_df %>%
    filter(padj < padj_threshold) %>%
    select(miRNA, parameters)
  
  # Convert to wide format (binary matrix)
  df_wide <- df_filtered %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = parameters, values_from = present, values_fill = 0)
  
  # Ensure all miRNAs are included, even if absent from all parameters
  df_complete <- expand.grid(miRNA = all_miRNAs) %>%
    left_join(df_wide, by = "miRNA") %>%
    replace(is.na(.), 0)  # Fill missing values with 0
  
  # Convert to data frame
  df_complete <- as.data.frame(df_complete)
  
  # Check if there is data to plot
  if (nrow(df_complete) == 0) {
    message("No significant miRNAs found for the given threshold.")
    return(NULL)
  }
  
  # Generate the UpSet plot
  p <- upset(
    df_complete,
    colnames(df_complete)[-1],  # Use all parameters columns
    name = "Parameter Overlap",
    width_ratio = 0.2,
    min_size = 1
  ) + ggtitle(paste("UpSet Plot by Parameters (padj <", padj_threshold, ")"))
  
  # Save plot to file
  ggsave(output_file, p, width = 12, height = 8, dpi = 300)
  
  # Return the plot object
  return(p)
}

create_upset_plot_by_parameters <- function(input_df, padj_threshold = 0.05, output_file = "upset_plot.png") {
  
  # Filter for significant miRNAs based on padj threshold
  df_filtered <- input_df %>%
    filter(padj < padj_threshold) %>%
    select(miRNA, parameters)
  
  # Count occurrences of each miRNA across different parameter sets
  miRNA_counts <- table(df_filtered$miRNA)
  
  # Determine the total number of unique parameters
  num_parameters <- length(unique(input_df$parameters))
  
  # Consensus miRNAs: those that appear in all parameter sets
  consensus_miRNAs <- names(miRNA_counts[miRNA_counts == num_parameters])
  
  # Convert df to wide format for UpSet plot
  df_wide <- df_filtered %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = parameters, values_from = present, values_fill = 0)
  
  # Ensure all miRNAs are included, even if absent from all parameters
  df_complete <- as.data.frame(df_wide)
  
  # Check if there is data to plot
  if (nrow(df_complete) == 0) {
    message("No significant miRNAs found for the given threshold.")
    return(NULL)
  }
  
  # Generate the UpSet plot
  p <- upset(
    df_complete,
    colnames(df_complete)[-1],  # Use all parameter columns
    name = "Parameter Overlap",
    width_ratio = 0.2,
    min_size = 1
  ) + ggtitle(paste("UpSet Plot by Parameters (padj <", padj_threshold, ")"))
  
  # Save plot to file
  ggsave(output_file, p, width = 8, height = 6, dpi = 300)
  
  # Return the plot and consensus miRNAs
  return(list(plot = p, consensus = consensus_miRNAs))
}

# Function to generate an UpSet plot using ComplexUpset with DESeq consensus
create_upset_plot_with_deseq <- function(original_df, consensus_miRNAs, padj_threshold = 0.05, output_file = "upset_plot.png") {
  
  # Generate a "fake" DESeq dataset
  fake_deseq_df <- data.frame(
    miRNA = consensus_miRNAs,  # Use the consensus miRNAs
    padj = runif(length(consensus_miRNAs), min = 0.0001, max = padj_threshold - 0.0001),  # Ensure padj < threshold
    Method = "DESeq"
  )
  
  # Merge fake DESeq data with original dataset
  combined_df <- bind_rows(original_df, fake_deseq_df)
  
  # Identify all unique miRNAs and methods
  all_miRNAs <- unique(combined_df$miRNA)
  all_methods <- unique(combined_df$Method)
  
  # Convert to a binary presence/absence table
  df_filtered <- combined_df %>%
    filter(padj < padj_threshold) %>%
    select(miRNA, Method)
  
  df_wide <- df_filtered %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = Method, values_from = present, values_fill = 0)
  
  # Ensure all miRNAs are included, even if absent from all methods
  df_complete <- expand.grid(miRNA = all_miRNAs) %>%
    left_join(df_wide, by = "miRNA") %>%
    replace(is.na(.), 0)  # Fill missing values with 0
  
  # Convert to data frame
  df_complete <- as.data.frame(df_complete)
  
  # Check if there is data to plot
  if (nrow(df_complete) == 0) {
    message("No significant miRNAs found for the given threshold.")
    return(NULL)
  }
  
  # Generate the UpSet plot
  p <- upset(
    df_complete,
    colnames(df_complete)[-1],  # Use all Method columns
    name = "Method Overlap",
    width_ratio = 0.2,
    min_size = 1
  ) + ggtitle(paste("UpSet Plot (padj <", padj_threshold, ")"))
  
  # Save plot to file
  ggsave(output_file, p, width = 8, height = 6, dpi = 300)
  
  # Find miRNAs detected by all methods (common_miRNAs)
  common_miRNAs <- df_complete %>%
    filter(rowSums(select(., -miRNA)) == length(all_methods)) %>%
    pull(miRNA)
  
  # Find miRNAs detected by no methods (unidentified_miRNAs)
  unidentified_miRNAs <- df_complete %>%
    filter(rowSums(select(., -miRNA)) == 0) %>%
    pull(miRNA)
  
  # Return the plot, common miRNAs, and unidentified miRNAs
  return(list(plot = p, common_miRNAs = common_miRNAs, unidentified_miRNAs = unidentified_miRNAs))
}

# Define the function to process and create ROC data using common ground truth
create_ROC_data_common <- function(comparison_df, negative_set, positive_set, padj_threshold = 0.05) {
  
  # Check if input dataframe is non-empty
  if (nrow(comparison_df) == 0) {
    stop("comparison_df is empty. Please provide valid data.")
  }
  
  # Check required columns
  required_columns <- c("Method", "miRNA", "padj", "pval")
  if (!all(required_columns %in% colnames(comparison_df))) {
    stop("comparison_df is missing required columns: Method, miRNA, padj, pval.")
  }
  
  methods <- unique(comparison_df$Method)
  roc_data <- list()
  
  for (method in methods) {
    # Subset the data for the current method
    test_data <- subset(comparison_df, Method == method)
    
    # Skip methods with no data
    if (nrow(test_data) == 0) {
      warning(paste("No data for method:", method))
      next
    }
    
    # Ensure negative and positive sets are properly defined
    true_DE <- intersect(positive_set, test_data$miRNA)
    true_no_DE <- intersect(negative_set, test_data$miRNA)
    
    # Sort unique p-values as thresholds
    thresholds <- sort(unique(test_data$pval))
    if (length(thresholds) == 0) {
      warning(paste("No thresholds found for method:", method))
      next
    }
    
    # Compute sensitivity and specificity
    for (threshold in thresholds) {
      predicted_DE <- test_data$miRNA[test_data$pval <= threshold]
      predicted_no_DE <- test_data$miRNA[test_data$pval > threshold]
      
      TP <- length(intersect(predicted_DE, true_DE))
      FP <- length(intersect(predicted_DE, true_no_DE))
      FN <- length(intersect(predicted_no_DE, true_DE))
      TN <- length(intersect(predicted_no_DE, true_no_DE))
      
      sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
      specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
      
      # Append results to roc_data list
      roc_data[[length(roc_data) + 1]] <- data.frame(
        Method = method,
        Threshold = threshold,
        Sensitivity = sensitivity,
        Specificity = specificity
      )
    }
  }
  
  # Combine all ROC data into a single dataframe
  if (length(roc_data) > 0) {
    roc_df <- do.call(rbind, roc_data)
    
    # Adjust specificity to avoid overlap at 1
    epsilon <- 1e-7
    specificity_1_indices <- which(roc_df$Specificity == 1)
    for (i in seq(length(specificity_1_indices) - 1, 1)) {
      roc_df$Specificity[specificity_1_indices[i]] <-
        roc_df$Specificity[specificity_1_indices[i]] - epsilon * (length(specificity_1_indices) - i)
    }
    return(roc_df)
  } else {
    warning("No ROC data generated. Check input data and thresholds.")
    return(NULL)
  }
}

# Define the function to process and create ROC data using monotonic increase/decrease as common ground truth
create_ROC_data_mono <- function(comparison_df, ground_truth_df, padj_threshold = 0.05) {
  
  # Check if input dataframe is non-empty
  if (nrow(comparison_df) == 0) {
    stop("comparison_df is empty. Please provide valid data.")
  }
  if (nrow(ground_truth_df) == 0) {
    stop("ground_truth_df is empty. Please provide valid data.")
  }
  
  # Check required columns in comparison_df
  required_columns_pred <- c("Method", "miRNA", "padj", "pval", "FC_direction")
  if (!all(required_columns_pred %in% colnames(comparison_df))) {
    stop("comparison_df is missing required columns: Method, miRNA, padj, pval, FC_direction.")
  }
  
  # Check required columns in ground_truth_df
  required_columns_truth <- c("miRNA", "FC_direction")
  if (!all(required_columns_truth %in% colnames(ground_truth_df))) {
    stop("ground_truth_df is missing required columns: miRNA, FC_direction.")
  }
  
  # Merge the two dataframes by miRNA.
  # Suffixes: _pred for comparison and _truth for ground truth.
  merged_df <- merge(comparison_df, ground_truth_df[, c("miRNA", "FC_direction")],
                     by = "miRNA", suffixes = c("_pred", "_truth"))
  if(nrow(merged_df) == 0){
    stop("No overlapping miRNAs between comparison_df and ground_truth_df.")
  }
  
  # Determine true class from ground truth:
  # DE if FC_direction is "up" or "down", non-DE if "neither"
  merged_df$true_class <- ifelse(merged_df$FC_direction_truth %in% c("up", "down"), "DE", "nonDE")
  # print merged_df
  # input()
  
  methods <- unique(merged_df$Method)
  roc_data <- list()
  
  for (method in methods) {
    print(method)
    # Subset the data for the current method
    test_data <- subset(merged_df, Method == method)
    test_data <- test_data[
      !is.na(test_data$log2FC) & !is.na(test_data$pval),
    ]
    print(dim(test_data))
    
    ## export test_data to TSV for method “NBSR”
       if (method == "NBSR") {
          write.table(
               test_data,
              file      = "NBSR_test_data.tsv",
             sep       = "\t",
              row.names = FALSE,
              quote     = FALSE
          )
       }
    # print (TP, TN, etc)
    if (nrow(test_data) == 0) {
      warning(paste("No data for method:", method))
      next
    }
    
    # Use unique sorted p-values from predictions as thresholds
    thresholds <- c(0,sort(unique(test_data$pval)))
    if (length(thresholds) == 0) {
      warning(paste("No thresholds found for method:", method))
      next
    }
    
    for (threshold in thresholds) {
      # Predicted as DE if pval <= threshold
      is_predicted_DE <- test_data$pval <= threshold
      
      # Among predicted DE, check if the predicted FC_direction matches ground truth.
      # A true positive occurs if:
      #   - The miRNA is truly DE (true_class == "DE")
      #   - It is predicted as DE (pval <= threshold)
      #   - Its predicted FC_direction equals the ground truth FC_direction.
      TP <- sum(is_predicted_DE & 
                  (test_data$true_class == "DE") & 
                  (test_data$FC_direction_pred == test_data$FC_direction_truth))
      
      # A false positive is counted when:
      #   - The miRNA is predicted as DE but is either truly nonDE or its FC_direction is mismatched.
      FP <- sum(is_predicted_DE & 
                  ((test_data$true_class == "nonDE") | 
                     (test_data$true_class == "DE" & test_data$FC_direction_pred != test_data$FC_direction_truth)))
      
      # False negatives: truly DE miRNAs that were not predicted as DE.
      FN <- sum((!is_predicted_DE) & (test_data$true_class == "DE"))
      
      # True negatives: miRNAs not predicted as DE that are truly nonDE.
      TN <- sum((!is_predicted_DE) & (test_data$true_class == "nonDE"))
      
      # Sensitivity and specificity calculations.
      # sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA
      # specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA
      print(method)
      print(TP)
      print(FN)
      print(TN)
      print(FP)
      TP <- ifelse(is.na(TP), 0, TP)
      FN <- ifelse(is.na(FN), 0, FN)
      TN <- ifelse(is.na(TN), 0, TN)
      FP <- ifelse(is.na(FP), 0, FP)
      
      sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else 0
      specificity <- if ((TN + FP) > 0) TN / (TN + FP) else 0
      
      roc_data[[length(roc_data) + 1]] <- data.frame(
        Method = method,
        Threshold = threshold,
        TP = TP,
        FP = FP,
        FN = FN,
        TN = TN,
        Sensitivity = sensitivity,
        Specificity = specificity
      )
    }
  }
  
  # Combine all ROC data into one dataframe.
  if (length(roc_data) > 0) {
    roc_df <- do.call(rbind, roc_data)
    
    # Adjust specificity to avoid duplicate values of 1.
    epsilon <- 1e-7
    specificity_1_indices <- which(roc_df$Specificity == 1)
    if(length(specificity_1_indices) > 1){
      for (i in seq(length(specificity_1_indices) - 1, 1)) {
        roc_df$Specificity[specificity_1_indices[i]] <-
          roc_df$Specificity[specificity_1_indices[i]] - epsilon * (length(specificity_1_indices) - i)
      }
    }
    return(roc_df)
  } else {
    warning("No ROC data generated. Check input data and thresholds.")
    return(NULL)
  }
}

# perform rlog

perform_rlog <- function(dds) {
  # Ensure the input is a DESeqDataSet object
  if (!inherits(dds, "DESeqDataSet")) {
    stop("Input must be a DESeqDataSet object")
  }
  
  # Define parameter combinations
  fit_types <- c("parametric", "local", "mean", "glmGamPoi")
  blinds <- c(TRUE, FALSE)
  
  # Initialize an empty list to store results
  rlog_results <- list()
  
  # Iterate over all combinations
  for (fit_type in fit_types) {
    for (blind in blinds) {
      # Generate a descriptive key name for this combination
      key_name <- paste0("fitType_", fit_type, "_blind_", blind)
      
      # Perform rlog transformation
      rlog_results[[key_name]] <- tryCatch({
        rlog(dds, blind = blind, fitType = fit_type)
      }, error = function(e) {
        message(paste("Error in combination:", key_name, "-", e$message))
        return(NULL)
      })
    }
  }
  
  # Return the list of results
  return(rlog_results)
}

# perform vst
perform_vst <- function(dds) {
  # Ensure the input is a DESeqDataSet object
  if (!inherits(dds, "DESeqDataSet")) {
    stop("Input must be a DESeqDataSet object")
  }
  
  # Define parameter combinations
  fit_types <- c("parametric", "local", "mean", "glmGamPoi")
  blinds <- c(TRUE, FALSE)
  
  # Initialize an empty list to store results
  vst_results <- list()
  
  # Iterate over all combinations
  for (fit_type in fit_types) {
    for (blind in blinds) {
      # Generate a descriptive key name for this combination
      key_name <- paste0("fitType_", fit_type, "_blind_", blind)
      
      # Perform VST
      vst_results[[key_name]] <- tryCatch({
        varianceStabilizingTransformation(dds, blind = blind, fitType = fit_type)
      }, error = function(e) {
        message(paste("Error in combination:", key_name, "-", e$message))
        return(NULL)
      })
    }
  }
  
  # Return the list of results
  return(vst_results)
}

# 
# if (TRUE){
#   groups <- list()
#   
#   for (x in c(20, 40, 60, 80)) {
#     current_group <- as.character(x)
#     print(current_group)
#     columns <- grep(paste0("^",current_group, "(?!.*_2$)"), colnames(input_df), perl = TRUE)
#     members <- colnames(input_df)[columns]
#     print(members)
#     groups[[paste0(current_group, "S")]] <- members
#     
#   }
#   columns <- grep("^[1-5]L(?!.*_.*_)", colnames(input_df), perl = TRUE)
#   groups[["0S"]] <- colnames(input_df)[columns]
#   
#   columns <- grep("^[1-5]S(?!.*_2$)", colnames(input_df), perl = TRUE)
#   groups[["100S"]] <- colnames(input_df)[columns]
# }