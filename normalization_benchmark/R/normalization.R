# library(DESeq2)
# change column names
fix_names <- function(input_cols){
  return(gsub("\\_R1\\|grp1", "", input_cols))
}

# define groups
define_groups <- function(input_df) {
  groups <- list()
  for (x in c(20, 40, 60, 80)) {
    current_group <- as.character(x)
    print(current_group)
    # columns <- grep(paste0("^",current_group), colnames(input_df))
    columns <- grep(paste0("^",current_group, "(?!.*_2$)"), colnames(input_df), perl = TRUE)
    members <- colnames(input_df)[columns]
    print(members)
    groups[[paste0(current_group, "S")]] <- members
    
  }
  columns <- grep("^[1-5]L(?!.*_.*_)", colnames(input_df), perl = TRUE)
  groups[["0S"]] <- colnames(input_df)[columns]
  columns <- grep("^[1-5]S(?!.*_2$)", colnames(input_df), perl = TRUE)
  groups[["100S"]] <- colnames(input_df)[columns]
  group_list <- names(groups)
  return(list(group_members = groups, group_list = group_list))
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

# function to aggregate per group

geometric_mean_by_group <- function(data, group_list) {
  # Initialize an empty list to store the results
  group_means <- list()
  
  # Iterate over each group in the list
  for (group_name in names(group_list)) {
    # Get the column names for the current group
    cols <- group_list[[group_name]]
    
    # Extract the relevant columns
    group_data <- data[, cols, drop = FALSE]
    
    # Calculate the geometric mean across rows
    geom_mean <- apply(group_data, 1, function(row) {
      exp(mean(log(row + 1), na.rm = TRUE)) # Calculate geometric mean for each row
    })
    
    # Store the result
    group_means[[group_name]] <- geom_mean
  }
  
  # Combine the results into a new dataframe
  result_df <- as.data.frame(group_means)
  colnames(result_df) <- names(group_list)
  desired_order <- c("0S", "20S", "40S", "60S", "80S", "100S")
  result_df <- result_df[, desired_order]
  
  # Return the result
  return(result_df)
}

# check if a row is monotonic

check_monotonicity <- function(x) {
  if (all(diff(x) > 0)) {
    return("increasing")
  } else if (all(diff(x) < 0)) {
    return("decreasing")
  } else {
    return("neither")
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

# calculate expected values (average)
# calculate expected values (per mouse)