library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)


source("R/normalization.R")
plot_dir <- "./plots/"
dir.create(plot_dir, recursive = TRUE)

# load RC matrix

RC_df <- read.delim("./data/mature_sense_minExpr0_RCadj.mat",
                    check.names = FALSE, row.names = 1)

# fix names
colnames(RC_df) <- fix_names(colnames(RC_df))

print(head(RC_df))

# separate groups
groups_object <- define_groups(RC_df)

# some mixture_group are missing

column_data <- data.frame(row.names = colnames(RC_df))
column_data$mixture_group <- sapply(rownames(column_data), 
                                    function(x) get_group(x,groups_object$group_members))

# log transform
logged_df <- log2(RC_df + 1)

# RPM total
RPMt_df <- read.delim("./data/mature_sense_minExpr0_RCadj_totalRPM.mat",
                      check.names = FALSE, row.names = 1)
colnames(RPMt_df) <- fix_names(colnames(RPMt_df))
# RPM lib

RPMl_df <- read.delim("./data/mature_sense_minExpr0_RCadj_libraryRPM.mat",
                      check.names = FALSE, row.names = 1)
colnames(RPMl_df) <- fix_names(colnames(RPMl_df))

list_of_transformations <- list(RC = RC_df, 
                                # logged = logged_df, 
                                RPM_total = RPMt_df
                                # RPM_lib = RPMl_df
                                )


# for blind False
# colData add another column (mixture group)
# design = ~ mixture_group)

dds <- DESeqDataSetFromMatrix(countData = RC_df, 
                              colData = data.frame(row.names = colnames(RC_df)), 
                              design = ~ 1)
                              
# DESEq transformations
if(TRUE){
# if(FALSE){
  # fitType = c("parametric", "local", "mean", "glmGamPoi")
  fit_types <- c("parametric", "local", "mean", "glmGamPoi")
  blind_options <- c(TRUE, FALSE)
  
  # Helper function to handle errors
  safe_transform <- function(transform_func, ...) {
    tryCatch(
      expr = transform_func(...),
      error = function(e) {
        message("Error in transformation: ", conditionMessage(e))
        return(NULL)  # Return NULL if an error occurs
      }
    )
  }
  
  # Add rlog transformations
  for (blind in blind_options) {
    if  (blind){
      dds <- DESeqDataSetFromMatrix(countData = RC_df, 
                                    colData = column_data, 
                                    design = ~ 1)
    }else{
      dds <- DESeqDataSetFromMatrix(countData = RC_df, 
                                    colData = column_data, 
                                    design = ~ mixture_group)
    }
    allowed_combinations <- c(
      "rlog_mean_blind_FALSE",
      "rlog_parametric_blind_TRUE",
      "vst_local_blind_TRUE"
    )
    for (fit in fit_types) {
      transformation_name <- paste0("rlog_", fit, "_blind_", blind)
      # print(paste0("rlog_", fit, "_blind_", blind))
      
      # print(head(result))
      
      if (transformation_name %in% allowed_combinations){
        result <- as.data.frame(assay(safe_transform(rlog, dds, fitType = fit, blind = blind)))
        list_of_transformations[[transformation_name]] <- 2^result  
      }
      
    }
  }
  fit_types <- c("parametric", "local", "mean")
  # Add VST transformations
  for (blind in blind_options) {
    if  (blind){
      dds <- DESeqDataSetFromMatrix(countData = RC_df, 
                                    colData = column_data, 
                                    design = ~ 1)
    }else{
      dds <- DESeqDataSetFromMatrix(countData = RC_df, 
                                    colData = column_data, 
                                    design = ~ mixture_group)
    }
    for (fit in fit_types) {
      # print(paste0("vst_", fit, "_blind_", blind))
      
      # print(head(result))
      transformation_name <- paste0("vst_", fit, "_blind_", blind)
      if (transformation_name %in% allowed_combinations){
        result <- as.data.frame(assay(safe_transform(varianceStabilizingTransformation, dds, fitType = fit, blind = blind)))
      list_of_transformations[[transformation_name]] <- 2^result
      }
    }
  }
}

# edgeR transformation TMM
dge <- DGEList(counts = RC_df)

# Calculate normalization factors using TMM
dge <- calcNormFactors(dge, method = "TMM")

# To extract the TMM-normalized values (on a per-million scale)
tmm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
list_of_transformations[["TMM"]] <- as.data.frame(tmm_counts)

# define groups
group_obj <- define_groups(RC_df)
groups <- group_obj$group_members
group_list <- group_obj$group_list

# make dataframe of reliable FC for monotonicity evaluation
if (TRUE){
  group_RPMs <- list()
  for (group in group_list) {
    # Get the samples corresponding to this group
    group_samples <- groups[[group]]
    
    # Subset the input_df to only the columns corresponding to this group
    group_RPM <- RPMl_df[, group_samples]
    
    # Store the RPM values for this group
    group_RPMs[[group]] <- group_RPM
  }
  
  # Initialize an empty matrix to store the average RPMs
  group_avg_RPM <- matrix(NA, nrow = nrow(RPMl_df), ncol = length(group_list))
  rownames(group_avg_RPM) <- rownames(RPMl_df)
  colnames(group_avg_RPM) <- group_list
  
  for (group in group_list) {
    # Extract the RPM data for this group
    group_RPM <- group_RPMs[[group]]
    
    # Calculate the geometric mean for each miRNA (log space, then exponentiate)
    log_RPM <- log2(group_RPM + 1)
    log_means <- rowMeans(log_RPM)
    group_avg_RPM[, group] <- log_means
  }
  
  group_avg_RPM <- 2 ^ group_avg_RPM
  
  group_extended <- group_avg_RPM
  group_extended <- as.data.frame(group_extended)
  
  # Add the average row value column
  group_extended$average <- rowMeans(group_extended)
  
  # Add the FC (fold change) column between 0S and 100S
  group_extended$FC <- group_extended$`100S` / (group_extended$`0S`) 
  
  group_extended$log2FC <- log2(group_extended$FC)
  group_extended$FC_direction <- ifelse(group_extended$log2FC > 0, "increasing", "decreasing")
  group_extended$abs_log2FC <- abs(group_extended$log2FC)
}

# 3 sets
# |logFC| > 0.1, .5 1 
# min 50 RPM average
# 
# increase_df_01 <- group_extended[abs(group_extended$log2FC) > 0.1, ]
# increase_df_05 <- group_extended[abs(group_extended$log2FC) > 0.5, ]
# increase_df_1 <- group_extended[abs(group_extended$log2FC) > 1, ]

## check monotonic trend for each method

if (TRUE){
  monotonic_list <- list()
  
  # for every key in transformations
  # calculate group geometric average
  # calculate monotonicity of each miRNA
  
  final_results <- data.frame(
    transformation = character(),
    top_n = integer(),
    percentage_agreement = numeric(),
    stringsAsFactors = FALSE
  )
  
  final_log2 <- data.frame(
    transformation = character(),
    top_n = integer(),
    percentage_agreement = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(list_of_transformations)) {
    # Access the current dataframe
    key_name <- names(list_of_transformations)[i]
    temp_df <- list_of_transformations[[i]]
    # make same row order (miRNAs)
    group_df <- geometric_mean_by_group(temp_df, groups)
    results_df <- data.frame(matrix(nrow = nrow(temp_df), ncol = 0))
    rownames(results_df) <- rownames(group_df)
    results_df$monotonicity <- apply(group_df, 1, check_monotonicity)
    # add putative column
    results_df$putative_direction <- group_extended$FC_direction
    # average expression column
    results_df$average_RPM <- group_extended$average
    # add log2FC column
    results_df$log2FC <- group_extended$abs_log2FC
    # calculate agreement
    results_df$agreement <- ifelse(results_df$monotonicity == results_df$putative_direction, 1, 0)
    filtered_df <- results_df[results_df$average_RPM >= 50, ]
    # Sort by average_RPM in descending order
    sorted_df <- filtered_df[order(-filtered_df$average_RPM), ]
    top_n_values <- c(1, 5, 10, 20, 50, seq(60, nrow(sorted_df), by = 10))
    
    for (top_n in top_n_values) {
      # Ensure top_n does not exceed the number of rows
      top_n <- min(top_n, nrow(sorted_df))
      
      # Subset the top N rows
      top_subset <- sorted_df[1:top_n, ]
      
      # Calculate the percentage of agreement
      percentage_agreement <- mean(top_subset$agreement) * 100
      
      # Add results to final_results
      final_results <- rbind(final_results, data.frame(
        transformation = key_name,
        top_n = top_n,
        percentage_agreement = percentage_agreement,
        stringsAsFactors = FALSE
      ))
    }
    
    sorted_df <- filtered_df[order(-filtered_df$log2FC), ]
    
    for (top_n in top_n_values) {
      # Ensure top_n does not exceed the number of rows
      top_n <- min(top_n, nrow(sorted_df))
      
      # Subset the top N rows
      top_subset <- sorted_df[1:top_n, ]
      
      # Calculate the percentage of agreement
      percentage_agreement <- mean(top_subset$agreement) * 100
      
      # Add results to final_results
      final_log2 <- rbind(final_log2, data.frame(
        transformation = key_name,
        top_n = top_n,
        percentage_agreement = percentage_agreement,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # print(head(final_results, 20))
  # print(head(final_log2, 20))
}

# Define a limited set of colors and shapes
# color_palette <- scales::hue_pal()(7)  # Use 7 distinct colors
# shape_palette <- c(16, 17, 18)         # Use 3 distinct shapes (circle, triangle, diamond)
# 
# # Map transformations to a combined color and shape grouping
# final_results$color_group <- factor(as.numeric(final_results$transformation) %% 7 + 1)  # 7 colors
# final_results$shape_group <- factor(as.numeric(final_results$transformation) %% 3 + 1)  # 3 shapes
# 
# final_log2$color_group <- factor(as.numeric(final_log2$transformation) %% 7 + 1)  # 7 colors
# final_log2$shape_group <- factor(as.numeric(final_log2$transformation) %% 3 + 1)  # 3 shapes


gg <- ggplot(final_results, aes(x = top_n, y = percentage_agreement, color = transformation)) +
  geom_line() +
  geom_point() +
  # scale_x_log10() +  # Optional: Log scale for top_n
  labs(x = "Top n miRNAs (expression)", y = "Percentage Monotonic (%)", 
       title = "percentage monotonic expression of the top expressed miRNAs") +
  theme_minimal()


pdf(paste0(plot_dir, "percentage_monotonic_expression.pdf"), width = 8, height = 6)  
print(gg)  
dev.off()  

gg <- ggplot(final_log2, aes(x = top_n, y = percentage_agreement, color = transformation)) +
  geom_line() +
  geom_point() +
  # scale_x_log10() +  # Optional: Log scale for top_n
  labs(x = "Top n miRNAs (absolute log2FC)", y = "Percentage Monotonic (%)", 
       title = "percentage of monotonic expression of the top DE miRNAs") +
  theme_minimal()

pdf(paste0(plot_dir, "percentage_monotonic_log2FC.pdf"), width = 8, height = 6)  
print(gg)  
dev.off() 

#######
##### repeat analysis but relative (using each transformation as internal reference)
#######

if (TRUE){
  final_results <- data.frame(
    transformation = character(),
    top_n = integer(),
    percentage_agreement = numeric(),
    stringsAsFactors = FALSE
  )
  
  final_log2 <- data.frame(
    transformation = character(),
    top_n = integer(),
    percentage_agreement = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(list_of_transformations)) {
    # Access the current dataframe
    key_name <- names(list_of_transformations)[i]
    temp_df <- list_of_transformations[[i]]
    # make same row order (miRNAs)
    group_df <- geometric_mean_by_group(temp_df, groups)
    results_df <- data.frame(matrix(nrow = nrow(temp_df), ncol = 0))
    rownames(results_df) <- rownames(group_df)
    results_df$monotonicity <- apply(group_df, 1, check_monotonicity)
    # calculate absolute log2FC
    results_df$FC <- group_df$`100S` / (group_df$`0S`)
    results_df$log2FC <- log2(results_df$FC)
    results_df$putative_direction <- ifelse(results_df$log2FC > 0, "increasing", "decreasing")
    results_df$abs_log2FC <- abs(results_df$log2FC)
    results_df$average <- rowMeans(group_df)
    results_df$average_RPM <- group_extended$average
    # calculate agreement
    results_df$agreement <- ifelse(results_df$monotonicity == results_df$putative_direction, 1, 0)
    # use RPM to filter so it's equivalent to all transformations
    filtered_df <- results_df[results_df$average_RPM >= 50, ]
    # Sort by average_RPM in descending order
    sorted_df <- filtered_df[order(-filtered_df$average), ]
    top_n_values <- c(1, 5, 10, 20, 50, seq(60, nrow(sorted_df), by = 10))
    
    for (top_n in top_n_values) {
      # Ensure top_n does not exceed the number of rows
      top_n <- min(top_n, nrow(sorted_df))
      
      # Subset the top N rows
      top_subset <- sorted_df[1:top_n, ]
      
      # Calculate the percentage of agreement
      percentage_agreement <- mean(top_subset$agreement) * 100
      
      # Add results to final_results
      parts <- unlist(strsplit(key_name, "_"))
      tgroup = parts[1]  # Return the second keyword
      final_results <- rbind(final_results, data.frame(
        transformation = key_name,
        transformation_group = tgroup,
        top_n = top_n,
        percentage_agreement = percentage_agreement,
        stringsAsFactors = FALSE
      ))
    }
    
    sorted_df <- filtered_df[order(-filtered_df$abs_log2FC), ]
    
    for (top_n in top_n_values) {
      # Ensure top_n does not exceed the number of rows
      top_n <- min(top_n, nrow(sorted_df))
      
      # Subset the top N rows
      top_subset <- sorted_df[1:top_n, ]
      
      # Calculate the percentage of agreement
      percentage_agreement <- mean(top_subset$agreement) * 100
      
      # Add results to final_results
      parts <- unlist(strsplit(key_name, "_"))
      tgroup = parts[1]  # Return the second keyword
      
      final_log2 <- rbind(final_log2, data.frame(
        transformation = key_name,
        transformation_group = tgroup,
        top_n = top_n,
        percentage_agreement = percentage_agreement,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  print(head(final_results, 20))
  print(head(final_log2, 20))
  
  old_group_extended <- group_extended
  
  
}

gg <- ggplot(final_results, aes(x = top_n, y = percentage_agreement, 
                                color = transformation, shape = transformation_group ) ) +
  geom_line() +
  geom_point() +
  # scale_x_log10() +  # Optional: Log scale for top_n
  labs(x = "Top n miRNAs (expression)", y = "Percentage Monotonic (%)", 
       title = "") +
  theme_minimal()

pdf(paste0(plot_dir, "percentage_monotonic_expres_relative.pdf"), width = 8, height = 6)  
print(gg)  
dev.off() 

gg <- ggplot(final_log2, aes(x = top_n, y = percentage_agreement, 
                             color = transformation, shape = transformation_group ) ) +
  geom_line() +
  geom_point() +
  # scale_x_log10() +  # Optional: Log scale for top_n
  labs(x = "Top n miRNAs (absolute log2FC)", y = "Percentage Monotonic (%)", 
       title = "percentage of monotonic expression of the top DE miRNAs internally calculated") +
  theme_minimal()

pdf(paste0(plot_dir, "percentage_monotonic_log2FC_relative.pdf"), width = 8, height = 6)
print(gg)  
dev.off()  

#
##### repeat analysis but using rlog_mean_blind_FALSE as reference
#

if (TRUE){
  final_results <- data.frame(
    transformation = character(),
    top_n = integer(),
    percentage_agreement = numeric(),
    stringsAsFactors = FALSE
  )
  
  final_log2 <- data.frame(
    transformation = character(),
    top_n = integer(),
    percentage_agreement = numeric(),
    stringsAsFactors = FALSE
  )
  
  temp_df <- list_of_transformations[["rlog_mean_blind_FALSE"]]
  group_extended <- geometric_mean_by_group(temp_df, groups)
  group_extended <- as.data.frame(group_extended)
  
  # Add the average row value column
  group_extended$average <- rowMeans(group_extended)
  
  # Add the FC (fold change) column between 0S and 100S
  group_extended$FC <- group_extended$`100S` / (group_extended$`0S`)  # Avoid division by zero
  
  # Print the extended data frame
  group_extended$log2FC <- log2(group_extended$FC)
  group_extended$FC_direction <- ifelse(group_extended$log2FC > 0, "increasing", "decreasing")
  group_extended$abs_log2FC <- abs(group_extended$log2FC)
  print("here")
  for (i in seq_along(list_of_transformations)) {
    # Access the current dataframe
    key_name <- names(list_of_transformations)[i]
    temp_df <- list_of_transformations[[i]]
    # make same row order (miRNAs)
    print(key_name)
    print("3")
    group_df <- geometric_mean_by_group(temp_df, groups)
    results_df <- data.frame(matrix(nrow = nrow(temp_df), ncol = 0))
    rownames(results_df) <- rownames(group_df)
    results_df$monotonicity <- apply(group_df, 1, check_monotonicity)
    print("4")
    # add putative column
    results_df$putative_direction <- group_extended$FC_direction
    # average expression column
    results_df$average_RPM <- old_group_extended$average
    # add log2FC column
    results_df$log2FC <- group_extended$abs_log2FC
    # calculate agreement
    results_df$agreement <- ifelse(results_df$monotonicity == results_df$putative_direction, 1, 0)
    filtered_df <- results_df[results_df$average_RPM >= 50, ]
    # Sort by average_RPM in descending order
    sorted_df <- filtered_df[order(-filtered_df$average_RPM), ]
    # print(sorted_df)
    print(dim(sorted_df))
    top_n_values <- c(1, 5, 10, 20, 50, seq(60, nrow(sorted_df), by = 10))
    
    for (top_n in top_n_values) {
      # Ensure top_n does not exceed the number of rows
      top_n <- min(top_n, nrow(sorted_df))
      
      # Subset the top N rows
      top_subset <- sorted_df[1:top_n, ]
      
      # Calculate the percentage of agreement
      percentage_agreement <- mean(top_subset$agreement) * 100
      
      # Add results to final_results
      parts <- unlist(strsplit(key_name, "_"))
      tgroup = parts[1]  # Return the second keyword
      
      # Add results to final_results
      final_results <- rbind(final_results, data.frame(
        transformation = key_name,
        transformation_group = tgroup,
        top_n = top_n,
        percentage_agreement = percentage_agreement,
        stringsAsFactors = FALSE
      ))
    }
    
    sorted_df <- filtered_df[order(-filtered_df$log2FC), ]
    
    for (top_n in top_n_values) {
      # Ensure top_n does not exceed the number of rows
      top_n <- min(top_n, nrow(sorted_df))
      
      # Subset the top N rows
      top_subset <- sorted_df[1:top_n, ]
      
      # Calculate the percentage of agreement
      percentage_agreement <- mean(top_subset$agreement) * 100
      # Add results to final_results
      parts <- unlist(strsplit(key_name, "_"))
      tgroup = parts[1]  # Return the second keyword
      # Add results to final_results
      final_log2 <- rbind(final_log2, data.frame(
        transformation = key_name,
        transformation_group = tgroup,
        top_n = top_n,
        percentage_agreement = percentage_agreement,
        stringsAsFactors = FALSE
      ))
    }
  }
}

gg <- ggplot(final_results, aes(x = top_n, y = percentage_agreement, 
                                color = transformation, shape = transformation_group ) ) +
  geom_line() +
  geom_point() +
  # scale_x_log10() +  # Optional: Log scale for top_n
  labs(x = "Top n miRNAs (expression)", y = "Percentage Monotonic (%)", 
       title = "") +
  
  theme_minimal()

pdf(paste0(plot_dir, "percentage_monotonic_expres_rlog.pdf"), width = 8, height = 6)
print(gg) 
dev.off()  

gg <- ggplot(final_log2, aes(x = top_n, y = percentage_agreement, 
                             color = transformation, shape = transformation_group ) ) +
  geom_line() +
  geom_point() +
  # scale_x_log10() +  # Optional: Log scale for top_n
  labs(x = "Top n miRNAs (absolute log2FC)", y = "Percentage Monotonic (%)", 
       title = "percentage of monotonic expression of the top DE miRNAs using rlog") +
  theme_minimal()

pdf(paste0(plot_dir, "percentage_monotonic_log2FC_rlog.pdf"), width = 8, height = 6)
print(gg)  
dev.off() 


# MA plot for each transformation (dot per feature)
if (TRUE){
  replicates1 <- c("1S", "20_5S_80_5L_1", "40_4S_60_4L_1", "60_3S_40_3L_1", "80_2S_20_2L_1")
  replicates2 <- c("1S_2", "20_5S_80_5L_2", "40_4S_60_4L_2", "60_3S_40_3L_2", "80_2S_20_2L_2")
  
  cum_result <- data.frame(
    A = numeric(),
    sd_M = numeric(),
    transformation = character(),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(list_of_transformations)) {
    # Access the current dataframe
    key_name <- names(list_of_transformations)[i]
    temp_df <- list_of_transformations[[i]]
    temp_df <- temp_df + 1
    rep1_df <- log2(temp_df[, replicates1])
    rep2_df <- log2(temp_df[, replicates2])
    M_df <- rep1_df - rep2_df
    A_df <- 1/2 * (rep1_df + rep2_df)
    M_vector <- as.vector(as.matrix(M_df))
    A_vector <- as.vector(as.matrix(A_df))
    MA_df <- data.frame(A=A_vector, M=M_vector)
    plot_df <- data.frame(A=A_vector, M=M_vector)
    
    MA_df <- MA_df %>% arrange(A)
    bin_rel_size <- 0.1
    bin_width <- bin_rel_size * diff(range(MA_df$A))
    total_bins <- 1000
    step_size <- 1/total_bins * (max(MA_df$A) - min(MA_df$A))
    # Define the breaks for overlapping bins
    breaks <- seq(min(MA_df$A), max(MA_df$A) - bin_width, by = step_size)
    assign_bins <- function(x) {
      lower_edges <- breaks
      upper_edges <- breaks + bin_width
      which(x >= lower_edges & x < upper_edges)
    }
    MA_df$Bins <- lapply(MA_df$A, assign_bins)
    # Explode the dataframe by bins
    binned_data <- tidyr::unnest(MA_df, cols = Bins)
    # Calculate mean A and mean M for each bin
    result <- binned_data %>%
      group_by(Bins) %>%
      summarise(
        A_mean = mean(A),
        M_mean = stats::sd(M),
        .groups = 'drop'
      )
    # Assign the mean A values as bin centers
    ####
    result$A_mean <- breaks[result$Bins] + bin_width / 2 ######### ASK MATT
    ####
    colnames(result) <- c("bin", "A", "sd_M")
    result$transformation <- key_name
    cum_result <- rbind(cum_result, result)
    
  }
  
  cum_result$transformation <- factor(cum_result$transformation, levels = names(list_of_transformations))
  line_types <- c("solid", "dashed")
  
  line_type_mapping <- rep(line_types, length.out = length(names(list_of_transformations)))
  names(line_type_mapping) <- names(list_of_transformations)
  
}

# TODO add levels with order

gg <- ggplot(cum_result, aes(x = A , y = sd_M, color = transformation, linetype = transformation)) +
  # gg <- ggplot(plot_df, aes(x = A, y = M)) +
  geom_line() +
  scale_linetype_manual(values = line_type_mapping) +
  theme_minimal() +  # Applies a clean theme
  labs(
    title = paste0("moving sd of M|A"),
    x = "A",
    y = "sd(M|A)"
  )
pdf(paste0(plot_dir, "sd_mean_plot.pdf"), width = 8, height = 6)
print(gg)  
dev.off()  


##
### O/E values  
##
##

RPM_threshold <- 50

if (TRUE){
  # for loop that iterates all transformations
  # calculate average value in pure samples
  # remove anything below RPM threshold (get list of miRNAs over threshold to subsect)
  over50_miRNAs <- rownames(old_group_extended[old_group_extended$average > RPM_threshold, ])
  
  group_lookup <- list()
  for (key in names(groups)) {
    for (value in groups[[key]]) {
      group_lookup[[value]] <- key
    }
  }
  
  OE_ratio_list <- list()
  
  for (i in seq_along(list_of_transformations)) {
    # Access the current dataframe
    key_name <- names(list_of_transformations)[i]
    temp_df <- list_of_transformations[[i]]
    group_df <- geometric_mean_by_group(temp_df, groups)
    temp_df <- temp_df[rownames(temp_df) %in% over50_miRNAs, ]
    group_df <- group_df[rownames(group_df) %in% over50_miRNAs, ]
    recalculated_df <- data.frame(matrix(nrow = nrow(temp_df), ncol = ncol(temp_df)))
    
    # print("")
    # print(key_name)
    # print(dim(temp_df))
    # print(dim(group_df))
    # print(dim(recalculated_df))
    # # print(head(temp_df))
    # # print(head(group_df))
    # # print(head(recalculated_df))
    # print("")
    
    colnames(recalculated_df) <- colnames(temp_df)
    columns_to_drop <- c()
    for (i in seq_along(colnames(temp_df))) {
      # Get the group name and extract the percentages
      col_name <- colnames(temp_df)[i]
      # get group name from group_assignment
      if (col_name %in% names(group_lookup)){
        group_name <- group_lookup[[col_name]]
        # print(col_name)
        # print(group_name)
        percentage <- as.numeric(gsub("S", "", group_name)) / 100
        # print(percentage)
        remaining_percentage <- 1 - percentage
        # Calculate the weighted value for each row using the columns 100S and 0S in summary_df
        recalculated_values <- (percentage * group_df[["100S"]]) + (remaining_percentage * group_df[["0S"]])
        # Assign the recalculated values to the corresponding column in recalculated_df
        recalculated_df[[col_name]] <- recalculated_values
      }else{
        # store columns not in list
        columns_to_drop <- c(columns_to_drop, col_name)
      }
      
    }
    # remove columns_to_drop from both dfs
    temp_df <- temp_df[, !(names(temp_df) %in% columns_to_drop)]
    recalculated_df <- recalculated_df[, !(names(recalculated_df) %in% columns_to_drop)]
    # print(columns_to_drop)
    # print(colnames(temp_df))
    # print(colnames(recalculated_df))
    print("")
    print(key_name)
    print(head(temp_df))
    print(head(recalculated_df))
    print("")
    OE_ratio <- temp_df/recalculated_df
    OE_ratio_list[[key_name]] <- OE_ratio
  }
  
  plot_df <- data.frame()
  
  for (i in seq_along(list_of_transformations)){
    temp_OE_ratio <- OE_ratio_list[[i]]
    key_name <- names(list_of_transformations)[i]
    print(key_name)
    print(head(temp_OE_ratio))
    temp_OE_ratio$miRNA <- rownames(temp_OE_ratio)
    long_ratio_data <- melt(temp_OE_ratio, id.vars = "miRNA", 
                            variable.name = "Sample", value.name = "ObservedExpectedRatio")
    long_ratio_data$transformation <- key_name
    # print(head(long_ratio_data))
    plot_df <- rbind(plot_df, long_ratio_data)
  }
  
}

gg <- ggplot(plot_df, aes(x = transformation, y = ObservedExpectedRatio)) +
  geom_boxplot() +  # Boxplot without displaying outliers
  # geom_jitter() +  
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue") +
  labs(x = "transformtion", y = "Observed/Expected Ratio") +
  ggtitle("Observed/Expected Ratio") +
  # scale_y_continuous(breaks = seq(0, max_y, by = 1)) +  # Add y-axis ticks at every unit
  # scale_y_continuous() +  # Add y-axis ticks at every unit
  scale_y_log10()+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(plot_dir, "log_OE.pdf"), plot = gg, width = 20, height = 11.25, device = "pdf")

# O/E values predicted per sample using source mouse data
# for loop that iterates all transformations
# calculate average value in pure samples

if (TRUE){
  OE_ratio_list <- list()
  relevant_columns <- colnames(recalculated_df)
  group_prefixes <- c("20", "40", "60", "80")
  percentages <- c(20, 40, 60, 80) / 100  # Convert to proportions
  remaining_percentages <- 1 - percentages  # Complementary percentages
  
  
  for (i in seq_along(list_of_transformations)) {
    # Access the current dataframe
    key_name <- names(list_of_transformations)[i]
    temp_df <- list_of_transformations[[i]]
    temp_df <- temp_df[rownames(temp_df) %in% over50_miRNAs, ]
    temp_df <- temp_df[, relevant_columns]
    
    mouse_recalculated_df <- data.frame(matrix(nrow = nrow(temp_df), ncol = 0))
    rownames(mouse_recalculated_df) <- rownames(temp_df)  # Set row names to match filtered_df
    # recalculated_df <- data.frame(matrix(nrow = nrow(temp_df), ncol = ncol(temp_df)))
    # colnames(recalculated_df) <- colnames(temp_df)
    print("step 1")
    prefix_columns <- c()
    for (i in seq_along(group_prefixes)) {
      prefix <- group_prefixes[i]
      percentage <- percentages[i]
      remaining_percentage <- remaining_percentages[i]
      print("step 2")
      # Select columns that start with the current prefix
      selected_columns <- grep(paste0("^", prefix), colnames(temp_df), value = TRUE)
      prefix_columns <- c(prefix_columns, selected_columns)
      # Loop through the selected columns and extract the corresponding mouse identifiers
      for (col in selected_columns) {
        # Extract the mouse identifier, which is located between the prefix and the "S" or "L" part
        mouse_id <- sub(paste0("^", prefix, "_(\\d+)[SL]_.*"), "\\1", col)  # Extracts the mouse number
        
        # Define the corresponding "S" and "L" columns using the extracted mouse ID
        S_column <- paste0(mouse_id, "S")
        L_column <- paste0(mouse_id, "L")
        print("step 3")
        # Check if both corresponding columns (e.g., 1S and 1L) exist in filtered_df
        if (S_column %in% colnames(temp_df) && L_column %in% colnames(temp_df)) {
          # Calculate the weighted value for each row
          print("step 4")
          recalculated_values <- (percentage * temp_df[[S_column]]) + (remaining_percentage * temp_df[[L_column]])
          print("step 5")
          # Add the recalculated values as a new column in mouse_recalculated_df
          mouse_recalculated_df[[col]] <- recalculated_values
        }
      }
    }
    temp_df <- temp_df[, prefix_columns]
    print(key_name)
    print(dim(temp_df))
    print(dim(mouse_recalculated_df))
    OE_ratio <- temp_df/mouse_recalculated_df
    OE_ratio_list[[key_name]] <- OE_ratio
  }
  
  # generate plot_df
  
  plot_df <- data.frame()
  
  for (i in seq_along(list_of_transformations)){
    temp_OE_ratio <- OE_ratio_list[[i]]
    key_name <- names(list_of_transformations)[i]
    temp_OE_ratio$miRNA <- rownames(temp_OE_ratio)
    long_ratio_data <- melt(temp_OE_ratio, id.vars = "miRNA", 
                            variable.name = "Sample", value.name = "ObservedExpectedRatio")
    long_ratio_data$transformation <- key_name
    # print(head(long_ratio_data))
    plot_df <- rbind(plot_df, long_ratio_data)
  }
  
}

gg <- ggplot(plot_df, aes(x = transformation, y = ObservedExpectedRatio)) +
  geom_boxplot() +  # Boxplot without displaying outliers
  # geom_jitter() +  
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue") +
  labs(x = "transformtion", y = "Observed/Expected Ratio") +
  ggtitle("Observed/Expected Ratio") +
  # scale_y_continuous(breaks = seq(0, max_y, by = 1)) +  # Add y-axis ticks at every unit
  # scale_y_continuous() +  # Add y-axis ticks at every unit
  scale_y_log10()+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(plot_dir, "log_mouse_OE.pdf"), plot = gg, width = 20, height = 11.25, device = "pdf")


### OLD STUFF


for (i in seq_along(group_prefixes)) {
  prefix <- group_prefixes[i]
  percentage <- percentages[i]
  remaining_percentage <- remaining_percentages[i]
  
  # Select columns that start with the current prefix
  selected_columns <- grep(paste0("^", prefix), colnames(recalculated_df), value = TRUE)
  
  # Loop through the selected columns and extract the corresponding mouse identifiers
  for (col in selected_columns) {
    # Extract the mouse identifier, which is located between the prefix and the "S" or "L" part
    mouse_id <- sub(paste0("^", prefix, "_(\\d+)[SL]_.*"), "\\1", col)  # Extracts the mouse number
    
    # Define the corresponding "S" and "L" columns using the extracted mouse ID
    S_column <- paste0(mouse_id, "S")
    L_column <- paste0(mouse_id, "L")
    
    # Check if both corresponding columns (e.g., 1S and 1L) exist in filtered_df
    if (S_column %in% colnames(temp_df) && L_column %in% colnames(temp_df)) {
      # Calculate the weighted value for each row
      recalculated_values <- (percentage * filtered_df[[S_column]]) + (remaining_percentage * filtered_df[[L_column]])
      
      # Add the recalculated values as a new column in mouse_recalculated_df
      mouse_recalculated_df[[col]] <- recalculated_values
    }
  }
}

if (FALSE){
  # O/E ratio
  OE_ratio <- filtered_df/recalculated_df
  miRNAs_over_50 <- names(average_RPM[average_RPM > 50])
  temp_OE_ratio <- OE_ratio
  temp_OE_ratio$miRNA <- rownames(temp_OE_ratio)
  temp_OE_ratio <- subset(temp_OE_ratio, miRNA %in% miRNAs_over_50)
  long_ratio_data <- melt(temp_OE_ratio, id.vars = "miRNA", variable.name = "Sample", value.name = "ObservedExpectedRatio")
  long_ratio_data$miRNA <- factor(long_ratio_data$miRNA, levels = miRNAs_over_50[order(average_RPM[miRNAs_over_50])])
  
  max_y <- ceiling(max(long_ratio_data$ObservedExpectedRatio))
  
  # Create the boxplot with custom hover text for individual points
  gg <- ggplot(long_ratio_data, aes(x = miRNA, y = ObservedExpectedRatio)) +
    geom_boxplot(outlier.shape = NA) +  # Boxplot without displaying outliers
    geom_jitter(aes(text = paste("miRNA:", miRNA, 
                                 "<br>Sample:", Sample, 
                                 "<br>Observed/Expected Ratio:", round(ObservedExpectedRatio, 2))), 
                width = 0.2, size = 0.5, color = "blue") +  # Smaller jittered points with hover text
    labs(x = "miRNA", y = "Observed/Expected Ratio") +
    ggtitle("Observed/Expected Ratio per miRNA (Filtered by RPM > 50)") +
    scale_y_continuous(breaks = seq(0, max_y, by = 1)) +  # Add y-axis ticks at every unit
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Create the boxplot with a log-scaled y-axis
  
  group_data_with_status <- as.data.frame(group_data_with_status)
  long_ratio_data$expected_direction <- group_data_with_status$monotonicity[match(long_ratio_data$miRNA, rownames(group_data_with_status))]
  
  
  gg <- ggplot(long_ratio_data, aes(x = miRNA, y = ObservedExpectedRatio, color = expected_direction)) +
    geom_boxplot(outlier.shape = NA) +  # Boxplot without displaying outliers
    geom_jitter(aes(text = paste("miRNA:", miRNA, 
                                 "<br>Sample:", Sample, 
                                 "<br>Observed/Expected Ratio:", round(ObservedExpectedRatio, 2))), 
                width = 0.2, size = 0.5, color = "blue") +  # Smaller jittered points with hover text
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +  # Add thick dashed line at y = 1
    labs(x = "miRNA", y = "Observed/Expected Ratio (log scale)") +
    ggtitle("Observed/Expected Ratio per miRNA (Filtered by RPM > 50)") +
    scale_y_log10(breaks = scales::log_breaks(base = 10)) +  # Apply log10 scale with base-10 breaks
    scale_color_manual(values = c("increasing" = "green", "decreasing" = "orange", "neither" = "black")) +  # Custom colors
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("/Users/ernesto/PycharmProjects/miRNA_reference/meetings/meetings_plots/observed_expected_ratio_plot.pdf", plot = gg, width = 20, height = 11.25, device = "pdf")
  
  # Convert to plotly with custom hover text
  ggplotly(gg, tooltip = "text")
  
  # boxplot each sample
  long_ratio_data$SampleGroup <- sapply(long_ratio_data$Sample, function(sample) {
    group_name <- names(groups)[sapply(groups, function(g) sample %in% g)]
    if (length(group_name) > 0) return(group_name) else return(NA)
  })
  
  gg <- ggplot(long_ratio_data, aes(x = Sample, y = ObservedExpectedRatio)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # Boxplot with reduced transparency
    geom_jitter(aes(text = paste("Sample:", Sample, 
                                 "<br>miRNA:", miRNA, 
                                 "<br>Observed/Expected Ratio:", round(ObservedExpectedRatio, 2))), 
                width = 0.2, size = 0.5, color = "blue") +  # Jittered points with custom hover info
    labs(x = "Sample", y = "Observed/Expected Ratio (Log Scale)") +
    ggtitle("Observed/Expected Ratio per Sample, Grouped by Sample Group") +
    scale_y_log10() +  # Set y-axis to logarithmic scale
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ SampleGroup, scales = "free_x")  # Facet by sample group
  
  ggsave("/Users/ernesto/PycharmProjects/miRNA_reference/meetings/meetings_plots/observed_expected_ratio_per_sample.pdf", plot = gg, width = 10, height = 5.65, device = "pdf")
  
  ggplotly(gg, tooltip = "text")
  
  
  gg <- ggplot(long_ratio_data, aes(x = Sample, y = ObservedExpectedRatio)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # Boxplot with reduced transparency
    geom_jitter(aes(text = paste("Sample:", Sample, 
                                 "<br>miRNA:", miRNA, 
                                 "<br>Observed/Expected Ratio:", round(ObservedExpectedRatio, 2))), 
                width = 0.2, size = 0.5, color = "blue") +  # Jittered points with custom hover info
    labs(x = "Sample", y = "Observed/Expected Ratio") +
    ggtitle("Observed/Expected Ratio per Sample") +
    scale_y_continuous(breaks = seq(0, max_y, by = 1)) +  # Add y-axis ticks at every unit
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

}

if (FALSE){
  # Get the relevant columns for the groups 0S to 100S
  columns_order <- c("0S", "20S", "40S", "60S", "80S", "100S")
  group_data <- group_avg_RPM[, columns_order]
  # monotonicity_status <- apply(group_data, 1, check_monotonicity)
  monotonicity_status <- apply(group_data, 1, check_monotonicity_lax)
  
  # Add the monotonicity status as a new column in the matrix
  # group_data <- cbind(group_data, monotonicity = monotonicity_status)
  group_data <- data.frame(group_data, monotonicity = monotonicity_status)
}

quit()

#### TO ADD

# constant miRNAs?
# iteratively remove miRNAs (especially high) and check if this norm improves


library(tidyverse)
library(reshape2)
library(ggplot2)

## decide ground truth (spleen vs liver adj p-val < 0.01)
# load edgeR and DESeq liver vs spleen matrix

DESeq_matrix <- read.csv("~/PycharmProjects/miRNA_reference/pairwise_DE_results/DESeq/pure_liver_VS_pure_spleen/DESeq.tsv", row.names=1)
edgeR_matrix <- read.delim("~/PycharmProjects/miRNA_reference/pairwise_DE_results/edgeR/pure_liver_VS_pure_spleen/edgeR.tsv", row.names=1)


# filter out everything below 50 avg RPM





### expected increase (strong evidence)
# logFC DESeq and edgeR > 1.5
DESeq_DE <- rownames(DESeq_matrix[!is.na(DESeq_matrix$log2FoldChange) & DESeq_matrix$log2FoldChange > 1, ])
edgeR_DE <- rownames(edgeR_matrix[!is.na(edgeR_matrix$logFC) & edgeR_matrix$logFC > 1, ])
strong_increase <- intersect(DESeq_DE, edgeR_DE)

# expected decrease (strong evidence)
# logFC DESeq and edgeR < -1.5
DESeq_DE <- rownames(DESeq_matrix[!is.na(DESeq_matrix$log2FoldChange) & DESeq_matrix$log2FoldChange < -1, ])
edgeR_DE <- rownames(edgeR_matrix[!is.na(edgeR_matrix$logFC) & edgeR_matrix$logFC < -1, ])
strong_decrease <- intersect(DESeq_DE, edgeR_DE)

# expected increase (weaker evidence)
# logFC DESeq or edge > 1.5 and the other > 0
DES_0 <- rownames(DESeq_matrix[!is.na(DESeq_matrix$log2FoldChange) & DESeq_matrix$log2FoldChange > 0, ])
edg_0 <- rownames(edgeR_matrix[!is.na(edgeR_matrix$logFC) & edgeR_matrix$logFC > 0, ])
DES_1 <- intersect(rownames(DESeq_matrix[!is.na(DESeq_matrix$log2FoldChange) & DESeq_matrix$log2FoldChange > 1, ]), edg_0)
edge_1 <- intersect(rownames(edgeR_matrix[!is.na(edgeR_matrix$logFC) & edgeR_matrix$logFC > 1, ]), DES_0)

weaker_increase <- union(DES_1, edge_1)

# expected decrease (weaker evidence)
# logFC DESeq or edge < -1.5 and the other < 0
DES_0 <- rownames(DESeq_matrix[!is.na(DESeq_matrix$log2FoldChange) & DESeq_matrix$log2FoldChange < 0, ])
edg_0 <- rownames(edgeR_matrix[!is.na(edgeR_matrix$logFC) & edgeR_matrix$logFC < 0, ])
DES_1 <- intersect(rownames(DESeq_matrix[!is.na(DESeq_matrix$log2FoldChange) & DESeq_matrix$log2FoldChange < 1, ]), edg_0)
edge_1 <- intersect(rownames(edgeR_matrix[!is.na(edgeR_matrix$logFC) & edgeR_matrix$logFC < 1, ]), DES_0)

weaker_decrease <- union(DES_1, edge_1)

# make matrix per group (geometric mean)
# load miRNA_data
input_df <- read.delim("/Users/ernesto/PycharmProjects/miRNA_reference/de/mature_sense_minExpr0_RCadj.mat", 
                       check.names = FALSE, row.names = 1)

colnames(input_df) <- gsub("\\_R1\\|grp1", "", colnames(input_df))



# define which samples belong to each group
if (TRUE){
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
  # columns <- grep("^[1-5]L", colnames(input_df))
  # columns <- columns[-5]
  groups[["pure_liver"]] <- colnames(input_df)[columns]
  
  columns <- grep("^[1-5]S(?!.*_2$)", colnames(input_df), perl = TRUE)
  groups[["pure_spleen"]] <- colnames(input_df)[columns]
}
group_list <- names(groups)

# generate a matrix per group
if (TRUE){
  # Assuming input_df is your miRNA raw count matrix with columns as sample names
  
  # Initialize an empty list to store RPM data per group
  group_RPMs <- list()
  
  # Function to normalize to RPM
  normalize_to_RPM <- function(counts_df) {
    total_counts <- colSums(counts_df)
    RPM <- sweep(counts_df, 2, total_counts, FUN = "/") * 1e6
    return(RPM)
  }
  
  # Add 1 to raw counts to avoid log(0)
  input_df <- input_df + 1
  
  # Normalize each group to RPM
  for (group in group_list) {
    # Get the samples corresponding to this group
    group_samples <- groups[[group]]
    
    # Subset the input_df to only the columns corresponding to this group
    group_df <- input_df[, group_samples]
    
    # Normalize to RPM
    group_RPM <- normalize_to_RPM(group_df)
    
    # Store the RPM values for this group
    group_RPMs[[group]] <- group_RPM
  }
  
  # Initialize an empty matrix to store the average RPMs
  group_avg_RPM <- matrix(NA, nrow = nrow(input_df), ncol = length(group_list))
  rownames(group_avg_RPM) <- rownames(input_df)
  colnames(group_avg_RPM) <- group_list
  
  # Calculate the geometric mean for each group and store in the matrix
  for (group in group_list) {
    # Extract the RPM data for this group
    group_RPM <- group_RPMs[[group]]
    
    # Calculate the geometric mean for each miRNA (log space, then exponentiate)
    log_RPM <- log2(group_RPM)
    log_means <- rowMeans(log_RPM)
    group_avg_RPM[, group] <- log_means
  }
  
  # Exponentiate to get the geometric means back in the original scale (if needed)
  group_avg_RPM <- 2 ^ group_avg_RPM
  
}

colnames(group_avg_RPM) [5:6] <- c("0S", "100S")

## monotonic function

check_monotonicity <- function(x) {
  if (all(diff(x) > 0)) {
    return("increasing")
  } else if (all(diff(x) < 0)) {
    return("decreasing")
  } else {
    return("neither")
  }
}

check_monotonicity_lax <- function(x) {
  if (all(diff(x) >= 0)) {
    return("increasing")
  } else if (all(diff(x) <= 0)) {
    return("decreasing")
  } else {
    return("neither")
  }
}

## check monotonic trend

if (TRUE){
  # Get the relevant columns for the groups 0S to 100S
  columns_order <- c("0S", "20S", "40S", "60S", "80S", "100S")
  group_data <- group_avg_RPM[, columns_order]
  # monotonicity_status <- apply(group_data, 1, check_monotonicity)
  monotonicity_status <- apply(group_data, 1, check_monotonicity_lax)
  
  # Add the monotonicity status as a new column in the matrix
  # group_data <- cbind(group_data, monotonicity = monotonicity_status)
  group_data <- data.frame(group_data, monotonicity = monotonicity_status)
}

## calculate percentages of correct recall in trend
if (TRUE){
  increase_status <- monotonicity_status[rownames(group_data) %in% strong_increase]
  decrease_status <- monotonicity_status[rownames(group_data) %in% strong_decrease]
  
  # 2. Calculate how many are correctly classified
  correct_increase <- sum(increase_status == "increasing")
  correct_decrease <- sum(decrease_status == "decreasing")
  
  # 3. Calculate the total number of miRNAs in each set
  total_increase <- length(strong_increase)
  total_decrease <- length(strong_decrease)
  
  # 4. Calculate percentages
  percent_increase_correct <- (correct_increase / total_increase) * 100
  percent_decrease_correct <- (correct_decrease / total_decrease) * 100
  print(percent_increase_correct)
  print(percent_decrease_correct)
}

## calculate percentages of correct recall in trend (weaker)
if (TRUE){
  increase_status <- monotonicity_status[rownames(group_data) %in% weaker_increase]
  decrease_status <- monotonicity_status[rownames(group_data) %in% weaker_decrease]
  
  # 2. Calculate how many are correctly classified
  correct_increase <- sum(increase_status == "increasing")
  correct_decrease <- sum(decrease_status == "decreasing")
  
  # 3. Calculate the total number of miRNAs in each set
  total_increase <- length(weaker_increase)
  total_decrease <- length(weaker_decrease)
  
  # 4. Calculate percentages
  percent_increase_correct <- (correct_increase / total_increase) * 100
  percent_decrease_correct <- (correct_decrease / total_decrease) * 100
}

### sequential dot plots

# Step 1: Create RPM matrix for all samples
# Normalize input_df to RPM, adding 1 to avoid log(0) issues
input_df <- input_df + 1
total_counts <- colSums(input_df)
all_samples_RPM <- sweep(input_df, 2, total_counts, FUN = "/") * 1e6

# Step 2: Filter out samples that are not in any group
# Flatten the list of groups to get all sample names within groups
grouped_samples <- unlist(groups)
filtered_RPM <- all_samples_RPM[, colnames(all_samples_RPM) %in% grouped_samples]

# Step 3: Melt the filtered RPM matrix for ggplot compatibility
# First, add rownames (miRNA) as a new column to preserve miRNA identifiers
filtered_RPM$miRNA <- rownames(filtered_RPM)
melted_RPM <- melt(filtered_RPM, id.vars = "miRNA", variable.name = "Sample", value.name = "RPM")

melted_RPM$Group <- sapply(melted_RPM$Sample, function(sample) {
  group <- names(groups)[sapply(groups, function(g) sample %in% g)]
  if (length(group) > 0) return(group) else return(NA)
})

plot_miRNA_dots <- function(miRNA_input, melted_RPM, group_order = c("0S", "20S", "40S", "60S", "80S", "100S")) {
  # Subset to the specified miRNA only
  miRNA_data <- subset(melted_RPM, miRNA == miRNA_input)
  
  # Ensure Group is a factor for ordering in ggplot
  miRNA_data$Group <- factor(miRNA_data$Group, levels = group_order)
  
  # Create the dot plot
  # ggplot(miRNA_data, aes(x = Group, y = RPM)) +
  ggplot(miRNA_data, aes(x = Group, y = RPM, text = paste("Sample:", Sample, "<br>RPM:", round(RPM, 2)))) +
    geom_point(aes(color = Group), position = position_jitter(width = 0.15, height = 0), size = 2) +
    scale_color_brewer(palette = "Dark2") +
    ggtitle(paste("RPM values for miRNA:", miRNA_input)) +
    xlab("Group") +
    ylab("RPM") +
    theme_minimal() +
    theme(legend.position = "none")
}

# plot in theory decreasing but not actually
# theoretical decrease strong_decrease
# labels different than decreasing 
failed_dec <- intersect(rownames(subset(group_data, group_data$monotonicity != "decreasing")), strong_decrease)

plot_miRNA_dots("Mmu-Mir-3090_3p", melted_RPM)
plot_miRNA_dots("Mmu-Mir-506-P18_3p", melted_RPM)
plot_miRNA_dots("Mmu-Mir-3073-as_3p*", melted_RPM)
plot_miRNA_dots("Mmu-Mir-3073_3p", melted_RPM)
plot_miRNA_dots("Mmu-Mir-3105_5p", melted_RPM)
plot_miRNA_dots("Mmu-Mir-192-P2_3p*", melted_RPM)







plot_miRNA_dots("Mmu-Mir-24-P3_5p*", melted_RPM)

plot_miRNA_dots("Mmu-Mir-1948_3p", melted_RPM)
plot_miRNA_dots("Mmu-Mir-459_3p", melted_RPM)
plot_miRNA_dots("Mmu-Mir-194-P2_5p", melted_RPM)
plot_miRNA_dots("Mmu-Mir-122_5p", melted_RPM)
plot_miRNA_dots("Mmu-Mir-192-P2_5p", melted_RPM)
plot_miRNA_dots("Mmu-Mir-148-P1_5p*", melted_RPM)



### FILTER only to miRNAs 50RPM or more
# discard all other miRNAs from all lists

#miRNAs with at least 50 RPM

average_RPM <- rowMeans(all_samples_RPM)
# Get the names of miRNAs with at least 50 RPM on average
miRNAs_above_50_RPM <- names(average_RPM[average_RPM >= 50])


strong_increase <- intersect(strong_increase, miRNAs_above_50_RPM)
strong_decrease <- intersect(strong_decrease, miRNAs_above_50_RPM)

strong_increase_mono <- intersect(strong_increase, rownames(subset(group_data, group_data$monotonicity == "increasing")))
strong_decrease_mono <- intersect(strong_decrease, rownames(subset(group_data, group_data$monotonicity == "decreasing")))

length(strong_decrease_mono)/length(strong_decrease)
length(strong_increase_mono)/length(strong_increase)

# failing cases
failed_dec <- intersect(rownames(subset(group_data, group_data$monotonicity != "decreasing")), strong_decrease)

gg <- plot_miRNA_dots("Mmu-Mir-193-P2b_3p", melted_RPM)

ggplotly(gg, tooltip = "text")

failed_inc <- intersect(rownames(subset(group_data, group_data$monotonicity != "increasing")), strong_increase)
gg <- plot_miRNA_dots("Mmu-Mir-223_3p", melted_RPM)
ggplotly(gg, tooltip = "text")

#####
#####
#####
########### MIRBASE #######
#####
#####
#####

# DE L vs S

input_df <- read.delim("/Users/ernesto/PycharmProjects/miRNA_reference/de_matrix/mirbase_MA/mature_sense_minExpr0_RCadj.mat", 
                       check.names = FALSE, row.names = 1)

colnames(input_df) <- gsub("\\_R1\\|grp1", "", colnames(input_df))

# LOAD DE MATRIX

DESeq_matrix <- read.csv("/Users/ernesto/PycharmProjects/miRNA_reference/DE_results/mirbase/DESeq/pure_liver_VS_pure_spleen/DESeq.tsv", row.names=1)
edgeR_matrix <- read.delim("/Users/ernesto/PycharmProjects/miRNA_reference/DE_results/mirbase/edgeR/pure_liver_VS_pure_spleen/edgeR.tsv", row.names=1)


# filter out everything below 50 avg RPM


### expected increase (strong evidence)
# logFC DESeq and edgeR > 1.5
DESeq_DE <- rownames(DESeq_matrix[!is.na(DESeq_matrix$log2FoldChange) & DESeq_matrix$log2FoldChange > 1, ])
edgeR_DE <- rownames(edgeR_matrix[!is.na(edgeR_matrix$logFC) & edgeR_matrix$logFC > 1, ])
strong_increase <- intersect(DESeq_DE, edgeR_DE)

# expected decrease (strong evidence)
# logFC DESeq and edgeR < -1.5
DESeq_DE <- rownames(DESeq_matrix[!is.na(DESeq_matrix$log2FoldChange) & DESeq_matrix$log2FoldChange < -1, ])
edgeR_DE <- rownames(edgeR_matrix[!is.na(edgeR_matrix$logFC) & edgeR_matrix$logFC < -1, ])
strong_decrease <- intersect(DESeq_DE, edgeR_DE)

# calculate RPM, reduce only to min 50 RPM average

input_df <- input_df + 1
total_counts <- colSums(input_df)
all_samples_RPM <- sweep(input_df, 2, total_counts, FUN = "/") * 1e6
average_RPM <- rowMeans(all_samples_RPM)
# Get the names of miRNAs with at least 50 RPM on average
miRNAs_above_50_RPM <- names(average_RPM[average_RPM >= 50])

# Step 2: Filter out samples that are not in any group
# Flatten the list of groups to get all sample names within groups
grouped_samples <- unlist(groups)
filtered_RPM <- all_samples_RPM[, colnames(all_samples_RPM) %in% grouped_samples]



## calculate monotonic
# make group matrix
if (TRUE){
  # Assuming input_df is your miRNA raw count matrix with columns as sample names
  
  # Initialize an empty list to store RPM data per group
  group_RPMs <- list()
  
  # Function to normalize to RPM
  normalize_to_RPM <- function(counts_df) {
    total_counts <- colSums(counts_df)
    RPM <- sweep(counts_df, 2, total_counts, FUN = "/") * 1e6
    return(RPM)
  }
  
  # Add 1 to raw counts to avoid log(0)
  input_df <- input_df + 1
  
  # Normalize each group to RPM
  for (group in group_list) {
    # Get the samples corresponding to this group
    group_samples <- groups[[group]]
    
    # Subset the input_df to only the columns corresponding to this group
    group_df <- input_df[, group_samples]
    
    # Normalize to RPM
    group_RPM <- normalize_to_RPM(group_df)
    
    # Store the RPM values for this group
    group_RPMs[[group]] <- group_RPM
  }
  
  # Initialize an empty matrix to store the average RPMs
  group_avg_RPM <- matrix(NA, nrow = nrow(input_df), ncol = length(group_list))
  rownames(group_avg_RPM) <- rownames(input_df)
  colnames(group_avg_RPM) <- group_list
  
  # Calculate the geometric mean for each group and store in the matrix
  for (group in group_list) {
    # Extract the RPM data for this group
    group_RPM <- group_RPMs[[group]]
    
    # Calculate the geometric mean for each miRNA (log space, then exponentiate)
    log_RPM <- log2(group_RPM)
    log_means <- rowMeans(log_RPM)
    group_avg_RPM[, group] <- log_means
  }
  
  # Exponentiate to get the geometric means back in the original scale (if needed)
  group_avg_RPM <- 2 ^ group_avg_RPM
  
}

colnames(group_avg_RPM) [5:6] <- c("0S", "100S")
columns_order <- c("0S", "20S", "40S", "60S", "80S", "100S")
group_data <- group_avg_RPM[, columns_order]
monotonicity_status <- apply(group_data, 1, check_monotonicity_lax)

# Add the monotonicity status as a new column in the matrix
# group_data <- cbind(group_data, monotonicity = monotonicity_status)
group_data <- data.frame(group_data, monotonicity = monotonicity_status)
colnames(group_data)[1:6] <- columns_order


# overlap between expected/observed monotonic increase

strong_increase <- intersect(strong_increase, miRNAs_above_50_RPM)
strong_decrease <- intersect(strong_decrease, miRNAs_above_50_RPM)

strong_increase_mono <- intersect(strong_increase, rownames(subset(group_data, group_data$monotonicity == "increasing")))
strong_decrease_mono <- intersect(strong_decrease, rownames(subset(group_data, group_data$monotonicity == "decreasing")))

length(strong_decrease_mono)/length(strong_decrease)
length(strong_increase_mono)/length(strong_increase)


failed_dec <- intersect(rownames(subset(group_data, group_data$monotonicity != "decreasing")), strong_decrease)
failed_inc <- intersect(rownames(subset(group_data, group_data$monotonicity != "increasing")), strong_increase)
# plots

filtered_RPM$miRNA <- rownames(filtered_RPM)
melted_RPM <- melt(filtered_RPM, id.vars = "miRNA", variable.name = "Sample", value.name = "RPM")

melted_RPM$Group <- sapply(melted_RPM$Sample, function(sample) {
  group <- names(groups)[sapply(groups, function(g) sample %in% g)]
  if (length(group) > 0) return(group) else return(NA)
})


gg <- plot_miRNA_dots("Mmu-Mir-193-P2b_3p", melted_RPM)

ggplotly(gg, tooltip = "text")


gg <- plot_miRNA_dots("Mmu-Mir-223_3p", melted_RPM)
ggplotly(gg, tooltip = "text")



## TODO
# try more robust filters

### iterate RPM thresholds
# build function that generates % increase and decrease for each threshold

get_monotonic_perc <- function(increase, decrease, mono_matrix, RPM_average, threshold){
  
  thresholded_miRNAs <-  names(RPM_average[RPM_average >= threshold])
  
  cons_increase <- intersect(increase, thresholded_miRNAs)
  cons_decrease <- intersect(decrease, thresholded_miRNAs)
  
  strong_increase_mono <- intersect(cons_increase, rownames(subset(mono_matrix, mono_matrix$monotonicity == "increasing")))
  strong_decrease_mono <- intersect(cons_decrease, rownames(subset(mono_matrix, mono_matrix$monotonicity == "decreasing")))
  
  perc_inc <- length(strong_increase_mono)/length(cons_increase)
  perc_dec <- length(strong_decrease_mono)/length(cons_decrease)
  
  return(c(perc_inc, perc_dec, length(strong_increase_mono), length(strong_decrease_mono)))
}

# input: strong increase, strong decrease, average RPMs, thresholds

# Define thresholds to iterate over
thresholds <- c(0:5, 10, 50, 100, 500, 1000)

# Initialize an empty data frame to store results
results_df <- data.frame()

# Loop through each threshold, using your get_monotonic_perc function
for (threshold in thresholds) {
  # Calculate the percentages for the current threshold
  result <- get_monotonic_perc(strong_increase, strong_decrease, group_data, average_RPM, threshold )
  
  # Add the results to the data frame in melted format
  results_df <- rbind(results_df, 
                      data.frame(threshold = threshold, direction = "increase", 
                                 percentage = result[1], total = result[3]),
                      data.frame(threshold = threshold, direction = "decrease", 
                                 percentage = result[2], total = result[4]))
}

# plot in bars or lines



### 1 average on each pure group

if (TRUE){
  # make average RPM matrix
  group_data$monotonicity <- NULL
  colnames(group_data) <- c("0S", "20S", "40S", "60S", "80S", "100S")
  group_ass <- groups
  names(group_ass) <- c( "20S", "40S", "60S", "80S", "0S", "100S")
  # make expected values matrix 
  
  filter_df_by_group <- function(group_list, df) {
    # Initialize a vector to store the group names for each column that is kept
    group_vector <- c()
    
    # Initialize a vector to store columns that should be kept
    columns_to_keep <- c()
    
    # Loop over each column in the dataframe
    for (col in colnames(df)) {
      # Find the group the column belongs to
      group_found <- NULL
      for (group_name in names(group_list)) {
        if (col %in% group_list[[group_name]]) {
          group_found <- group_name
          break
        }
      }
      
      # If the column is in the list, add it to the columns to keep
      # and add the group name to the group vector
      if (!is.null(group_found)) {
        columns_to_keep <- c(columns_to_keep, col)
        group_vector <- c(group_vector, group_found)
      }
    }
    
    # Subset the dataframe to only include the columns in columns_to_keep
    filtered_df <- df[, columns_to_keep, drop = FALSE]
    
    # Return the filtered dataframe and the group vector
    list(filtered_df = filtered_df, group_vector = group_vector)
  }
  
  result <- filter_df_by_group(group_ass, all_samples_RPM)
  
  # Extract filtered_df and group_vector from the result list
  filtered_df <- result$filtered_df
  group_vector <- result$group_vector
  
  # Check column names in summary_df to ensure they are correct
  print("Column names in summary_df:")
  print(colnames(group_avg_RPM_df))
  
  # Initialize an empty data frame to store recalculated values
  recalculated_df <- data.frame(matrix(nrow = nrow(filtered_df), ncol = ncol(filtered_df)))
  colnames(recalculated_df) <- colnames(filtered_df)  # Set column names to match filtered_df
  
  # Loop through each column in filtered_df and corresponding group in group_vector
  for (i in seq_along(group_vector)) {
    # Get the group name and extract the percentages
    group_name <- group_vector[i]
    percentage <- as.numeric(gsub("S", "", group_name)) / 100
    remaining_percentage <- 1 - percentage
    
    # Print to inspect each calculated percentage and remaining percentage
    cat("Group:", group_name, "| Percentage:", percentage * 100, "| Remaining:", remaining_percentage * 100, "\n")
    
    # Calculate the weighted value for each row using the columns 100S and 0S in summary_df
    recalculated_values <- (percentage * group_avg_RPM_df[["100S"]]) + (remaining_percentage * group_avg_RPM_df[["0S"]])
    
    # Assign the recalculated values to the corresponding column in recalculated_df
    recalculated_df[[i]] <- recalculated_values
  }
  
  # View the recalculated dataframe to inspect the final output
  print("Recalculated Data Frame:")
  print(head(recalculated_df))
  rownames(recalculated_df) <- rownames(filtered_df)
  mean_expected_df <- recalculated_df
  
}

# O/E ratio
OE_ratio <- filtered_df/recalculated_df
miRNAs_over_50 <- names(average_RPM[average_RPM > 50])
temp_OE_ratio <- OE_ratio
temp_OE_ratio$miRNA <- rownames(temp_OE_ratio)
temp_OE_ratio <- subset(temp_OE_ratio, miRNA %in% miRNAs_over_50)
long_ratio_data <- melt(temp_OE_ratio, id.vars = "miRNA", variable.name = "Sample", value.name = "ObservedExpectedRatio")
long_ratio_data$miRNA <- factor(long_ratio_data$miRNA, levels = miRNAs_over_50[order(average_RPM[miRNAs_over_50])])

max_y <- ceiling(max(long_ratio_data$ObservedExpectedRatio))

# Create the boxplot with custom hover text for individual points
gg <- ggplot(long_ratio_data, aes(x = miRNA, y = ObservedExpectedRatio)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without displaying outliers
  geom_jitter(aes(text = paste("miRNA:", miRNA, 
                               "<br>Sample:", Sample, 
                               "<br>Observed/Expected Ratio:", round(ObservedExpectedRatio, 2))), 
              width = 0.2, size = 0.5, color = "blue") +  # Smaller jittered points with hover text
  labs(x = "miRNA", y = "Observed/Expected Ratio") +
  ggtitle("Observed/Expected Ratio per miRNA (Filtered by RPM > 50)") +
  scale_y_continuous(breaks = seq(0, max_y, by = 1)) +  # Add y-axis ticks at every unit
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create the boxplot with a log-scaled y-axis

group_data_with_status <- as.data.frame(group_data_with_status)
long_ratio_data$expected_direction <- group_data_with_status$monotonicity[match(long_ratio_data$miRNA, rownames(group_data_with_status))]

# boxplot each sample
long_ratio_data$SampleGroup <- sapply(long_ratio_data$Sample, function(sample) {
  group_name <- names(groups)[sapply(groups, function(g) sample %in% g)]
  if (length(group_name) > 0) return(group_name) else return(NA)
})


### 2 only corresponding mouse
# calculate expected RPM

group_prefixes <- c("20", "40", "60", "80")
percentages <- c(20, 40, 60, 80) / 100  # Convert to proportions
remaining_percentages <- 1 - percentages  # Complementary percentages

# Initialize an empty data frame to store recalculated values for each mouse/sample
mouse_recalculated_df <- data.frame(matrix(nrow = nrow(filtered_df), ncol = 0))
rownames(mouse_recalculated_df) <- rownames(filtered_df)  # Set row names to match filtered_df

for (i in seq_along(group_prefixes)) {
  prefix <- group_prefixes[i]
  percentage <- percentages[i]
  remaining_percentage <- remaining_percentages[i]
  
  # Select columns that start with the current prefix
  selected_columns <- grep(paste0("^", prefix), colnames(filtered_df), value = TRUE)
  
  # Loop through the selected columns and extract the corresponding mouse identifiers
  for (col in selected_columns) {
    # Extract the mouse identifier, which is located between the prefix and the "S" or "L" part
    mouse_id <- sub(paste0("^", prefix, "_(\\d+)[SL]_.*"), "\\1", col)  # Extracts the mouse number
    
    # Define the corresponding "S" and "L" columns using the extracted mouse ID
    S_column <- paste0(mouse_id, "S")
    L_column <- paste0(mouse_id, "L")
    
    # Check if both corresponding columns (e.g., 1S and 1L) exist in filtered_df
    if (S_column %in% colnames(filtered_df) && L_column %in% colnames(filtered_df)) {
      # Calculate the weighted value for each row
      recalculated_values <- (percentage * filtered_df[[S_column]]) + (remaining_percentage * filtered_df[[L_column]])
      
      # Add the recalculated values as a new column in mouse_recalculated_df
      mouse_recalculated_df[[col]] <- recalculated_values
    }
  }
}

# O/E ratio
mouse_filtered_df <- filtered_df[rownames(mouse_recalculated_df), colnames(mouse_recalculated_df)]
OE_ratio <- mouse_filtered_df/mouse_recalculated_df
miRNAs_over_50 <- names(average_RPM[average_RPM > 50])
temp_OE_ratio <- OE_ratio
temp_OE_ratio$miRNA <- rownames(temp_OE_ratio)
temp_OE_ratio <- subset(temp_OE_ratio, miRNA %in% miRNAs_over_50)
long_ratio_data <- melt(temp_OE_ratio, id.vars = "miRNA", variable.name = "Sample", value.name = "ObservedExpectedRatio")
long_ratio_data$miRNA <- factor(long_ratio_data$miRNA, levels = miRNAs_over_50[order(average_RPM[miRNAs_over_50])])

long_ratio_data$SampleGroup <- sapply(long_ratio_data$Sample, function(sample) {
  group_name <- names(groups)[sapply(groups, function(g) sample %in% g)]
  if (length(group_name) > 0) return(group_name) else return(NA)
})

long_ratio_data$expected_direction <- group_data_with_status$monotonicity[match(long_ratio_data$miRNA, rownames(group_data_with_status))]


# rlog, vst, TMM


# expected matrix (group and per mouse)
# 


## assessment
# the sum of absolute errors
# O/E ratio
# log O/E ratios



# start with "exclusive" miRNAs
# plot bunch
# check by mouse
# median




