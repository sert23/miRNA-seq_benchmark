library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(htmlwidgets)
library(plotly)
library(tidyr)
library(ggnewscale)
library(compositions)
library(zCompositions)


source("R/normalization.R")
plot_dir <- "./plots/"
dir.create(plot_dir, recursive = TRUE)
minRPM <- 50

# name map
name_map <- c(
  "RPM_total"                = "RPM (total)",
  "RPM_lib"                = "RPM",
  "TMM"                      = "edgeR - TMM",
  "rlog_parametric_blind_TRUE"  = "DESeq2 - rlog (parametric)",
  "rlog_mean_blind_FALSE"       = "DESeq2 - rlog (mean)",
  "vst_local_blind_TRUE"        = "DESeq2 - VST (local)",
  "RC"                       = "Raw read count",
  "MAP" = "MAP"
)

# make monotonicity example figure first
if (TRUE){
  set.seed(42)  # make the jitter reproducible
  
  p_df <- data.frame(
    Sample = factor(
      c("0S (pure liver)", "20S", "40S", "60S", "80S", "100S (pure spleen)"),
      levels = c("0S (pure liver)", "20S", "40S", "60S", "80S", "100S (pure spleen)")
    ),
    Counts = c(5, 18, 35, 55, 75, 95)
  )
  
  # compute the decreasing values, with jitter on 60S
  p_df2 <- p_df %>%
    mutate(
      Dec = 100 - Counts,
      Dec = ifelse(
        Sample == "60S",
        Dec + sample(-10:10, 1),
        Dec
      )
    )
  
  # choose a palette
  inc_col <- "#1f78b4"  # blue = increasing
  dec_col <- "#e31a1c"  # red  = decreasing
  
  gg <- ggplot(p_df2, aes(x = Sample)) +
    # blue = increasing
    geom_col(
      aes(y    = Counts, fill = "Monotonic increase"),
      position     = position_nudge(x = -0.2),
      width        = 0.4
    ) +
    # red = decreasing
    geom_col(
      aes(y    = Dec, fill = "Monotonic decrease"),
      position     = position_nudge(x =  0.2),
      width        = 0.4
    ) +
    scale_fill_manual(
      name   = NULL,
      values = c(
        "Monotonic increase" = inc_col,
        "Monotonic decrease" = dec_col
      )
    ) +
    labs(
      title = "Monotonic trends in the expression of a miRNA",
      x     = NULL,
      y     = "miRNA-X counts"
    ) +
    theme_minimal() +
    theme_minimal(base_size = 14) +        # increase base font size
    theme(
      text                = element_text(size = 14),       # all text
      axis.title          = element_text(size = 14),       # axis titles
      axis.text.x         = element_text(size = 14, angle = 45, hjust = 1),
      axis.text.y         = element_text(size = 14),
      plot.title          = element_text(size = 16, hjust = 0.5),
      legend.title        = element_text(size = 14),
      legend.text         = element_text(size = 14),
      legend.position     = "bottom",
      legend.justification= "center",
      legend.direction    = "horizontal"
    )
  
  pdf(paste0(plot_dir, "example_monotonicity.pdf"), width = 8, height = 6)  
  print(gg)  
  dev.off()  
}

#### create resampled RC matrix
if (TRUE){
  # read original RC matrix
  RC_df <- read.delim("./data/mature_sense_minExpr0_RCadj.mat",
                      check.names = FALSE, row.names = 1)
  # count min amount of total reads
  col_sums <- colSums(RC_df)
  min_am  <- min(col_sums) # 12214826, let's sample 10M
  
  # calculate probabilities
  prob_df <- sweep(RC_df, 
                   MARGIN = 2, 
                   STATS  = colSums(RC_df), 
                   FUN    = "/")
  # sample 10M for each column
  n_feats   <- nrow(prob_df)
  n_draws   <- 1e7  # 10 million
  
  sampled_mat <- sapply(seq_len(ncol(prob_df)), function(i) {
    sample_counts_col(prob_df[[i]])
  })
  
  # Convert to data.frame and restore row/col names
  sampled_df <- as.data.frame(sampled_mat,
                              row.names = rownames(prob_df))
  colnames(sampled_df) <- colnames(prob_df)
  
  RC_df <- sampled_df
  
  
}

# load RC matrix

# RC_df <- read.delim("./data/mature_sense_minExpr0_RCadj.mat",
#                     check.names = FALSE, row.names = 1)
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

# RPM lib (recalculated from resampled 10M RC)
RPMl_df <- RC_df / colSums(RC_df) * 1e6

# loading sRNAbench matrix
# RPMl_df <- read.delim("./data/mature_sense_minExpr0_RCadj_libraryRPM.mat",
#                       check.names = FALSE, row.names = 1)
# colnames(RPMl_df) <- fix_names(colnames(RPMl_df))

list_of_transformations <- list(RC = RC_df, 
                                # logged = logged_df, 
                                RPM_total = RPMt_df,
                                RPM_lib = RPMl_df
                                )


# add Matt's function for MAP here
if (TRUE){
  # 1. Compute library sizes and counts
  lib_sizes <- colSums(RC_df)
  I <- nrow(RC_df)
  alpha <- 1               # uniform prior concentration
  n_j <- lib_sizes
  alpha_j <- rep(alpha, length(n_j))
  
  # 2. Compute MLE proportions
  pi_mle <- sweep(RC_df, 2, lib_sizes, FUN = "/")
  
  # 3. Apply MAP shrinkage column‐wise
  pi_map <- sapply(seq_len(ncol(pi_mle)), function(j) {
    map_estimator(
      pi_mle = pi_mle[, j],
      I      = I,
      alpha_j= alpha_j[j],
      n_j    = n_j[j]
    )
  })
  # restore dimnames
  rownames(pi_map) <- rownames(pi_mle)
  colnames(pi_map) <- colnames(pi_mle)
}
list_of_transformations$MAP <- pi_map

# CLR
if(TRUE){
  # pi_mle proportions
  prop_df <- pi_mle
  prop_df_filt <- cmultRepl(prop_df, method = "CZM", output = "prop")
  clr_df <- clr(acomp(prop_df_filt))
}
# list_of_transformations$CLR <- clr_df

# check for any NAs
any_na <- any(is.na(clr_df))
na_positions <- which(is.na(clr_df), arr.ind = TRUE)

cat("Any NAs?", any_na, "\n")
print(na_positions)

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
if(TRUE){dge <- DGEList(counts = RC_df)

# Calculate normalization factors using TMM
dge <- calcNormFactors(dge, method = "TMM")

# To extract the TMM-normalized values (on a per-million scale)
tmm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
list_of_transformations[["TMM"]] <- as.data.frame(tmm_counts)
}

# define groups
group_obj <- define_groups(RC_df)
groups <- group_obj$group_members
group_list <- group_obj$group_list

# group RPM
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


#######
##### evaluate monotonic trend for each method
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
    print(key_name)
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
    # results_df$agreement <- ifelse(results_df$monotonicity == results_df$putative_direction, 1, 0)
    results_df$agreement <- ifelse(
      results_df$monotonicity %in% c("increasing", "decreasing"),
      1,
      0
    )
    # use RPM to filter so it's equivalent to all transformations
    filtered_df <- results_df[results_df$average_RPM >= minRPM, ]
    # Sort by average_RPM in descending order
    sorted_df <- filtered_df[order(-filtered_df$average), ]
    top_n_values <- c(5, 10, 20, 50, seq(60, nrow(sorted_df), by = 10))
    
    for (top_n in top_n_values) {
      # Ensure top_n does not exceed the number of rows
      top_n <- min(top_n, nrow(sorted_df))
      
      # Subset the top N rows
      top_subset <- sorted_df[1:top_n, ]
      
      # Calculate the percentage of agreement
      # percentage_agreement <- mean(top_subset$agreement) * 100
      # Calculate the percentage of monotonic
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
  
  methods <- unique(final_results$transformation)
  pal_hex <- unname(grDevices::palette.colors())
  pal_hex <- pal_hex[pal_hex != "#000000"]
  cols_use <- pal_hex[seq_along(methods)]
  colour_map <- setNames(cols_use, methods)
  
  
  ## fix levels to mimic order in legend
  
  order_levels <- final_results %>%
    filter(top_n == max(top_n)) %>%               # pick the last x‐value
    arrange(desc(percentage_agreement)) %>%        # sort so highest % is first
    pull(transformation)                           # extract the transformation names
  
  final_results <- final_results %>%
    mutate(transformation = factor(transformation, levels = order_levels))
  
}


# lineplot % monotonic for top n (expression)
gg <- ggplot(final_results, aes(x = top_n, y = percentage_agreement, 
                                color = transformation) ) +
  geom_line() +
  geom_point() +
  # scale_x_log10() +  # Optional: Log scale for top_n
  labs(x = "Top n miRNAs (sorted by expression)", y = "Percentage with monotonic trend (%)", 
       title = "Proportion of top-expressed miRNAs displaying a monotonic trend", 
       color = "normalization method") +
  scale_color_manual(
    name   = "normalization method",
    values = colour_map,
    labels = name_map,
    breaks = order_levels
  ) +
  theme_minimal()

pdf(paste0(plot_dir, "percentage_monotonic_expres_relative.pdf"), width = 8, height = 6)  
print(gg)  
dev.off() 

# lineplot % monotonic for top n (absolute FC)
gg <- ggplot(final_log2, aes(x = top_n, y = percentage_agreement, 
                             color = transformation) ) +
  geom_line() +
  geom_point() +
  # scale_x_log10() +  # Optional: Log scale for top_n
  labs(x = "Top n miRNAs (sorted by  |log2FC|)", y = "Percentage with monotonic trend (%)", 
       title = "Proportion of top DE miRNAs displaying a monotonic trend", 
       color = "normalization method") +
  scale_color_manual(
    name   = "normalization method",
    values = colour_map,
    labels = name_map,
    breaks = order_levels
  ) +
  theme_minimal()

pdf(paste0(plot_dir, "percentage_monotonic_log2FC_relative.pdf"), width = 8, height = 6)
print(gg)  
dev.off()

# heatmaps
if (TRUE){
  rpm_key <- grep("RPM", names(list_of_transformations), value = TRUE)[1]
  rpm_key <- "RPM_lib"
  # rpm_key <- grep("RPM", names(list_of_transformations), value = TRUE)[1]
  
  # Compute geometric mean by group for RPM
  rpm_grp       <- geometric_mean_by_group(list_of_transformations[[rpm_key]], groups)
  average_RPM   <- rowMeans(rpm_grp)
  
  # Get the top 100 miRNAs by that average
  topN <- 100
  top100 <- names(sort(average_RPM, decreasing = TRUE))[1:topN]
  
  # ADDING FC COLUMN TO HEATMAP
  if (FALSE){
    
    # generate RAW_DF here
    RAW_DF <- rpm_grp  
    print(colnames(RAW_DF))
    
    
    # Compute log2FC from two raw columns (placeholders RAW_DF, COL_A, COL_B)
    # Ensure miRNA names align to top100 order
    fc_vec <- log2(
      RAW_DF[top100, "100S", drop = TRUE] /
        RAW_DF[top100, "0S" , drop = TRUE]
    )
    print(fc_vec)
    
    # Build a data frame for the extra column
    fc_df <- data.frame(
      miRNA          = top100,
      transformation = "log2FC",
      log2FC         = fc_vec,
      stringsAsFactors = FALSE
    )
  
  
    
    # 4. Build a long table of monotonicity for those 100
    mono_list <- lapply(names(list_of_transformations), function(key) {
      grp  <- geometric_mean_by_group(list_of_transformations[[key]], groups)
      grp100 <- grp[top100, , drop = FALSE]                 # only top100 rows
      mono <- apply(grp100, 1, check_monotonicity)
      data.frame(
        miRNA          = rownames(grp100),
        transformation = key,
        monotonicity   = mono,
        stringsAsFactors = FALSE
      )
    })
    heat_df <- bind_rows(mono_list)
    
    # 5. Factor levels to fix the x-axis order
    heat_df <- heat_df %>%
      mutate(
        miRNA          = factor(miRNA, levels = top100),
        transformation = factor(transformation, levels = names(colour_map))
      )
    # Make sure factor levels include the new “log2FC” method
    heat_df <- heat_df %>%
      bind_rows(fc_df %>%
                  mutate(
                    miRNA          = factor(miRNA, levels = top100),
                    transformation = factor(transformation,
                                            levels = c(names(colour_map), "log2FC"))
                  ))
    gg <- ggplot() +
      # 1) existing monotonicity tiles
      geom_tile(
        data = filter(heat_df, transformation != "log2FC"),
        aes(x = transformation, y = miRNA, fill = monotonicity),
        color = "white", linewidth = 0.2
      ) +
      # 2) new log2FC tiles
      geom_tile(
        data = filter(heat_df, transformation == "log2FC"),
        aes(x = transformation, y = miRNA, fill = log2FC),
        color = "white", linewidth = 0.2
      ) +
      # scales for monotonicity
      scale_fill_manual(
        name   = "Trend",
        values = c(
          increasing = "#1b9e77",
          decreasing = "#d95f02",
          none       = "#FFFF00"
        ),
        labels = c(
          increasing = "Increasing",
          decreasing = "Decreasing",
          none       = "Non-monotonic"
        ),
        na.value = "white",
        guide = guide_legend(order = 1)
      ) +
      # continuous scale for the log2FC column only
      scale_fill_gradient2(
        name     = "log₂FC",
        low      = "red",
        mid      = "white",
        high     = "green",
        midpoint = 0,
        guide    = guide_colorbar(order = 2)
      ) +
      scale_x_discrete(
        labels = c(name_map, log2FC = "log₂FC"),
        expand = c(0, 0)
      ) +
      labs(
        x     = "Normalization method",
        y     = paste0("miRNA (top ", topN, " by RPM)"),
        title = "Per-miRNA monotonic trend + log₂FC"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        panel.grid   = element_blank()
      )
    
    pdf(paste0(plot_dir, "top100_wFC.pdf"), width = 8, height = 10)
    print(gg)  
    dev.off()
    
    quit()
      
  }
  
  
  # Build a long table of monotonicity for those 100
  mono_list <- lapply(names(list_of_transformations), function(key) {
    grp  <- geometric_mean_by_group(list_of_transformations[[key]], groups)
    grp100 <- grp[top100, , drop = FALSE]                 # only top100 rows
    mono <- apply(grp100, 1, check_monotonicity)
    data.frame(
      miRNA          = rownames(grp100),
      transformation = key,
      monotonicity   = mono,
      stringsAsFactors = FALSE
    )
  })
  heat_df <- bind_rows(mono_list)
  
  # Factor levels to fix the x-axis order
  heat_df <- heat_df %>%
    mutate(
      miRNA          = factor(miRNA, levels = top100),
      transformation = factor(transformation, levels = names(colour_map))
    )
  
  # Draw the heatmap
  gg <- ggplot(heat_df,
               aes(x = transformation, y = miRNA, fill = monotonicity)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_manual(
      name   = "trend",
      values = c(increasing = "#1b9e77",
                 decreasing = "#d95f02",
                 none       = "#FFFF00"),
      labels = c(increasing = "Increasing",
                 decreasing = "Decreasing",
                 none       = "Non-monotonic"),
      na.value = "white"
    ) +
    scale_x_discrete(
      labels = name_map,
      expand = c(0, 0)
    ) +
    labs(
      x     = "Normalization method",
      y     = paste0("miRNA (top ", topN, " by RPM)"),
      title = "Per-miRNA monotonic trend"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6),
      panel.grid   = element_blank()
    )
  
  pdf(paste0(plot_dir, "top100.pdf"), width = 8, height = 10)
  print(gg)  
  dev.off()
  
  # interactive heatmap
  # Recompute the top 100 miRNAs by RPM
  rpm_key     <- grep("RPM", names(list_of_transformations), value = TRUE)[1]
  rpm_grp     <- geometric_mean_by_group(list_of_transformations[[rpm_key]], groups)
  average_RPM <- rowMeans(rpm_grp)
  top100      <- names(sort(average_RPM, decreasing = TRUE))[1:100]
  
  # Build a long df with monotonicity + expression text
  mono_list <- lapply(names(list_of_transformations), function(key) {
    grp     <- geometric_mean_by_group(list_of_transformations[[key]], groups)
    grp100  <- grp[top100, , drop = FALSE]
    mono    <- apply(grp100, 1, check_monotonicity)
    # pack all group‐means into one tooltip string per miRNA
    expr_txt <- apply(grp100, 1, function(vals) {
      paste0(names(vals), ": ", round(vals, 2), collapse = "\n")
    })
    data.frame(
      miRNA          = rownames(grp100),
      transformation = key,
      monotonicity   = mono,
      expr_text      = expr_txt,
      stringsAsFactors = FALSE
    )
  })
  heat_df <- bind_rows(mono_list) %>%
    mutate(
      miRNA          = factor(miRNA, levels = top100),
      transformation = factor(transformation, levels = names(colour_map)),
      # build the full hover text
      tooltip = paste0(
        "miRNA: ", miRNA, "\n",
        "Method: ", name_map[transformation], "\n",
        "Trend: ", monotonicity, "\n\n",
        expr_text
      )
    )
  
  # Static ggplot with text aesthetic
  p <- ggplot(heat_df,
              aes(x = transformation,
                  y = miRNA,
                  fill = monotonicity,
                  text = tooltip)) +
    geom_tile(color = "white", size = 0.2) +
    scale_fill_manual(
      name   = "Trend",
      values = c(increasing = "#1b9e77",
                 decreasing = "#d95f02",
                 none       = "#FFFFFF"),
      labels = c(increasing = "Increasing",
                 decreasing = "Decreasing",
                 none       = "Non-monotonic")
    ) +
    scale_x_discrete(labels = name_map, expand = c(0, 0)) +
    labs(x = "Normalization method",
         y = "miRNA (top 100 by RPM)",
         title = "Per-miRNA monotonic trend") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6),
      panel.grid   = element_blank()
    )
  
  # Turn it interactive and save as HTML
  interactive_heatmap <- ggplotly(p, tooltip = "text")
  saveWidget(interactive_heatmap,
             file          =  paste0(plot_dir, "interactive_heatmap.html"),
             selfcontained = TRUE)
  
}


if (FALSE) {
    RAW_DF     <- rpm_grp
    
    # pick top100 by RPM
    topN   <- 100
    top100 <- names(sort(average_RPM, decreasing = TRUE))[1:topN]
    
    # monotonicity for each normalization
    mono_list <- lapply(names(list_of_transformations), function(key) {
      grp100 <- geometric_mean_by_group(list_of_transformations[[key]], groups)[ top100, , drop = FALSE ]
      data.frame(
        miRNA          = rownames(grp100),
        transformation = key,
        monotonicity   = apply(grp100, 1, check_monotonicity),
        stringsAsFactors = FALSE
      )
    })
    heat_df <- bind_rows(mono_list)
    
    # compute log2FC and turn into the same categorical trend
    idx    <- match(top100, rownames(RAW_DF))
    fc_vec <- log2(RAW_DF[idx, "100S"] / RAW_DF[idx, "0S"])
    fc_df  <- data.frame(
      miRNA          = top100,
      transformation = "log2FC",
      monotonicity   = ifelse(
        fc_vec >  0, "increasing",
        ifelse(fc_vec <  0, "decreasing", "none")
      ),
      stringsAsFactors = FALSE
    )
    heat_df <- bind_rows(heat_df, fc_df)
    
    # set factor levels (so log2FC appears last)
    heat_df <- heat_df %>%
      mutate(
        miRNA          = factor(miRNA, levels = top100),
        transformation = factor(
          transformation,
          levels = c(names(colour_map), "log2FC")
        )
      )
    
    # single‐scale heatmap
    gg <- ggplot(heat_df, aes(transformation, miRNA, fill = monotonicity)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_manual(
        name   = "Trend",
        values = c(
          increasing = "#1b9e77",
          decreasing = "#d95f02",
          none       = "#FFFF00"
        ),
        labels   = c(
          increasing = "Increasing",
          decreasing = "Decreasing",
          none       = "Non-monotonic"
        ),
        na.value = "white"
      ) +
      scale_x_discrete(
        labels = c(name_map, log2FC = "log₂FC"),
        expand = c(0, 0)
      ) +
      labs(
        x     = "Normalization method",
        y     = paste0("miRNA (top ", topN, " by RPM)"),
        title = "Per-miRNA monotonic trend + log₂FC"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        panel.grid.major = element_line(color = "black")
      )
    
    
    pdf(paste0(plot_dir, "top100_wFC.pdf"), width = 8, height = 10)
    print(gg)
    dev.off()
    
  }
  
# heatmap with FC between pure groups
generate_monotonic_heatmap(
  list_of_transformations,
  groups,
  rpm_grp,
  name_map,
  plot_dir
)

tissue_specific <- c ("mmu-miR-192-5p", "mmu-miR-142a-5p", "mmu-miR-146a-5p")

for (mir in tissue_specific){
  generate_correlation_heatmap(
    list_of_transformations,
    groups,
    rpm_grp,
    name_map,
    plot_dir,
    mir
  )
}



quit()

# AFTER THIS QUIT, IT'S OLD STUFF OR TESTS


# heatmap of correlations to tissue specific miRNAs
# markers <- c("mmu-miR-143-3p", "mmu-miR-192-5p", "mmu-miR-122-5p",
# "mmu-miR-142a-3p", "mmu-miR-142a-5p", "mmu-miR-194-5p",
# "mmu-miR-101b-3p", "mmu-miR-146a-5p", "mmu-miR-150-5p") 
tissue_specific <- c ("mmu-miR-192-5p", "mmu-miR-142a-5p", "mmu-miR-146a-5p")





quit()


if (TRUE) {
  # 1. Identify RPM and pick top 100 miRNAs
  rpm_key     <- grep("RPM", names(list_of_transformations), value = TRUE)[1]
  rpm_grp     <- geometric_mean_by_group(list_of_transformations[[rpm_key]], groups)
  average_RPM <- rowMeans(rpm_grp)
  topN        <- 100
  top100      <- names(sort(average_RPM, decreasing = TRUE))[1:topN]
  
  # 2. Build monotonicity for each method
  mono_list <- lapply(names(list_of_transformations), function(key) {
    grp     <- geometric_mean_by_group(list_of_transformations[[key]], groups)
    grp100  <- grp[top100, , drop = FALSE]
    mono    <- apply(grp100, 1, check_monotonicity)
    data.frame(
      miRNA          = rownames(grp100),
      transformation = key,
      monotonicity   = mono,
      stringsAsFactors = FALSE
    )
  })
  heat_df <- bind_rows(mono_list)
  
  # 3. Compute log2FC between two raw‐counts columns and categorise
  #    — placeholders RAW_DF, "COL_A", "COL_B"
  
  RAW_DF <- rpm_grp  
  # Compute log2FC from two raw columns (placeholders RAW_DF, COL_A, COL_B)
  # Ensure miRNA names align to top100 order
  fc_vec <- log2(
    RAW_DF[top100, "100S", drop = TRUE] /
      RAW_DF[top100, "0S" , drop = TRUE]
  )
  
  
  fc_df <- data.frame(
    miRNA          = top100,
    transformation = "log2FC",
    monotonicity   = ifelse(
      fc_vec > 0, "increasing",
      ifelse(fc_vec < 0, "decreasing", "none")
    ),
    stringsAsFactors = FALSE
  )
  heat_df <- bind_rows(heat_df, fc_df)
  
  # 4. Fix factor levels so the x-axis order includes log2FC last
  heat_df <- heat_df %>%
    mutate(
      miRNA          = factor(miRNA, levels = top100),
      transformation = factor(
        transformation,
        levels = c(names(colour_map), "log2FC")
      )
    )
  
  # 5. Plot single‐scale heatmap
  gg <- ggplot(heat_df, aes(x = transformation, y = miRNA, fill = monotonicity)) +
    geom_tile() +
    # geom_tile(color = "black", linewidth = 0.2) +
    # geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_manual(
      name   = "Trend",
      values = c(
        increasing = "#1b9e77",
        decreasing = "#d95f02",
        none       = "#FFFF00"
      ),
      labels   = c(
        increasing = "Increasing",
        decreasing = "Decreasing",
        none       = "Non-monotonic"
      ),
      na.value = "white"
    ) +
    scale_x_discrete(
      labels = c(name_map, log2FC = "log₂FC"),
      expand = c(0, 0)
    ) +
    labs(
      x     = "Normalization method",
      y     = paste0("miRNA (top ", topN, " by RPM)"),
      title = "Per-miRNA monotonic trend + log₂FC"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6),
      # panel.grid   = element_blank()
    )
  
  # 6. Save
  pdf(paste0(plot_dir, "top100_with_fc.pdf"), width = 8, height = 10)
  print(gg)
  dev.off()
}

# interactive heatmap
if (FALSE){
  library(conflicted)
  conflicts_prefer(plotly::layout)
  
  # 0. libraries
  library(dplyr)
  library(tidyr)
  library(plotly)
  library(crosstalk)
  library(htmlwidgets)
  
  # 1. define human labels and colours (adjust to your palette)
  name_map <- c(
    RC                      = "Raw read count",
    RPM                     = "RPM",
    `edgeR-TMM`             = "edgeR TMM",
    `DESeq2-rlog_mean`      = "DESeq2 rlog (mean)",
    `DESeq2-rlog_parametric`= "DESeq2 rlog (parametric)",
    `DESeq2-vst_local`      = "DESeq2 VST (local)"
  )
  colour_map <- c(
    RC                      = "#E69F00",
    RPM                     = "#56B4E9",
    `edgeR-TMM`             = "#009E73",
    `DESeq2-rlog_mean`      = "#F0E442",
    `DESeq2-rlog_parametric`= "#0072B2",
    `DESeq2-vst_local`      = "#D55E00"
  )
  
  # 2. pick top100 by average RPM
  rpm_key   <- grep("RPM", names(list_of_transformations), value=TRUE)[1]
  rpm_grp   <- geometric_mean_by_group(list_of_transformations[[rpm_key]], groups)
  top100    <- names(sort(rowMeans(rpm_grp), decreasing=TRUE))[1:100]
  
  # 3. build heatmap data frame
  heat_df <- lapply(names(list_of_transformations), function(key) {
    grp100 <- geometric_mean_by_group(list_of_transformations[[key]], groups)[top100, , drop=FALSE]
    mono   <- apply(grp100, 1, check_monotonicity)
    expr   <- apply(grp100, 1, function(v) paste0(names(v), ": ", round(v,2), collapse="\n"))
    data.frame(
      miRNA          = rownames(grp100),
      transformation = key,
      monotonicity   = mono,
      tooltip        = paste0(
        "miRNA: ", rownames(grp100), "\n",
        "Method: ", name_map[key], "\n",
        "Trend: ", mono, "\n\n",
        expr
      ),
      cell_id = paste(rownames(grp100), key, sep=" | "),
      stringsAsFactors=FALSE
    )
  }) %>% bind_rows() %>%
    mutate(
      miRNA          = factor(miRNA, levels=top100),
      transformation = factor(transformation, levels=names(name_map))
    )
  
  # 4. build long data for dot-plot
  long_df <- lapply(names(list_of_transformations), function(key) {
    grp100 <- geometric_mean_by_group(list_of_transformations[[key]], groups)[top100, , drop=FALSE]
    df     <- as.data.frame(grp100)
    df$miRNA          <- rownames(df)
    df$transformation <- key
    pivot_longer(df,
                 cols = -c(miRNA, transformation),
                 names_to  = "group",
                 values_to = "value"
    ) %>% mutate(cell_id = paste(miRNA, transformation, sep=" | "))
  }) %>% bind_rows()
  
  # 5. wrap in SharedData
  sd <- SharedData$new(long_df, key=~cell_id)
  
  # 6. heatmap as a grid of squares
  p_heat <- plot_ly(
    sd,
    x        = ~transformation,
    y        = ~miRNA,
    color    = ~monotonicity,
    colors   = colour_map,
    type     = "scatter",
    mode     = "markers",
    key      = ~cell_id,
    hoverinfo= "text",
    text     = ~tooltip,
    marker   = list(
      symbol = "square",
      size   = 20,
      line   = list(color="white", width=1)
    )
  ) %>%
    layout(
      xaxis = list(
        title     = "Normalization method",
        tickangle = 45,
        tickmode  = "array",
        tickvals  = levels(heat_df$transformation),
        ticktext  = name_map
      ),
      yaxis = list(
        title    = "miRNA (top 100 by RPM)",
        autorange= "reversed"
      )
    ) %>%
    highlight(on        = "plotly_click",
              off       = "plotly_doubleclick",
              selectize = FALSE)
  
  # 7. dot-plot that auto-filters by the same key
  p_dot <- plot_ly(
    sd,
    x        = ~group,
    y        = ~value,
    type     = "scatter",
    mode     = "markers+lines",
    marker   = list(size=8),
    line     = list(shape="linear"),
    hoverinfo= "text",
    text     = ~paste0("Group: ", group, "<br>Value: ", round(value,2))
  ) %>%
    layout(
      xaxis = list(title="Group"),
      yaxis = list(title="Mean expression")
    )
  
  # 8. combine and save
  fig <- subplot(p_heat, p_dot, widths=c(0.7,0.3), shareY=FALSE) %>%
    plotly::layout(title="Heatmap ↔ Linked dot-plot")
  
  saveWidget(fig,
             file          =  paste0(plot_dir, "interactive_heatmap_with_dotplot.html"),
             selfcontained = TRUE)
  
}

# facet scatter
if (TRUE){
  top20 <- top100
  sample_vec <- unlist(groups, use.names = FALSE)
  group_vec  <- rep(names(groups), lengths(groups))
  # Now:
  #   sample_vec[i] is a sample name
  #   group_vec[i] is the corresponding group label (the list‐element name)
  
  sample_to_group_df <- data.frame(
    sample = sample_vec,
    group  = group_vec,
    stringsAsFactors = FALSE
  )
  
  sample_to_group_df$mouse <- as.numeric(
    sub(".*?([0-9]+)(S|L).*", "\\1", sample_to_group_df$sample)
  ) 
  
  # ──────── STEP 2: Build a “long” table of sample‐level values for the top20 miRNAs ────────
  # We'll loop over each normalization method, subset to top20, pivot longer, and join group info.
  
  long_list <- lapply(names(list_of_transformations), function(method_key) {
    mat <- list_of_transformations[[method_key]]
    
    # Subset to the top‐20 miRNAs
    mat20 <- mat[top20, , drop = FALSE]
    
    # Convert to a data.frame and pivot longer
    df <- as.data.frame(mat20)
    df$miRNA <- rownames(df)
    
    long_df <- df %>%
      pivot_longer(
        cols = -miRNA,
        names_to  = "sample",
        values_to = "value"
      ) %>%
      # Join on sample_to_group_df to assign group
      left_join(sample_to_group_df, by = "sample") %>%
      mutate(
        method = method_key
      )
    
    long_df
  })
  
  expr_long <- bind_rows(long_list) %>%
    # Drop any sample that didn’t match a group (i.e. group == NA)
    filter(!is.na(group)) %>%
    mutate(
      miRNA  = factor(miRNA, levels = top20),
      method = factor(method, levels = names(list_of_transformations)),
      # Enforce the exact order you want:
      group  = factor(group,
                      levels = c("0S", "20S", "40S", "60S", "80S", "100S"))
    )
  expr_long <- expr_long %>%
    mutate(
      mouse = factor(mouse)   # *** ADDED ***
    )
  
  # ──────── STEP 3: Open a multi‐page PDF device ────────
  pdf(
    # file   = paste0(plot_dir, "top20_scatter_by_method.pdf"),
    file   = paste0(plot_dir, "top100_scatter_by_method.pdf"),
    width  = 8,
    height = 6
  )
  
  # ──────── STEP 4: Loop over each miRNA and draw the facetted plot ────────
  
  for (mir in top20) {
    df_sub <- filter(expr_long, miRNA == mir)
    
    summ_df <- df_sub %>%
      group_by(method, group) %>%
      summarise(
        mean_val = mean(value),                    # arithmetic mean
        geo_mean = exp(mean(log(value))),          # geometric mean
        geo_trim = {                               # geo mean without furthest point
          vals    <- value
          drop_i  <- which.max(abs(vals - mean(vals)))
          exp(mean(log(vals[-drop_i])))
        },
        .groups = "drop"
      )
    
    p <- ggplot(df_sub, aes(x = group,  group = mouse, y = value, color = mouse)
) +
    # p <- ggplot(df_sub, aes(x = group, y = value)) +
      # geom_boxplot(outlier.shape = NA) +
      geom_line(alpha = 0.8, size = 0.7) +
      
      # ─── summary “step” markers ───
      geom_point(
        data  = summ_df,
        aes(x = group, y = mean_val),
        inherit.aes = FALSE,
        shape  = 45,    # horizontal bar
        size   = 6,
        color  = "black"
      ) +  # *** ADDED: arithmetic mean marker ***
      geom_point(
        data  = summ_df,
        aes(x = group, y = geo_mean),
        shape  = 95,
        size   = 6,
        inherit.aes = FALSE,
        color  = "darkgreen"
      ) +  # *** ADDED: geometric mean marker ***
      geom_point(
        data  = summ_df,
        aes(x = group, y = geo_trim),
        shape  = 95,
        size   = 6,
        inherit.aes = FALSE,
        color  = "darkorange"
      ) +  # *** ADDED: trimmed geometric mean marker ***
      
      # geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
      scale_color_brewer(palette = "Set1") +
      geom_point(size = 2) +
      facet_wrap(~ method, ncol = 3, 
                 scales = "free_y") +
      # scale_x_discrete(drop = FALSE) +   # ensure every named group appears on each facet
      labs(
        title = paste0(mir, "– expression by method"),
        x     = "Transformation",
        y     = "Expression value"
      ) +
      theme_minimal() +
      theme(
        strip.text       = element_text(size = 8),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    print(p)
  }
  
  dev.off()
  
  
}


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

order_levels <- cum_result %>%
  group_by(transformation) %>%
  slice_max(order_by = A, n = 1, with_ties = FALSE) %>%
  arrange(desc(sd_M)) %>%
  pull(transformation)

# print(length(order_levels))
# print(length(name_map))
# quit()
gg <- ggplot(cum_result, aes(x = A , y = sd_M, color = transformation, linetype = transformation)) +
  # gg <- ggplot(plot_df, aes(x = A, y = M)) +
  geom_line() +
  scale_linetype_manual(values = line_type_mapping,
                        labels = name_map,
                        breaks = order_levels
                        ) +
  scale_colour_manual(
    name   = "Normalization method",
    values = colour_map,
    breaks = order_levels,
    labels = name_map
  ) +
  scale_linetype_manual(
    name   = "Normalization method",
    values = line_type_mapping,
    breaks = order_levels,
    labels = name_map
  ) +
  theme_minimal() + 
  labs(
    title = paste0(""),
    x = "Average expression",
    y = "sd(log2FC)"
  )

pdf(paste0(plot_dir, "sd_mean_plot.pdf"), width = 8, height = 6)
print(gg)  
dev.off()  

##
### O/E values  
##
##

RPM_threshold <- minRPM

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

ggsave(paste0(plot_dir, "log_OE.pdf"), plot = gg, width = 8, height = 6, device = "pdf")

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

gg <- ggplot(plot_df,
             aes(x = transformation,
                 y = ObservedExpectedRatio,
                 fill = transformation)) +
  geom_boxplot(alpha = 0.8) +                # show.legend defaults to TRUE
  guides(alpha = "none") +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "blue") +
  scale_y_log10() +
  
  # apply your human-readable names on the x-axis
  scale_x_discrete(
    labels = name_map
  ) +
  
  # manual fill scale with your fixed colours and labels
  scale_fill_manual(
    name   = "normalization method",
    values = colour_map,
    labels = name_map
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = "Observed/Expected Ratio"
  )

ggsave(paste0(plot_dir, "log_mouse_OE.pdf"), plot = gg, width = 8, height = 6, device = "pdf")

quit()


#### TO ADD

# constant miRNAs?
# iteratively remove miRNAs (especially high) and check if this norm improves