library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
# Load necessary libraries
library(DESeq2)
library(edgeR)
library(tidyverse)
library(limma)
source("../source/utils.R")

# output folder
benchmark_folder <- "./pairwise_DE"

# load matrix and change labels

# load miRNA_data
input_df <- read.delim("input_data/mature_sense_minExpr0_RCadj.mat", 
                       check.names = FALSE, row.names = 1)
colnames(input_df) <- fix_names(colnames(input_df))

# define which samples belong to each group
groups <- define_groups(input_df)
group_list <- names(groups)

# generate all possible comparisons
pairwise_comparisons <- combn(group_list, 2, simplify = FALSE)

#
## limma-voom 
# columns average, log2FC, pval, padj, comparison, method, parameters
#

method_name <- "limma"
method_dir <- file.path(benchmark_folder, method_name)
combo_name <- "default"
combo_dir <- file.path(method_dir, combo_name)
dir.create(combo_dir, recursive = TRUE)

for (pair in pairwise_comparisons) {
  comparison_name <- paste(c(pair), collapse = "_VS_")
  comparison_dir <- file.path(combo_dir, comparison_name)
  if (!dir.exists(comparison_dir)) {
    dir.create(comparison_dir, recursive = TRUE)
  }
  
  tryCatch({
    # Subset the expression matrix to include only the samples present in the sample_info
    group1_samples <- groups[[pair[1]]]
    group2_samples <- groups[[pair[2]]]
    group1_matrix <- input_df[, group1_samples]
    group2_matrix <- input_df[, group2_samples]
    expression_matrix <- cbind(group1_matrix, group2_matrix)
    column_names <- colnames(expression_matrix)
    group <- factor(c(rep(pair[1], ncol(group1_matrix)), rep(pair[2], ncol(group2_matrix))))
    
    dge <- DGEList(counts = expression_matrix, group = group)
    
    # Normalize library sizes (TMM normalization)
    dge <- calcNormFactors(dge)
    
    # Apply voom transformation
    design <- model.matrix(~ group)  # Design matrix for two groups
    voom_data <- voom(dge, design)
    
    # Fit a linear model
    fit <- lmFit(voom_data, design)
    fit <- eBayes(fit)
    
    # Extract results
    limma_results <- topTable(fit, coef = 2, number = Inf, sort.by = "p")
    
    # Add log2FC and adjusted p-values
    limma_results$log2FC <- limma_results$logFC
    limma_results$padj <- limma_results$adj.P.Val
    
    # Calculate average expression for each miRNA
    limma_results$average <- rowMeans(voom_data$E)  # Voom-transformed average expression
    
    # Add metadata columns
    limma_results$comparison <- comparison_name
    limma_results$DE_method <- method_name
    limma_results$parameters <- combo_name
    
    # Reformat output to include required columns
    limma_summary <- data.frame(
      miRNA = rownames(limma_results),
      average = limma_results$average,
      log2FC = limma_results$log2FC,
      pval = limma_results$P.Value,
      padj = limma_results$padj,
      comparison = limma_results$comparison,
      DE_method = limma_results$DE_method,
      parameters = limma_results$parameters
    )
    
    write.table(limma_summary,
                file = file.path(comparison_dir, "limma.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
  }, error = function(e) {
    # Print an error message and remove the directory in case of failure
    print(paste("Error in limma comparison", comparison_name, combo_name))
    print(e)
    # Remove the subdirectory for the failed comparison
    unlink(comparison_dir, recursive = TRUE)
  })
  # Remove empty directory if no results are generated
  if (length(list.files(combo_dir)) == 0) {
    unlink(combo_dir, recursive = TRUE)
  }
}

#
## t-test on RPM

method_name <- "t-test"
method_dir <- file.path(benchmark_folder, method_name)
combo_name <- "default"
combo_dir <- file.path(method_dir, combo_name)
dir.create(combo_dir, recursive = TRUE)

for (pair in pairwise_comparisons) {
  comparison_name <- paste(c(pair), collapse = "_VS_")
  comparison_dir <- file.path(combo_dir, comparison_name)
  if (!dir.exists(comparison_dir)) {
    dir.create(comparison_dir, recursive = TRUE)
  }
  
  tryCatch({
    # Subset the expression matrix to include only the samples present in the sample_info
    group1_samples <- groups[[pair[1]]]
    group2_samples <- groups[[pair[2]]]
    group1_matrix <- input_df[, group1_samples]
    group2_matrix <- input_df[, group2_samples]
    expression_matrix <- cbind(group1_matrix, group2_matrix)
    column_names <- colnames(expression_matrix)
    group <- factor(c(rep(pair[1], ncol(group1_matrix)), rep(pair[2], ncol(group2_matrix))))
    
    # Normalize to RPM
    total_counts <- colSums(expression_matrix) # Total counts per sample
    rpm_matrix <- sweep(expression_matrix, 2, total_counts, FUN = "/") * 1e6
    
    # Calculate log2 fold change (log2FC)
    group1_indices <- group == levels(group)[1] # Indices for group 1
    group2_indices <- group == levels(group)[2] # Indices for group 2
    log2fc <- apply(rpm_matrix, 1, function(row) {
      mean_group1 <- mean(row[group1_indices])
      mean_group2 <- mean(row[group2_indices])
      log2(mean_group2 + 1) - log2(mean_group1 + 1) # Add 1 to avoid log2(0)
    })
    
    # Calculate average expression for each miRNA
    average_expression <- rowMeans(rpm_matrix)
    
    # Run t-test
    t_test_results <- apply(rpm_matrix, 1, function(row) {
      t.test(row ~ group) # Perform t-test for each miRNA
    })
    
    # Create summary table
    t_test_summary <- data.frame(
      miRNA = rownames(rpm_matrix),
      average = average_expression,
      log2FC = log2fc,
      pval = sapply(t_test_results, function(x) x$p.value),
      padj = p.adjust(sapply(t_test_results, function(x) x$p.value), method = "BH"),
      comparison = comparison_name,
      DE_method = method_name,
      parameters = combo_name
    )
    
    # Write the results to a file
    write.table(t_test_summary,
                file = file.path(comparison_dir, "t-test.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
  }, error = function(e) {
    # Print an error message and remove the directory in case of failure
    print(paste("Error in comparison", comparison_name, combo_name))
    print(e)
    # Remove the subdirectory for the failed comparison
    unlink(comparison_dir, recursive = TRUE)
  })
  # Remove empty directory if no results are generated
  if (length(list.files(combo_dir)) == 0) {
    unlink(combo_dir, recursive = TRUE)
  }
}



#
## wilcoxon on RPM

method_name <- "Wilcoxon"
method_dir <- file.path(benchmark_folder, method_name)
combo_name <- "default"
combo_dir <- file.path(method_dir, combo_name)
dir.create(combo_dir, recursive = TRUE)

for (pair in pairwise_comparisons) {
  comparison_name <- paste(c(pair), collapse = "_VS_")
  comparison_dir <- file.path(combo_dir, comparison_name)
  if (!dir.exists(comparison_dir)) {
    dir.create(comparison_dir, recursive = TRUE)
  }
  
  tryCatch({
    # Subset the expression matrix to include only the samples present in the sample_info
    group1_samples <- groups[[pair[1]]]
    group2_samples <- groups[[pair[2]]]
    group1_matrix <- input_df[, group1_samples]
    group2_matrix <- input_df[, group2_samples]
    expression_matrix <- cbind(group1_matrix, group2_matrix)
    column_names <- colnames(expression_matrix)
    group <- factor(c(rep(pair[1], ncol(group1_matrix)), rep(pair[2], ncol(group2_matrix))))
    
    # Normalize to RPM
    total_counts <- colSums(expression_matrix) # Total counts per sample
    rpm_matrix <- sweep(expression_matrix, 2, total_counts, FUN = "/") * 1e6
    
    # Calculate log2 fold change (log2FC)
    group1_indices <- group == levels(group)[1] # Indices for group 1
    group2_indices <- group == levels(group)[2] # Indices for group 2
    log2fc <- apply(rpm_matrix, 1, function(row) {
      mean_group1 <- mean(row[group1_indices])
      mean_group2 <- mean(row[group2_indices])
      log2(mean_group2 + 1) - log2(mean_group1 + 1) # Add 1 to avoid log2(0)
    })
    
    # Calculate average expression for each miRNA
    average_expression <- rowMeans(rpm_matrix)
    
    # Perform Wilcoxon rank-sum test
    wilcox_results <- apply(rpm_matrix, 1, function(row) {
      wilcox.test(row[group1_indices], row[group2_indices]) # Wilcoxon test
    })
    
    # Create summary table
    wilcox_summary <- data.frame(
      miRNA = rownames(rpm_matrix),
      average = average_expression,
      log2FC = log2fc,
      pval = sapply(wilcox_results, function(x) x$p.value),
      padj = p.adjust(sapply(wilcox_results, function(x) x$p.value), method = "BH"),
      comparison = comparison_name,
      DE_method = method_name,
      parameters = combo_name
    )
    
    # Write the results to a file
    write.table(wilcox_summary,
                file = file.path(comparison_dir, "wilcoxon.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
  }, error = function(e) {
    # Print an error message and remove the directory in case of failure
    print(paste("Error in comparison", comparison_name, combo_name))
    print(e)
    # Remove the subdirectory for the failed comparison
    unlink(comparison_dir, recursive = TRUE)
  })
  # Remove empty directory if no results are generated
  if (length(list.files(combo_dir)) == 0) {
    unlink(combo_dir, recursive = TRUE)
  }
}

#
## edgeR
#

method_name <- "edgeR"
method_dir <- file.path(benchmark_folder, method_name)
combo_name <- "default"
combo_dir <- file.path(method_dir, combo_name)
dir.create(combo_dir, recursive = TRUE)

for (pair in pairwise_comparisons) {
  comparison_name <- paste(c(pair), collapse = "_VS_")
  comparison_dir <- file.path(combo_dir, comparison_name)
  if (!dir.exists(comparison_dir)) {
    dir.create(comparison_dir, recursive = TRUE)
  }
  
  tryCatch({
    # Subset the expression matrix to include only the samples present in the sample_info
    group1_samples <- groups[[pair[1]]]
    group2_samples <- groups[[pair[2]]]
    group1_matrix <- input_df[, group1_samples]
    group2_matrix <- input_df[, group2_samples]
    expression_matrix <- cbind(group1_matrix, group2_matrix)
    column_names <- colnames(expression_matrix)
    group <- factor(c(rep(pair[1], ncol(group1_matrix)), rep(pair[2], ncol(group2_matrix))))
    
    # Create DGEList object
    dge <- DGEList(counts = as.matrix(expression_matrix), group = group)
    
    # Estimate common dispersion
    dge <- estimateCommonDisp(dge)
    
    # Exact test for differences in the mean between two groups of negative-binomial counts
    et <- exactTest(dge)
    
    # Extract differentially expressed genes
    topTags_obj <- topTags(et, n = Inf)  # Extract all genes
    res_table <- topTags_obj$table
    
    # Calculate average expression (log2 CPM)
    cpm_matrix <- cpm(dge, log = TRUE)  # Log2-transformed CPM
    average_expression <- rowMeans(cpm_matrix)  # Average log2 CPM
    
    # Add metadata columns
    res_table$average <- average_expression[rownames(res_table)]
    res_table$comparison <- comparison_name
    res_table$DE_method <- method_name
    res_table$parameters <- combo_name
    
    # Reformat output to include required columns
    edgeR_summary <- data.frame(
      miRNA = rownames(res_table),
      average = res_table$average,
      log2FC = res_table$logFC,         # Log2 fold change
      pval = res_table$PValue,          # Raw p-value
      padj = p.adjust(res_table$PValue, method = "BH"),  # Adjusted p-value
      comparison = res_table$comparison,
      DE_method = res_table$DE_method,
      parameters = res_table$parameters
    )
    
    # Save the results to a file
    write.table(
      edgeR_summary,
      file = file.path(comparison_dir, "edgeR.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    
  }, error = function(e) {
    # Print an error message and remove the directory in case of failure
    print(paste("Error in edgeR comparison", comparison_name, combo_name))
    print(e)
    # Remove the subdirectory for the failed comparison
    unlink(comparison_dir, recursive = TRUE)
  })
  
  # Remove empty directory if no results are generated
  if (length(list.files(combo_dir)) == 0) {
    unlink(combo_dir, recursive = TRUE)
  }
}

#
## DESeq
#
# check there are 15 comparisons per parameter

method_name <- "DESeq2"
method_dir <- file.path(benchmark_folder, method_name)

# Define the parameter combinations
test_options <- c("Wald", "LRT")
fitType_options <- c("parametric", "local", "mean", "glmGamPoi")
sfType_options <- c("ratio", "poscounts", "iterate")

# Iterate over each combination of test, fitType, and sfType
for (test_type in test_options) {
  for (fit_type in fitType_options) {
    for (sf_type in sfType_options) {
      
      # Create a unique directory name based on the parameters
      combo_name <- paste("test", test_type, "fit", fit_type, "sf", sf_type, sep = "_")
      combo_dir <- file.path(method_dir, combo_name)
      dir.create(combo_dir, recursive = TRUE)
      
      # Iterate over all pairwise comparisons
      for (pair in pairwise_comparisons) {
        comparison_name <- paste(c(pair), collapse = "_VS_")
        comparison_dir <- file.path(combo_dir, comparison_name)
        if (!dir.exists(comparison_dir)) {
          dir.create(comparison_dir, recursive = TRUE)
        }
        
        tryCatch({
          # Subset the expression matrix to include only the samples present in the sample_info
          group1_samples <- groups[[pair[1]]]
          group2_samples <- groups[[pair[2]]]
          group1_matrix <- input_df[, group1_samples]
          group2_matrix <- input_df[, group2_samples]
          expression_matrix <- cbind(group1_matrix, group2_matrix)
          column_names <- colnames(expression_matrix)
          group_info <- c(rep(pair[1], ncol(group1_matrix)), rep(pair[2], ncol(group2_matrix)))
          sample_info <- data.frame(Column = column_names, Group = group_info, row.names = column_names)
          
          # Define control parameters for sfType = "iterate"
          par_list <- list()
          if (sf_type == "iterate"){
            par_list[["maxit"]] <- 1000
          }
          
          # Create DESeq2 dataset
          dds <- DESeqDataSetFromMatrix(
            countData = expression_matrix,
            colData = sample_info,
            design = ~ Group
          )
          
          # Run DESeq2 with selected parameters
          if (test_type == "LRT") {
            dds <- DESeq(dds, test = test_type, fitType = fit_type, sfType = sf_type, reduced = ~ 1)  
          } else {
            dds <- DESeq(dds, test = test_type, fitType = fit_type, sfType = sf_type)
          }
          
          # Get DE results
          res <- results(dds)
          
          # Compute average expression (log2-transformed normalized counts)
          normalized_counts <- counts(dds, normalized = TRUE)
          average_expression <- rowMeans(log2(normalized_counts + 1))
          
          # Format results to match limma & edgeR output
          deseq_summary <- data.frame(
            miRNA = rownames(res),
            average = average_expression[rownames(res)], # Log2-transformed average expression
            log2FC = res$log2FoldChange,
            pval = res$pvalue,
            padj = res$padj,
            comparison = comparison_name,
            DE_method = method_name,
            parameters = combo_name,
            stringsAsFactors = FALSE
          )
          
          # Save the results to a file
          write.table(
            deseq_summary,
            file = file.path(comparison_dir, "DESeq2.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE
          )
          
        }, error = function(e) {
          # Print an error message and remove the directory in case of failure
          print(paste("Error in comparison", comparison_name, "with parameters:", 
                      "test =", test_type, "fitType =", fit_type, "sfType =", sf_type))
          print(e)
          
          # Remove the subdirectory for the failed comparison
          unlink(comparison_dir, recursive = TRUE)
        })
        
        # Remove empty directory if no results were saved
        if (length(list.files(combo_dir)) == 0) {
          unlink(combo_dir, recursive = TRUE)
        }
      }
    }
  }
}

# Run cleanup for each method 
# (removing parameter combinations where some comparisons fail)
cleanup_directories(method_dir)


