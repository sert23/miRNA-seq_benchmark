library(ggplot2)
library(ComplexUpset)
library(dplyr)
library(tidyr)
library(grid)

# Marc comments about housekeeping miRNAs
# From a paper "Different authors analyzed the suitability of other miRNAs such as miR-16, miR-192, miR-374a, miR-374b and let-7d or combinations of miRNAs (such as miR-146a, miR-16, miR-195, miR-30e and miR-744) for normalization of circulating miRNAs levels.
# From another paper: HKG = housekeeping genes. "The potential HKGs were Let-7a, Let-7d, Let-7g, miR-16, RNU6, RNU48, miR-191, miR-223, miR-484, and miR-520d-5p." - Of these, miR-223 is neutrophil/macrophage, so that won't work for you.  The others might
# source("/shared/results/MNM_miRNA_reference/normalizations/R/normalization.R")
source("../source/utils.R")

benchmark_folder <- "./pairwise_DE"

RPMl_df <- read.delim("./input_data/mature_sense_minExpr0_RCadj_libraryRPM.mat",
                      check.names = FALSE, row.names = 1)
colnames(RPMl_df) <- fix_names(colnames(RPMl_df))

# read monotonic
mono_df <- read.csv("./input_data/monotonic.csv")

# setting RPM threshold (currently 20 RPM)
RPMl_df <- RPMl_df[, colnames(RPMl_df) %in% c("1L", "2L", "3L", "4L",  "5L", 
                                              "1S", "2S", "3S", "4S",  "5S")]
RPMl_df <- RPMl_df[rowMeans(RPMl_df) >= 20, ]
miRNAs_over_thres <- rownames(RPMl_df)

### ROC curves per method (internal)
# 0S_VS_100S is ground truth for each method

aggregated_0vs100S <- aggregate_tsv_files(benchmark_folder, "0S_VS_100S")

# remove entries below RPM threshold
aggregated_0vs100S <- aggregated_0vs100S[aggregated_0vs100S$miRNA %in% miRNAs_over_thres, ]

aggregated_20vs80S <- aggregate_tsv_files(benchmark_folder, "20S_VS_80S")
# remove entries below RPM threshold
aggregated_20vs80S <- aggregated_20vs80S[aggregated_20vs80S$miRNA %in% miRNAs_over_thres, ]

aggregated_40vs80S <- aggregate_tsv_files(benchmark_folder, "40S_VS_80S")
# remove entries below RPM threshold
aggregated_40vs80S <- aggregated_40vs80S[aggregated_40vs80S$miRNA %in% miRNAs_over_thres, ]

aggregated_40vs60S <- aggregate_tsv_files(benchmark_folder, "40S_VS_60S")
# remove entries below RPM threshold
aggregated_40vs60S <- aggregated_40vs60S[aggregated_40vs60S$miRNA %in% miRNAs_over_thres, ]


# aggregated_0vs100S <- reverse_fc(aggregated_0vs100S, c("NBSR", "miRglmm" ))
# aggregated_20vs80S <- reverse_fc(aggregated_20vs80S, c("NBSR", "miRglmm" ))
# aggregated_40vs80S <- reverse_fc(aggregated_40vs80S, c("NBSR", "miRglmm" ))
# aggregated_40vs60S <- reverse_fc(aggregated_40vs60S, c("NBSR", "miRglmm" ))

# create plot folder
plot_dir <- "./plots/"
dir.create(plot_dir, recursive = TRUE)

roc_data <- create_ROC_data(aggregated_0vs100S, aggregated_20vs80S)
roc_plot <- plot_ROC(roc_data, paste0(plot_dir, "20vs80_internal.pdf"), "20Svs80S internal consistency")
                     # "/shared/results/MNM_miRNA_reference/DE_benchmark/20vs80test.pdf")

roc_data <- create_ROC_data(aggregated_0vs100S, aggregated_40vs80S)
roc_plot <- plot_ROC(roc_data, paste0(plot_dir,  "40vs80_internal.pdf"), "40Svs80S internal consistency")

roc_data <- create_ROC_data(aggregated_0vs100S, aggregated_40vs60S)
roc_plot <- plot_ROC(roc_data, paste0(plot_dir,  "40vs60_internal.pdf"), "40Svs60S internal consistency")

## overlap among methods for internal ground truth
# removing DESeq2 (too many parameter combinations)

aggregated_0vs100S_minus_DESeq <- subset(aggregated_0vs100S, Method != "DESeq2")

plot_dest <- file.path(paste0(plot_dir,  "upset_0vs100_up_down.pdf"))

# formatting monotonic miRNAs data as the rest
mono_df$Method <- mono_df$DE_method
mono_df$Parameter <- mono_df$parameters
mono_df$Comparison <- mono_df$comparison
mono_df <- mono_df[c("miRNA",	"average",	"log2FC",	"pval",	"padj",	"comparison",
                     "DE_method",	"parameters", "Method", "Parameter", "Comparison",
                     "FC_direction")]

# print(head(mono_df))
aggregated_0vs100S_minus_DESeq <- rbind(aggregated_0vs100S_minus_DESeq, mono_df)
test_LRT_fit_local_sf_poscounts <- subset(aggregated_0vs100S, parameters == "test_LRT_fit_local_sf_poscounts")
agg <- rbind(aggregated_0vs100S_minus_DESeq, test_LRT_fit_local_sf_poscounts)
agg <- agg[agg$miRNA %in% miRNAs_over_thres, ]
create_upset_plot(agg, 0.05, 2, plot_dest)

# Upset DESeq2 parameter combinations

aggregated_0vs100S_DESeq_only <- subset(aggregated_0vs100S, Method == "DESeq2")
aggregated_0vs100S_DESeq_only <- aggregated_0vs100S_DESeq_only[aggregated_0vs100S_DESeq_only$miRNA 
                                                                %in% miRNAs_over_thres, ]
plot_dest <- file.path(paste0(plot_dir,  "upset_DESeq_param.pdf"))
par_res<- create_upset_plot_by_parameters(aggregated_0vs100S_DESeq_only, 0.05, plot_dest)

## make ROC curves with common ground truth ("absolute" ground truth)
# TODO they need to be fixed because the negative ground truth is empty
# we have to decide a different way of declaring a negative set

# roc_results <- create_ROC_data_common(aggregated_20vs80S, negative_set, positive_set, padj_threshold = 0.05)
# # roc_plot <- plot_ROC(roc_results, 
# #                      file.path(paste0(plot_dir, "20vs80_absolute.pdf")))
# 
# roc_results <- create_ROC_data_common(aggregated_40vs80S, negative_set, positive_set, padj_threshold = 0.05)
# # roc_plot <- plot_ROC(roc_results, 
# #                      file.path(paste0(plot_dir, "40vs80_absolute.pdf")))
# 
# roc_results <- create_ROC_data_common(aggregated_40vs60S, negative_set, positive_set, padj_threshold = 0.05)
# # roc_plot <- plot_ROC(roc_results, 
# #                      file.path(paste0(plot_dir, "40vs60_absolute.pdf")))
# 
# 

### Monotonic increase/decrease ground truth ROC
# load data
aggregated_0vs100S <- aggregate_tsv_files(benchmark_folder, "0S_VS_100S")
aggregated_0vs100S <- aggregated_0vs100S[aggregated_0vs100S$miRNA %in% miRNAs_over_thres, ]

aggregated_20vs80S <- aggregate_tsv_files(benchmark_folder, "20S_VS_80S")
aggregated_20vs80S <- aggregated_20vs80S[aggregated_20vs80S$miRNA %in% miRNAs_over_thres, ]

aggregated_40vs80S <- aggregate_tsv_files(benchmark_folder, "40S_VS_80S")
aggregated_40vs80S <- aggregated_40vs80S[aggregated_40vs80S$miRNA %in% miRNAs_over_thres, ]

aggregated_40vs60S <- aggregate_tsv_files(benchmark_folder, "40S_VS_60S")
aggregated_40vs60S <- aggregated_40vs60S[aggregated_40vs60S$miRNA %in% miRNAs_over_thres, ]

# aggregated_0vs100S <- reverse_fc(aggregated_0vs100S, c("NBSR", "miRglmm" ))
# aggregated_20vs80S <- reverse_fc(aggregated_20vs80S, c("NBSR", "miRglmm" ))
# aggregated_40vs80S <- reverse_fc(aggregated_40vs80S, c("NBSR", "miRglmm" ))
# aggregated_40vs60S <- reverse_fc(aggregated_40vs60S, c("NBSR", "miRglmm" ))

# load monotonic and filter miRNAs below threshold
mono_df <- read.csv("./input_data/monotonic.csv")
mono_df <- mono_df[mono_df$miRNA %in% miRNAs_over_thres, ]

# ROC DESeq only
# 
DESeq_only <- aggregated_0vs100S[aggregated_0vs100S$DE_method == "DESeq2", ]
DESeq_only$Method <- DESeq_only$parameters
roc_results <- create_ROC_data_mono(DESeq_only, mono_df, padj_threshold = 0.05)
roc_plot <- plot_ROC(roc_results, 
                     file.path(paste0(plot_dir, "mono_ROC_DESeq.pdf")))

# keep only DESEq test_LRT_fit_local_sf_poscounts
aggregated_0vs100S_minus_DESeq <- subset(aggregated_0vs100S, Method != "DESeq2")
test_LRT_fit_local_sf_poscounts <- subset(aggregated_0vs100S, parameters == "test_LRT_fit_local_sf_poscounts")
agg1 <- rbind(aggregated_0vs100S_minus_DESeq, test_LRT_fit_local_sf_poscounts)

roc_results <- create_ROC_data_mono(agg1, mono_df, padj_threshold = 0.05)
# print("after")
# wilcoxon_only <- subset(roc_results, Method == "Wilcoxon")
# print(wilcoxon_only)
roc_plot <- plot_ROC(roc_results, 
                     file.path(paste0(plot_dir, "mono_ROC_0vs100.pdf")))
                     
# aggregated_0vs100S, aggregated_20vs80S, aggregated_40vs80S, aggregated_40vs60S
comparisons <- list(c_0vs100S = aggregated_0vs100S, 
                    c_20vs80S = aggregated_20vs80S,
                    c_40vs80S = aggregated_40vs80S,
                    c_40vs60S = aggregated_40vs60S
                    )

# generate plots for all comparisons
for (name in names(comparisons)){
  df <- comparisons[[name]]
  print(name)
  df_minus_DESeq <- subset(df, Method != "DESeq2")
  test_LRT_fit_local_sf_poscounts <- subset(df, parameters == "test_LRT_fit_local_sf_poscounts")
  agg <- rbind(df_minus_DESeq, test_LRT_fit_local_sf_poscounts)
  roc_results <- create_ROC_data_mono(agg, mono_df, padj_threshold = 0.05)
  dest_plot <- file.path(paste0(plot_dir, "mono_ROC_", 
                   name, ".pdf"))
  roc_plot <- plot_ROC(roc_results, dest_plot, substring(name, 3))
}








