# miRNA-seq_benchmark

## step 0: install dependencies
I developed this using R version 4.4.2 (2024-10-31)
R packages needed are in R_packages_snapshot.rds
To install go in R and run:
```{r}
pkgs_snapshot <- readRDS("R_packages_snapshot.rds")
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(pkgs_snapshot)
```
## instructions to run normalization benchmark
go in normalization_benchmark:
Rscript normalization_benchmark.R

## instructions to run differential expression benchmark
go in DE_benchmark and:
Rscript 1_generate_DE_results.R
Rscript 2_benchmark_DE.R
