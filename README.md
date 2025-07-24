# HAUDI

## Introduction

HAUDI constructs polygenic scores for individuals with recent admixture.
It does so by using local ancestry information to explicitly model
ancestry-specific effect sizes.

HAUDI takes a similar approach to [GAUDI](https://github.com/quansun98/GAUDI),
which uses a fused lasso approach to encourage sparsity while penalizing
differences in ancestry-specific effects. HAUDI re-parameterizes this model
to a standard LASSO problem, allowing for efficient PGS estimation
using the `bigstatsr` package.

## Installation

```{r}
remotes::install_github("frankp-0/HAUDI")
```

## Instructions

### 1. Prepare input

HAUDI requires as input a set of phased plink2 pgen files, with corresponding
local ancestry files in the .lanc format described by [Admix-kit](https://kangchenghou.github.io/admix-kit/prepare-dataset.html).
We provide a helper function `convert_to_lanc` for converting
[FLARE](https://github.com/browning-lab/flare) and [RFMix](https://github.com/slowkoni/rfmix)
local ancestry files to .lanc format.

### 2. Creating file-backed matrices

HAUDI uses file-backed matrices, implemented in the `bigstatsr`
package to store genotype/ancestry data. To create this file-backed
matrix, use the function `make_fbm`, which takes plink2 and local
ancestry files as input. An example may look like:

```{r}
lanc_files <- system.file(paste0("extdata/chr", 20:22, ".lanc"), package = "HAUDI")
plink_prefixes <- system.file(
  paste0("extdata/chr", 20:22, ".pgen"),
  package = "HAUDI"
) |>
  gsub(pattern = ".pgen", replacement = "")

input <- HAUDI::make_fbm(
  lanc_files = lanc_files,
  plink_prefixes = plink_prefixes,
  fbm_prefix = "toy_data")
```

This command would return a list containing:

1) an object of class FBM.code256 containing genotype/ancestry data
2) a data frame containing SNP information (chromosome, position, etc.)

Users may save these objects with: `saveRDS(input$fbm, file="input.rds")`
and `write.table(input$info, file="input_info.txt")`

### 3. Run HAUDI/GAUDI models

Users may run HAUDI with the `haudi` function.
An example may look like:

```{r}
pheno_file <- system.file("extdata/toy.pheno", package = "HAUDI")
y <- read.csv(pheno_file, sep = "\t")$phenotype
result <- haudi(
  fbm = input$fbm, fbm_info = input$info, y_train = y,
  gamma_vec = seq(1, 2, 0.2), family = "gaussian")
```

To obtain a data frame with ancestry-specific effect estimates, use the
helper function `get_beta_haudi`.

The `HAUDI` package also provides provides a function to run GAUDI,
using the same input data as for HAUDI.
