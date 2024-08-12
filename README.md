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

### 1. Prepare input VCF file

HAUDI requires as input a tabix-indexed vcf file, with AN1 and AN2 subfields
for haplotype local ancestry (as is produced by [flare](https://github.com/browning-lab/flare)).
If you use flare to estimate local ancestry, you can annotate your original vcf
file using a command like: `bcftools annotate -c FORMAT -a flare.anc.vcf.gz
target.vcf.gz -Oz -o target.anc.vcf.gz`.

This will produce a VCF file where the AN1, AN2 fields are missing for variants
not included in the reference panel. To interpolate local ancestry in
the following step,ensure that the VCF file is sorted.

### 2. Creating file-backed matrices

HAUDI uses file-backed matrices, implemented in the `bigstatsr`
package to store genotype/ancestry data. To create this file-backed
matrix, use the function `make_fbm`, which takes the VCF file prepared
in the previous step as input. An example may be:

```{r}
fbm_result <- make_fbm(
  vcf_file = "target.anc.vcf.gz",
  fbm_pref = "target",
  chunk_size = 400,
  min_ac = 10,
  geno_format = "GT",
  anc_names = c("Pop_01", "Pop_02", "Pop_03") 
)
```

This command would return a list containing:

1) an object of class FBM.code256 containing genotype/ancestry data
2) a data frame containing SNP information (chromosome, position, etc.)

Users may save these objects with: `saveRDS(fbm_result$FBM, file="target.rds")`
and `write.table(fbm_result$info, file="target_info.txt")`

### 3. Run HAUDI/GAUDI models

Users may run HAUDI with the `haudi` function.
An example may look like:

```{r}
haudi_model <- haudi(
  fbm_obj = fbm_result$FBM,
  fbm_info = fbm_result$info,
  y = y,
  gamma = 2,
  ind_train = NULL,
  family = "gaussian",
  snps = NULL,
  K = 10
)
```

To obtain a data frame with ancestry-specific effect estimates, use the
helper function `get_beta_haudi`.

The `HAUDI` package also provides provides a function to run GAUDI,
using the same input data as for HAUDI.
