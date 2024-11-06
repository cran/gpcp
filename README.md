# Genomic Prediction of Cross Performance (`gpcp`)

## Overview

`gpcp` is an R package that performs genomic prediction of cross performance using both genotype and phenotype data. The package supports diploid and polyploid species and processes data in several steps including loading necessary software, converting genotype data, processing phenotype data, fitting mixed models, and predicting cross performance based on weighted marker effects.

The package uses the `sommer`, `dplyr`, and `AGHmatrix` R packages for mixed model analysis and genomic data processing.

## Installation

To install `gpcp` directly from GitHub, use the `devtools` package:

1. Install the `devtools` package if you don't have it:

   ```r
   install.packages("devtools")
   devtools::install_github("cmn92/gpcp")
2. Load the package:
   ```r
   library(gpcp)

## Usage

Here is an example of how to use the runGPCP function for genomic prediction of cross performance:

  ```r
  # Load phenotype data from a CSV file
phenotypeFile <- read.csv("~/Documents/GCPC_input_files/2020_TDr_PHENO (1).csv")

# Specify the genotype file path (VCF or HapMap format)
genotypeFile <- "~/Documents/GCPC_input_files/genotypeFile.vcf"

# Define necessary inputs
genotypes <- "Accession"  # Column name for genotype IDs in phenotype data
traits <- c("rAUDPC_YMV", "YIELD", "DMC")  # List of traits to predict
weights <- c(0.2, 3, 1)  # Weights corresponding to traits
userFixed <- c("LOC", "REP")  # Fixed effects variables
Ploidy <- 2  # Ploidy level of the organism
NCrosses <- 150  # Number of top crosses to output

# Run genomic prediction of cross performance
finalcrosses <- runGPCP(
    phenotypeFile = phenotypeFile,
    genotypeFile = genotypeFile,
    genotypes = genotypes,
    traits = paste(traits, collapse = ","),
    weights = weights,
    userFixed = paste(userFixed, collapse = ","),
    Ploidy = Ploidy,
    NCrosses = NCrosses
)

# View the predicted crosses
print(finalcrosses)
```
## Input Arguments
phenotypeFile: A data frame containing phenotypic data, typically read from a CSV file.

genotypeFile: A file path to the genotypic data, either in VCF format or as a HapMap.

genotypes: A character string representing the column name in the phenotype file that corresponds to the genotype IDs.

traits: A string of comma-separated trait names from the phenotype file, which will be used for genomic prediction.

weights: A numeric vector specifying the weights for the traits. The order of weights should correspond to the order of traits.

userFixed: A string of comma-separated fixed effect variables from the phenotype file.

Ploidy: An integer representing the ploidy level of the organism.

NCrosses: An integer specifying the number of top crosses to output.


## Output
The runGPCP function returns a data frame with predicted crosses, including:

Parent1: The first parent genotype ID.

Parent2: The second parent genotype ID.

CrossPredictedMerit: The predicted merit of the cross.

## Dependencies
sommer: Mixed models for genome-wide prediction and association studies.

dplyr: Data manipulation package.

AGHmatrix: Implements additive and dominance genomic relationship matrices.

You can install these dependencies with:

```r
install.packages(c("sommer", "dplyr", "AGHmatrix"))
```
## References
Xiang, J., et al. (2016). "Mixed Model Methods for Genomic Prediction." Nature Genetics.

Batista, L., et al. (2021). "Genetic Prediction and Relationship Matrices." Theoretical and Applied Genetics.

## License
This project is licensed under the MIT License - see the LICENSE file for details.











   

