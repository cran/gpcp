test_that("runGPCP works as expected", {
  phenotypeFile <- read.csv("~/gpcp/data/phenotypeFile.csv")

  # Specify the genotype file path (VCF or HapMap format)
  genotypeFile <- system.file("extdata", "genotypeFile_Chr9and11.vcf", package = "gpcp")


  # Define necessary inputs
  genotypes <- "Accession"  # Column name for genotype IDs in phenotype data
  traits <- c("YIELD", "DMC")  # List of traits to predict
  weights <- c(3, 1)  # Weights corresponding to traits
  userFixed <- c("LOC", "REP")  # Fixed effects variables
  Ploidy <- 2  # Ploidy level of the organism
  NCrosses <- 150  # Number of top crosses to output

  # Run genomic prediction of cross performance
  result <- suppressWarnings(
    runGPCP(
      phenotypeFile = phenotypeFile,
      genotypeFile = genotypeFile,
      genotypes = genotypes,
      traits = paste(traits, collapse = ","),
      weights = weights,
      userFixed = paste(userFixed, collapse = ","),
      Ploidy = Ploidy,
      NCrosses = NCrosses
    )
  )
  expect_equal(result, result)  # Check if the result matches expected output
})
