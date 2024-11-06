################################################################################
# Genomic prediction of cross performance for YamBase
################################################################################

#Authors: Marlee Labroo, Christine Nyaga, Lukas Mueller

# There are ten main steps to this protocol:
# 1. Load the software needed.
# 2. Declare user-supplied variables.
# 3. Read in the genotype data and convert to numeric allele counts.
# 4. Get the genetic predictors needed.
# 5. Process the phenotypic data.
# 6. Fit the mixed models in sommer.
# 7. Backsolve from individual estimates to marker effect estimates / GBLUP -> RR-BLUP
# 8. Weight the marker effects and add them together to form an index of merit.
# 9. Predict the crosses.
# 10. Format the information needed for output.

# Declare global variables to avoid R CMD check NOTE
utils::globalVariables(c("P1Sex", "P2Sex"))

#' @title Genomic Prediction of Cross Performance
#' This function performs genomic prediction of cross performance using genotype and phenotype data.
#'
#' @param phenotypeFile A data frame containing phenotypic data, typically read from a CSV file.
#' @param genotypeFile Path to the genotypic data, either in VCF or HapMap format.
#' @param genotypes A character string representing the column name in the phenotype file for the genotype IDs.
#' @param traits A string of comma-separated trait names from the phenotype file.
#' @param weights A numeric vector specifying weights for the traits.
#' @param userSexes A string representing the column name corresponding to the individuals' sexes.
#' @param userFixed A string of comma-separated fixed effect variables.
#' @param userRandom A string of comma-separated random effect variables.
#' @param Ploidy An integer representing the ploidy level of the organism.
#' @param NCrosses An integer specifying the number of top crosses to output.
#' @return A data frame containing predicted cross performance.
#' @export
#' @useDynLib gpcp, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom tools file_ext
#' @importFrom magrittr %>%
#' @importFrom VariantAnnotation readVcf
#' @importFrom VariantAnnotation genotypeToSnpMatrix
#' @importFrom methods as
#' @importFrom stats as.formula filter na.omit sd
#' @importFrom utils combn head read.delim
#' @examples
#' # Load phenotype data from CSV
#' phenotypeFile <- read.csv(system.file("extdata", "phenotypeFile.csv", package = "gpcp"))
#' genotypeFile <- system.file("extdata", "genotypeFile_Chr9and11.vcf", package = "gpcp")
#' finalcrosses <- runGPCP(
#'     phenotypeFile = phenotypeFile,
#'     genotypeFile = genotypeFile,
#'     genotypes = "Accession",
#'     traits = "YIELD,DMC",
#'     weights = c(3, 1),
#'     userFixed = "LOC,REP",
#'     Ploidy = 2,
#'     NCrosses = 150
#' )
#' message(finalcrosses)
runGPCP = function(phenotypeFile, genotypeFile, genotypes, traits,
                   weights = NA, userSexes = "", userFixed = NA, userRandom = NA, Ploidy = NA, NCrosses = NA) {

  # Rcpp::sourceCpp("~/gpcp/R/QuantGenResources/CalcCrossMeans.cpp") # this is called CalcCrossMean.cpp on Github

  if (!requireNamespace("snpStats", quietly = TRUE)) {
    stop("The 'snpStats' package is required but not installed. Please install it using BiocManager::install('snpStats').")
  }
  ################################################################################
  # 2. Declare user-supplied variables
  ################################################################################

  # a. Define path with internal YamBase instructions such that the object 'userGeno'
  #    is defined as a VCF file of genotypes.

  userGeno <- genotypeFile


  # b. Define path2 with internal YamBase instructions such that the object 'userPheno'
  #    is defined as the phenotype file.
  userPheno <- phenotypeFile

  # c. The user should be able to select their fixed variables from a menu
  #    of the column names of the userPheno object. The possible interaction terms
  #    also need to be shown somehow. Then, those strings should be passed
  #    to this vector, 'userFixed'. Please set userFixed to NA if no fixed effects
  #    besides f are requested.
  #    f is automatically included as a fixed effect- a note to the user would be good.

  # userFixed <- c()
  # userFixed <- c("studyYear") # for testing only
  userFixed <- unlist(strsplit(userFixed, split = ",", fixed = T))


  # d. The user should be able to select their random variables from a menu
  #    of the column names of the userPheno object. The possible interaction terms
  #    also need to be shown somehow. Then, those strings should be passed
  #    to this vector, 'userRandom'.

  # userRandom <- c()
  # userRandom <- "blockNumber" # for testing only
  if(is.na(userRandom)){
    userRandom=NA
  } else {
    userRandom <- unlist(strsplit(userRandom, split = ",", fixed = T))
  }
  # e. The user should be able to indicate which of the userPheno column names
  #    represents individual genotypes identically as they are represented in the VCF
  #    column names. No check to ensure matching at this stage. This single string
  #    should be passed to this vector, userID.

  userID <- genotypes
  #userID <- "germplasmName" # for testing only


  # f. The user must indicate the ploidy level of their organism, and the integer
  #    provided should be passed to the vector 'userPloidy'. CalcCrossMeans.cpp
  #    currently supports ploidy = {2, 4, 6}. Ideally, the user could select
  #    their ploidy from a drop-down to avoid errors here, and there would be a note
  #    that other ploidies are not currently supported. If not, a possible error is
  #    provided.

  userPloidy <- Ploidy
  userPloidy <- 2 # for testing only

  # if(userPloidy %in% c(2, 4, 6) != TRUE){
  #   stop("Only ploidies of 2, 4, and 6 are supported currently. \n
  #        Please confirm your ploidy level is supported.")
  # }


  # g. The user should be able to select their response variables from a drop-down menu
  #    of the column names of the userPheno object. Then, those strings should be passed
  #    to this vector, 'userResponse'.

  # userResponse <- c()
  # userResponse <- c("YIELD", "DMC", "OXBI") # for testing only
  userResponse <- unlist(strsplit(traits, split = ",", fixed = T))

  # h. The user must indicate weights for each response. The order of the vector
  #    of response weights must match the order of the responses in userResponse.

  userWeights <- weights
  # userWeights <- c(1, 0.8, 0.2) # for YIELD, DMC, and OXBI respectively; for testing only
  # userWeights <- as.numeric(unlist(strsplit(weights, split = ",", fixed = T)))

  # i. The user can indicate the number of crosses they wish to output.
  #    The maximum possible is a full diallel.

  # userNCrosses <- c()
  userNCrosses <- NCrosses # for testing only


  # j. The user can (optionally) input the individuals' sexes and indicate the column
  #    name of the userPheno object which corresponds to sex. The column name
  #   string should be passed to the 'userSexes' object. If the user does not wish
  #   to remove crosses with incompatible sexes (e.g. because the information is not available),
  #   then userSexes should be set to NA.


  # userSexes <- c()

  # userSexes <- "Sex" # for testing only
  # userPheno$Sex <- sample(c("M", "F"), size = nrow(userPheno), replace = TRUE, prob = c(0.7, 0.3)) # for testing only
  # Please note that for the test above, sex is sampled randomly for each entry, so the same accession can have
  # different sexes. This does not matter for the code or testing.

  ################################################################################
  # 3. Read in the genotype data and convert to numeric allele counts.
  ################################################################################

  # a. The VCF file object 'userGeno' needs to be converted to a numeric matrix
  #    of allele counts in whic:
  #    Rownames represent the individual genotype IDs
  #    Colnames represent the site IDs
  #    A cell within a given row and column represents the row individual's
  #    genotype at the site in the column.

  #   The individual's genotype should be an integer from 0... ploidy to represent
  #   counts of the alternate allele at the site. Diploid example:
  #    0 = homozygous reference
  #    1 = heterozygous
  #    2 = homozygous alternate

  #    The genotypes must not contain monomorphic or non-biallelic sites.
  #    Users need to pre-process their VCF to remove these (e.g. in TASSEL or R)
  #    I can put an error message into this script if a user tries to input
  #    monomorphic or biallelic sites which could be communicated through the GUI.
  #    It's also possible to filter them here.

  if (file_ext(genotypeFile) == "vcf") {
    message("READING VARIANT FILE ")
    #  Import VCF with VariantAnnotation package and extract matrix of dosages
    myVCF <- VariantAnnotation::readVcf(genotypeFile)
    # G <- t(geno(myVCF)$DS) # Individual in row, genotype in column
    mat <- VariantAnnotation::genotypeToSnpMatrix(myVCF)
    # G <- t(geno(myVCF)$DS) # Individual in row, genotype in column
    G <- methods::as(mat$genotypes, "numeric")
    G <- G[, colSums(is.na(G)) < nrow(G)]

    #   TEST temporarily import the genotypes via HapMap:
    # source("R/hapMap2numeric.R") # replace and delete
    # G <- hapMap2numeric(genotypeFile) # replace and delete
  } else {
    # accession_names     abc      abc2    abc3
    # marker1                   0      0        2
    # marker2                   1      0        0
    # marker3                   0      0        0

    message("READING DOSAGE FILE ")
    GF <- utils::read.delim(genotypeFile)
    GD <- GF[, -1]
    GM <- as.matrix(GD)
    G <- t(GM)
  }

  # message("G Matrix start --------")
  # message(G[1:5, 1:5])
  # message("G Matrix end =========")

  ################################################################################
  # 4. Get the genetic predictors needed.
  ################################################################################

  message("GENETIC PREDICTIONS...")
  # 4a. Get the inbreeding coefficent, f, as described by Xiang et al., 2016
  # The following constructs f as the average heterozygosity of the individual
  # The coefficient of f estimated later then needs to be divided by the number of markers
  # in the matrix D before adding it to the estimated dominance marker effects
  # One unit of change in f represents changing all loci from homozygous to heterozygous

  ### GC <- G - (userPloidy/2) #this centers G
  GC <- G * (userPloidy - G) * (2 / userPloidy)^2 # center at G
  f <- rowSums(GC, na.rm = TRUE) / apply(GC, 1, function(x) sum(!is.na(x)))

  # Another alternate way to construct f is the total number of heterozygous loci in the individual
  # The coefficient of this construction of f does not need to be divided by the number of markers
  # It is simply added to each marker dominance effect
  # The coefficient of this construction of f represents the average dominance effect of a marker
  # One unit of change in f represents changing one locus from homozygous to heterozygous
  # f <- rowSums(D, na.rm = TRUE)


  message("DISTANCE MATRIX...")
  # 4b. Get the additive and dominance relationship matrices following Batista et al., 2021
  # https://doi.org/10.1007/s00122-021-03994-w

  # Additive: this gives a different result than AGHmatrix VanRaden's Gmatrix
  # AGHmatrix: Weights are implemented for "VanRaden" method as described in Liu (2020)?
  allele_freq <- colSums(G) / (userPloidy * nrow(G))
  W <- t(G) - userPloidy * allele_freq
  WWt <- crossprod(W)
  denom <- sum(userPloidy * allele_freq * (1 - allele_freq))
  A <- WWt / denom

  # Check with paper equation:
  # w <- G - (userPloidy/2)
  # num <- w %*% t(w)
  # denom = sum(userPloidy * allele_freq * (1 - allele_freq))
  # A2 <- num/denom
  # table(A == A2)
  # cor(as.vector(A), as.vector(A2)) # 0.9996...


  # Dominance or digenic dominance
  if (userPloidy == 2) {
    D <- AGHmatrix::Gmatrix(G, method = "Su", ploidy = userPloidy, missingValue = NA)
  }

  if (userPloidy > 2) {
    # Digenic dominance
    C_matrix <- matrix(length(combn(userPloidy, 2)) / 2,
                       nrow = nrow(t(G)),
                       ncol = ncol(t(G))
    )

    Ploidy_matrix <- matrix(userPloidy,
                            nrow = nrow(t(G)),
                            ncol = ncol(t(G))
    )

    Q <- (allele_freq^2 * C_matrix) -
      (Ploidy_matrix - 1) * allele_freq * t(G) +
      0.5 * t(G) * (t(G) - 1)

    Dnum <- crossprod(Q)
    denomDom <- sum(C_matrix[, 1] * allele_freq^2 * (1 - allele_freq)^2)
    D <- Dnum / denomDom
  }


  ################################################################################
  # 5. Process the phenotypic data.
  ################################################################################

  # write(summary(userPheno), stderr())

  # a. Paste f into the phenotype dataframe
  message("processing phenotypic data...")
  userPheno$f <- f[as.character(userPheno[, userID])]


  # b. Scale the response variables.
  for (i in 1:length(userResponse)) {
    userPheno[, userResponse[i]] <- (userPheno[, userResponse[i]] - mean(userPheno[, userResponse[i]], na.rm = TRUE)) / sd(userPheno[, userResponse[i]], na.rm = TRUE)
  }

  # c. Paste in a second ID column for the dominance effects.
  dominanceEffectCol <- paste(userID, "2", sep = "")

  userPheno[, dominanceEffectCol] <- userPheno[, userID]

  uniq <- length(sapply(lapply(userPheno, unique), length))


  # Additional steps could be added here to remove outliers etc.

  ################################################################################
  # 6. Fit the mixed models in sommer.
  ################################################################################

  message("Fitting mixed model in sommer")
  # 6a. Make a list to save the models.

  userModels <- list()

  for (i in 1:length(userResponse)) {
    message(paste("User response: ", userResponse[i]))
    # check if fixed effects besides f are requested, then paste together
    # response variable and fixed effects
    if (!is.na(userFixed[1])) {
      fixedEff <- paste(userFixed, collapse = " + ")
      fixedEff <- paste(fixedEff, "f", sep = " + ")
      fixedArg <- paste(userResponse[i], " ~ ", fixedEff, sep = "")
    }
    if (is.na(userFixed[1])) {
      fixedArg <- paste(userResponse[i], " ~ ", "f")
    }


    # check if random effects besides genotypic additive and dominance effects
    # are requested, then paste together the formula

    message("Generating formula...")

    if (!is.na(userRandom[1])) {
      randEff <- paste(userRandom, collapse = " + ")
      ID2 <- paste(userID, 2, sep = "")
      randEff2 <- paste("~sommer::vsr(", userID, ", Gu = A) + sommer::vsr(", ID2, ", Gu = D)", sep = "")
      randArg <- paste(randEff2, randEff, sep = " + ")
    }
    if (is.na(userRandom[1])) {
      ID2 <- paste(userID, 2, sep = "")
      randArg <- paste("~sommer::vsr(", userID, ", Gu = A) + sommer::vsr(", ID2, ", Gu = D)", sep = "")
    }

    message(paste("Fit mixed GBLUP model...", randArg))

    #  write(paste("USER PHENO:", userPheno), stderr())
    #  write(paste("COLNAMES: ", colnames(userPheno)), stderr())
    # fit the mixed GBLUP model
    myMod <- sommer::mmer(
      fixed = stats::as.formula(fixedArg),
      random = stats::as.formula(randArg),
      rcov = ~units,
      getPEV = FALSE,
      data = userPheno
    )


    # save the fit model


    userModels[[i]] <- myMod
  }


  ######################################################################################
  # 7. Backsolve from individual estimates to marker effect estimates / GBLUP -> RR-BLUP
  ######################################################################################

  # a. Get the matrices and inverses needed
  #    This is not correct for polyploids yet.
  A.G <- G - (userPloidy / 2) # this is the additive genotype matrix (coded -1 0 1 for diploids)
  D.G <- 1 - abs(A.G) # this is the dominance genotype matrix (coded 0 1 0 for diploids)


  A.T <- A.G %*% t(A.G) ## additive genotype matrix
  # inverse; may cause an error sometimes, if so, add a small amount to the diag
  epsilon <- 1e-8  # A small value to add to the diagonal
  # Try to invert the matrix
  A.Tinv <- tryCatch({
    solve(A.T)  # Try solving normally
  }, error = function(e) {
    # If there is an error (like singular matrix), add epsilon to diagonal and retry
    warning("Matrix is singular; adding small value to diagonal and retrying inversion.")
    A.T.reg <- A.T + diag(epsilon, nrow(A.T))
    solve(A.T.reg)  # Solve the regularized matrix
  })
  A.TTinv <- t(A.G) %*% A.Tinv # M'%*% (M'M)-

  D.T <- D.G %*% t(D.G) ## dominance genotype matrix
  ## inverse
  D.Tinv <- tryCatch({
    solve(D.T)  # Try solving normally
  }, error = function(e) {
    # If there is an error (like singular matrix), add epsilon to diagonal and retry
    warning("Matrix is singular; adding small value to diagonal and retrying inversion.")
    D.T.reg <- D.T + diag(epsilon, nrow(D.T))
    solve(D.T.reg)  # Solve the regularized matrix
  })
  D.TTinv <- t(D.G) %*% D.Tinv # M'%*% (M'M)-


  # b. Loop through and backsolve to marker effects.

  userAddEff <- list() # save them in order
  userDomEff <- list() # save them in order

  for (i in 1:length(userModels)) {
    myMod <- userModels[[i]]

    # get the additive and dominance effects out of the sommer list
    subMod <- myMod$U
    subModA <- subMod[[1]]
    subModA <- subModA[[1]]
    subModD <- subMod[[2]]
    subModD <- subModD[[1]]

    # backsolve
    addEff <- A.TTinv %*% matrix(subModA[colnames(A.TTinv)], ncol = 1) # these must be reordered to match A.TTinv
    domEff <- D.TTinv %*% matrix(subModD[colnames(D.TTinv)], ncol = 1) # these must be reordered to match D.TTinv

    # add f coefficient back into the dominance effects
    subModf <- myMod$Beta
    fCoef <- subModf[subModf$Effect == "f", "Estimate"] # raw f coefficient
    fCoefScal <- fCoef / ncol(G) # divides f coefficient by number of markers
    dirDomEff <- domEff + fCoefScal

    # save
    userAddEff[[i]] <- addEff
    userDomEff[[i]] <- dirDomEff
  }






  ################################################################################
  # 8. Weight the marker effects and add them together to form an index of merit.
  ################################################################################

  ai <- 0
  di <- 0
  for (i in 1:length(userWeights)) {

    ai <- ai + userAddEff[[i]] * userWeights[i]
    di <- di + userDomEff[[i]] * userWeights[i]
  }






  ################################################################################
  # 9. Predict the crosses.
  ################################################################################

  # If the genotype matrix provides information about individuals for which
  # cross prediction is not desired, then the genotype matrix must be subset
  # for use in calcCrossMean(). calcCrossMean will return predicted cross
  # values for all individuals in the genotype file otherwise.

  message("Predict crosses...")

  GP <- G[rownames(G) %in% userPheno[, userID], ]

  message("GP:")
  # message(head(GP))

  crossPlan <- calcCrossMean(
    GP,
    ai,
    di,
    userPloidy
  )


  message("Done with calcCrossMean!!!!!!")


  ################################################################################
  # 10. Format the information needed for output.
  ################################################################################

  # Add option to remove crosses with incompatible sexes.


  # hash <- new.env(hash = TRUE, parent = emptyenv(), size = 100L)

  # assign_hash(userPheno$germplasmName, userPheno$userSexes, hash)

  if (userSexes != "") { # "plant sex estimation 0-4"
    # !is.na(userSexes)  && !is.na(sd(userPheno[, userSexes]))

    # Reformat the cross plan
    crossPlan <- as.data.frame(crossPlan)
    crossPlan <- crossPlan[order(crossPlan[, 3], decreasing = TRUE), ] # orders the plan by predicted merit
    crossPlan[, 1] <- rownames(GP)[crossPlan[, 1]] # replaces internal ID with genotye file ID
    crossPlan[, 2] <- rownames(GP)[crossPlan[, 2]] # replaces internal ID with genotye file ID
    colnames(crossPlan) <- c("Parent1", "Parent2", "CrossPredictedMerit")

    # Look up the parent sexes and subset
    crossPlan$P1Sex <- userPheno[match(crossPlan$Parent1, userPheno$germplasmName), userSexes] # get sexes ordered by Parent1

    crossPlan$P2Sex <- userPheno[match(crossPlan$Parent2, userPheno$germplasmName), userSexes] # get sexes ordered by Parent2

    col_repl <- c("P1Sex", "P2Sex")
    crossPlan %>% dplyr::filter(P1Sex == 0 | P2Sex == 0) # remove the 0s
    crossPlan %>% dplyr::filter(P1Sex == 1 & P2Sex == 1) # remove same sex crosses with score of 1
    crossPlan %>% dplyr::filter(P1Sex == 2 & P2Sex == 2) # remove same sex crosses with score of 2
    # crossPlan <- crossPlan[crossPlan$P1Sex != crossPlan$P2Sex, ] # remove crosses with same-sex parents

    ## replace plant sex numbers to male, female etc

    crossPlan[col_repl] <- sapply(crossPlan[col_repl], function(x) replace(x, x %in% "NA", "NA"))
    crossPlan[col_repl] <- sapply(crossPlan[col_repl], function(x) replace(x, x %in% 1, "Male"))
    crossPlan[col_repl] <- sapply(crossPlan[col_repl], function(x) replace(x, x %in% 2, "Female"))
    crossPlan[col_repl] <- sapply(crossPlan[col_repl], function(x) replace(x, x %in% 3, "Monoecious male (m>f)"))
    crossPlan[col_repl] <- sapply(crossPlan[col_repl], function(x) replace(x, x %in% 4, "Monoecious female(f>m)"))




    # subset the number of crosses the user wishes to output
    if (nrow(crossPlan)<100) {
      finalcrosses = crossPlan
    } else {
      crossPlan[1:userNCrosses, ]
      finalcrosses=crossPlan[1:userNCrosses, ]
    }

  } else {
    # only subset the number of crosses the user wishes to output
    crossPlan <- as.data.frame(crossPlan)
    crossPlan <- na.omit(crossPlan)
    crossPlan <- crossPlan[order(crossPlan[, 3], decreasing = TRUE), ] # orders the plan by predicted merit
    crossPlan[, 1] <- rownames(GP)[crossPlan[, 1]] # replaces internal ID with genotye file ID
    crossPlan[, 2] <- rownames(GP)[crossPlan[, 2]] # replaces internal ID with genotye file ID
    colnames(crossPlan) <- c("Parent1", "Parent2", "CrossPredictedMerit")


    ## save the best 100 predictions
    if (nrow(crossPlan)<100) {
      finalcrosses = crossPlan
    } else {
      crossPlan[1:userNCrosses, ]
      finalcrosses=crossPlan[1:userNCrosses, ]
    }
  }

  return(finalcrosses)
}
