# Copyright 2020, Georgetown University and University of Pittsburgh.  All Rights Reserved.
#
#============================================================================
# simpleM_Ex.R
# from http://simplem.sourceforge.net/
#============================================================================
# License:  GPL version 2 or newer.
# NO warranty.

#============================================================================
# citation:
#
# Gao X, Starmer J and Martin ER (2008) A Multiple Testing Correction Method for
# Genetic Association Studies Using Correlated Single Nucleotide Polymorphisms.
# Genetic Epidemiology 32:361-369
#
# Gao X, Becker LC, Becker DM, Starmer J, Province MA (2009) Avoiding the high
# Bonferroni penalty in genome-wide association studies. Genetic Epidemiology
# (Epub ahead of print)

#============================================================================
# readme:
# example SNP file format:
# row => SNPs
# column => Unrelated individuals

# The data file should contain only POLYMORPHIC SNPs.

# Missing values should be imputed.
# There should be NO missing values in the SNP data file.
# SNPs are coded as 0, 1 and 2 for the number of reference alleles.
# SNPs are separated by one-character spaces.

# You may need to change file path (search for "fn_In" variable)
# depending on where your snp file is stored at.

#============================================================================
# Meff through the PCA approach
# use a part of the eigen values according to how much percent they contribute
# to the total variation
Meff_PCA <- function(eigenValues, percentCut){
  totalEigenValues <- sum(eigenValues)
  myCut <- percentCut*totalEigenValues
  num_Eigens <- length(eigenValues)
  myEigenSum <- 0
  index_Eigen <- 0

  for(i in 1:num_Eigens){
    if(myEigenSum <= myCut){
      myEigenSum <- myEigenSum + eigenValues[i]
      index_Eigen <- i
    }
    else{
      break
    }
  }
  return(index_Eigen)
}

#============================================================================
# infer the cutoff => Meff
#' @importFrom stats cor
inferCutoff <- function(dt_My){
  CLD <- cor(dt_My)
  eigen_My <- eigen(CLD)

  # PCA approach
  eigenValues_dt <- abs(eigen_My$values)
  Meff_PCA_gao <- Meff_PCA(eigenValues_dt, PCA_cutoff)
  return(Meff_PCA_gao)
}

#============================================================================
PCA_cutoff <- 0.995

#============================================================================
# fix length, simpleM
#fn_In <- "D:/simpleM_Ex/snpSample.txt"				# <---- change path here!!!
#mySNP_nonmissing <- read.table(fn_In, colClasses="integer")
#' M_GAO
#'
#' Compute the effective number of independent SNPs in a region.
#'
#' @param SNP_mx A matrix of polymorphic SNPs (coded 0, 1, 2) with SNPs in rows and individuals in columns.
#'
#' @export
#'
#' @references
#' Gao X, Starmer J and Martin ER (2008) A Multiple Testing Correction Method for Genetic Association Studies Using Correlated Single Nucleotide Polymorphisms. Genetic Epidemiology 32:361-369
#'
#' @source http://simplem.sourceforge.net/
#'
M_GAO = function(SNP_mx){

  mySNP_nonmissing <- t(SNP_mx)
  numLoci <- length(mySNP_nonmissing[, 1])

  simpleMeff <- NULL
  fixLength <- 133
  i <- 1
  myStart <- 1
  myStop <- 1
  while(myStop < numLoci){
    myDiff <- numLoci - myStop
    if(myDiff <= fixLength) break

    myStop <- myStart + i*fixLength - 1
    snpInBlk <- t(mySNP_nonmissing[myStart:myStop, ])
    MeffBlk <- inferCutoff(snpInBlk)
    simpleMeff <- c(simpleMeff, MeffBlk)
    myStart <- myStop+1
  }
  snpInBlk <- t(mySNP_nonmissing[myStart:numLoci, ])
  MeffBlk <- inferCutoff(snpInBlk)
  simpleMeff <- c(simpleMeff, MeffBlk)

  #cat("Total number of SNPs is: ", numLoci, "\n")
  #cat("Inferred Meff is: ", sum(simpleMeff), "\n")

  return(sum(simpleMeff)) # Inferred Meff


}
#============================================================================
# end

#' dRule
#'
#' This function applies the dynamic rule to determine the number of basis functions to use
#'
#' @param geno.only The input matrix of SNP genotypes, coded 0, 1, 2.
#'
#' @export
#'
dRule = function(geno.only){
  # SNP matrix
  cond <- colSums(geno.only) > 0
  SNP_mx = geno.only[ ,cond]

  # dynamic rule for number of basis function use M_GAO #added by CY, Nov 2019
  M_gao = M_GAO(SNP_mx)

  if (M_gao < 19)             {betabasis_Bsp = genobasis_Bsp = 6; betabasis_Fsp = genobasis_Fsp = 7}
  if (M_gao > 18 & M_gao <25) {betabasis_Bsp = genobasis_Bsp = 8; betabasis_Fsp = genobasis_Fsp = 9}
  if (M_gao > 24 & M_gao <31) {betabasis_Bsp = genobasis_Bsp = 10; betabasis_Fsp = genobasis_Fsp = 11}
  if (M_gao > 30 & M_gao <37) {betabasis_Bsp = genobasis_Bsp = 12; betabasis_Fsp = genobasis_Fsp = 13}
  if (M_gao > 36 & M_gao <43) {betabasis_Bsp = genobasis_Bsp = 14; betabasis_Fsp = genobasis_Fsp = 15}
  if (M_gao > 42)             {betabasis_Bsp = genobasis_Bsp = 16; betabasis_Fsp = genobasis_Fsp = 17}

  out = list(betabasis_Bsp=betabasis_Bsp, genobasis_Bsp=genobasis_Bsp, betabasis_Fsp=betabasis_Fsp, genobasis_Fsp=genobasis_Fsp, M_gao=M_gao)
  return(out)

}
