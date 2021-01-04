# Copyright 2020, Georgetown University and University of Pittsburgh.  All Rights Reserved.
#
#' Example data for the PedGLMM package
#'
#' @name exampleData
#'
#' @docType data
#'
#' @usage data(exampleData)
#'
#' @format Four data frames: Ped, geno, cov, snpPos
#'
#' @return None
#'
#' @seealso \code{\link{Ped}}, \code{\link{geno}}, \code{\link{cov}}, \code{\link{snpPos}}
#'
#' @examples
#' data(exampleData)
#' dim(Ped)
#' head(Ped)
#' dim(geno)
#' head(geno[,1:10])
#' dim(cov)
#' head(cov)
#' dim(snpPos)
#' head(snpPos)
NULL
#'
#' Ped
#'
#' Example pedigree data frame, available via data(exampleData)
#'
#' A data frame containing the pedigree information with the following columns:
#' \describe{
#'   \item{ID}{Person ID}
#'   \item{ped}{pedigree ID,character or numeric allowed.}
#'   \item{person}{person ID, a unique ID within each pedigree, numeric or character allowed.}
#'   \item{father}{father ID, NA if no father.}
#'   \item{mother}{mother ID,  NA if no mother.}
#'   \item{sex}{sex, coded as 1 for male, 2 for female.}
#'   \item{trait}{trait phenotype, either case-control status coded as 1 for affected and 0 for unaffected. Subjects with missing (NA) will be removed from the analysis.}
#' }
#'
#' @seealso \code{\link{exampleData}}, \code{\link{geno}}, \code{\link{cov}}, \code{\link{snpPos}}

"Ped"
#'
#' geno
#'
#' Example genotype data frame, available via data(exampleData).
#'
#' A data frame containing the genotype information.  This is a matrix with genotypes for subjects (rows) at each variant position (columns). The first two columns are required to be named “ped” and “person”, which are used to match subjects to their data in the pedigree data.frame. The genotypes are coded as 0, 1, 2 for autosomal markers (typically a count of the number of the minor alleles).
#'
#'
#' @seealso \code{\link{Ped}}, \code{\link{exampleData}}, \code{\link{cov}}, \code{\link{snpPos}}
"geno"
#'
#' cov
#'
#' Example covariate data frame, available via data(exampleData).
#'
#' A data frame containing the covariate information. The first two columns are required to be named “ped” and “person”, which are used to match subjects to their data in the pedigree data frame.
#'
#' @seealso \code{\link{Ped}}, \code{\link{geno}}, \code{\link{exampleData}}, \code{\link{snpPos}}
"cov"
#'
#' snpPos
#'
#' Example marker position data frame, available via data(exampleData).
#'
#' This data frame provides marker positions for each SNP. The first column, chr, contains
#' the chromosome number, the second column, snp, contains the SNP name, and the third
#' column, pos, contains the position of the SNP in base pairs.
#'
#' @seealso \code{\link{Ped}}, \code{\link{geno}}, \code{\link{cov}}, \code{\link{exampleData}}

"snpPos"
