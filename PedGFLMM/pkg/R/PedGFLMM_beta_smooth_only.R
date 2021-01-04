# Copyright 2020, Georgetown University and University of Pittsburgh.  All Rights Reserved.
#
### Updated by Chi-Yang on 03/19/2019: block distribution = "binomial", make the routine work for no covariate
### Written by Chi-Yang Chiu and Ruzong Fan, May 2016
### This function use r function "pedigreemm" to fit GFLMM for pedigree data
### need to provide ped_list that can be created by r function "pedigree" from library "pedigreemm"

# library(fda)
# library(MASS)
# library(Matrix)
# library(nlme)
# library(pedigreemm)



#' PedGFLMM_beta_smooth_only
#'
#' Computes the PedGFLMM statistics under the beta smooth only model.
#'
#' @param ped A data frame containing the pedigree information with the following columns:
#' \describe{
#'   \item{ID}{Person ID}
#'   \item{ped}{pedigree ID,character or numeric allowed.}
#'   \item{person}{person ID, a unique ID within each pedigree, numeric or character allowed.}
#'   \item{father}{father ID, NA if no father.}
#'   \item{mother}{mother ID,  NA if no mother.}
#'   \item{sex}{sex, coded as 1 for male, 2 for female.}
#'   \item{trait}{trait phenotype, case-control status coded as 1 for affected and 0 for unaffected. Subjects with missing (NA) will be removed from the analysis.}
#' }
#' @param geno A data frame containing the genotype information.  This is a matrix with genotypes for subjects (rows) at each variant position (columns). The first two columns are required to be named “ped” and “person”, which are used to match subjects to their data in the pedigree data.frame. The genotypes are coded as 0, 1, 2 for autosomal markers (typically a count of the number of the minor alleles).
#' @param covariate A data frame containing the covariate information. The first two columns are required to be named “ped” and “person”, which are used to match subjects to their data in the pedigree data frame. This is optional and the default "covariate = NULL" is for the case when the covariate matrix is not provided.
#' @param pos Position of the markers in base pairs.
#' @param order The order used to generate the B-spline basis.
#' @param beta_basis The number of basis functions used to estimate the genetic effect function.
#' @param base Can be either 'bspline' or 'fspline'.
#' @param optimizer Optimizer to use (default = "bobyqa").
#' @param Wald If Wald is set to true, return the Wald p-value in addition to the LRT p-value (Default: Wald = FALSE).
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{LRT}{The p-value based on a likelihood ratio test}
#'   \item{Wald}{The p-value based on a Wald test, returned if 'Wald' is TRUE}
#'   \item{nbetabasis}{The number of basis functions used to estimate the genetic effect function}
#'   \item{M_gao}{The effective number of variants in the region, as computed by M_GAO function}
#' }
#'
#' @importFrom stats anova pchisq vcov
#' @import fda
#' @import MASS
#' @import Matrix
#' @import nlme
#' @importFrom lme4 glmerControl
#' @import pedigreemm
#'
#' @export
#'
#' @seealso \code{\link{PedGLMM_additive_effect_model}}, \code{\link{PedGFLMM_fixed_model}}, \code{\link{exampleData}}
#'
#' @references
#' Chiu CY, Yuan F, Zhang BS, Yuan A, Li X, Fang HB, Lange K, Weeks DE, Wilson AF, Bailey-Wilson JE, Lakhal-Chaieb ML, Cook RJ, McMahon FJ, Amos CI, Xiong MM, and Fan RZ (2019) Pedigree-based linear mixed models for association analysis of quantitative traits with next-generation sequencing data. Genetic Epidemiology 43(2):189-206.
#'
#' Fan RZ, Wang YF, Mills JL, Wilson AF, Bailey-Wilson JE, and Xiong MM (2013) Functional linear models for association analysis of quantitative traits. Genetic Epidemiology 37 (7):726- 742.
#'
#' Fan RZ, Wang YF, Mills JL, Carter TC, Lobach I, Wilson AF, Bailey-Wilson JE, Weeks DE, and Xiong MM (2014) Generalized functional linear models for case-control association studies. Genetic Epidemiology 38 (7):622-637.
#'
#' Jiang YD, Chiu CY, Yan Q, Chen W, Gorin MB, Conley YP, Lakhal-Chaieb ML, Cook RJ, Amos CI, Wilson AF, Bailey-Wilson JE, McMahon FJ, Vazquez AI, Yuan A, Zhong XG, Xiong MM, Weeks DE, and Fan RZ (2020) Gene-based association testing of dichotomous traits with generalized linear mixed models for family data.
#'
#' Schaid DJ, McDonnell SK, Sinnwell JP, and Thibodeau SN (2013) Multiple genetic variant association testing by collapsing and kernel methods with pedigree or population structured data. Genetic Epidemiology 37:409-418.
#'
#' @examples
#' data(exampleData)
#'
#' order  =   4
#'
#' bsmooth_bsp=PedGFLMM_beta_smooth_only(ped = Ped, geno = as.matrix(geno),
#'    pos = snpPos$pos, order = order, covariate = as.matrix(cov),
#'    base = "bspline")
#' bsmooth_bsp
#'
#' bsmooth_fsp=PedGFLMM_beta_smooth_only(ped = Ped, geno = as.matrix(geno),
#'    pos = snpPos$pos, order = order, covariate = as.matrix(cov),
#'    base = "fspline")
#' bsmooth_fsp
#'
#' bsmooth_bsp_no_cov=PedGFLMM_beta_smooth_only(ped = Ped, geno = as.matrix(geno),
#'    pos = snpPos$pos, order = order, covariate = NULL,
#'    base = "bspline")
#' bsmooth_bsp_no_cov
#'
#' bsmooth_fsp_no_cov=PedGFLMM_beta_smooth_only(ped = Ped, geno = as.matrix(geno),
#'    pos = snpPos$pos, order = order, covariate = NULL,
#'    base = "fspline")
#' bsmooth_fsp_no_cov
PedGFLMM_beta_smooth_only = function(ped, geno, covariate=NULL, pos, order, beta_basis=NULL,
                                     base = "bspline", optimizer = "bobyqa", Wald = FALSE)
{
  ### GENERATING "pedigree_list" FOR USING R FUNCTION "pedigremm"
  dad = ped$father
  dad[dad == 0] = NA
  mom = ped$mother
  mom[mom == 0] = NA

  dat1 = ped
  dat1$dad_ID = dat1$mom_ID = NA

  fam = unique(dat1$ped)
  for (j in 1:length(fam))
  {
    dat=dat1[dat1$ped==fam[j],]

    for (i in 1:nrow(dat))
    {
      if (dat$father[i]!=0)
      {
        dad_id_in_family=dat$father[i]
        dat$dad_ID[i]=dat[dat$per==dad_id_in_family,]$ID
      }
      if (dat$mother[i]!=0)
      {
        mom_id_in_family=dat$mother[i]
        dat$mom_ID[i]=dat[dat$per==mom_id_in_family,]$ID
      }
    }

    dat1[dat1$ped==fam[j],]$dad_ID= dat$dad_ID
    dat1[dat1$ped==fam[j],]$mom_ID= dat$mom_ID
  }

  dad=dat1$dad_ID
  mom=dat1$mom_ID

  ### PEDIGREE LIST
  pede <- editPed( sire = as.integer(dad),
                   dam  = as.integer(mom),
                   label = ped$ID)
  ped_list = pedigreemm:: pedigree(label = pede$label, sire = pede$sire, dam = pede$dam)

  ### EXTRACT GENO
  geno = geno[,-c(1:2), drop = FALSE]  # REMOVE THE FIRST TWO COLUMNS WHICH INCLUDE PEDIGREE ID and PERSON ID WITHIN FAMILY

  ### number basis functions use dynamic rule
  if(is.null(beta_basis)) {NB = dRule(geno.only = geno)}

  ### Create basis
  if (base ==  "bspline")
  {
    nbeta_basis = ifelse(is.null(beta_basis), NB$betabasis_Bsp, beta_basis)
    betabasis  = create.bspline.basis(norder = order, nbasis = nbeta_basis)
  } else if (base == "fspline")
  {
    nbeta_basis = ifelse(is.null(beta_basis), NB$betabasis_Fsp, beta_basis)
    betabasis  = create.fourier.basis(c(0,1), nbasis = nbeta_basis)
  }else { }

  #### CASE WTIH NO MISSING GENO
  if (sum(is.na(geno))==0)
  {
    #dqr     = qr(geno)
    #index   = dqr$pivot[1:dqr$rank]
    #geno    = as.matrix(geno[, index])
    #pos     = pos[index]

    cond <- colSums(geno) > 0
    geno = geno[ ,cond]
    pos  = pos[cond]

    ### Normalize the region to [0,1] if needed
    if (max(pos) > 1)
    {
      pos = (pos - min(pos)) / (max(pos) - min(pos))
    } else{ }

    B  = eval.basis(pos, betabasis)
    UJ = as.matrix(geno) %*% B
  }

  #### CASE WITH MISSING GENO
  if (sum(is.na(geno))>0)
  {
    ### Normalize the region to [0,1] if needed
    if (max(pos) > 1)
    {
      pos = (pos - min(pos)) / (max(pos) - min(pos))
    } else{ }

    B  = eval.basis(pos, betabasis)
    UJ = NULL
    for (i in 1:nrow(geno))
    {
      idx = which(is.na(geno[i,]))
      if (length(idx) == 0 || length(idx) == ncol(geno))
      {
        gi = geno[i,]; Bi  = B
      } else
      {
        gi = geno[i,-idx]; Bi = B[-idx,]
      }

      gi_m = matrix(gi, nrow = 1)
      tmp =  unlist(gi_m) %*% Bi
      UJ  = rbind(UJ,tmp)
    }
  }

  #### SET TRAIT TO NA IF WHOLE GENO IS MISSING, SO THAT "fitnull" BELOW WILL USE DATA WITHOUT SUBJECT WITH WHOLE GENO MISSING
  if (sum(is.na(UJ))>0)
  {
    tmp1 = is.na(UJ)
    tmp2 = apply (tmp1, 1, sum)
    ped[tmp2>0,]$trait = NA
  }

  ### Model
  if(is.null(covariate))
  {
    dimcov = 0 # for Wald's test statistics

    fit     = pedigreemm( trait ~ UJ + (1|ID), pedigree = list(ID = ped_list), data=ped,
                          family = 'binomial',
                          control = glmerControl(optimizer, optCtrl = list(maxfun = 100000)))
    fitnull = pedigreemm( trait ~    + (1|ID), pedigree = list(ID = ped_list), data=ped,
                          family = 'binomial',
                          control = glmerControl(optimizer, optCtrl = list(maxfun = 100000)))
  } else
  {
    dimcov = ncol(covariate) - 2 # for Wald's test statistics
    covariate = covariate[,-c(1:2), drop = FALSE]  # REMOVE THE FIRST TWO COLUMNS WHICH INCLUDE PEDIGREE ID and PERSON ID WITHIN FAMILY
    #covariate[is.na(covariate)] = 0
    covariate <- as.matrix(covariate)

    fit     = pedigreemm( trait ~ covariate + UJ + (1|ID), pedigree = list(ID = ped_list), data=ped,
                          family = 'binomial',
                          control = glmerControl(optimizer, optCtrl = list(maxfun = 100000)))
    fitnull = pedigreemm( trait ~ covariate      + (1|ID), pedigree = list(ID = ped_list),
                          data=ped, family = 'binomial',
                          control = glmerControl(optimizer, optCtrl = list(maxfun = 100000)))
  }

  out = list()
  ### LRT
  nested_model_comparison = anova(fit,fitnull)
  out$LRT = nested_model_comparison[2,8]

  ### Wald test
  if (Wald)
  {
    coef_summary = summary(fit)$coefficients  # Sumary of coefficients
    coef         = coef_summary[,1]           # Obtain MLE of all parameters from the first column
    coef         = as.vector(coef)

    n_UJ    = length(coef)-(1+dimcov)                                # Number of coefficients for basis
    coef_UJ = coef[(2+dimcov):length(coef)]                          # Remove first 1+dimcov rows for intercept and covariate
    covm    = vcov(fit)                                              # Covariance matrix of all parameter's MLE
    covm_UJ = covm[(2+dimcov):length(coef),(2+dimcov):length(coef)]  # Covariance of coefficient MLE for geno
    covm_UJ = as.matrix(covm_UJ)

    Wald      = t(coef_UJ) %*% ginv(covm_UJ) %*% coef_UJ                # Wald's test statistic
    out$Wald = pchisq(as.numeric(Wald), df=ncol(UJ), lower.tail=FALSE) # df = number of restrictions i.e number of coef estimates for basis function ###
  }

  ###
  out$nbetabasis = nbeta_basis
  out$M_gao = ngeno_basis = ifelse(is.null(beta_basis), NB$M_gao, NA)

  out
}



#' PedGFLMM package
#'
#' @description This package implements family-based additive generalized linear mixed
#' models (GLMM) and generalized functional linear mixed models (GFLMM) for
#' gene-based association testing of dichotomous traits (Jiang et al, 2020).
#'
#' @author Yingda Jiang, Chi-Yang Chiu, Daniel E. Weeks, Ruzong Fan
#' @docType package
#' @name PedGFLMM
NULL
