#### THIS PROGRAM IS FOR ILLUSTRATION OF R FUNCTIONS "PedGFLMM"  
#### CREATED BY CHI-YANG, MAY, 2016
### Updated by Chi-Yang on 01/08/2020: add control to change optimizer
### Updated by Chi-Yang on 12/13/2019: block Fan's QR, use DW's polymorphic variant matrix
### Updated by Chi-Yang on 03/19/2019: block distribution = "binomial", make the routine work for no covariate 
### Written by Chi-Yang Chiu and Ruzong Fan, May 2016
### This function use r function "pedigreemm" to fit GFLMM for pedigree data
### need to provide ped_list that can be created by r function "pedigree" from library "pedigreemm"
### optimizer = "bobyqa" or "Nelder_Mead"

rm(list=ls())

## change working directory if necessary
#setwd("C:/Users/cchiu1/Desktop/JASA_Supplementary/PedGFLMM")

## Set the working directory to the location of this file when it is sourced:
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

source("MGAO.R")
source("dynamic_rule.R")

source("PedGLMM_additive_effect_model.R")
source("PedGFLMM_beta_smooth_only.R")
source("PedGFLMM_fixed_model.R")

library(fda)
###library(MASS)
###library(Matrix)
###library(nlme)
###require("pedgene")

library(pedigreemm) 

###
Ped     = read.csv("data/Ped.csv")
geno    = read.csv("data/geno.csv")
cov     = read.csv("data/covariate.csv")
snpPos  = read.csv("data/snpPos.csv")

order  =   4

#### Additive model
add=PedGLMM_additive_effect_model(ped=Ped, geno = as.matrix(geno), covariate = as.matrix(cov), optimizer = "bobyqa", Wald = TRUE)  
add

#### use dynamic rules to determine number of basis function
## fixed model
fixed_bsp1=PedGFLMM_fixed_model(ped = Ped, geno = as.matrix(geno), pos = snpPos$pos, order = order, beta_basis=NULL, geno_basis = NULL,covariate = as.matrix(cov), base = "bspline", optimizer = "bobyqa", Wald = TRUE)  
fixed_fsp1=PedGFLMM_fixed_model(ped = Ped, geno = as.matrix(geno), pos = snpPos$pos, order = order, beta_basis=NULL, geno_basis = NULL,covariate = as.matrix(cov), base = "fspline", optimizer = "bobyqa", Wald = TRUE)  
fixed_bsp1
fixed_fsp1

## beta-smoott only
bsmooth_bsp1=PedGFLMM_beta_smooth_only(ped = Ped, geno = as.matrix(geno), pos = snpPos$pos, order = order, beta_basis=NULL, covariate = as.matrix(cov), base = "bspline", optimizer = "bobyqa", Wald = TRUE)  
bsmooth_fsp1=PedGFLMM_beta_smooth_only(ped = Ped, geno = as.matrix(geno), pos = snpPos$pos, order = order, beta_basis=NULL, covariate = as.matrix(cov), base = "fspline", optimizer = "bobyqa", Wald = TRUE)  
bsmooth_bsp1
bsmooth_fsp1

#### use predetermined number of basis function
betabasis_Bsp = genobasis_Bsp = 10
betabasis_Fsp = genobasis_Fsp = 11

## fixed model
fixed_bsp2=PedGFLMM_fixed_model(ped = Ped, geno = as.matrix(geno), pos = snpPos$pos, order = order, beta_basis=betabasis_Bsp, geno_basis = genobasis_Bsp ,covariate = as.matrix(cov), base = "bspline", optimizer = "bobyqa", Wald = TRUE)  
fixed_fsp2=PedGFLMM_fixed_model(ped = Ped, geno = as.matrix(geno), pos = snpPos$pos, order = order, beta_basis=betabasis_Fsp, geno_basis = genobasis_Fsp,covariate = as.matrix(cov), base = "fspline", optimizer = "bobyqa", Wald = TRUE)  
fixed_bsp2
fixed_fsp2

## beta-smoott only
bsmooth_bsp2=PedGFLMM_beta_smooth_only(ped = Ped, geno = as.matrix(geno), pos = snpPos$pos, order = order, beta_basis=betabasis_Bsp, covariate = as.matrix(cov), base = "bspline", optimizer = "bobyqa", Wald = TRUE)  
bsmooth_fsp2=PedGFLMM_beta_smooth_only(ped = Ped, geno = as.matrix(geno), pos = snpPos$pos, order = order, beta_basis=betabasis_Fsp, covariate = as.matrix(cov), base = "fspline", optimizer = "bobyqa", Wald = TRUE)  
bsmooth_bsp2
bsmooth_fsp2
