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

setwd("C:/NICHD/Research/software/Fan/PedGFLMM")

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

betabasis_Bsp = 10
genobasis_Bsp = 10

betabasis_Fsp = 11
genobasis_Fsp = 11
order  =   4

add=PedGLMM_additive_effect_model(ped=Ped, geno = as.matrix(geno), covariate = as.matrix(cov), optimizer = "bobyqa", Wald = TRUE)  
add

fixed_bsp=PedGFLMM_fixed_model(ped = Ped, geno = as.matrix(geno), pos = snpPos$pos, order = order, beta_basis=betabasis_Bsp, geno_basis = genobasis_Bsp,covariate = as.matrix(cov), base = "bspline", optimizer = "bobyqa", Wald = TRUE)  
fixed_fsp=PedGFLMM_fixed_model(ped = Ped, geno = as.matrix(geno), pos = snpPos$pos, order = order, beta_basis=betabasis_Fsp, geno_basis = genobasis_Fsp,covariate = as.matrix(cov), base = "fspline", optimizer = "bobyqa", Wald = TRUE)  
fixed_bsp
fixed_fsp

bsmooth_bsp=PedGFLMM_beta_smooth_only(ped = Ped, geno = as.matrix(geno), pos = snpPos$pos, order = order, beta_basis=betabasis_Bsp, covariate = as.matrix(cov), base = "bspline", optimizer = "bobyqa", Wald = TRUE)  
bsmooth_fsp=PedGFLMM_beta_smooth_only(ped = Ped, geno = as.matrix(geno), pos = snpPos$pos, order = order, beta_basis=betabasis_Fsp, covariate = as.matrix(cov), base = "fspline", optimizer = "bobyqa", Wald = TRUE)  
bsmooth_bsp
bsmooth_fsp
