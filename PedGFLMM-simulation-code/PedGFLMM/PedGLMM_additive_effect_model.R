### Updated by Chi-Yang on 01/08/2020: add control to change optimizer
### Updated by Chi-Yang on 03/19/2019: block distribution = "binomial", make the routine work for no covariate  
### Written by Chi-Yang Chiu and Ruzong Fan, May 2016
### This function use r function "pedigreemm" to fit GFLMM for pedigree data
### need to provide ped_list that can be created by r function "pedigree" from library "pedigreemm"
### optimizer = "bobyqa" or "Nelder_Mead"

library(fda)
library(MASS)
library(Matrix)
library(nlme)
library(pedigreemm)

PedGLMM_additive_effect_model = function(ped, geno, covariate = NULL, optimizer= "bobyqa", Wald = FALSE)
   {
   length_index = dim(geno)[2]  ### added by Bingsong Zhang ###
   #### GENERATING "pedigree_list" FOR USING R FUNCTION "pedigremm"
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

   ###
   geno = geno[,-c(1:2)]  # REMOVE THE FIRST TWO COLUMNS WHICH INCLUDE PEDIGREE ID and PERSON ID WITHIN FAMILY
   geno[is.na(geno)]=0

   #dqr     = qr(geno)
   #index   = dqr$pivot[1:dqr$rank]
   if (length_index>3)                    ### added by Bingsong Zhang ###
      {
      #geno= as.matrix(geno[, index, drop = FALSE])
      cond <- colSums(geno) > 0
      geno = geno[ ,cond]
      } else
         {
         geno= as.matrix(geno)
         }

   ### Model
   if(is.null(covariate))
      {
      dimcov = 0 # for Wald's test statistics 

      fit     = pedigreemm( trait ~ geno + (1|ID), pedigree = list(ID = ped_list), data=ped, family = 'binomial', control = glmerControl(optimizer, optCtrl = list(maxfun = 100000)))
      fitnull = pedigreemm( trait ~        (1|ID), pedigree = list(ID = ped_list), data=ped, family = 'binomial', control = glmerControl(optimizer, optCtrl = list(maxfun = 100000)))
      }else
         {
         dimcov = ncol(covariate) - 2 # for Wald's test statistics 
         covariate = covariate[,-c(1:2), drop = FALSE]  # REMOVE THE FIRST TWO COLUMNS WHICH INCLUDE PEDIGREE ID and PERSON ID WITHIN FAMILY
                                                        #covariate[is.na(covariate)] = 0
 
         fit     = pedigreemm( trait ~ covariate + geno + (1|ID), pedigree = list(ID = ped_list), data=ped, family = 'binomial',control = glmerControl(optimizer, optCtrl = list(maxfun = 100000)))
         fitnull = pedigreemm( trait ~ covariate + (1|ID),        pedigree = list(ID = ped_list), data=ped, family = 'binomial',control = glmerControl(optimizer, optCtrl = list(maxfun = 100000)))
         }

   pval = list()
   ### LRT
   nested_model_comparison = anova(fit, fitnull)
   pval$LRT = nested_model_comparison[2,8]

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
      pval$Wald = pchisq(as.numeric(Wald), df=ncol(geno), lower.tail=FALSE) # df = number of restrictions i.e number of coef estimates for basis function ###
      }

   pval
   }
