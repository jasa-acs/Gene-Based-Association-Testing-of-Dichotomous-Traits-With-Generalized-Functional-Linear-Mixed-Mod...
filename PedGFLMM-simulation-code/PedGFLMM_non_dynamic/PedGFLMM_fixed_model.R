### Updated by Chi-Yang on 01/08/2020: add control to change optimizer
### Updated by Chi-Yang on 12/13/2019: block Fan's QR, use DW's polymorphic variant matrix
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

PedGFLMM_fixed_model = function(ped, geno, covariate=NULL, pos, order, beta_basis, geno_basis, base = "bspline", optimizer= "bobyqa", Wald = FALSE)
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

   ### Create basis
   if (base ==  "bspline")
      {
      betabasis  = create.bspline.basis(norder = order, nbasis = beta_basis)
      genobasis  = create.bspline.basis(norder = order, nbasis = geno_basis)
      } else if (base == "fspline"){
      betabasis  = create.fourier.basis(c(0,1), nbasis = beta_basis)
      genobasis  = create.fourier.basis(c(0,1), nbasis = geno_basis)
      }else { }

   ### EXTRACT GENO
   geno = geno[,-c(1:2), drop = FALSE]  # REMOVE THE FIRST TWO COLUMNS WHICH INCLUDE PEDIGREE ID and PERSON ID WITHIN FAMILY

   #### DEALING WITH MISSING GENO
   idx.na = is.na(geno)

   #### CASE WTIH NO MISSING GENO
   if (sum(idx.na)==0)
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

      B = eval.basis(pos, genobasis)

      to_mul = ginv(t(B) %*% B) %*% t(B)
      U      = as.matrix(geno) %*% t( to_mul )
      J      = inprod(genobasis, betabasis)
      UJ = matrix( U %*% J, ncol = ncol(J) )
      }

   #### CASE WITH MISSING GENO
   if (sum(idx.na)>0)
      {
      ### Normalize the region to [0,1] if needed
      if (max(pos) > 1) 
         {
         pos = (pos - min(pos)) / (max(pos) - min(pos))
         } else{ }

      B  = eval.basis(pos, genobasis)
      J  = inprod(genobasis, betabasis)
      UJ = NULL

      for (i in 1:nrow(geno))
         {
         idx = which(is.na(geno[i,]))
         if (length(idx)== ncol(geno))
            {
            tmp =  rep(NA, ncol(B))
            }else if (length(idx)==  0)
               {
               gi  = geno[i,]
               Bi  = B
               to_mul = ginv(t(Bi) %*% Bi) %*% t(Bi)
               gi_m   = matrix(gi,nrow = 1)
               U      = unlist(gi_m) %*% t( to_mul )
               tmp    =  U %*% J
               }else
                  {
                  gi  = geno[i,-idx]
                  Bi  = B[-idx,]
                  to_mul = ginv(t(Bi) %*% Bi) %*% t(Bi)
                  gi_m   = matrix(gi,nrow = 1)
                  U      = unlist(gi_m) %*% t( to_mul )
                  tmp    =  U %*% J
                  }
         UJ = rbind(UJ,tmp)
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

      fit     = pedigreemm( trait ~ UJ + (1|ID), pedigree = list(ID = ped_list), data=ped, family = 'binomial',control = glmerControl(optimizer, optCtrl = list(maxfun = 100000)))
      fitnull = pedigreemm( trait ~    + (1|ID), pedigree = list(ID = ped_list), data=ped, family = 'binomial',control = glmerControl(optimizer, optCtrl = list(maxfun = 100000)))
      } else
         {
         dimcov = ncol(covariate) - 2 # for Wald's test statistics 
         covariate = covariate[,-c(1:2), drop = FALSE]  # REMOVE THE FIRST TWO COLUMNS WHICH INCLUDE PEDIGREE ID and PERSON ID WITHIN FAMILY
                                                        #covariate[is.na(covariate)] = 0
 
         fit     = pedigreemm( trait ~ covariate + UJ + (1|ID), pedigree = list(ID = ped_list), data=ped, family = 'binomial',control = glmerControl(optimizer, optCtrl = list(maxfun = 100000)))
         fitnull = pedigreemm( trait ~ covariate      + (1|ID), pedigree = list(ID = ped_list), data=ped, family = 'binomial',control = glmerControl(optimizer, optCtrl = list(maxfun = 100000)))
         }

   pval = list()
   ### LRT
   nested_model_comparison = anova(fit,fitnull)
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
      pval$Wald = pchisq(as.numeric(Wald), df=ncol(UJ), lower.tail=FALSE) # df = number of restrictions i.e number of coef estimates for basis function ###
      }

   pval
   }
