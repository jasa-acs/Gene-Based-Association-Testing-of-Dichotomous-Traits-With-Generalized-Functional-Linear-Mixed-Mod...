### UPDATED BY CY, NOV 2019, USE EFFECTIVE MNUMBER OF VARIANTS, GAO et. al 2008 Gen Epi 
### 1st MODIFIED BY Fan, Feb 2019 ###

###################################################################
# Make batch file with different index
# idx1: causal percent, look at Parameters.r for
#       Causal_percent.A <- c(0.05, 0.10, 0.15)
#       Causal_percent   <- Causal_percent.A[idx1]
# So: idx1 = 1 implies causal percent = 0.05
#     idx1 = 2 implies causal percent = 0.10
#     idx1 = 3 implies causal percent = 0.15
############################################
# idx2 : region size,  look at Parameters.r for
# Range.A  <- c(6000, 9000, 12000, 15000, 18000, 21000)
# Region_L <- Range.A[idx2]
# So: idx2 = 1 implies Region_L = 3k 
###########################################
# idx3 : sample size, in loop
#######################################################
# idx4 : % of negative betas, look at Parameters.r for
# Beta1.Sign.A <- c(0, 0.2, 0.5)  # % of betas with negative sign
# Beta1.Sign <- Beta1.Sign.A[idx4]
# So: idx4 = 1 implies Beta1.Sign = 0
#     idx4 = 2 implies Beta1.Sign = 0.2
#     idx4 = 3 implies Beta1.Sign = 0.5
###################################################################
# idx5: upper limit of rare variant freq, look at Parameters.r for
#    CutOff.A <- c(0.03, 0.04, 0.05, 0.06) 
#    CutOff   <- CutOff.A [idx5]
# So idx5 = 1 implies CutOff = 0.03 used in the papers of Lin et al.
#    idx5 = 2 implies CutOff = 0.04
#    idx5 = 3 implies CutOff = 0.05
#    idx5 = 4 implies CutOff = 0.06
####################################################################
#setwd("/data/chiuc2/NICHD/Research/paper/PedGFLMM/JASA_Supplementary")

## Set the working directory to the location of this file when it is sourced:
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

#source("PedGFLMM/PedGFLMM_beta_smooth_only.R")
#source("PedGFLMM/PedGFLMM_fixed_model.R")

source("../../PedGFLMM/PedGFLMM_beta_smooth_only.R")
source("../../PedGFLMM/PedGFLMM_fixed_model.R")

library(fda)
library(MASS)
library(Matrix)
library(nlme)
library(glmm)
require("pedgene")
library(pedigreemm)
library(tidyverse)

#source("simulations/Function.r")
#source("simulations/Parameters.r")

source("../Function.r")
source("../Parameters.r")

# custom tryCatch to return result and warnings -- http://stackoverflow.com/a/24569739/2271856
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

######################################################
#   Main Simulation Part
#DLL_Path  <- "simulations/Util.so"
#BIN_FILE  <- "simulations/SKAT/Simulation_Code_Example/out.50.hap-1.bin"
#DATA_FILE <- "simulations/SKAT/Simulation_Code_Example/out.50.hap-1"
#Info_FILE <- "simulations/SKAT/Simulation_Code_Example/out.50.pos-1"

DLL_Path  <- "../Util.so"
BIN_FILE  <- "../SKAT/Simulation_Code_Example/out.50.hap-1.bin"
DATA_FILE <- "../SKAT/Simulation_Code_Example/out.50.hap-1"
Info_FILE <- "../SKAT/Simulation_Code_Example/out.50.pos-1"
if (!file.exists(DLL_Path)) {
  stop("ERROR: The required file 'Util.so' is missing.")
}

# Number of markers
N_Marker       <- 20142
# Number of haplotypes
N              <- 10000
# Total region length
Total_Region_L <- 1000000 # Total 1 Mb

dyn.load(DLL_Path)
INFO    <- read.table(Info_FILE, header=TRUE)
snpPos  <- as.vector( INFO$CHROM_POS )

Init(BIN_FILE, 100)

###################################################
### Read in the pedigree data and kinCoef files ### 
###################################################
#kinCoef <- get(load("simulations/Input/kinCoef.RData"))
#pedstr  <- get(load("simulations/Input/pedstr.RData"))

kinCoef <- get(load("../Input/kinCoef.RData"))
pedstr  <- get(load("../Input/pedstr.RData"))

### Second copy of the pedigree ###
pedstr2 = pedstr
pedstr2$Ped = pedstr2$Ped + 100
pedstr2$ID  = pedstr2$ID + 228

### Combine the two pedigree files ###
pedstr = rbind(pedstr, pedstr2)
###

###  Make kinship Coefficients for the new pedigree ###
kin1 = cbind(kinCoef, 0*kinCoef)
kin2 = cbind(0*kinCoef, kinCoef)
kinCoef =rbind(kin1, kin2)

#################################
### parameters for simulation ###
#################################

POLYGEN = 0.2  # polygenic effect
PREV  <- 0.01   # disease prevalence
beta0 <- log(PREV/(1 - PREV))   # intercept adjusting population prevalence

N.SIM = 3000  # Number of replicates
stop = 999 # stop rule

#### number of basis functions
#betabasis_Bsp = genobasis_Bsp = 16
#betabasis_Fsp = genobasis_Fsp = 17
order  =   4

ROW_NAME =  c("0.05", "0.01", "1E-3", "1E-4", "1E-5", "1E-6")

power_func = function (idx1 = 1, idx2 = 2, idx4 =1, idx5 = 1, causal_variant = "rare_and_common", test_C, Delta)
   {
   ranseed = 100
   set.seed(ranseed + idx1 + idx2 + idx4 + idx5)
   Region_L   <- Range.A[idx2]   # Region length
   
   Causal_percent <- Causal_percent.A[idx1]
   Beta1.Sign <- Beta1.Sign.A[idx4]
   CutOff     <- CutOff.A[idx5]
   print(CutOff)
   
   ###
   bcnt_LRT_fixed <- rep(0, 6)
   fcnt_LRT_fixed <- rep(0, 6)
   bcnt_LRT_beta  <- rep(0, 6)
   fcnt_LRT_beta  <- rep(0, 6)  
   
   bpval_LRT_fixed_total = 0
   fpval_LRT_fixed_total = 0
   bpval_LRT_beta_total  = 0  
   fpval_LRT_beta_total  = 0

   cnt_kernel_BT = cnt_kernel_MB = cnt_kernel_UW <- rep(0, 6)
   cnt_burden_BT = cnt_burden_MB = cnt_burden_UW <- rep(0, 6)
   
   kernel_BT_total = kernel_MB_total = kernel_UW_total = 0 
   burden_BT_total = burden_MB_total = burden_UW_total = 0 
   
   # Get Rare SNPs
   Marker.Rare <- which(INFO$FREQ1 <= CutOff)

   sumTable <- data.frame(type=NA,rep=NA,LRT=NA,warning=NA,error=NA,warnings=NA)

   total_run = 0
   for (i in 1:N.SIM)
      {
      total_run = total_run + 1
      # Generate Region
      Region.Start <- runif(1)*(Total_Region_L - Region_L)
      Region.End   <- Region.Start + Region_L
  
      Marker.Idx1  <- which(INFO$CHROM_POS >= Region.Start)
      Marker.Idx2  <- which(INFO$CHROM_POS <= Region.End)
      Marker.Idx3  <- sort(intersect(Marker.Idx1,Marker.Idx2))
      IDX_IN       <- which(INFO$FREQ1[Marker.Idx3] -  INFO$FREQ2[Marker.Idx3] < 0)
      Marker.Idx   <- Marker.Idx3[IDX_IN]  
      Marker.N     <- length(Marker.Idx)
  
      # Select SNPs to use: rare variants or rare & common variants        
      if (causal_variant == "rare") 
         {
         Marker.Causal1  <- intersect(Marker.Rare, Marker.Idx)      
         }else if (causal_variant == "rare_and_common")
            {
            Marker.Causal1  <- Marker.Idx
            }else {print("Error: causal variants is not properly defined!")}
       
      Marker.Causal   <- sort(sample(Marker.Causal1, length(Marker.Causal1) *Causal_percent))      
      Marker.Causal.N <- length(Marker.Causal)
		
      Marker.Causal.Index <- NULL
      for (k in Marker.Causal)
            {
            temp <- which(Marker.Idx == k)
            Marker.Causal.Index <- c(Marker.Causal.Index, temp)
            }

      # MAF of Marker_Causal
      Marker.Causal.MAF <- INFO$FREQ1[Marker.Causal]
    
      # If Beta1 < 0
      if (Delta == 0)
         {Beta1.VEC <- Get_BETA_C(idx1, idx2, Marker.Causal.MAF)
         }else
            {Beta1.VEC <- Get_BETA_C_Delta(idx1, idx2, Marker.Causal.MAF, Delta)}
      
      # Beta1.Sign
      if (Beta1.Sign > 0)
            {
            temp.n   <- Marker.Causal.N * Beta1.Sign
            temp.idx <- sample(1:Marker.Causal.N, temp.n)
            Beta1.VEC[temp.idx] <- - Beta1.VEC[temp.idx]
            } 

      out_geno_idx <- Simu_Control(10000)
      hap <- Get_Haplotype(Marker.Idx, out_geno_idx)

      NRhap <- nrow(hap)
      NChap <- ncol(hap)

      ############################################################### 
      ## Simulate and acertain 25 pedigrees with at least 2 cases.
      ## Transmit haplotypes to offspring.
      ## Convert haplotypes to genotypes, and assign phenotypes.
      ###############################################################

      ## define the pedigree structure ##
      pedstr$trait <- -1    # preassign the trait #
      pedstr$cov1  <- -1    # preassign the covariate 1 #
      pedstr$cov2  <- -1    # preassign the covariate 2 #
      family <- unique(pedstr$Ped)
      pedstr$m2 <- pedstr$m1 <- 0

      ## assign haplotypes and genotypes ##
      powerPed  <- data.frame()
      powerGeno <- data.frame()

      for (f in family) 
         {
         substr    <- pedstr[pedstr$Ped == f, ]
         founderID <- substr[substr$Father == 0, "ID"]
         Nfounder  <- length(founderID)
         nonFounderID <- substr[substr$Father != 0, "ID"]
    
         index <- 0    # ascertainment index # 
  
         while (index == 0) 
            { ## sample haplotypes for founders ##
            substr[substr$ID %in% founderID, "m1"] <- sample(c(1:NRhap), Nfounder, replace= TRUE)
            substr[substr$ID %in% founderID, "m2"] <- sample(c(1:NRhap), Nfounder, replace= TRUE)
    
            ## sample haplotypes and convert to genotypes ##
            for (j in nonFounderID) 
               {
               fatherPer  <- substr[substr$ID == j, "Father"]
               fatherPool <- as.vector(substr[substr$Per == fatherPer, c("m1","m2")])
               substr[substr$ID == j, "m1"] <- sample(fatherPool, 1, replace= TRUE)
               motherPer <- substr[substr$ID == j, "Mother"]
               motherPool <- as.vector(substr[substr$Per == motherPer, c("m1","m2")])
               substr[substr$ID == j, "m2"] <- sample(motherPool, 1, replace= TRUE)
               }
    
            NRsubstr <- nrow(substr)
            geno <- hap[substr$m1, ] + hap[substr$m2, ]
            colnames(geno) <- colnames(hap)
 
            geno_causal <- geno[, Marker.Causal.Index]
            g = geno_causal %*% Beta1.VEC / test_C
 
            ## simulate polygenic effect ##
            muMean   <- rep(0, NRsubstr)
            kinStart <- substr[1,"ID"]
            kinEnd   <- substr[NRsubstr,"ID"]
            corSigma <- kinCoef[kinStart:kinEnd,kinStart:kinEnd]
            a <- mvrnorm(n= 1, mu= muMean, Sigma= corSigma, tol= 1e-16, empirical= FALSE, EISPACK= FALSE)    

            ## compute the Bernoulli probability of affected status conditional on genotype ##
            ## based on Schaid's proposed model ##
            for (k in 1:NRsubstr) 
               {
               cov1 = rbinom(1, size=1, 0.5)
               cov2 = rnorm(1, mean = 0, sd = 1)     
               A <- beta0 + cov1 + cov2 + g[k] + POLYGEN*a[k]
               mu <- exp(A) / (1 + exp(A))
               substr[k, "trait"] <- sample(c(1,0), 1, replace= FALSE, prob = c(mu,1-mu))
               substr[k, "cov1"]  <- cov1
               substr[k, "cov2"]  <- cov2
               }
            ## Ascertain at least 2 affected siblings in the pedigree ##
            affectMomID <- substr[substr$trait == 1, "Mother"]
            dup <- affectMomID[duplicated(affectMomID)]
            if (sum(dup) != 0) 
               {
               powerPed  <- rbind(powerPed, substr)
               powerGeno <- rbind(powerGeno, geno)
               index <- 1
               }    
            }
         }  # end of for loop #
  
      # Covariates
      ped       = data.frame(ID=powerPed$ID, ped=powerPed$Ped, person=powerPed$Per, father=powerPed$Father, mother=powerPed$Mother, sex=powerPed$Sex, trait=powerPed$trait)
      geno      = data.frame(ped=powerPed$Ped, person=powerPed$Per,powerGeno)
      covariate = data.frame(ped=powerPed$Ped, person=powerPed$Per, X1 = powerPed$cov1, X2 = powerPed$cov2)

      geno.only = geno[,-c(1:2)]  # THE FIRST TWO COLUMNS ARE PEDIGREE ID AND PERSON ID WITHIN FAMILY
      map       = data.frame(chr=rep(1,ncol(geno.only)), snp=colnames(geno.only) , pos = snpPos[Marker.Idx])
         
      # SNP matrix
      cond <- colSums(geno.only) > 0
      SNP_mx = geno.only[ ,cond]

      # Number of polymorphic variants
      n.poly <- ncol(SNP_mx)  ## BY DW
      # Number of variants
      n.vars <- ncol(geno.only)
 
      tryCatch({
      ### fixed GFLMM ###
      bpval_LRT_fixed = NA
      msgs <- NA
      msgs <- capture.output(myWarnErr <- myTryCatch(PedGFLMM_fixed_model(ped = ped, 
              geno = as.matrix(geno), covariate = as.matrix(covariate), pos = map$pos, 
              order = order, beta_basis = NULL, geno_basis = NULL, 
              base = "bspline", optimizer = "bobyqa", Wald = FALSE)),type="message")

     if (is.null(myWarnErr$value)) {
        myWarnErr$value$LRT <- NA
      }      
      bpval_PedGFLMM_fixed <- myWarnErr$value
      if (is.null(myWarnErr$warning) == TRUE) {
        msg.w <- NA
      } else {
        msg.w <- conditionMessage(myWarnErr$warning)
      }
      if (is.null(myWarnErr$error) == TRUE) {
        msg.e <- NA
      } else {
        msg.e <- conditionMessage(myWarnErr$error)
      }

      
      sumTable <- bind_rows(sumTable,data.frame(type="Fixed",rep=i,n.poly=n.poly,n.vars=n.vars,M_gao=bpval_PedGFLMM_fixed$M_gao,betabasis=bpval_PedGFLMM_fixed$nbetabasis, genobasis=bpval_PedGFLMM_fixed$ngenobasis,LRT=bpval_PedGFLMM_fixed$LRT, 
                                                warning=msg.w,
                                                error=msg.e,warnings=paste0(msgs,";",collapse=""),stringsAsFactors = FALSE))
      
      bpval_LRT_fixed   = bpval_PedGFLMM_fixed$LRT
      wrs = paste0(msgs,";",collapse="")      

      if (is.null(myWarnErr$value) != TRUE & is.null(myWarnErr$warning) == TRUE & is.null(myWarnErr$error) == TRUE & wrs == ";")
         {
         bpval_LRT_fixed_total = bpval_LRT_fixed_total + 1
         if (bpval_LRT_fixed < 0.000001) {bcnt_LRT_fixed[6] = bcnt_LRT_fixed[6] + 1}
         if (bpval_LRT_fixed < 0.00001)  {bcnt_LRT_fixed[5] = bcnt_LRT_fixed[5] + 1}
         if (bpval_LRT_fixed < 0.0001)   {bcnt_LRT_fixed[4] = bcnt_LRT_fixed[4] + 1}
         if (bpval_LRT_fixed < 0.001)    {bcnt_LRT_fixed[3] = bcnt_LRT_fixed[3] + 1}
         if (bpval_LRT_fixed < 0.01)     {bcnt_LRT_fixed[2] = bcnt_LRT_fixed[2] + 1}
         if (bpval_LRT_fixed < 0.05)     {bcnt_LRT_fixed[1] = bcnt_LRT_fixed[1] + 1}
         }

      ###fspline fixed
      fpval_LRT_fixed = NA
      msgs <- NA

      msgs <- capture.output(myWarnErr <- myTryCatch(PedGFLMM_fixed_model(ped = ped, 
                 geno = as.matrix(geno), covariate = as.matrix(covariate), pos = map$pos, 
                 order = order, beta_basis = NULL, geno_basis = NULL, 
                 base = "fspline", optimizer = "bobyqa", Wald = FALSE)), type="message")
      
      if (is.null(myWarnErr$value)) {
        myWarnErr$value$LRT <- NA
      }
      fpval_PedGFLMM_fixed <- myWarnErr$value
      if (is.null(myWarnErr$warning) == TRUE) {
        msg.w <- NA
      } else {
        msg.w <- conditionMessage(myWarnErr$warning)
      }
      if (is.null(myWarnErr$error) == TRUE) {
        msg.e <- NA
      } else {
        msg.e <- conditionMessage(myWarnErr$error)
      }
       
      
      sumTable <- bind_rows(sumTable,data.frame(type="fFixed",rep=i,n.poly=n.poly,n.vars=n.vars,M_gao=fpval_PedGFLMM_fixed$M_gao, betabasis=fpval_PedGFLMM_fixed$nbetabasis, genobasis=fpval_PedGFLMM_fixed$ngenobasis,LRT=fpval_PedGFLMM_fixed$LRT, 
                                                warning=msg.w,
                                                error=msg.e,warnings=paste0(msgs,";",collapse=""),stringsAsFactors = FALSE))
      
      fpval_LRT_fixed   = fpval_PedGFLMM_fixed$LRT
      wrs = paste0(msgs,";",collapse="")

      if (is.null(myWarnErr$value) != TRUE & is.null(myWarnErr$warning) == TRUE & is.null(myWarnErr$error) == TRUE & wrs == ";")
         {
         fpval_LRT_fixed_total =  fpval_LRT_fixed_total + 1
         if (fpval_LRT_fixed < 0.000001) {fcnt_LRT_fixed[6] = fcnt_LRT_fixed[6] + 1}
         if (fpval_LRT_fixed < 0.00001)  {fcnt_LRT_fixed[5] = fcnt_LRT_fixed[5] + 1}
         if (fpval_LRT_fixed < 0.0001)   {fcnt_LRT_fixed[4] = fcnt_LRT_fixed[4] + 1}
         if (fpval_LRT_fixed < 0.001)    {fcnt_LRT_fixed[3] = fcnt_LRT_fixed[3] + 1}
         if (fpval_LRT_fixed < 0.01)     {fcnt_LRT_fixed[2] = fcnt_LRT_fixed[2] + 1}
         if (fpval_LRT_fixed < 0.05)     {fcnt_LRT_fixed[1] = fcnt_LRT_fixed[1] + 1}
         }
     
      #### beta-smooth only
      ###bspline beta smooth only
      bpval_LRT_beta = NA
      msgs <- NA
      msgs <- capture.output(myWarnErr <- myTryCatch(PedGFLMM_beta_smooth_only(ped = ped, geno = as.matrix(geno), 
                                   covariate = as.matrix(covariate), pos = map$pos, order = order, 
                                   beta_basis = NULL, base = "bspline", optimizer = "bobyqa",
                                   Wald = FALSE)), type="message")

      if (is.null(myWarnErr$value)) {
        myWarnErr$value$LRT <- NA
      }      
      bpval_PedGFLMM_beta_smooth_only <- myWarnErr$value
      if (is.null(myWarnErr$warning) == TRUE) {
        msg.w <- NA
      } else {
        msg.w <- conditionMessage(myWarnErr$warning)
      }
      if (is.null(myWarnErr$error) == TRUE) {
        msg.e <- NA
      } else {
        msg.e <- conditionMessage(myWarnErr$error)
      }
  
      
      sumTable <- bind_rows(sumTable,data.frame(type="bBetaSmooth",rep=i,n.poly=n.poly,n.vars=n.vars,M_gao=bpval_PedGFLMM_beta_smooth_only$M_gao,betabasis=bpval_PedGFLMM_beta_smooth_only$nbetabasis, genobasis=NA,LRT=bpval_PedGFLMM_beta_smooth_only$LRT, 
                                                warning=msg.w,
                                                error=msg.e,warnings=paste0(msgs,";",collapse=""),stringsAsFactors = FALSE))
      
      bpval_LRT_beta  = bpval_PedGFLMM_beta_smooth_only$LRT
      wrs = paste0(msgs,";",collapse="")

      if (is.null(myWarnErr$value) != TRUE & is.null(myWarnErr$warning) == TRUE & is.null(myWarnErr$error) == TRUE & wrs == ";")
         {
         bpval_LRT_beta_total = bpval_LRT_beta_total + 1
         if (bpval_LRT_beta < 0.000001) {bcnt_LRT_beta[6] = bcnt_LRT_beta[6] + 1}
         if (bpval_LRT_beta < 0.00001)  {bcnt_LRT_beta[5] = bcnt_LRT_beta[5] + 1}
         if (bpval_LRT_beta < 0.0001)   {bcnt_LRT_beta[4] = bcnt_LRT_beta[4] + 1}
         if (bpval_LRT_beta < 0.001)    {bcnt_LRT_beta[3] = bcnt_LRT_beta[3] + 1}
         if (bpval_LRT_beta < 0.01)     {bcnt_LRT_beta[2] = bcnt_LRT_beta[2] + 1}
         if (bpval_LRT_beta < 0.05)     {bcnt_LRT_beta[1] = bcnt_LRT_beta[1] + 1}
         }

      ###fspline beta smooth only
      fpval_LRT_beta = NA
      msgs <- NA
      msgs <- capture.output(myWarnErr <- myTryCatch(PedGFLMM_beta_smooth_only(ped = ped, 
                 geno = as.matrix(geno), covariate = as.matrix(covariate), pos = map$pos, 
                 order = order, beta_basis = NULL, base = "fspline", optimizer = "bobyqa",
                 Wald = FALSE)), type="message")

      if (is.null(myWarnErr$value)) {
        myWarnErr$value$LRT <- NA
      }
      fpval_PedGFLMM_beta_smooth_only <- myWarnErr$value
      if (is.null(myWarnErr$warning) == TRUE) {
        msg.w <- NA
      } else {
        msg.w <- conditionMessage(myWarnErr$warning)
      }
      if (is.null(myWarnErr$error) == TRUE) {
        msg.e <- NA
      } else {
        msg.e <- conditionMessage(myWarnErr$error)
      }
       
      
      sumTable <- bind_rows(sumTable,data.frame(type="fBetaSmooth",rep=i,n.poly=n.poly,n.vars=n.vars,M_gao=fpval_PedGFLMM_beta_smooth_only$M_gao,betabasis=fpval_PedGFLMM_beta_smooth_only$nbetabasis, genobasis=NA,LRT=fpval_PedGFLMM_beta_smooth_only$LRT, 
                                                warning=msg.w,
                                                error=msg.e,warnings=paste0(msgs,";",collapse=""),stringsAsFactors = FALSE))
         
      fpval_LRT_beta  = fpval_PedGFLMM_beta_smooth_only$LRT
      wrs = paste0(msgs,";",collapse="")

      if (is.null(myWarnErr$value) != TRUE & is.null(myWarnErr$warning) == TRUE & is.null(myWarnErr$error) == TRUE & wrs == ";")
         {
         fpval_LRT_beta_total =  fpval_LRT_beta_total + 1
         if (fpval_LRT_beta < 0.000001) {fcnt_LRT_beta[6] = fcnt_LRT_beta[6] + 1}
         if (fpval_LRT_beta < 0.00001)  {fcnt_LRT_beta[5] = fcnt_LRT_beta[5] + 1}
         if (fpval_LRT_beta < 0.0001)   {fcnt_LRT_beta[4] = fcnt_LRT_beta[4] + 1}
         if (fpval_LRT_beta < 0.001)    {fcnt_LRT_beta[3] = fcnt_LRT_beta[3] + 1}
         if (fpval_LRT_beta < 0.01)     {fcnt_LRT_beta[2] = fcnt_LRT_beta[2] + 1}
         if (fpval_LRT_beta < 0.05)     {fcnt_LRT_beta[1] = fcnt_LRT_beta[1] + 1}
         }
         
         ### pedgene results ###
         n_snp = ncol(geno)-2 # NUMBER OF SNPS, TAKING AWAY 2 BECAUSE THE FIRST TWO COLUMNS ARE PEDIGRRE ID AND PERSON ID WITHIN FAMILY
         weight1 <- rep(1, n_snp)
         schaidPed <- ped
         schaidPed$ID <- NULL   #REMOVE ID FROM DATA SET
         names(schaidPed)[names(schaidPed)=="ped"] = "famid"

         cov = covariate[,-c(1:2)]
         cov = as.matrix(cov)
         binfit <- glm(ped$trait ~ cov, family="binomial")
         schaidPed$trait.adjusted[!is.na(ped$trait)] <- fitted.values(binfit)
        
         # ADD ADJUSTED TRAIT
         schaidPowerG <- geno
         names(schaidPowerG)[names(schaidPowerG)=="ped"] = "famid"

         pKernelPowerBT = pBurdenPowerBT = NA
         schaid_powerBT <- pedgene(schaidPed, schaidPowerG, male.dose= 2, checkpeds= FALSE, weights= NULL, weights.mb= FALSE, method= "kounen", acc.davies= 1e-9) 
         pKernelPowerBT <- schaid_powerBT$pgdf$pval.kernel
         pBurdenPowerBT <- schaid_powerBT$pgdf$pval.burden
         if (is.na(pKernelPowerBT) != TRUE)
            { 
            kernel_BT_total = kernel_BT_total + 1  

            if (pKernelPowerBT < 0.000001) {cnt_kernel_BT[6] = cnt_kernel_BT[6] + 1}
            if (pKernelPowerBT < 0.00001)  {cnt_kernel_BT[5] = cnt_kernel_BT[5] + 1}
            if (pKernelPowerBT < 0.0001)   {cnt_kernel_BT[4] = cnt_kernel_BT[4] + 1}
            if (pKernelPowerBT < 0.001)    {cnt_kernel_BT[3] = cnt_kernel_BT[3] + 1}
            if (pKernelPowerBT < 0.01)     {cnt_kernel_BT[2] = cnt_kernel_BT[2] + 1}
            if (pKernelPowerBT < 0.05)     {cnt_kernel_BT[1] = cnt_kernel_BT[1] + 1}
            }
            
         if (is.na(pBurdenPowerBT) != TRUE)
            {
            burden_BT_total = burden_BT_total+1          

            if (pBurdenPowerBT < 0.000001) {cnt_burden_BT[6] = cnt_burden_BT[6] + 1}
            if (pBurdenPowerBT < 0.00001)  {cnt_burden_BT[5] = cnt_burden_BT[5] + 1}
            if (pBurdenPowerBT < 0.0001)   {cnt_burden_BT[4] = cnt_burden_BT[4] + 1}
            if (pBurdenPowerBT < 0.001)    {cnt_burden_BT[3] = cnt_burden_BT[3] + 1}
            if (pBurdenPowerBT < 0.01)     {cnt_burden_BT[2] = cnt_burden_BT[2] + 1}
            if (pBurdenPowerBT < 0.05)     {cnt_burden_BT[1] = cnt_burden_BT[1] + 1}    
            }
            
         ###
         pKernelPowerMB = pBurdenPowerMB = NA
         schaid_powerMB <- pedgene(schaidPed, schaidPowerG, male.dose= 2, checkpeds= FALSE, weights= NULL, weights.mb= TRUE, method= "kounen", acc.davies= 1e-9) 
         pKernelPowerMB <- schaid_powerMB$pgdf$pval.kernel
         pBurdenPowerMB <- schaid_powerMB$pgdf$pval.burden
         if (is.na(pKernelPowerMB) != TRUE)
            { 
            kernel_MB_total = kernel_MB_total + 1  

            if (pKernelPowerMB < 0.000001) {cnt_kernel_MB[6] = cnt_kernel_MB[6] + 1}
            if (pKernelPowerMB < 0.00001)  {cnt_kernel_MB[5] = cnt_kernel_MB[5] + 1}
            if (pKernelPowerMB < 0.0001)   {cnt_kernel_MB[4] = cnt_kernel_MB[4] + 1}
            if (pKernelPowerMB < 0.001)    {cnt_kernel_MB[3] = cnt_kernel_MB[3] + 1}
            if (pKernelPowerMB < 0.01)     {cnt_kernel_MB[2] = cnt_kernel_MB[2] + 1}
            if (pKernelPowerMB < 0.05)     {cnt_kernel_MB[1] = cnt_kernel_MB[1] + 1}
            }
            
         if (is.na(pBurdenPowerMB) != TRUE)
            {
            burden_MB_total = burden_MB_total+1          

            if (pBurdenPowerMB < 0.000001) {cnt_burden_MB[6] = cnt_burden_MB[6] + 1}
            if (pBurdenPowerMB < 0.00001)  {cnt_burden_MB[5] = cnt_burden_MB[5] + 1}
            if (pBurdenPowerMB < 0.0001)   {cnt_burden_MB[4] = cnt_burden_MB[4] + 1}
            if (pBurdenPowerMB < 0.001)    {cnt_burden_MB[3] = cnt_burden_MB[3] + 1}
            if (pBurdenPowerMB < 0.01)     {cnt_burden_MB[2] = cnt_burden_MB[2] + 1}
            if (pBurdenPowerMB < 0.05)     {cnt_burden_MB[1] = cnt_burden_MB[1] + 1}    
            }
         
         ###
         pKernelPowerUW = pBurdenPowerUW = NA
         schaid_powerUW <- pedgene(schaidPed, schaidPowerG, male.dose= 2, checkpeds= FALSE, weights= weight1, weights.mb= NULL, method= "kounen", acc.davies= 1e-9) 
         pKernelPowerUW <- schaid_powerUW$pgdf$pval.kernel
         pBurdenPowerUW <- schaid_powerUW$pgdf$pval.burden
         if (is.na(pKernelPowerUW) != TRUE)
            { 
            kernel_UW_total = kernel_UW_total + 1  

            if (pKernelPowerUW < 0.000001) {cnt_kernel_UW[6] = cnt_kernel_UW[6] + 1}
            if (pKernelPowerUW < 0.00001)  {cnt_kernel_UW[5] = cnt_kernel_UW[5] + 1}
            if (pKernelPowerUW < 0.0001)   {cnt_kernel_UW[4] = cnt_kernel_UW[4] + 1}
            if (pKernelPowerUW < 0.001)    {cnt_kernel_UW[3] = cnt_kernel_UW[3] + 1}
            if (pKernelPowerUW < 0.01)     {cnt_kernel_UW[2] = cnt_kernel_UW[2] + 1}
            if (pKernelPowerUW < 0.05)     {cnt_kernel_UW[1] = cnt_kernel_UW[1] + 1}
            }
            
         if (is.na(pBurdenPowerUW) != TRUE)
            {
            burden_UW_total = burden_UW_total+1          

            if (pBurdenPowerUW < 0.000001) {cnt_burden_UW[6] = cnt_burden_UW[6] + 1}
            if (pBurdenPowerUW < 0.00001)  {cnt_burden_UW[5] = cnt_burden_UW[5] + 1}
            if (pBurdenPowerUW < 0.0001)   {cnt_burden_UW[4] = cnt_burden_UW[4] + 1}
            if (pBurdenPowerUW < 0.001)    {cnt_burden_UW[3] = cnt_burden_UW[3] + 1}
            if (pBurdenPowerUW < 0.01)     {cnt_burden_UW[2] = cnt_burden_UW[2] + 1}
            if (pBurdenPowerUW < 0.05)     {cnt_burden_UW[1] = cnt_burden_UW[1] + 1}    
            }
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})   
      
      if (i %% 100 == 0) print(i)
 
       #### break for-loop if the effective number of replicates is obtained
       if (fpval_LRT_fixed_total > stop & bpval_LRT_beta_total > stop & fpval_LRT_beta_total > stop){
       print(paste("Coming out from for loop Where i =  ", i))
       break
       }

      }

   ###
   fixed_LRT_bpval  = fixed_LRT_fpval  = beta_LRT_bpval  = beta_LRT_fpval  = rep(NA,6)
 
   if(bpval_LRT_fixed_total != 0){ fixed_LRT_bpval = bcnt_LRT_fixed/bpval_LRT_fixed_total}
   if(fpval_LRT_fixed_total != 0){ fixed_LRT_fpval = fcnt_LRT_fixed/fpval_LRT_fixed_total}
   if(bpval_LRT_beta_total  != 0){  beta_LRT_bpval = bcnt_LRT_beta/bpval_LRT_beta_total}
   if(fpval_LRT_beta_total  != 0){  beta_LRT_fpval = fcnt_LRT_beta/fpval_LRT_beta_total}
   
   kernel_BT_pval    = cnt_kernel_BT/kernel_BT_total
   kernel_MB_pval    = cnt_kernel_MB/kernel_MB_total
   kernel_UW_pval    = cnt_kernel_UW/kernel_UW_total
   
   burden_BT_pval    = cnt_burden_BT/burden_BT_total
   burden_MB_pval    = cnt_burden_MB/burden_MB_total
   burden_UW_pval    = cnt_burden_UW/burden_UW_total

   pval = data.frame( fixed_LRT_bpval,  fixed_LRT_fpval,  beta_LRT_bpval,  beta_LRT_fpval, 
                      kernel_BT_pval, kernel_MB_pval, kernel_UW_pval, burden_BT_pval, burden_MB_pval, burden_UW_pval,
                      bpval_LRT_fixed_total,  fpval_LRT_fixed_total,  bpval_LRT_beta_total,  fpval_LRT_beta_total,
                      bcnt_LRT_fixed,         fcnt_LRT_fixed,         bcnt_LRT_beta,         fcnt_LRT_beta, total_run)   
  
   if (causal_variant == "rare" && Delta != 0) 
      {
      File_Pval     <- sprintf("./power_test_C=%d_Delta=%.2f/Power_causal_rare/Power_rare_region_size_%d_Neg_beta_pct_%d_Causal_pct_%d_Nsim_%d.csv", test_C, Delta, Region_L, Beta1.Sign * 100, Causal_percent * 100, N.SIM)
      File_Summary  <- sprintf("./power_test_C=%d_Delta=%.2f/Power_causal_rare/Rare_SumTable_Power_region_size_%d_Neg_beta_pct_%d_Causal_pct_%d_Nsim_%d.csv.csv", test_C, Delta, Region_L, Beta1.Sign * 100, Causal_percent * 100, N.SIM)                                                                            
       }else if (causal_variant == "rare_and_common" && Delta != 0)
         {
         File_Pval     <- sprintf("./power_test_C=%d_Delta=%.2f/Power_causal_rare_and_common/Power_rare_and_common_region_size_%d_Neg_beta_pct_%d_Causal_pct_%d_Nsim_%d.csv", test_C, Delta, Region_L, Beta1.Sign * 100, Causal_percent * 100, N.SIM)
         File_Summary  <- sprintf("./power_test_C=%d_Delta=%.2f/Power_causal_rare_and_common/Rare_and_Common_SumTable_Power_region_size_%d_Neg_beta_pct_%d_Causal_pct_%d_Nsim_%d.csv.csv", test_C, Delta, Region_L, Beta1.Sign * 100, Causal_percent * 100, N.SIM)                                                                            
              }
  
   if (causal_variant == "rare" && Delta == 0) 
      {
      File_Pval     <- sprintf("./power_test_C=%d/Power_causal_rare/Power_rare_region_size_%d_Neg_beta_pct_%d_Causal_pct_%d_Nsim_%d.csv", test_C, Region_L, Beta1.Sign * 100, Causal_percent * 100, N.SIM)
      File_Summary  <- sprintf("./power_test_C=%d/Power_causal_rare/Rare_SumTable_Power_region_size_%d_Neg_beta_pct_%d_Causal_pct_%d_Nsim_%d.csv.csv", test_C, Region_L, Beta1.Sign * 100, Causal_percent * 100, N.SIM)                                                                            
          }else if (causal_variant == "rare_and_common" && Delta == 0)
         {
         File_Pval     <- sprintf("./power_test_C=%d/Power_causal_rare_and_common/Power_rare_and_common_region_size_%d_Neg_beta_pct_%d_Causal_pct_%d_Nsim_%d.csv", test_C, Region_L, Beta1.Sign * 100, Causal_percent * 100, N.SIM)
         File_Summary  <- sprintf("./power_test_C=%d/Power_causal_rare_and_common/Rare_and_Common_SumTable_Power_region_size_%d_Neg_beta_pct_%d_Causal_pct_%d_Nsim_%d.csv.csv", test_C, Delta, Region_L, Beta1.Sign * 100, Causal_percent * 100, N.SIM)                                                                                
           }        
   write.csv(pval, file = File_Pval, row.names = ROW_NAME)   
   # Remove row 1
   sumTable <- sumTable[-1,]
   # Reset the row names
   rownames(sumTable) <- NULL
   write.csv(sumTable, file = File_Summary)

   }

