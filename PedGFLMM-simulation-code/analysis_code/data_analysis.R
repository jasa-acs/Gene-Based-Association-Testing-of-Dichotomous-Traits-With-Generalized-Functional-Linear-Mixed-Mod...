#### Rare only, BY CY,03/10/2020
rm(list=ls())

v_type = "R"

memory.limit(size = 40650000)
library(tidyverse)
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

#### READ IN FUNCTIONS AND DADA
DIR <- "C:/Users/cchiu1/Desktop/GitHub/JASA_Supplementary"
#DIR <-  "/data/chiuc2/NICHD/Research/paper/PedGFLMM/JASA_Supplementary"

setwd(DIR)
require(pedgene)
require(pedigreemm)

source("PedGFLMM/dynamic_rule.R")
source("PedGFLMM/MGAO.R")
source("PedGFLMM/PedGFLMM_beta_smooth_only.R")
source("PedGFLMM/PedGFLMM_fixed_model.R")

#####
genoData <- get(load("AMD/data/exomeAMD_cleaned_noMissing.RData"))
refGene <- read.table("AMD/data/refseq_genes.txt", header= TRUE, stringsAsFactors= FALSE)
snpAnnot <- get(load("AMD/data/snpAnnotDataV5.RData"))
snpAnnot <- snpAnnot[snpAnnot$filterReason == "okay", ]

## generate pedigree structure ##
schaidPed <- genoData[ ,c("FID","IID","PAT","MAT","SEX","PHENOTYPE")]
colnames(schaidPed) <- c("ped","person","father","mother","sex","trait")

trait_tmp = schaidPed$trait
trait_tmp[trait_tmp==-9]=NA     # IN RAW DATA "-9" MEANS "NA"
schaidPed$trait = trait_tmp - 1 # TRANSFORM TRAIT FROM 1,2 TO 0,1

pedPer <- schaidPed[ ,c("ped","person")]
genoData$PAT <- genoData$MAT <- genoData$SEX <- genoData$PHENOTYPE <- NULL
genoData$FID <- genoData$IID <- NULL
colnames(genoData) <- sub("_.*","",colnames(genoData))

#### SETTING FOR FUNCTIONAL MODELS
flmPed = schaidPed
flmPed$ID = 1:dim(schaidPed)[[1]] 
covariate = data.frame(sped=schaidPed$ped,person=schaidPed$person,sex=as.factor(schaidPed$sex))

if (v_type =="R"){##extract rare variants for analysis ##
rareMarker <- snpAnnot[snpAnnot$MAF <0.05, "marker"]
rareMarker <- intersect(rareMarker,colnames(genoData))
genoData <- genoData[ ,rareMarker]
}

refGene$name <- refGene$cdsStart <-refGene$cdsEnd <- NULL
refGene <- refGene[nchar(refGene$chrom) <= 5, ]
refGene <- refGene[refGene$chrom != "chrX" & refGene$chrom != "chrY", ]
refGene <- refGene[!duplicated(refGene), ]
nGene <- nrow(refGene)

posAnnotData <- snpAnnot[ ,c("marker","rsid_Hum","pos_Hum","chrom")]

data_analysis = function(LB=1, UB= nGene, cf1 = 2, cf2 = 2){
for (z in LB:UB) {
  tryCatch({  
  testGene <- refGene[z, ]
  
  ## get positions of variants within a gene ##
  posAnnot <- posAnnotData[posAnnotData$pos_Hum >= testGene$txStart & posAnnotData$pos_Hum <= testGene$txEnd
			   & posAnnotData$chrom %in% testGene$chrom, ]
  selectMarker <- intersect(colnames(genoData),posAnnot$marker)

  if (length(selectMarker) > cf1){
    geno <- genoData[ ,selectMarker]
 
    snpNumber <- ncol(geno)
    for (k in 1:snpNumber) {
   	marker <- geno[ ,k]
    	genotype0 <- sum(marker == 0, na.rm= TRUE)
    	genotype2 <- sum(marker == 2, na.rm= TRUE)
    	if (genotype0 < genotype2) {
    		geno[ ,k] <- 2-geno[ ,k]
    	}    # IMPORTANT: rare allele coded as 2 #
    }

    geno[is.na(geno)] <- 0
    pos <- posAnnot[posAnnot$marker %in% selectMarker, "pos_Hum"]
  
    maf <- colMeans(geno)
    pos <- pos[maf > 0]          
    if (length(pos) >= cf2) {     # cf2 = minimum number of variants #
      geno <- geno[ ,maf > 0]     # remove nonpolymorphic variants #
      nsnp    <- ncol(geno)
      weight <- rep(1, ncol(geno))
     
      pedgeno <- cbind(pedPer, geno)
  
      ## READOUT THE RESULTS##
      chr   <- testGene$chrom
      gene  <- testGene$name2
      start <- testGene$txStart
      end   <- testGene$txEnd

      sumTable = NULL
      
      #### use dynamic rule to determine order
      order = dRule(geno)$order
      
      #### FUNCTIONAL MODELS
      ## bspline + fixed
      msgs <- NA
      msgs <- capture.output(myWarnErr <- myTryCatch(PedGFLMM_fixed_model(ped = flmPed, 
              geno = pedgeno, covariate = as.matrix(covariate), pos = pos, 
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

      sumTable <- bind_rows(sumTable,data.frame(
      chr=chr, gene=gene, start=start, end=end, geneID=z,
      type="bFixed",n.poly=nsnp,n.vars=snpNumber,M_gao=bpval_PedGFLMM_fixed$M_gao,
      betabasis=bpval_PedGFLMM_fixed$nbetabasis, 
      genobasis=bpval_PedGFLMM_fixed$ngenobasis,
      LRT=bpval_PedGFLMM_fixed$LRT, 
      warning=msg.w,
      error=msg.e,
      warnings=paste0(msgs,";",collapse=""),stringsAsFactors = FALSE))

      ## fspline + fixed
      msgs <- NA
      msgs <- capture.output(myWarnErr <- myTryCatch(PedGFLMM_fixed_model(ped = flmPed, 
              geno = pedgeno, covariate = as.matrix(covariate), pos = pos, 
              order = order, beta_basis = NULL, geno_basis = NULL, 
              base = "fspline", optimizer = "bobyqa", Wald = FALSE)),type="message")

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

      sumTable <- bind_rows(sumTable,data.frame(
      chr=chr, gene=gene, start=start, end=end, geneID=z,
      type="fFixed",n.poly=nsnp,n.vars=snpNumber,M_gao=fpval_PedGFLMM_fixed$M_gao,
      betabasis=fpval_PedGFLMM_fixed$nbetabasis, 
      genobasis=fpval_PedGFLMM_fixed$ngenobasis,
      LRT=fpval_PedGFLMM_fixed$LRT, 
      warning=msg.w,
      error=msg.e,
      warnings=paste0(msgs,";",collapse=""),stringsAsFactors = FALSE))

      ## bspline + beta
      msgs <- NA
      msgs <- capture.output(myWarnErr <- myTryCatch(PedGFLMM_beta_smooth_only(ped = flmPed, 
              geno = pedgeno, covariate = as.matrix(covariate), pos = pos, 
              order = order, beta_basis = NULL, 
              base = "bspline", optimizer = "bobyqa", Wald = FALSE)),type="message")

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

      sumTable <- bind_rows(sumTable,data.frame(
      chr=chr, gene=gene, start=start, end=end, geneID=z,
      type="bBetaSmooth",n.poly=nsnp,n.vars=snpNumber,M_gao=bpval_PedGFLMM_beta_smooth_only$M_gao,
      betabasis=bpval_PedGFLMM_beta_smooth_only$nbetabasis, 
      genobasis=NA,
      LRT=bpval_PedGFLMM_beta_smooth_only$LRT, 
      warning=msg.w,
      error=msg.e,
      warnings=paste0(msgs,";",collapse=""),stringsAsFactors = FALSE))
            
      ## fspline + beta
      msgs <- NA
      msgs <- capture.output(myWarnErr <- myTryCatch(PedGFLMM_beta_smooth_only(ped = flmPed, 
              geno = pedgeno, covariate = as.matrix(covariate), pos = pos, 
              order = order, beta_basis = NULL, 
              base = "fspline", optimizer = "bobyqa", Wald = FALSE)),type="message")

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

      sumTable <- bind_rows(sumTable,data.frame(
      chr=chr, gene=gene, start=start, end=end, geneID=z,
      type="fBetaSmooth",n.poly=nsnp,n.vars=snpNumber,M_gao=fpval_PedGFLMM_beta_smooth_only$M_gao,
      betabasis=fpval_PedGFLMM_beta_smooth_only$nbetabasis, 
      genobasis=NA,
      LRT=fpval_PedGFLMM_beta_smooth_only$LRT, 
      warning=msg.w,
      error=msg.e,
      warnings=paste0(msgs,";",collapse=""),stringsAsFactors = FALSE))
             
      FILE = sprintf("AMD/analysis/results/%s/gene_%d.csv", v_type, z)
      write.csv(sumTable, file= FILE, row.names= FALSE)     
      } 
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})  
     print(z)
}     
}

#### run the analysis
data_analysis(LB=1, UB=5, cf1 = 1, cf2 = 2)
