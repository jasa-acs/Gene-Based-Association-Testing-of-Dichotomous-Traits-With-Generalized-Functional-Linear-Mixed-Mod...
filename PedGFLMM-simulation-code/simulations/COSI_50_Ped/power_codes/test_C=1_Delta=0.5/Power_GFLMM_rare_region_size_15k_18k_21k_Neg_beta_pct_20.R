# The output generated from running this file will be written in the folder:
#   simulations/COSI_June_50_Ped/power_test_C=1_Delta=0.50/Power_causal_rare

# Set the working directory to the location of this file
# when it is sourced:
this.dir <- dirname(parent.frame(2)$ofile)
this.dir
setwd(this.dir)


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
  
# source("C:/NICHD/Research/paper/PedGFLMM/JASA_Supplementary/simulations/COSI_June_50_Ped/Power_func.R")
source("../../Power_func.R")

power_idx1_EQ_1_idx4_EQ_1 = power_func(idx1 = 1, idx2 = 5, idx4 = 2, idx5 = 1, causal_variant = "rare", test_C = 1, Delta = 0.5)
power_idx1_EQ_2_idx4_EQ_1 = power_func(idx1 = 2, idx2 = 5, idx4 = 2, idx5 = 1, causal_variant = "rare", test_C = 1, Delta = 0.5)
power_idx1_EQ_3_idx4_EQ_1 = power_func(idx1 = 3, idx2 = 5, idx4 = 2, idx5 = 1, causal_variant = "rare", test_C = 1, Delta = 0.5)

###
power_idx1_EQ_1_idx4_EQ_1 = power_func(idx1 = 1, idx2 = 6, idx4 = 2, idx5 = 1, causal_variant = "rare", test_C = 1, Delta = 0.5)
power_idx1_EQ_2_idx4_EQ_1 = power_func(idx1 = 2, idx2 = 6, idx4 = 2, idx5 = 1, causal_variant = "rare", test_C = 1, Delta = 0.5)
power_idx1_EQ_3_idx4_EQ_1 = power_func(idx1 = 3, idx2 = 6, idx4 = 2, idx5 = 1, causal_variant = "rare", test_C = 1, Delta = 0.5)

power_idx1_EQ_1_idx4_EQ_1 = power_func(idx1 = 1, idx2 = 7, idx4 = 2, idx5 = 1, causal_variant = "rare", test_C = 1, Delta = 0.5)
power_idx1_EQ_2_idx4_EQ_1 = power_func(idx1 = 2, idx2 = 7, idx4 = 2, idx5 = 1, causal_variant = "rare", test_C = 1, Delta = 0.5)
power_idx1_EQ_3_idx4_EQ_1 = power_func(idx1 = 3, idx2 = 7, idx4 = 2, idx5 = 1, causal_variant = "rare", test_C = 1, Delta = 0.5)
