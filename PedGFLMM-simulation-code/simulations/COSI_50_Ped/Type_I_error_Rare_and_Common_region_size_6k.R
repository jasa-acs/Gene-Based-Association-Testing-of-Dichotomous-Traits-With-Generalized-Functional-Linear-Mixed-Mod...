# This will output a summary file into the folder:
# simulations/COSI_50_Ped/Type_I_Error_Rare_and_Common

# Set the working dirctory to the location of this file
# when it is sourced:
this.dir <- dirname(parent.frame(2)$ofile)
this.dir
setwd(this.dir)

# source("C:/NICHD/Research/paper/PedGFLMM/JASA_Supplementary/simulations/COSI_June_50_Ped/Type_I_error_func.R")
source("Type_I_error_func.R")


DUMX=5  # Number of simulations (we used 2500)
DUMZ=1  # Seed number (we used 1200)

# To time how long this takes:
ptm <- proc.time()
print(paste0("Start time: ",Sys.time()))

Type_I_error_rates = Type_I_error_func (idx2 = 2, N.SIM = DUMX, causal_variant = "rare_and_common", seed =DUMZ)

print(paste0("End time: ",Sys.time()))
# Print out the timing results
print(proc.time() - ptm)