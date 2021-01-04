# This needs to be run with the working directory set to
#   simulations/Type_I_error_rate_aggregate
# But after it runs, the working directory will be
#   simulations

# Set the working dirctory to the location of this file
# when it is sourced:
this.dir <- dirname(parent.frame(2)$ofile)
this.dir
setwd(this.dir)

# With the current settings, the output will be generated in the two sub-folders
#   Type_I_Error_Rare
# and
#   Type_I_Error_Rare_and_Common
# of simulations/COSI_50_Ped/Type_I_Error_combined_Results
source("50_Ped_function_Type_I_aggregate.R")


date="COSI_50_Ped"
#### 1200 seeds
variant_type=c("Rare", "Rare_and_Common")
region_size=3000*(2:7)

Nseed = 1200
Nsim  = 2500 

for (i in 1:length(date)){
   for(j in 1:length(variant_type)){
      for (k in 1:length(region_size)){

        tryCatch({
         Type_I_aggregate(date=date[i], variant_type=variant_type[j], region_size=region_size[k], reps=Nseed, sims = Nsim)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

                            }
                         }
                      }