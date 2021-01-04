#UPDATED BY Fan 02/20/2019
#CREATED BY CHIYANG 04/16/2017
#FOR COSI_June_50_Ped

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
# of simulations/COSI_June_50_Ped/Type_I_Error_combined_Results

source("50_Ped_function_Type_I_aggregate.R")

# setwd("C:/NICHD/Research/paper/PedGFLMM/JASA_Supplementary/simulations")
setwd("..")


# source("C:/NICHD/Research/paper/PedGFLMM/JASA_Supplementary/simulations/Type_I_error_rate_aggregate/50_Ped_function_Type_I_aggregate.R")

date="COSI_June_50_Ped"
variant_type=c("Rare", "Rare_and_Common")
region_size=c(6000)

for (i in 1:length(date)){
   for(j in 1:length(variant_type)){
      for (k in 1:length(region_size)){

        tryCatch({
         Type_I_aggregate(date=date[i], variant_type=variant_type[j], region_size=region_size[k], reps1=400, sims1 = 2500, reps2=20, sims2 = 500)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

                            }
                         }
                      }

