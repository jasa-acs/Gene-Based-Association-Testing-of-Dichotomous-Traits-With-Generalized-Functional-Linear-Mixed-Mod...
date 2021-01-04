# This function needs to be run from within the same folder
# that contains the "bar_dynamic_function_plot.R" file.

# Set the working dirctory to the location of this file
# when it is sourced:
this.dir <- dirname(parent.frame(2)$ofile)
this.dir
setwd(this.dir)

source("bar_dynamic_function_plot.R")

vtype=c("rare_and_common", "rare")
path=c("COSI_50_Ped/power_test_C=1_Delta=0.50")

# This function reads its input from csv files in two sub-folders of the folder
# defined by the 'path' variable.  These two sub-folders are named:
#
# Power_causal_rare
# Power_causal_rare_and_common
#
# and this function expects to write its plots to the sub-folder 'Fig_Pow' of
# the folder defined by the 'path' variable.

alpha=0.0001
region_size=c(6000,9000,12000,15000,18000,21000)

cat("Writing bar plots to the folder:\n")
cat(paste0(path,"/Fig_pow/\n"))


for (k in 1:length(path))
   {
   for (i in 1:length(vtype))
      {
      for (j in 1:length(region_size))
         {
         tryCatch({
         bar_plot(vtype=vtype[i], path=path[k], alpha=alpha, region_size=region_size[j])
                       },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
         }
      }
   }

