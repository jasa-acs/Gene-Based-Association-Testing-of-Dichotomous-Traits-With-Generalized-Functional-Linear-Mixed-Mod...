# Set the working dirctory to the location of this file
# when it is sourced:
this.dir <- dirname(parent.frame(2)$ofile)
this.dir
setwd(this.dir)

# setwd("C:/NICHD/Research/paper/PedGFLMM/simulations")
source("function_aggregate_sum_table.R")

Date   = c("COSI_25_Ped", "COSI_25_Ped_basis=6_7")
Vnt    = c("Rare","Rare_and_Common")
Length = 3000*(2:7)

Nseed = 250
Nsim  = 4000

for (i in 1:length(Date)){
   for (j in 1:length(Vnt)){
      for (k in 1:length(Length)){
   sum_table(date= Date[i], vnt=Vnt[j], length=Length[k], Nseed = Nseed, Nsim=Nsim)
                             }
                          }
                       }
