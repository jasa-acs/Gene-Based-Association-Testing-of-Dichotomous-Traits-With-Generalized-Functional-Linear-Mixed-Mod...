#UPDATED BY CHIYANG 08/26/2019, fOR COSI_50_Ped
#UPDATED BY CHIYANG 08/19/2019
#UPDATED BY CHIYANG 05/02/2017
#CREATED BY CHIYANG 04/16/2017

Type_I_aggregate=function(date, variant_type,region_size, reps=250, sims = 4000)
   {
   f_list1=list()
   dir1=sprintf("%s/Type_I_Error_%s/Type_I_error_region_size_%d_Nsim_%d",date,variant_type,region_size, sims)

   for(i in 1:reps)
         {
         f1 = NULL
         oname1 = paste(dir1,"_seed_",i, sep="")

         tryCatch({
            f1=assign(oname1, read.csv(paste(oname1, ".csv", sep="")))
            f1=as.matrix(f1)

            vars <- colnames(f1) %in% c("X", "bcnt_LRT_fixed", "fcnt_LRT_fixed", "bcnt_LRT_beta", "fcnt_LRT_beta",
                            "bpval_LRT_fixed_total", "fpval_LRT_fixed_total", "bpval_LRT_beta_total", "fpval_LRT_beta_total")            
            f1 =f1[,vars]
               }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      f_list1[[i]] = f1
      }

   out1=Reduce("+", f_list1)  # add by CY 08/26/2019
   out1 = as.data.frame(out1)

   out = NULL
   out$alpha = out1[,1]/reps

   # Empirical Type I error rate  
   out$fixed_LRT_bpval =  out1$bcnt_LRT_fixed/out1$bpval_LRT_fixed_total       
   out$fixed_LRT_fpval =  out1$fcnt_LRT_fixed/out1$fpval_LRT_fixed_total      
   out$beta_LRT_bpval  =  out1$bcnt_LRT_beta/out1$bpval_LRT_beta_total
   out$beta_LRT_fpval  =  out1$fcnt_LRT_beta/out1$fpval_LRT_beta_total 
   
   # number of available p-values          
   out$bpval_LRT_fixed_total =  out1$bpval_LRT_fixed_total       
   out$fpval_LRT_fixed_total =  out1$fpval_LRT_fixed_total      
   out$bpval_LRT_beta_total  =  out1$bpval_LRT_beta_total
   out$fpval_LRT_beta_total  =  out1$fpval_LRT_beta_total 

   # Total number of simulations
   out$nsim = reps*sims 
  
   dir=sprintf("%s/Type_I_Error_combined_Results/Type_I_Error_%s/Type_I_error_region_size_%d_Nsim_%d",date,variant_type,region_size, reps*sims )
   dir_out=paste(dir,".csv",sep="")
   write.csv(out,dir_out,row.names=FALSE)

      }