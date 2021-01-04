#### To combine SumTable

# Set the working dirctory to the location of this file
# when it is sourced:
this.dir <- dirname(parent.frame(2)$ofile)
this.dir
setwd(this.dir)

sum_table = function(date= "COSI_50_Ped", vnt="Rare", length=6000, Nseed = 100, Nsim=100){

   Add = Bfixed = Ffixed = Bbeta = Fbeta = NULL
   length_1=0
   for (i in 1:Nseed){
      if(vnt =="Rare"){
         file_name = sprintf("%s/Type_I_Error_%s/Summary/%s_SumTable_Type_I_region_size_%d_Nsim_%d", date, vnt, vnt, length, Nsim)
                      }
      if(vnt =="Rare_and_Common"){
         file_name = sprintf("%s/Type_I_Error_%s/Summary/%s_SumTable_TypeI_region_size_%d_Nsim_%d", date, vnt, vnt, length, Nsim)
      }

      oname1 = paste(file_name,"_seed_",i, sep="")
      tryCatch({
         tmp = assign(oname1, read.csv(paste(oname1, ".csv", sep="")))
         tmp$rep = tmp$rep + (i-1)*Nsim
         Add     = rbind(Add,    tmp[which(tmp$type=="Additive"),]   )
         Bfixed  = rbind(Bfixed, tmp[which(tmp$type=="Fixed"),]      )
         Ffixed  = rbind(Ffixed, tmp[which(tmp$type=="fFixed"),]     ) 
         Bbeta   = rbind(Bbeta,  tmp[which(tmp$type=="bBetaSmooth"),])
         Fbeta   = rbind(Fbeta,  tmp[which(tmp$type=="fBetaSmooth"),])
         tmp = as.matrix(tmp)
         length_2 = !is.null(tmp)
         length_1 = length_1+length_2
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
     
   }
   
   Add_name   = sprintf("%s/Type_I_Error_combined_Results/Type_I_Error_%s/sum_table_Add_Type_I_error_region_size_%d_Nsim_%d.csv",    date, vnt, length, Nsim*length_1)
   Bfixed_name= sprintf("%s/Type_I_Error_combined_Results/Type_I_Error_%s/sum_table_Bfixed_Type_I_error_region_size_%d_Nsim_%d.csv", date, vnt, length, Nsim*length_1)
   Ffixed_name= sprintf("%s/Type_I_Error_combined_Results/Type_I_Error_%s/sum_table_Ffixed_Type_I_error_region_size_%d_Nsim_%d.csv", date, vnt, length, Nsim*length_1)
   Bbeta_name = sprintf("%s/Type_I_Error_combined_Results/Type_I_Error_%s/sum_table_Bbeta_Type_I_error_region_size_%d_Nsim_%d.csv",  date, vnt, length, Nsim*length_1)
   Fbeta_name = sprintf("%s/Type_I_Error_combined_Results/Type_I_Error_%s/sum_table_Fbeta_Type_I_error_region_size_%d_Nsim_%d.csv",  date, vnt, length, Nsim*length_1)
  
   write.csv(Add,       Add_name,row.names=FALSE)
   write.csv(Bfixed, Bfixed_name,row.names=FALSE)
   write.csv(Ffixed, Ffixed_name,row.names=FALSE)
   write.csv(Bbeta,   Bbeta_name,row.names=FALSE)
   write.csv(Fbeta,   Fbeta_name,row.names=FALSE)
}






