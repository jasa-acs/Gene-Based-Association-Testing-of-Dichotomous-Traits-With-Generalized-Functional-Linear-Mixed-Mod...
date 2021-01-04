#UPDATED BY CHIYANG 05/02/2017
#CREATED BY CHIYANG 04/16/2017
#FOR COSI_June_50_Ped


Type_I_aggregate=function(date, variant_type,region_size, reps1=400, sims1 = 2500, reps2=18, sims2 = 500)
   {
   #### FIRST READ THOSE COMPLETED RESULTS WITH ORINGINAL SETTING "reps1" AND "sims1"
   f_list1=list()
   dir1=sprintf("%s/Type_I_Error_%s/Type_I_error_region_size_%d_Nsim_%d",date,variant_type,region_size, sims1)

   for(i in 1:reps1)
         {
         f1 = NULL
         oname1 = paste(dir1,"_seed_",i, sep="")

         tryCatch({
            f1=assign(oname1, read.csv(paste(oname1, ".csv", sep="")))
            f1=as.matrix(f1)
               }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      f_list1[[i]] = f1
      }

   f_list_out1 = f_list1[!sapply(f_list1, is.null)]
   out1=Reduce("+", f_list_out1)/ length(f_list_out1)
   colnames(out1)[1]="alpha"
   NSIM_TOTAL1 = sims1*length(f_list_out1)
   out = cbind(out1, NSIM = NSIM_TOTAL1)
   Nsim_total = NSIM_TOTAL1

   #### THIS PART IS FOR MAKE-UP RESULTS
   if (length(f_list_out1)!=reps1){
      f_list2=list()
      dir2=sprintf("%s/Type_I_Error_%s/Type_I_error_region_size_%d_Nsim_%d",date,variant_type,region_size, sims2)
      idx = 0
      for(i in 1:reps2)
         {
         f2 = NULL
         oname2 = paste(dir2,"_seed_",(reps1+i), sep="")

         tryCatch({
            f2=assign(oname2, read.csv(paste(oname2, ".csv", sep="")))
            f2=as.matrix(f2)
               }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

      f_list2[[i]] = f2
      if(is.null(f2)!=TRUE) idx = idx+1
      if ((NSIM_TOTAL1+ idx*sims2)== reps1*sims1) break
      }

      f_list_out2 = f_list2[!sapply(f_list2, is.null)]
      out2=Reduce("+", f_list_out2)/ length(f_list_out2)
      colnames(out2)[1]="alpha"
      NSIM_TOTAL2 = sims2*length(f_list_out2)
      out = (out2*NSIM_TOTAL2 + out1*NSIM_TOTAL1)/(NSIM_TOTAL2+NSIM_TOTAL1)
      out = cbind(out, NSIM=NSIM_TOTAL2+NSIM_TOTAL1)
      Nsim_total = NSIM_TOTAL2 + NSIM_TOTAL1
   }
   dir3=sprintf("%s/Type_I_Error_combined_Results/Type_I_Error_%s/Type_I_error_region_size_%d_Nsim_%d",date,variant_type,region_size, Nsim_total)
   dir_out=paste(dir3,".csv",sep="")
   write.csv(out,dir_out,row.names=FALSE)

      }
