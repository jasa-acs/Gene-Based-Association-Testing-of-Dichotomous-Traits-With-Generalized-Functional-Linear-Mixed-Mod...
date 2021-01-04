#############################################
# Make batch file with different index
# idx1: causal percent, look at Parameters.r for
#       Causal_percent.A <- c(0.05, 0.10, 0.15)
#       Causal_percent   <- Causal_percent.A[idx1]
# So: idx1 = 1 implies causal percent = 0.05
#     idx1 = 2 implies causal percent = 0.10
#     idx1 = 3 implies causal percent = 0.15
############################################
# idx2 : region size,  look at Parameters.r for
# Range.A  <- c(3000, 6000, 9000, 12000, 15000, 18000, 21000, 24000, 27000, 30000)
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

Beta1.Sign.A <- c(0, 0.2, 0.5)  # % of betas with negative sign
CutOff.A     <- c(0.03, 0.04, 0.05, 0.06)

Causal_percent.A <- c(0.05, 0.10, 0.15)
Range.A          <- c(3000, 6000, 9000, 12000, 15000, 18000, 21000, 24000, 27000, 30000)
#N.Sample.A       <- c(250, 500, 1000)
Is.null.a        <- c(0, 1)
Beta1.IDX        <- 2
CutOff           <- 0.03

#Causal_percent <- Causal_percent.A[idx1]

#Region_L   <- Range.A[idx2]   # Region length, total 4
#Beta1.Sign <- Beta1.Sign.A[idx4]
Is.null    <- 0

Weight_Type <- 1 # beta weights

###############################
# Beta function

Get_LogRR_MAF <- function(MAF, p_min = 0.15, theta_max =4 ,s=4)
   {
   re <- 1 + (theta_max -1) * exp( -s*MAF/p_min)
   
   return(log(re))
   }

Generate_Covariates <- function(n)
   {	
	 X1 <- rnorm(n)
	 X2 <- rbinom(n,1,0.5)
	
	 return(cbind(X1,X2))
   }

Get_Alpha0_C <- function(X)
   {	
	 n  <- dim(X)[1]
	 re <- X[,1] * 0.5 / sqrt(2) + Beta0 /2
	 
   return(re)
   }	

Get_Alpha0_Q <- function(X)
   {
	 re <- X[,1] * 0.5 + X[,2] * 0.5 
	 return(re)
   }	

Get_Weights <- function(MAF, IDX)
   {
	 if(IDX == 1)
      {
		  re <- dbeta(MAF, shape1 = 1, shape2 = 25)
	    } 
      else if (IDX == 2)
         {		
		     re <- 1/sqrt(MAF*(1 - MAF))
         }

	 return(re)
   }

Get_Beta_Q_MAF <- function(MAF)
   {
	 re <- abs(0.5 * log10(MAF))
	 return(re)
   }

Get_BETA_Q<-function( IDX, MAF)
   {
	 n <- length(MAF)

	 if(IDX == 1)
      {
		  re<-Get_Beta_Q_MAF(MAF) * log(7)
	    } 
      else if (IDX == 2)
         {
		     re <- Get_Beta_Q_MAF(MAF) * log(5)
	       } 
         else if (IDX == 3)
            {
		        re <- Get_Beta_Q_MAF(MAF) * log(2.0)
	          } 

	 return(re)
   }

Get_BETA_C<-function(idx1, idx2, MAF)
   {
	 n <- length(MAF)
   if (idx2 == 1)
      {
	    if (idx1 == 1)
         {
		     re <- Get_Beta_Q_MAF(MAF) * log(90)
	       } else if (idx1 == 2)
            {
		        re <- Get_Beta_Q_MAF(MAF) * log(70)
            } else if (idx1 == 3)
               {
		           re <- Get_Beta_Q_MAF(MAF) * log(50) 
	             }
	    } else if (idx2 == 2)
         {
	       if (idx1 == 1)
            {
		        re <- Get_Beta_Q_MAF(MAF) * log(90)/1.25
	          } else if (idx1 == 2)
               {
		           re <- Get_Beta_Q_MAF(MAF) * log(70)/1.25
               } else if (idx1 == 3)
                  {
		              re <- Get_Beta_Q_MAF(MAF) * log(50)/1.25 
	                }
	       } else if (idx2 == 3)
            {
	          if (idx1 == 1)
               {
		           re <- Get_Beta_Q_MAF(MAF) * log(90)/1.5
	             } else if (idx1 == 2)
                  {
		              re <- Get_Beta_Q_MAF(MAF) * log(70)/1.5
                  } else if (idx1 == 3)
                     {
		                 re <- Get_Beta_Q_MAF(MAF) * log(50)/1.5 
	                   }
	          } else if (idx2 == 4)
               {
	             if (idx1 == 1)
                  {
		              re <- Get_Beta_Q_MAF(MAF) * log(90)/1.75
	                } else if (idx1 == 2)
                     {
		                 re <- Get_Beta_Q_MAF(MAF) * log(70)/1.75
                     } else if (idx1 == 3)
                        {
		                    re <- Get_Beta_Q_MAF(MAF) * log(50)/1.75 
	                      }
	             } else if (idx2 == 5)
                  {
	                if (idx1 == 1)
                     {
		                 re <- Get_Beta_Q_MAF(MAF) * log(90)/2.0
	                   } else if (idx1 == 2)
                        {
		                    re <- Get_Beta_Q_MAF(MAF) * log(70)/2.0
                        } else if (idx1 == 3)
                           {
		                       re <- Get_Beta_Q_MAF(MAF) * log(50)/2.0 
	                         }
	                } else if (idx2 == 6)
                     {
	                   if (idx1 == 1)
                        {
		                    re <- Get_Beta_Q_MAF(MAF) * log(90)/2.25
	                      } else if (idx1 == 2)
                           {
		                       re <- Get_Beta_Q_MAF(MAF) * log(70)/2.25
                           } else if (idx1 == 3)
                              {
		                          re <- Get_Beta_Q_MAF(MAF) * log(50)/2.25 
                              }
                     } else if (idx2 == 7)
                        {
	                      if (idx1 == 1)
                           {
		                       re <- Get_Beta_Q_MAF(MAF) * log(90)/2.5
	                         } else if (idx1 == 2)
                              {
                              re <- Get_Beta_Q_MAF(MAF) * log(70)/2.5
                              } else if (idx1 == 3)
                                 {
		                             re <- Get_Beta_Q_MAF(MAF) * log(50)/2.5 
                                 } 
                            } else if (idx2 == 8)
                        {
	                      if (idx1 == 1)
                           {
		                       re <- Get_Beta_Q_MAF(MAF) * log(90)/2.75
	                         } else if (idx1 == 2)
                              {
                              re <- Get_Beta_Q_MAF(MAF) * log(70)/2.75
                              } else if (idx1 == 3)
                                 {
		                             re <- Get_Beta_Q_MAF(MAF) * log(50)/2.75 
                                 }
                                  }else if (idx2 == 9)
                        {
	                      if (idx1 == 1)
                           {
		                       re <- Get_Beta_Q_MAF(MAF) * log(90)/3
	                         } else if (idx1 == 2)
                              {
                              re <- Get_Beta_Q_MAF(MAF) * log(70)/3
                              } else if (idx1 == 3)
                                 {
		                             re <- Get_Beta_Q_MAF(MAF) * log(50)/3 
                                 }
                                 }else if (idx2 == 10)
                        {
	                      if (idx1 == 1)
                           {
		                       re <- Get_Beta_Q_MAF(MAF) * log(90)/3.25 
	                         } else if (idx1 == 2)
                              {
                              re <- Get_Beta_Q_MAF(MAF) * log(70)/3.25 
                              } else if (idx1 == 3)
                                 {
		                             re <- Get_Beta_Q_MAF(MAF) * log(50)/3.25 
                                 }
                        }else { } 
   return(re)
   }

###
Get_BETA_C_Delta<-function(idx1, idx2, MAF, Delta)
   {
	 n <- length(MAF)
   if (idx2 == 1)
      {
	    if (idx1 == 1)
         {
		     re <- Get_Beta_Q_MAF(MAF) * log(90)/(1.0 + Delta)
	       } else if (idx1 == 2)
            {
		        re <- Get_Beta_Q_MAF(MAF) * log(70) /(1.0 + Delta)
            } else if (idx1 == 3)
               {
		           re <- Get_Beta_Q_MAF(MAF) * log(50) /(1.0 + Delta)
	             }
	    } else if (idx2 == 2)
         {
	       if (idx1 == 1)
            {
		        re <- Get_Beta_Q_MAF(MAF) * log(90)/(1.25 + Delta)
	          } else if (idx1 == 2)
               {
		           re <- Get_Beta_Q_MAF(MAF) * log(70)/(1.25 + Delta)
               } else if (idx1 == 3)
                  {
		              re <- Get_Beta_Q_MAF(MAF) * log(50)/(1.25 + Delta) 
	                }
	       } else if (idx2 == 3)
            {
	          if (idx1 == 1)
               {
		           re <- Get_Beta_Q_MAF(MAF) * log(90)/(1.5 + Delta)
	             } else if (idx1 == 2)
                  {
		              re <- Get_Beta_Q_MAF(MAF) * log(70)/(1.5 + Delta)
                  } else if (idx1 == 3)
                     {
		                 re <- Get_Beta_Q_MAF(MAF) * log(50)/(1.5 + Delta) 
	                   }
	          } else if (idx2 == 4)
               {
	             if (idx1 == 1)
                  {
		              re <- Get_Beta_Q_MAF(MAF) * log(90)/(1.75 + Delta)
	                } else if (idx1 == 2)
                     {
		                 re <- Get_Beta_Q_MAF(MAF) * log(70)/(1.75 + Delta)
                     } else if (idx1 == 3)
                        {
		                    re <- Get_Beta_Q_MAF(MAF) * log(50)/(1.75 + Delta) 
	                      }
	             } else if (idx2 == 5)
                  {
	                if (idx1 == 1)
                     {
		                 re <- Get_Beta_Q_MAF(MAF) * log(90)/(2.0 + Delta)
	                   } else if (idx1 == 2)
                        {
		                    re <- Get_Beta_Q_MAF(MAF) * log(70)/(2.0 + Delta)
                        } else if (idx1 == 3)
                           {
		                       re <- Get_Beta_Q_MAF(MAF) * log(50)/(2.0 + Delta) 
	                         }
	                } else if (idx2 == 6)
                     {
	                   if (idx1 == 1)
                        {
		                    re <- Get_Beta_Q_MAF(MAF) * log(90)/(2.25 + Delta)
	                      } else if (idx1 == 2)
                           {
		                       re <- Get_Beta_Q_MAF(MAF) * log(70)/(2.25 + Delta)
                           } else if (idx1 == 3)
                              {
		                          re <- Get_Beta_Q_MAF(MAF) * log(50)/(2.25 + Delta) 
                              }
                     } else if (idx2 == 7)
                        {
	                      if (idx1 == 1)
                           {
		                       re <- Get_Beta_Q_MAF(MAF) * log(90)/(2.5 + Delta)
	                         } else if (idx1 == 2)
                              {
                              re <- Get_Beta_Q_MAF(MAF) * log(70)/(2.5 + Delta)
                              } else if (idx1 == 3)
                                 {
		                             re <- Get_Beta_Q_MAF(MAF) * log(50)/(2.5 + Delta) 
                                 } 
                            } else if (idx2 == 8)
                        {
	                      if (idx1 == 1)
                           {
		                       re <- Get_Beta_Q_MAF(MAF) * log(90)/(2.75 + Delta)
	                         } else if (idx1 == 2)
                              {
                              re <- Get_Beta_Q_MAF(MAF) * log(70)/(2.75 + Delta)
                              } else if (idx1 == 3)
                                 {
		                             re <- Get_Beta_Q_MAF(MAF) * log(50)/(2.75 + Delta) 
                                 }
                                  }else if (idx2 == 9)
                        {
	                      if (idx1 == 1)
                           {
		                       re <- Get_Beta_Q_MAF(MAF) * log(90)/(3 + Delta)
	                         } else if (idx1 == 2)
                              {
                              re <- Get_Beta_Q_MAF(MAF) * log(70)/(3 + Delta)
                              } else if (idx1 == 3)
                                 {
		                             re <- Get_Beta_Q_MAF(MAF) * log(50)/(3 + Delta) 
                                 }
                                 }else if (idx2 == 10)
                        {
	                      if (idx1 == 1)
                           {
		                       re <- Get_Beta_Q_MAF(MAF) * log(90)/(3.25 + Delta)
	                         } else if (idx1 == 2)
                              {
                              re <- Get_Beta_Q_MAF(MAF) * log(70)/(3.25 + Delta)
                              } else if (idx1 == 3)
                                 {
		                             re <- Get_Beta_Q_MAF(MAF) * log(50)/(3.25 + Delta) 
                                 }
                        }else { } 
   return(re)
   }

