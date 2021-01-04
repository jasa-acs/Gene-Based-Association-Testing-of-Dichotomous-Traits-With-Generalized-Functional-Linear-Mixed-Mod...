
######################
### Init Function  ###
######################

Init <- function(FILE, seed=-1)
   {
   nMarker<-0
   nSample<-0
  
   if (seed < 0)
      {
      seed = as.integer(runif(1)*10000)
      }
  
   temp <- .C("Init",as.integer(nMarker) ,as.integer(nSample), as.character(BIN_FILE),as.integer(seed))

   nMarker <- temp[[1]]
   nSample <- temp[[2]]

   c(nMarker,nSample)
   }

#
# Get ALL Data (Sample)
#
Get_DATA_ALL <- function(Marker_Idx, N)
   {
   n <- length(Marker_Idx)
   Genotype <- rep(0, n*N);
   temp <- .C("Get_DATA_ALL", as.integer(n), as.integer(Marker_Idx), as.integer(Genotype))

   out <- list()
   out[[1]] <- matrix(temp[[3]], byrow=TRUE, nrow=N);
   out[[2]] <- temp[[3]]

   out
   }

Get_Haplotype<-function(Marker_Idx, Sample_Idx)
   {
   p <- length(Marker_Idx)
   n <- dim(Sample_Idx)[1]

   Haplotype <- rep(0, n*p);
   temp <- .C("Get_Haplotype", as.integer(p), as.integer(Marker_Idx)
   , as.integer(n), as.integer(Sample_Idx[,1])
   , as.integer(Haplotype))

   out <- temp[[5]]

   matrix(out, byrow=TRUE, ncol=p)
   }

#
# Set effect
#
Set_Beta <- function(Marker_Idx, beta1, N)
   {
   n <- length(Marker_Idx)
   OUT_Effect <- rep(0,N);
   temp <- .C("Set_Beta", as.integer(n), as.integer(Marker_Idx), as.double(beta1), as.double(OUT_Effect))

   temp[[4]]
   }

#
# When each sample has different beta0 due to other covariates. 
#
Set_Beta_wBeta0 <- function(Marker_Idx, beta1, beta0, N)
   {
   n <- length(Marker_Idx)
   OUT_Effect <- rep(0,N)
   temp <- .C("Set_Beta_wBeta0", as.integer(n), as.integer(Marker_Idx), as.double(beta1), as.double(beta0), as.double(OUT_Effect))

   temp[[5]]
   }

# Simulate Case
#
Simu_Case <- function(nCase)
   {
   Idx1 <- rep(-1,nCase)
   Idx2 <- rep(-1,nCase)
   temp <- .C("Simu_Case", as.integer(nCase), as.integer(Idx1), as.integer(Idx2))

   out  <- cbind(temp[[2]] + 1, temp[[3]] + 1)
   out
   }

# Simulate Control
#
Simu_Control <- function(nControl)
   {
   Idx1 <- rep(-1, nControl)
   Idx2 <- rep(-1,nControl)
   temp <- .C("Simu_Control", as.integer(nControl), as.integer(Idx1), as.integer(Idx2))

   out  <- cbind(temp[[2]] + 1, temp[[3]] + 1)
   out
   }

# Simulate Case
#
Simu_Case_wCov<-function(nCase, beta_cov){

  Idx1<-rep(-1,nCase);
  Idx2<-rep(-1,nCase);
  COV<-rep(-1,nCase);
  temp<-.C("Simu_Case_wCov",as.integer(nCase),as.double(beta_cov),as.integer(Idx1), as.integer(Idx2)
  ,as.integer(COV))

  out<-cbind(temp[[3]] +1,temp[[4]] +1,temp[[5]])
  out
  }

# Simulate Control
#
Simu_Control_wCov <- function(nControl,beta_cov)
   {
   Idx1 <- rep(-1,nControl)
   Idx2 <- rep(-1,nControl)
   COV  <- rep(-1,nControl)
   temp <- .C("Simu_Control_wCov", as.integer(nControl), as.double(beta_cov), as.integer(Idx1)
   , as.integer(Idx2), as.integer(COV))

   out  <- cbind(temp[[3]] + 1,temp[[4]] + 1,temp[[5]])
   out
   }
