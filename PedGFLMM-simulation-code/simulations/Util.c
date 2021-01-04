#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>


/*#define BIN_FILE "~/Project/Kernel_Machine/DATA/out.50.hap-1.bin" */
#define N_Sample 10000
#define N_Marker 20142

char OUT[N_Sample * N_Marker ];   /* N_Sample * N_Marker Matrix  */
double  Effect[N_Sample];


double runif_new(){

 double val ;
 val = ((double) rand() ) / RAND_MAX ; 
 return val ;
 
}

void  setseed(int * seed){
  
/*Rprintf("Seed Number %d\n",*seed); */
  srand ( *seed );

}


int GetData(int sample, int marker){
  
  int re;
  
  if(sample >= N_Sample || marker >= N_Marker){
    Rprintf("Error: Index [%d][%d]!\n",sample,marker);
    return 0;
  }
  
  re = (int) OUT[ sample*N_Marker + marker]  ;

  return re ;
}

void  Init(int * nMarker, int * nSample, char ** FileName, int * seed){

  int temp;
  FILE *fpin;
  char * file = FileName[0];

  fpin = fopen(file,"r");
  if(fpin == NULL){
    Rprintf("Error: file open [%s]!\n",file);
    exit(1);
  }

  temp = fread ( OUT, sizeof(char) , N_Sample * N_Marker, fpin );
  if(temp !=   N_Sample * N_Marker){
    Rprintf("Error: file read size [%d][%d]!\n",temp,N_Sample * N_Marker);
    exit(1);
  }

  Rprintf("File Read \n");
  *nMarker = N_Marker;
  *nSample = N_Sample ;
  
  setseed(seed);
  }


/*
  Haplotype : Sample_N * Marker_N Matrix, row-wise matrix
*/
void Get_Haplotype(int *Marker_N, int * Marker_Idx, int *Sample_N , int *Sample_Idx1, int * Haplotype)
   {
   int i, j, temp_m, temp_s1, temp_s2, p, n;
   p = *Marker_N;
   n = *Sample_N;
   for(i=0;i< p;i++)
      {
      temp_m = Marker_Idx[i] -1;
      for(j=0;j < n ;j++)
         {
         temp_s1 = Sample_Idx1[j] -1;      
         Haplotype[j*p + i] = GetData(temp_s1, temp_m); /*** modified by Fan on July 17 2015 ***/
         }
      }
   }

/*
  Genotype : N * Marker_N Matrix, row-wide
*/
void  Get_DATA_ALL(int *Marker_N, int * Marker_Idx, int * Genotype){

  int i,j,temp,p;
  p = *Marker_N;
  for(i=0;i< p;i++){
    temp = Marker_Idx[i] -1;
    for(j=0;j < N_Sample ;j++){
      Genotype[j*p + i] = GetData(j, temp); 
    }
  }
}

void Generate_Sample_Idx(int * i1, int *i2){

  int i, j ;
  i = (int)(runif_new()* N_Sample);
  if(i >= N_Sample){
	i = N_Sample -1;
  }

  do {
    j = (int)(runif_new()* N_Sample);
	if(j >= N_Sample){
		j = N_Sample -1;
	}
  }while( i == j);
  
  *i1 = i;
  *i2 = j;
}

/*
# beta[0] : For intercept
*/

void  Set_Beta(int * Marker_N, int * Marker_Idx, double * beta, double * OUT_Effect){

  int i,j,temp, temp1,p;
  double beta_0 ;
  p = *Marker_N;
  
  for(j=0;j < N_Sample ;j++){
    beta_0 = beta[0];
    for(i=0;i< p;i++){
      temp = Marker_Idx[i] -1;
      temp1 = GetData(j, temp);
      beta_0 += ((double) temp1) * beta[i+1];
    }
    Effect[j] = beta_0;
    OUT_Effect[j] = beta_0;
  }

}


void  Set_Beta_wBeta0(int * Marker_N, int * Marker_Idx, double * beta, double * beta0, double * OUT_Effect){

  int i,j,temp, temp1,p;
  double beta_0 ;
  p = *Marker_N;
  
  for(j=0;j < N_Sample ;j++){
    beta_0 = beta0[j];
    for(i=0;i< p;i++){
      temp = Marker_Idx[i] -1;
      temp1 = GetData(j, temp);
      beta_0 += ((double) temp1) * beta[i];
    }
    Effect[j] = beta_0;
    OUT_Effect[j] = beta_0;
  }

}



void  Simu_Case(int *n, int * Index1, int * Index2){

  int i1,i2, temp_n, n1,count;
  double beta1, p1;
  
  count = 0;
  do {
    Generate_Sample_Idx(&i1,&i2);
    beta1 =  Effect[i1] + Effect[i2];
    p1 = 1 - 1/(1 + exp(beta1));
    n1 = rbinom(1,p1);
    
    if(n1 == 1){
      Index1[count] = i1;
      Index2[count] = i2;
      count++;
      /* Rprintf("[%f][%f] - [%d][%d][%d]\n",p1,beta1,count,i1,i2);     */
    }
    
  }while(count < *n);

}


void  Simu_Control(int *n, int * Index1, int * Index2){

  int i1,i2, temp_n, n1,count;
  double beta1, p1;

  count = 0;
  do {
    Generate_Sample_Idx(&i1,&i2);
    beta1 =  Effect[i1] + Effect[i2];
    p1 = 1 - 1/(1 + exp(beta1));
    n1 = rbinom(1,p1);

    if(n1 == 0){
      Index1[count] = i1;
      Index2[count] = i2;
      count++;
    }

  }while(count < *n);

}



void  Simu_Case_wCov(int *n,double * beta_cov, int * Index1, int * Index2, int * Cov){

  int i1,i2, temp_n, n1,count, cov1;
  double beta1, p1 ;
  
  count = 0;
  do {
    Generate_Sample_Idx(&i1,&i2);
    cov1 = rbinom(1,0.5);
    beta1 =  Effect[i1] + Effect[i2] + ((double)cov1) * (*beta_cov);
    p1 = 1 - 1/(1 + exp(beta1));
    n1 = rbinom(1,p1);
    
    if(n1 == 1){
      Index1[count] = i1;
      Index2[count] = i2;
      Cov[count] = cov1;
      count++;
      /* Rprintf("[%f][%f] - [%d][%d][%d]\n",p1,beta1,count,i1,i2);     */
    }
    
  }while(count < *n);

}


void  Simu_Control_wCov(int *n,double * beta_cov, int * Index1, int * Index2, int * Cov){

  int i1,i2, temp_n, n1,count,cov1;
  double beta1, p1;

  count = 0;
  do {
    Generate_Sample_Idx(&i1,&i2);
    cov1 = rbinom(1,0.5);
    beta1 =  Effect[i1] + Effect[i2]  + ((double)cov1) * (*beta_cov);
    p1 = 1 - 1/(1 + exp(beta1));
    n1 = rbinom(1,p1);

    if(n1 == 0){
      Index1[count] = i1;
      Index2[count] = i2;
      Cov[count] = cov1;
      count++;
    }

  }while(count < *n);

}

