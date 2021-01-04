#### Chi-Yang Chiu 05/16/2016

#alpha="05","01","001","0001","00001","000001"
# Neg_beta_pct1=0
# Neg pct2 =  0, Neg pct3 = 0, 20, 50
# Neg pct2 = 20, Neg pct3 = 0, 20, 50
# Neg pct2 = 50, Neg pct3 = 0, 20, 50

#path: COSI_50_Ped
#vtype: rare_and_common/rare
#Neg_beta_pct1=0, 20, 50

#setwd("C:/NICHD/Research/paper/PedGFLMM/simulations")

# Set the working dirctory to the location of this file
# when it is sourced:
this.dir <- dirname(parent.frame(2)$ofile)
this.dir
setwd(this.dir)

bar_plot=function(vtype, path, alpha, region_size)
   {
   #vtype="rare_and_common";alpha=0.01; path="COSI_50_Ped"; region_size=6000

   f1 <- sprintf("%s/Power_causal_%s/Power_%s_region_size_%d_Neg_beta_pct_0_Causal_pct_5_Nsim_3000.csv", path, vtype, vtype, region_size)
   f2 <- sprintf("%s/Power_causal_%s/Power_%s_region_size_%d_Neg_beta_pct_0_Causal_pct_10_Nsim_3000.csv", path, vtype, vtype, region_size)
   f3 <- sprintf("%s/Power_causal_%s/Power_%s_region_size_%d_Neg_beta_pct_0_Causal_pct_15_Nsim_3000.csv", path, vtype, vtype, region_size)
   f4 <- sprintf("%s/Power_causal_%s/Power_%s_region_size_%d_Neg_beta_pct_20_Causal_pct_5_Nsim_3000.csv", path, vtype, vtype, region_size)
   f5 <- sprintf("%s/Power_causal_%s/Power_%s_region_size_%d_Neg_beta_pct_20_Causal_pct_10_Nsim_3000.csv", path, vtype, vtype, region_size)
   f6 <- sprintf("%s/Power_causal_%s/Power_%s_region_size_%d_Neg_beta_pct_20_Causal_pct_15_Nsim_3000.csv", path, vtype, vtype, region_size)
   f7 <- sprintf("%s/Power_causal_%s/Power_%s_region_size_%d_Neg_beta_pct_50_Causal_pct_5_Nsim_3000.csv", path, vtype, vtype, region_size)
   f8 <- sprintf("%s/Power_causal_%s/Power_%s_region_size_%d_Neg_beta_pct_50_Causal_pct_10_Nsim_3000.csv", path, vtype, vtype, region_size)
   f9 <- sprintf("%s/Power_causal_%s/Power_%s_region_size_%d_Neg_beta_pct_50_Causal_pct_15_Nsim_3000.csv", path, vtype, vtype, region_size)

   d1 = read.csv(f1); d1$X=c(5e-2,1e-2,1e-3,1e-4,1e-5,1e-6); colnames(d1)[1]="level"
   d2 = read.csv(f2); d2$X=c(5e-2,1e-2,1e-3,1e-4,1e-5,1e-6); colnames(d2)[1]="level"
   d3 = read.csv(f3); d3$X=c(5e-2,1e-2,1e-3,1e-4,1e-5,1e-6); colnames(d3)[1]="level"
   d4 = read.csv(f4); d4$X=c(5e-2,1e-2,1e-3,1e-4,1e-5,1e-6); colnames(d4)[1]="level"
   d5 = read.csv(f5); d5$X=c(5e-2,1e-2,1e-3,1e-4,1e-5,1e-6); colnames(d5)[1]="level"
   d6 = read.csv(f6); d6$X=c(5e-2,1e-2,1e-3,1e-4,1e-5,1e-6); colnames(d6)[1]="level"
   d7 = read.csv(f7); d7$X=c(5e-2,1e-2,1e-3,1e-4,1e-5,1e-6); colnames(d7)[1]="level"
   d8 = read.csv(f8); d8$X=c(5e-2,1e-2,1e-3,1e-4,1e-5,1e-6); colnames(d8)[1]="level"
   d9 = read.csv(f9); d9$X=c(5e-2,1e-2,1e-3,1e-4,1e-5,1e-6); colnames(d9)[1]="level"

   column= c(#"fixed_LRT_bpval", 
             "fixed_LRT_fpval", "beta_LRT_bpval", "beta_LRT_fpval",
             "kernel_BT_pval",  "kernel_MB_pval", "kernel_UW_pval",
             "burden_BT_pval",  "burden_MB_pval",  "burden_UW_pval")

   a1=d1[d1$level==alpha,names(d1)%in%column]
   a2=d2[d2$level==alpha,names(d2)%in%column]
   a3=d3[d3$level==alpha,names(d3)%in%column]
   b1=d4[d4$level==alpha,names(d4)%in%column]
   b2=d5[d5$level==alpha,names(d5)%in%column]
   b3=d6[d6$level==alpha,names(d6)%in%column]
   c1=d7[d7$level==alpha,names(d7)%in%column]  #Fan
   c2=d8[d8$level==alpha,names(d8)%in%column]  #Fan
   c3=d9[d9$level==alpha,names(d9)%in%column]  #Fan

   filename <- sprintf("%s/Fig_Pow/bar_%s_region_size_%d_Nsim_3000.ps", path, vtype, region_size)
   postscript(file = filename, height=8, width=10, horizontal=F)

   par(mfrow=c(3,3))
   ncolor = length(column)
   color=c("cornflowerblue","cadetblue1",terrain.colors(ncolor-2))

   barplot(t(a1), main="(a1) Neg_beta_pct=0,  Causal_pct=5",  xlab="", ylab="Empirical Power", cex.axis=0.8, cex.names=0.68, ylim=c(0,1), beside=T, col=color, names.arg="", cex.main=0.8)
   barplot(t(a2), main="(a2) Neg_beta_pct=0,  Causal_pct=10", xlab="", ylab="Empirical Power", cex.axis=0.8, cex.names=0.68, ylim=c(0,1), beside=T, col=color, names.arg="", cex.main=0.8)
   barplot(t(a3), main="(a3) Neg_beta_pct=0,  Causal_pct=15", xlab="", ylab="Empirical Power", cex.axis=0.8, cex.names=0.68, ylim=c(0,1), beside=T, col=color, names.arg="", cex.main=0.8)

   barplot(t(b1), main="(b1) Neg_beta_pct=20, Causal_pct=5",  xlab="", ylab="Empirical Power", cex.axis=0.8, cex.names=0.68, ylim=c(0,1), beside=T, col=color, names.arg="", cex.main=0.8)
   barplot(t(b2), main="(b2) Neg_beta_pct=20, Causal_pct=10", xlab="", ylab="Empirical Power", cex.axis=0.8, cex.names=0.68, ylim=c(0,1), beside=T, col=color, names.arg="", cex.main=0.8)
   barplot(t(b3), main="(b3) Neg_beta_pct=20, Causal_pct=15", xlab="", ylab="Empirical Power", cex.axis=0.8, cex.names=0.68, ylim=c(0,1), beside=T, col=color, names.arg="", cex.main=0.8)

   barplot(t(c1), main="(c1) Neg_beta_pct=50, Causal_pct=5",  xlab="", ylab="Empirical Power", cex.axis=0.8, cex.names=0.68, ylim=c(0,1), beside=T, col=color, names.arg="", cex.main=0.8)
   leg=c(#"LRT_GFLMM_BS", 
   "LRT_GFLMM_FR", "LRT_beta_BS", "LRT_beta_FR", "Kernel_BT", "Kernel_MB", "Kernel_UW", "Burden_BT", "Burden_MB", "Burden_UW")
   legend("topright", legend=leg, fill=color, bty="n", cex=0.568)

   barplot(t(c2), main="(c2) Neg_beta_pct=50, Causal_pct=10", xlab="", ylab="Empirical Power", cex.axis=0.8, cex.names=0.68, ylim=c(0,1), beside=T, col=color, names.arg="", cex.main=0.8)
   barplot(t(c3), main="(c3) Neg_beta_pct=50, Causal_pct=15", xlab="", ylab="Empirical Power", cex.axis=0.8, cex.names=0.68, ylim=c(0,1), beside=T, col=color, names.arg="", cex.main=0.8)

   graphics.off()
   }

