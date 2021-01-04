#### dynamic rule

# source("MGAO.R")

dRule = function(geno.only){
      # SNP matrix
      cond <- colSums(geno.only) > 0
      SNP_mx = geno.only[ ,cond]

      # dynamic rule for number of basis function use M_GAO #added by CY, Nov 2019
      M_gao = M_GAO(SNP_mx)

      if (M_gao < 19)             {betabasis_Bsp = genobasis_Bsp = 6; betabasis_Fsp = genobasis_Fsp = 7}
      if (M_gao > 18 & M_gao <25) {betabasis_Bsp = genobasis_Bsp = 8; betabasis_Fsp = genobasis_Fsp = 9}
      if (M_gao > 24 & M_gao <31) {betabasis_Bsp = genobasis_Bsp = 10; betabasis_Fsp = genobasis_Fsp = 11}
      if (M_gao > 30 & M_gao <37) {betabasis_Bsp = genobasis_Bsp = 12; betabasis_Fsp = genobasis_Fsp = 13}
      if (M_gao > 36 & M_gao <43) {betabasis_Bsp = genobasis_Bsp = 14; betabasis_Fsp = genobasis_Fsp = 15}
      if (M_gao > 42)             {betabasis_Bsp = genobasis_Bsp = 16; betabasis_Fsp = genobasis_Fsp = 17}

     out = list(betabasis_Bsp=betabasis_Bsp, genobasis_Bsp=genobasis_Bsp, betabasis_Fsp=betabasis_Fsp, genobasis_Fsp=genobasis_Fsp, M_gao=M_gao)
     return(out)

}
