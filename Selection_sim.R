# This file contains codes for numerical results in Table 2 for the simulations of the model selections. 
#------------------------------------------------------------------------------------------------------------------------------------
setwd("/disk2/home/yanhan/Sim/SLR_Code")
rm(list = ls())
source('Selection.R')

library(doParallel)
library(foreach)

p.x = 3; p.z = 2
p = p.x + p.z
d.z = p.z+1; d.x = p.x+1
U.gm = c(10,1, 10)
L.gm = c(-10.1,1,-10)
U.bt = rep(20, d.x)
L.bt = rep(-10, d.x)

Grid = Grid_gen(U = U.gm, L = L.gm, p = p.z,s=20)


Rep = 500
for(n in c(200,400,800,1600)){
  for(m in c("1,2,3,4")){
    
    cl <- makePSOCKcluster(70)
    registerDoParallel(cl)
    
    
    Output_A_1 = foreach(i = 1 : Rep,.combine='rbind') %dopar%{
      #setwd("D:/Han/Rwork/SLR/FRLM/")
      #rm(list = ls())
      #source('lib_FRLM.R')
      set.seed(i+500)
      
      U.gm = c(10,1, 10)
      L.gm = c(-10.1,1,-10)
      U.bt = rep(20, d.x)
      L.bt = rep(-10, d.x)
      
      Data = Model_DGP(model = m, n = n)
      
      try = Submodels(Data = Data, beta_true = Data$beta, Regime_true = Data$regime, Grid = Grid, 
                      U.gm = U.gm, L.gm = L.gm, L.bt = L.bt, U.bt = U.bt,tau1 = 0.01 )
      
      output = c(try$SSR4, try$SSR3,try$SSR2, try$SSR1, try$mis4, try$mis3, try$mis2, try$mis1,
                 try$beta_e4, try$beta_e3, try$beta_e2, try$beta_e1)
      
    }  
    write.csv(Output_A_1, sprintf("/disk2/home/yanhan/Sim/SLR_Code/Selection/Model_%s_N_%s.csv",m,n))
    stopImplicitCluster()
    stopCluster(cl)
    
  }
}