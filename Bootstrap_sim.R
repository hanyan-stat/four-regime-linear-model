# This file contains codes for numerical results in Table 3 for simulations of the Bootstrap. 
#------------------------------------------------------------------------------------------------------------------------------------
setwd("/disk2/home/yanhan/Sim/SLR_Code")
rm(list = ls())
source('Bootstrap.R')

library(doParallel)
library(foreach)

p.x = 3; p.z = 2
p = p.x + p.z
d.z = p.z+1; d.x = p.x+1
U.gm = c(10,1, 10)
L.gm = c(-10.1,1,-10)
U.bt = rep(20, d.x)
L.bt = rep(-10, d.x)

Grid = Grid_gen(U = U.gm, L = L.gm, p = p.z,s=25)


Rep = 500
B = 500
for(n in c(200,400,800,1600)){
    
    cl <- makePSOCKcluster(70)
    registerDoParallel(cl)
    
    Output_A_1 = foreach(i = 1 : Rep,.combine='rbind') %dopar%{
      #setwd("D:/Han/Rwork/SLR/FRLM/")
      #rm(list = ls())
      #source('lib_FRLM.R')
      set.seed(i)
      U.gm = c(10,1, 10)
      L.gm = c(-10.1,1,-10)
  
      
      
      Data = Model_DGP(model = "4", n = n)
      bt.1 = step1_grid(y = Data$y, x = Data$x, z = Data$z, grid = Grid)
      bt.1 = bt.1$bt.hat
      
      const = iter_constraint_wm(L.gm = L.gm, U.gm = U.gm, d.z = d.z, d.x = d.x, n.obs = n, tau1 = 0)
      
      est = MIQP_iter_fast(x = Data$x, y = Data$y, z=  Data$z, gs = T, bt = bt.1, Grid = Grid,
                           const = const, U.gm = U.gm, L.gm = L.gm, max.iter = 10)
      
      boot_result = replicate(n = B, boot_est(Data = Data, est = est, method = "residual", bt = bt.1, const = const))
      output = c(boot_result)
      
    }  
    write.csv(Output_A_1, sprintf("/disk2/home/yanhan/Sim/SLR_Code/Bootstrap/N_%s.csv",n))
    stopImplicitCluster()
    stopCluster(cl)
  
}