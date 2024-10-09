# This file contains codes for numerical results in Table 1 for the  simulations of estimation. 
#------------------------------------------------------------------------------------------------------------------------------------
setwd("/disk2/home/yanhan/Sim/SLR_Code")
rm(list = ls())
source('MIQP.R')

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


Settings = read.csv("Temporal_settings.csv")
Rep = 500

for(s in c(22)){
  N = Settings$N[s]; ar = Settings$AR[s]; ma = Settings$MA[s]; level = Settings$Level[s]
  
  Start_time = Sys.time()
  
  n = N; sg = 1
  
  noise = rnorm(n, mean = 0, sd = sg)
  p.x = 3; p.z = 2
  p = p.x + p.z
  d.z = p.z+1; d.x = p.x+1; n.obs = N
  
  
  cl <- makePSOCKcluster(64)
  registerDoParallel(cl)
  
  
  Output_A_1 = foreach(i = 1 : Rep,.combine='rbind') %dopar%{
    #setwd("D:/Han/Rwork/SLR/FRLM/")
    #rm(list = ls())
    #source('lib_FRLM.R')
    set.seed(i+1111)
    
    data = matrix(0.0, nrow = n, ncol = p)
    # for(i in 1:p){data[,i] = runif(n = N, min = -5, max = 5)}
    for(c in 1:p.x){data[,c] = DGP(ar = ar, ma = ma, level = level, n = N)}
    for(c in (p.x+1):p){data[,c] = DGP(ar = ar, ma = ma, level = level, n = N); 
    data[,c] = data[,c]/max(abs(data[,c]))}
    
    data = as.data.frame(data)
    names(data) = c(sprintf("x%s", 1:p.x),sprintf("z%s", 1:p.z))
    
    gamma.0 = cbind(c(0,1,-1),c(0,1,1))
    beta.0 = cbind(c(1,1,1,1), c(-3,2,-1,0),
                   c(0,1,3,-1), c(2,-1,0,2))
    
    x = as.matrix(data[,1:p.x])
    x = cbind(rep(1, n), x)
    
    z = as.matrix(data[,(1+p.x):p])
    z = cbind(rep(1, n), z)
    regime = regime_class(Z=z, gamma = gamma.0)
    
    y = rowSums((x%*% beta.0) * regime) + noise
    
    U.gm = c(10,1, 10)
    L.gm = c(-10.1,1,-10)
    U.bt = rep(20, d.x)
    L.bt = rep(-10, d.x)
    
    
    iter_result = MIQP_iter(x = x, y = y, z=  z, Grid = Grid, U.gm = U.gm, L.gm = L.gm, L.bt = L.bt, U.bt = U.bt, max.iter = 5)
    
    iter_gamma = iter_result$gamma; iter_beta = iter_result$beta
    iter_dist = dist_v_H(v = gamma.0[,1], H = iter_gamma)
    if(iter_dist$ii == 2){
      iter_gamma[,1] = iter_result$gamma[,2]; iter_gamma[,2] = iter_result$gamma[,1]
      iter_beta[,2] = iter_result$beta[,4]; iter_beta[,4] = iter_result$beta[,2] 
    }
    
    
    iter_gamma_bias = sqrt(sum((iter_gamma - gamma.0)^2))
    iter_beta_bias = sqrt(sum((iter_beta - beta.0)^2))
    
    iter_regime = regime_class(Z=z, gamma = iter_gamma)
    iter_regime_bias = sum((iter_regime - regime)^2)/(2*N)
    
    output = c(N, ar, ma, level, iter_gamma_bias, iter_beta_bias, iter_regime_bias, 
               c(iter_gamma), c(iter_beta), iter_result$time)
  }
  
  
  End_time = Sys.time()
  write.csv(Output_A_1, sprintf("/disk2/home/yanhan/Sim/SLR_Code/Temporal/ar_%s_ma_%s_level_%sN_%s.csv",ar,ma,level, N))
  stopImplicitCluster()
  stopCluster(cl)
  
  time = (End_time - Start_time)
  cat("setting=",s,"is done",'\n')
  print(time)
  print(colMeans(Output_A_1[,1:7]))
}