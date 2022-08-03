# File Name: sim1.R
# AR settings
#
# Latest Update: 
#   
#   2022-08-03  Original Code
# 

#--------------------------------------------------

setwd("D:/Han/Rwork/SLR/FRLM/")
rm(list = ls())
source('lib_FRLM2.R')

library(doParallel)
library(foreach)

dist_v_H = function(v, H){
  bias = H-v
  mse = sapply(1:ncol(H), function(x){sum(bias[,x]^2)})
  ii = which.min(mse)
  bias = bias[,ii]
  return(list(ii = ii, bias = bias))
}

regime_class = function(Z, gamma){
  z.prod = Z %*% gamma
  z.1 = (z.prod[,1]>=0)*(z.prod[,2]>=0)
  z.2 = (z.prod[,1]<0)*(z.prod[,2]>=0)
  z.3 = (z.prod[,1]<0)*(z.prod[,2]<0)
  z.4 = (z.prod[,1]>=0)*(z.prod[,2]<0)
  z.class = cbind(z.1, z.2, z.3, z.4)
  return(z.class)
}

#---------------------------------------------------------------------------------------------------------------------------------
#Setting 1. Independent data

#---------------------------------------------------------------------------------------------------------------------------------
n = 1600;sg = 1
noise = rnorm(n, mean = 0, sd = sg)
p.x = 3; p.z = 2
p = p.x + p.z
d.z = p.z+1; d.x = p.x+1; n.obs = n
data = matrix(0.0, nrow = n, ncol = p)
# for(i in 1:p){data[,i] = runif(n = N, min = -5, max = 5)}
for(i in 1:p.x){data[,i] = arima.sim(model=list(ar=0.5), n=n)}
for(i in (p.x+1):p){data[,i] = arima.sim(model=list(ar=0.5), n=n)}
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

U.gm = c(5,1, 5)
L.gm = c(-5,1,-5)
U.bt = rep(5, d.x)
L.bt = rep(-5, d.x)

Grid = Grid_gen(U = U.gm, L = L.gm, p = p.z,s=30)


step1.out = step1_grid(y=y,x=x,z = z,grid=Grid)
rm(Grid); gc()
bt.hat.step1 = step1.out$bt.hat

for(N in c(400)){
  Start_time = Sys.time()
  Rep = 500
  cl <- makePSOCKcluster(30)
  registerDoParallel(cl)
  Output_A_1 = foreach(i = 1 : Rep,.combine='rbind') %dopar%{
    #setwd("D:/Han/Rwork/SLR/FRLM/")
    #rm(list = ls())
    #source('lib_FRLM.R')
    set.seed(i)
    
    n = N; sg = 1
    
    noise = rnorm(n, mean = 0, sd = sg)
    p.x = 3; p.z = 2
    p = p.x + p.z
    d.z = p.z+1; d.x = p.x+1; n.obs = N
    data = matrix(0.0, nrow = n, ncol = p)
    for(i in 1:p.x){data[,i] = arima.sim(model=list(ar=0.5), n=N)}
    for(i in (p.x+1):p){data[,i] = arima.sim(model=list(ar=0.5), n=N)}
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
    
    U.gm = c(5,1, 5)
    L.gm = c(-5,1,-5)
    U.bt = rep(5, d.x)
    L.bt = rep(-5, d.x)
    

    try = MIQP_iter(x = x, y = y, z=  z, bt.hat.step1 = bt.hat.step1, U.gm = U.gm, L.gm = L.gm, L.bt = L.bt, U.bt = U.bt, max.iter = 50)
    
    gamma.hat = try$gamma; beta.hat = try$beta
    dist_1 = dist_v_H(v = gamma.0[,1], H = try$gamma)
    if(dist_1$ii == 2){
      gamma.hat[,1] = try$gamma[,2]; gamma.hat[,2] = try$gamma[,1]
      beta.hat[,2] = try$beta[,4]; beta.hat[,4] = try$beta[,2] 
    }
    
    regime.hat = regime_class(Z=z, gamma = gamma.hat)
    
    gamma.bias = sum((gamma.hat - gamma.0)^2)
    beta.bias = sum((beta.hat - beta.0)^2)
    regime.bias = sum((regime.hat - regime)^2)/(2*N)
    
    output = c(gamma.bias, beta.bias, regime.bias, c(gamma.hat), c(beta.hat))
    
  }
  
  End_time = Sys.time()
  write.csv(Output_A_1, sprintf("D:/Han/Rwork/SLR/FRLM/sim_2_outputs/Setting_1_N_%s.csv", N))
  stopImplicitCluster()
  stopCluster(cl); gc()
  
  time = (End_time - Start_time)
  cat("N=",N,"is done",'\n')
  print(time)
  print(colMeans(Output_A_1[,1:3]))
  
}

set_1_N_200 = read.csv("sim_2_outputs/Setting_1_N_200.csv")

means = colMeans(set_1_N_200[,-1])

#---------------------------------------------------------------------------------------------------------------------------------
#Setting 2. 3-regime mode-1

#---------------------------------------------------------------------------------------------------------------------------------



for(N in c(200, 400, 800,1600)){
  Start_time = Sys.time()
  
  Rep = 50
  cl <- makePSOCKcluster(30)
  registerDoParallel(cl)
  Output_A_2 = foreach(i = 1 : Rep,.combine='rbind') %dopar%{
    #setwd("D:/Han/Rwork/SLR/FRLM/")
    #rm(list = ls())
    #source('lib_FRLM.R')
    
    n = N; sg = 1
    set.seed(i)
    
    noise = rnorm(n, mean = 0, sd = sg)
    p.x = 3; p.z = 2
    p = p.x + p.z
    d.z = p.z+1; d.x = p.x+1; n.obs = N
    data = matrix(0.0, nrow = n, ncol = p)
    for(i in 1:p.x){data[,i] = arima.sim(model=list(ar=0.5), n=N)}
    for(i in (p.x+1):p){data[,i] = arima.sim(model=list(ar=0.5), n=N)}
    data = as.data.frame(data)
    names(data) = c(sprintf("x%s", 1:p.x),sprintf("z%s", 1:p.z))
    
    gamma.0 = cbind(c(-1,1,0),c(1,1,0))
    beta.0 = cbind(c(1,1,1,1), c(-3,2,-1,0),
                   c(0,1,3,-1), c(2,-1,0,2))
    
    
    x = as.matrix(data[,1:p.x])
    x = cbind(rep(1, n), x)
    
    z = as.matrix(data[,(1+p.x):p])
    z = cbind(rep(1, n), z)
    regime = regime_class(Z=z, gamma = gamma.0)
    
    y = rowSums((x%*% beta.0) * regime) + noise
    
    U.gm = c(5,1, 5)
    L.gm = c(-5,1,-5)
    U.bt = rep(5, d.x)
    L.bt = rep(-5, d.x)
    
    Grid = Grid_gen(U = U.gm, L = L.gm, p = p.z)
    
    try = MIQP_iter(x = x, y = y, z=  z, Grid = Grid, U.gm = U.gm, L.gm = L.gm, L.bt = L.bt, U.bt = U.bt, max.iter = 50)
    
    gamma.hat = try$gamma; beta.hat = try$beta
    dist_1 = dist_v_H(v = gamma.0[,1], H = try$gamma)
    if(dist_1$ii == 2){
      gamma.hat[,1] = try$gamma[,2]; gamma.hat[,2] = try$gamma[,1]
      beta.hat[,2] = try$beta[,4]; beta.hat[,4] = try$beta[,2] 
    }
    
    
    gamma.bias = sum((gamma.hat - gamma.0)^2)
    beta.bias = 0
    for(bi in 1:3){
      dist_b = dist_v_H(v = beta.0[,bi], H = beta.hat)
      beta.bias = beta.bias+sum(dist_b$bias^2)
    }
    
    output = c(gamma.bias, beta.bias,  c(gamma.hat), c(beta.hat))
    
  }
  
  End_time = Sys.time()
  write.csv(Output_A_2, sprintf("D:/Han/Rwork/SLR/FRLM/sim_2_outputs/Setting_2_N_%s.csv", N))
  stopImplicitCluster()
  stopCluster(cl)
  
  time = (End_time - Start_time)
  cat("N=",N,"is done",'\n')
  print(time)
  print(colMeans(Output_A_2[,1:3]))
}


#---------------------------------------------------------------------------------------------------------------------------------
#Setting 3. 3-regime mode-2

#---------------------------------------------------------------------------------------------------------------------------------



for(N in c(200, 400, 800,1600)){
  Start_time = Sys.time()
  
  Rep = 50
  cl <- makePSOCKcluster(30)
  registerDoParallel(cl)
  Output_A_3 = foreach(i = 1 : Rep,.combine='rbind') %dopar%{
    #setwd("D:/Han/Rwork/SLR/FRLM/")
    #rm(list = ls())
    #source('lib_FRLM.R')
    
    n = N; sg = 1
    set.seed(i)
    
    noise = rnorm(n, mean = 0, sd = sg)
    p.x = 3; p.z = 2
    p = p.x + p.z
    d.z = p.z+1; d.x = p.x+1; n.obs = N
    data = matrix(0.0, nrow = n, ncol = p)
    for(i in 1:p.x){data[,i] = arima.sim(model=list(ar=0.5), n=N)}
    for(i in (p.x+1):p){data[,i] = arima.sim(model=list(ar=0.5), n=N)}
    data = as.data.frame(data)
    names(data) = c(sprintf("x%s", 1:p.x),sprintf("z%s", 1:p.z))
    
    gamma.0 = cbind(c(0,1,-1),c(0,1,1))
    beta.0 = cbind(c(1,1,1,1), c(-3,2,-1,0),
                   c(0,1,3,-1))
    
    x = as.matrix(data[,1:p.x])
    x = cbind(rep(1, n), x)
    
    z = as.matrix(data[,(1+p.x):p])
    z = cbind(rep(1, n), z)
    regime.o = regime_class(Z=z, gamma = gamma.0)
    regime = matrix(0.0, nrow = nrow(regime.o), ncol = 3)
    regime[,1:2] = regime.o[,1:2]; regime[,3] = regime.o[,3]+ regime.o[,4]
    
    y = rowSums((x%*% beta.0) * regime) + noise
    
    U.gm = c(5,1, 5)
    L.gm = c(-5,1,-5)
    U.bt = rep(5, d.x)
    L.bt = rep(-5, d.x)
    
    Grid = Grid_gen(U = U.gm, L = L.gm, p = p.z)
    
    try = MIQP_iter(x = x, y = y, z=  z, Grid = Grid, U.gm = U.gm, L.gm = L.gm, L.bt = L.bt, U.bt = U.bt, max.iter = 50)
    
    gamma.hat = try$gamma; beta.hat = try$beta
    dist_1 = dist_v_H(v = gamma.0[,1], H = try$gamma)
    if(dist_1$ii == 2){
      gamma.hat[,1] = try$gamma[,2]; gamma.hat[,2] = try$gamma[,1]
      beta.hat[,2] = try$beta[,4]; beta.hat[,4] = try$beta[,2] 
    }
    
    
    gamma.bias = sum((gamma.hat - gamma.0)^2)
    beta.bias = 0
    for(bi in 1:3){
      dist_b = dist_v_H(v = beta.0[,bi], H = beta.hat)
      beta.bias = beta.bias+sum(dist_b$bias^2)
    }
    
    output = c(gamma.bias, beta.bias,  c(gamma.hat), c(beta.hat))
    
  }
  
  End_time = Sys.time()
  write.csv(Output_A_3, sprintf("D:/Han/Rwork/SLR/FRLM/sim_2_outputs/Setting_3_N_%s.csv", N))
  stopImplicitCluster()
  stopCluster(cl)
  
  time = (End_time - Start_time)
  cat("N=",N,"is done",'\n')
  print(time)
  print(colMeans(Output_A_3[,1:3]))
  
}

## Oracle 
for(N in c(200, 400, 800,1600)){
  Rep = 500
  cl <- makePSOCKcluster(3)
  clusterEvalQ(cl, .libPaths("D:/R/library" ))
  source('E:/Rwork/SLR/FRLM/ib_FRLM2.R')
  
  registerDoParallel(cl)
  Output_A_3 = foreach(i = 1 : Rep,.combine='rbind') %dopar%{
  n = N; sg = 1
  set.seed(i)
  
  noise = rnorm(n, mean = 0, sd = sg)
  p.x = 3; p.z = 2
  p = p.x + p.z
  d.z = p.z+1; d.x = p.x+1; n.obs = N
  data = matrix(0.0, nrow = n, ncol = p)
  for(i in 1:p.x){data[,i] = arima.sim(model=list(ar=0.5), n=N)}
  for(i in (p.x+1):p){data[,i] = arima.sim(model=list(ar=0.5), n=N)}
  data = as.data.frame(data)
  names(data) = c(sprintf("x%s", 1:p.x),sprintf("z%s", 1:p.z))
  
  gamma.0 = cbind(c(0,1,-1),c(0,1,1))
  beta.0 = cbind(c(1,1,1,1), c(-3,2,-1,0),
                 c(0,1,3,-1))
  
  x = as.matrix(data[,1:p.x])
  x = cbind(rep(1, n), x)
  
  z = as.matrix(data[,(1+p.x):p])
  z = cbind(rep(1, n), z)
  regime.o = regime_class(Z=z, gamma = gamma.0)
  regime = matrix(0.0, nrow = nrow(regime.o), ncol = 3)
  regime[,1:2] = regime.o[,1:2]; regime[,3] = regime.o[,3]+ regime.o[,4]
  
  y = rowSums((x%*% beta.0) * regime) + noise
  
  z1 = regime[,1]; z2= regime[,2]; z3 = regime[,3]
  x1 = x; x1[z2==1,]=0; x1[z3==1,]=0
  x2 = x; x2[z1==1,]=0; x2[z3==1,]=0
  x3 = x; x3[z2==1,]=0; x3[z1==1,]=0
  
  xx = cbind(x1,x2,x3)
  lm0 = lm(y~xx+0)
  beta.hat.0 = coef(lm0)
  beta.hat.0 = matrix(beta.hat.0,ncol = 3, byrow = F)
  
  beta.bias = 0
  for(bi in 1:3){
    dist_b = dist_v_H(v = beta.0[,bi], H = beta.hat.0)
    beta.bias = beta.bias+sum(dist_b$bias^2)
  }
  output = c(beta.bias, c(beta.hat.0))
  }
  
  End_time = Sys.time()
  write.csv(Output_A_3, sprintf("E:/Rwork/SLR/FRLM/sim_2_outputs/Setting_3_N_%s_0.csv", N))
  stopImplicitCluster()
  stopCluster(cl)

}

#---------------------------------------------------------------------------------------------------------------------------------
#Setting 4. 2-regime mode-2

#---------------------------------------------------------------------------------------------------------------------------------

n = 1600; sg = 1
set.seed(i)

noise = rnorm(n, mean = 0, sd = sg)
p.x = 3; p.z = 2
p = p.x + p.z
d.z = p.z+1; d.x = p.x+1; n.obs = N
data = matrix(0.0, nrow = n, ncol = p)
for(i in 1:p){data[,i] = runif(n, min = 0, max = 20)}
for(i in 1:p.x){data[,i] = arima.sim(model=list(ar=0.5), n=N)}  
data = as.data.frame(data)
names(data) = c(sprintf("x%s", 1:p.x),sprintf("z%s", 1:p.z))

gamma.0 = cbind(c(0,1,-1),c(0,1,1))
beta.0 = cbind(c(1,1,1,1), c(-3,2,-1,0),
               c(0,1,3,-1))

x = as.matrix(data[,1:p.x])
x = cbind(rep(1, n), x)

z = as.matrix(data[,(1+p.x):p])
z = cbind(rep(1, n), z)
regime.o = regime_class(Z=z, gamma = gamma.0)
regime = matrix(0.0, nrow = nrow(regime.o), ncol = 3)
regime[,1:2] = regime.o[,1:2]; regime[,3] = regime.o[,3]+ regime.o[,4]

y = rowSums((x%*% beta.0) * regime) + noise

U.gm = c(5,1, 5)
L.gm = c(-5,1,-5)
U.bt = rep(5, d.x)
L.bt = rep(-5, d.x)

Grid = Grid_gen(U = U.gm, L = L.gm, p = p.z,s=30)


step1.out = step1_grid(y=y,x=x,z = z,grid=Grid)
rm(Grid); gc()
bt.hat.step1 = step1.out$bt.hat



for(N in c(1600)){
  Start_time = Sys.time()
  
  Rep = 40
  cl <- makePSOCKcluster(30)
  registerDoParallel(cl)
  Output_A_4 = foreach(i = 1 : Rep,.combine='rbind') %dopar%{
    
    n = N; sg = 1
    set.seed(i)
    
    noise = rnorm(n, mean = 0, sd = sg)
    p.x = 3; p.z = 2
    p = p.x + p.z
    d.z = p.z+1; d.x = p.x+1; n.obs = N
    data = matrix(0.0, nrow = n, ncol = p)
    for(i in 1:p){data[,i] = runif(n, min = 0, max = 20)}
    for(i in 1:p.x){data[,i] = arima.sim(model=list(ar=0.5), n=N)}
    for(i in (p.x+1):p){data[,i] = arima.sim(model=list(ar=0.5), n=N)}
    data = as.data.frame(data)
    names(data) = c(sprintf("x%s", 1:p.x),sprintf("z%s", 1:p.z))
    
    gamma.0 = c(0,1,-1)
    beta.0 = cbind(c(1,1,1,1), c(-3,2,-1,0))
    
    x = as.matrix(data[,1:p.x])
    x = cbind(rep(1, n), x)
    
    z = as.matrix(data[,(1+p.x):p])
    z = cbind(rep(1, n), z)
    # regime.o = regime_class(Z=z, gamma = gamma.0)
    
    z.prod = z %*% gamma.0
    z.1 = (z.prod>=0)
    z.2 = (z.prod<0)
    regime = cbind(z.1, z.2)
    
    
    y = rowSums((x%*% beta.0) * regime) + noise
    
    U.gm = c(5,1, 5)
    L.gm = c(-5,1,-5)
    U.bt = rep(5, d.x)
    L.bt = rep(-5, d.x)
    
    try = MIQP_iter(x = x, y = y, z=  z, bt.hat.step1 = bt.hat.step1, U.gm = U.gm, L.gm = L.gm,
                    L.bt = L.bt, U.bt = U.bt, tau1= 0, tau2 = 0.9, max.iter = 50)
    
    gamma.hat = try$gamma; beta.hat = try$beta
    dist_gm = dist_v_H(v = gamma.0, H = try$gamma)
    if(dist_gm$ii == 2){
      gamma.hat[,1] = try$gamma[,2]; gamma.hat[,2] = try$gamma[,1]
      beta.hat[,2] = try$beta[,4]; beta.hat[,4] = try$beta[,2] 
    }
    
    gamma.bias = sum((gamma.hat[,1] - gamma.0)^2)
    beta.bias = 0
    for(bi in 1:2){
      dist_b = dist_v_H(v = beta.0[,bi], H = beta.hat)
      beta.bias = beta.bias+sum(dist_b$bias^2)
    }
    
    output = c(gamma.bias, beta.bias,  c(gamma.hat), c(beta.hat))
    
  }
  
  End_time = Sys.time()
  write.csv(Output_A_4, sprintf("D:/Han/Rwork/SLR/FRLM/sim_2_outputs/Setting_4_N_%s.csv", N))
  stopImplicitCluster()
  stopCluster(cl)
  
  time = (End_time - Start_time)
  cat("N=",N,"is done",'\n')
  print(time)
  print(colMeans(Output_A_4[,1:2]))
  
}




#---------------------------------------------------------------------------------------------------------------------------------
#Setting 5. 1-regime model

#---------------------------------------------------------------------------------------------------------------------------------

for(N in c(200,400,800)){
  Start_time = Sys.time()
  
  Rep = 40
  cl <- makePSOCKcluster(30)
  registerDoParallel(cl)
  Output_A_5 = foreach(i = 1 : Rep,.combine='rbind') %dopar%{
    
    n = N; sg = 1
    
    noise = rnorm(n, mean = 0, sd = sg)
    p.x = 3; p.z = 2
    p = p.x + p.z
    d.z = p.z+1; d.x = p.x+1; n.obs = N
    data = matrix(0.0, nrow = n, ncol = p)
    for(i in 1:p){data[,i] = runif(n, min = 0, max = 20)}
    data = as.data.frame(data)
    names(data) = c(sprintf("x%s", 1:p.x),sprintf("z%s", 1:p.z))
    
    gamma.0 = c(0,1,-1)
    beta.0 = cbind(c(1,1,1,1))
    
    x = as.matrix(data[,1:p.x])
    x = cbind(rep(1, n), x)
    
    z = as.matrix(data[,(1+p.x):p])
    z = cbind(rep(1, n), z)
    # regime.o = regime_class(Z=z, gamma = gamma.0)
    
    
    y = (x%*% beta.0) + noise
    
    U.gm = c(5,1, 5)
    L.gm = c(-5,1,-5)
    U.bt = rep(5, d.x)
    L.bt = rep(-5, d.x)
    
    Grid = Grid_gen(U = U.gm, L = L.gm, p = p.z)
    
    try = MIQP_iter(x = x, y = y, z=  z, Grid = Grid, U.gm = U.gm, L.gm = L.gm, L.bt = L.bt, U.bt = U.bt, max.iter = 50)
    
    gamma.hat = try$gamma; beta.hat = try$beta
    
    dist_b = dist_v_H(v = beta.0, H = beta.hat)
    
    
    output = c(gamma.bias, beta.bias,  c(gamma.hat), c(beta.hat))
    
  }
  
  End_time = Sys.time()
  write.csv(Output_A_5, sprintf("D:/Han/Rwork/SLR/FRLM/sim_2_outputs/Setting_5_N_%s.csv", N))
  stopImplicitCluster()
  stopCluster(cl)
  
  time = (End_time - Start_time)
  cat("N=",N,"is done",'\n')
  print(time)
  print(colMeans(Output_A_5[,2]))
}

