# This file contains codes for conducting the Bootstrap procedure

#------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
setwd("/disk2/home/yanhan/SLR")
rm(list = ls())
source('MIQP.R')


library(ks)



#-------------------------------------------
#     Fitting with the estimated FRM        #
#-------------------------------------------

fit_frlm = function(Data, est){
  x = Data$x; z = Data$z; y = Data$y
  gm = est$gamma; bt = est$beta
  
  Regime4 =  regime_class(Z =  Data$z, gamma =  gm)
  fit = rowSums((Data$x%*% est$beta) * Regime4)
  res = y - fit
  
  return(list(fit = fit, res = res))
  
}



#-------------------------------------------
# Smoothed and wild Bootstrap sampling     #
#-------------------------------------------


boot_sample = function(Data, est, method){
  library(simukde)
  x = Data$x; z = Data$z; y = Data$y
  m = nrow(x)
  tz = z[,-1]
  fit = fit_frlm(Data = Data, est = est)
  res = fit$res; res_center = res - mean(res)  
  
  if(method == "smooth"){
    fhat = kde(x = tz)
    tz_boot = rkde(fhat = fhat, n = m)
    z_boot = cbind(1, tz_boot)   
    res_boot = sample(res_center, size = m, replace = T)
    Regime_b = regime_class(Z = z_boot, gamma = est$gamma)
    fit_boot = rowSums((Data$x%*% est$beta) * Regime_b)
    
    y_boot = fit_boot + res_boot
  }else if(method == "residual"){
    z_boot = z
    wild_weight = sample(c((1-sqrt(5))/2,(1+sqrt(5))/2), size = m, 
                         prob = c((1+sqrt(5))/(2*sqrt(5)),(sqrt(5)-1)/(2*sqrt(5))), replace = T)
    res_boot = res_center * wild_weight
    Regime_b = regime_class(Z = z_boot, gamma = est$gamma)
    fit_boot = rowSums((Data$x%*% est$beta) * Regime_b)
    
    y_boot = fit_boot + res_boot
  }
  
  
  return(list(x = Data$x, z = z_boot, y = y_boot))  
}



#-------------------------------------------
#           Bootstrap estimate             #
#-------------------------------------------


boot_est = function(Data, est,method, bt, const){
  b_sample = boot_sample(Data = Data, est = est,method = method)
  
  b_x = b_sample$x; b_z = b_sample$z; b_y = b_sample$y
  
  boot_result = MIQP_iter_fast(x = b_x, y = b_y, z=  b_z, gs = T, bt = bt, const = const,
                               Grid = Grid, U.gm = U.gm, L.gm = L.gm, max.iter = 10)
  return(c(boot_result$gamma[c(1,3),]))
}