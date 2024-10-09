# This file contains all codes for the MIQP algorithm, the estimation for gamma and beta, the evaluation for the estimation, and the data generation process

#------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())

library(dplyr)
library(Matrix)
library(gurobi)

#-------------------------------------------
#           Auxillary functions            #
#-------------------------------------------

get_m <-function(z,A,b,flag.print=FALSE){
  
  # Call necessary libraries
  library('gurobi')   # Gurobi library for Linear Programming
  
  # Declare parameters
  dim.gm = length(as.numeric(z))    # The dimension of gamma
  dim.b = length(b)                 # The number of constraint for parameter space
  params <- list(OutputFlag=0)      # Parameters for Gurobi
  
  # Declare a model for positive side
  model.pos <- list()
  
  A.1 = rbind( A, -z )
  model.pos$A          = A.1
  model.pos$obj        = z
  model.pos$modelsense = "max"
  model.pos$rhs        = c(b,0)
  model.pos$sense      = c(rep('<=',dim.b+1))
  model.pos$lb         = c(rep(-10^5,dim.gm))   # w/o lb Gurobi automatically impose 0 as a lower bound. 
  
  result.pos <- gurobi(model.pos, params)
  
  # Declare a model for negative side
  model.neg <- list()
  
  A.2 = rbind( A, z )
  model.neg$A          = A.2
  model.neg$obj        = -z
  model.neg$modelsense = "max"
  model.neg$rhs        = c( b , 0)
  model.neg$sense      = c(rep('<=',dim.b+1))
  model.neg$lb         = c(rep(-10^5,dim.gm))    # Lower bound for gamma is set to be -10^10 to allow negative values. 
  
  result.neg <- gurobi(model.neg, params)
  
  if (flag.print){
    cat('---------------------','\n')
    cat('Solutions','\n')
    cat('---------------------','\n')
    cat('objective value (positive):', result.pos$objval,'\n')
    cat('maximizer (positive):', result.pos$x,'\n')
    cat('objective value (negative):', result.neg$objval,'\n')
    cat('maximizer (negative):', result.neg$x,'\n')
  }
  
  # Check if the solution exists; then relocate the value by add/subtract the constant value f[1]
  if (!is.numeric(result.pos$objval)){
    result.pos$objval = 0
  }
  if (!is.numeric(result.neg$objval)){
    result.neg$objval = 0
  }
  
  if ( result.neg$objval > result.pos$objval ) {
    result = result.neg
  } else {
    result = result.pos
  }
  
  # Raise a flag when the result is not optmal
  if (result$status!="OPTIMAL") {
    cat('Warning: the result is not optimal!','\n')
  }
  
  return(list(m=result$objval, sol=result$x))
}


block_four = function(M){
  M.1 = cbind(M,M,M,M)
  M.2 = rbind(M.1,M.1,M.1,M.1)
  return(M.2)
}



get_Q = function(x, L.bt, d.z){
  library('Matrix')
  
  T = nrow(as.matrix(x))
  d.x = ncol(as.matrix(x))
  
  outer.x = array(NA, dim=c(d.x,d.x,T)) # d.x x d.x x T arrays to store column outer products
  
  for (i in c(1:T)){
    outer.x[,,i] = x[i,] %*% t(x[i,])
  }
  
  x11 = outer.x[,,1]
  Q.row = block_four(x11)
  Q.row.L = x11 %*% L.bt
  Q.diag.LL = block_four(x11)
  Q.diag.LI = block_four(x11 %*% L.bt)
  Q.diag.II = block_four(t(L.bt) %*% x11 %*% L.bt)
  
  for (i in c(2:T)){
    x.ii = outer.x[,,i]
    LL.i = block_four(x.ii)
    LI.i = block_four(x.ii %*% L.bt)
    II.i = block_four(t(L.bt) %*% x.ii %*% L.bt)
    
    Q.diag.LL = bdiag(Q.diag.LL, LL.i)
    Q.diag.LI = bdiag(Q.diag.LI, LI.i)
    Q.diag.II = bdiag(Q.diag.II, II.i)
  }
  
  Q1 = cbind(Q.diag.LL, Q.diag.LI)   # The first d.x x (T+1)d.x matrix
  Q2 = cbind(t(Q.diag.LI), Q.diag.II)  # Next rows of Q
  Q = rbind(Q1,Q2) / T
  
  
  # Add additional 0 columns
  Q = cbind(Q, matrix(0, nrow=nrow(Q), ncol=(4*d.x + 2*d.z+2*T)))
  # Add additional 0 rows
  Q = rbind(Q, matrix(0, nrow=(4*d.x + 2*d.z+2*T), ncol=ncol(Q)))
  
  return(Q)
}


get_L = function(x, y, L.bt, d.z){
  n.obs = length(y)
  y = as.matrix(y,n.obs,1)
  d.x = ncol(as.matrix(x))
  
  
  A.1 = y[1]%*%rep(x[1,],4)
  A.2 = t(y[1]*rep((x[1,]%*%L.bt),4))
  for (i in c(2:n.obs)){
    A.1 = cbind(A.1, y[i]%*%rep(x[i,],4))
    A.2 = cbind((A.2), t(y[i]*rep((x[i,]%*%L.bt),4)))
  }
  
  A = cbind(A.1,A.2)
  A = (-2/n.obs) * A
  L = cbind(A, matrix(0,1,(4*d.x + 2*d.z+2*n.obs))) 
  
  return(L)
}



#-------------------------------------------
#          MIQP--Joint optimization        #
#-------------------------------------------

joint_constraint = function(A, b, d.x, d.z, n.obs,M, z, eta = 1e-6){
  ncol.A = 4*d.x*n.obs + 4*n.obs + 2*n.obs + 2*d.z + 4*d.x
  nrow.A = 8*d.x + 8*n.obs*d.x + 4*n.obs + 8*n.obs +8*n.obs + 12*n.obs + 8 + 4*d.z
  A.const = matrix(0, nrow = nrow.A, ncol = ncol.A)
  I.dx = diag(d.x); I.dx.all = diag(4*d.x)
  I.df = diag(d.z)
  
  # const #3:  Left inequality of z_t'gamma
  A3.a = diag(c(M + 2*eta, M+2*eta)) 
  A3.b = diag(2)%x% (-z)
  A3.ab= cbind(matrix(0, nrow = nrow(A3.a), ncol = 4*d.x*n.obs + 4*n.obs), A3.a, A3.b, 
               matrix(0, nrow = nrow(A3.a), ncol =  4*d.x)) # 2*n.obs
  rm(A3.a, A3.b)#; gc()
  
  
  # const #3.5:  Right inequality of z_t'gamma
  # nrow.c = nrow.f + 1; nrow.f = nrow.f + 2*n.obs
  A3.c = -diag(c(M,M))
  A3.d = diag(2)%x% (z)
  A3.cd = cbind(matrix(0, nrow = nrow(A3.c), ncol = 4*d.x*n.obs + 4*n.obs), A3.c, A3.d, 
                matrix(0, nrow = nrow(A3.c), ncol =  4*d.x)) # 2*n.obs
  #  rm(A3.c, A3.d)#; gc()
  
  A.const = rbind(A, A3.ab, A3.cd)
  
  b.3 = c(M+eta, M+eta)
  #b.const = c(b.const, b.3)
  # b.const #3.5  Right inequality of z_t'gamma
  b.3.5 = rep(0,2*n.obs)
  #  b.const = c(b.const, b.3.5)
  
  b.const = c(b, b.3, b.3.5)
  return(list(A.const=A.const,b.const=b.const))
  
}



build_constraint = function(L.bt, U.bt, L.gm, U.gm, M, d.x, d.z, n.obs, z, tau1, tau2, eta=1e-6){
  
  # dim of parameters
  ncol.A = 4*d.x*n.obs + 4*n.obs + 2*n.obs + 2*d.z + 4*d.x
  nrow.A = 8*d.x + 8*n.obs*d.x + 4*n.obs + 8*n.obs +8*n.obs + 12*n.obs + 8 + 4*d.z
  A.const = matrix(0, nrow = nrow.A, ncol = ncol.A)
  I.dx = diag(d.x); I.dx.all = diag(4*d.x)
  I.df = diag(d.z)
  
  # const #1: beta tilde upper and bounds
  nrow.c = 1; nrow.f = 8*d.x
  A1 = rbind(I.dx.all, -I.dx.all) 
  A.const[nrow.c:nrow.f ,] = cbind(matrix(0, nrow=nrow(A1), ncol=ncol.A-4*d.x), A1) # 8*d.x
  
  # const #2:  0 <= l.tilde <= beta.tilde  
  # LHS is imposed by lower bounds later
  # RHS is l.tilde - beta.tilde <= 0 
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4*n.obs*d.x
  A2.a = I.dx.all %x% diag(n.obs) 
  A2.b = -rep(1,n.obs) %x% I.dx.all
  A.const[nrow.c:nrow.f ,] = cbind(A2.a, matrix(0, nrow = nrow(A2.a), ncol = ncol.A - 4*d.x*(n.obs+1)),A2.b)#4*n.obs*d.x
  rm(A2.a, A2.b)#; gc()
  
  # const #2.5: 0 <= l.tilde <= beta.tilde  
  # LHS inequality: -l.tilde <= 0
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4*n.obs*d.x
  A2.c = -I.dx.all %x% diag(n.obs)
  
  A.const[nrow.c:nrow.f ,] = cbind(A2.c, matrix(0, nrow = nrow(A2.c), ncol = ncol.A - 4*d.x*n.obs)) #4*n.obs*d.x
  rm(A2.c) 
  #  gc()
  
  # const #3:  Left inequality of z_t'gamma
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 2*n.obs
  A3.a = diag(c(M + 2*eta, M+2*eta)) 
  A3.b = diag(2)%x% (-z)
  A.const[nrow.c:nrow.f ,]= cbind(matrix(0, nrow = nrow(A3.a), ncol = 4*d.x*n.obs + 4*n.obs), A3.a, A3.b, 
                                  matrix(0, nrow = nrow(A3.a), ncol =  4*d.x)) # 2*n.obs
  rm(A3.a, A3.b)#; gc()
  
  
  # const #3.5:  Right inequality of z_t'gamma
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 2*n.obs
  A3.c = -diag(c(M,M))
  A3.d = diag(2)%x% (z)
  A.const[nrow.c:nrow.f ,] = cbind(matrix(0, nrow = nrow(A3.c), ncol = 4*d.x*n.obs + 4*n.obs), A3.c, A3.d, 
                                   matrix(0, nrow = nrow(A3.c), ncol =  4*d.x)) # 2*n.obs
  rm(A3.c, A3.d)#; gc()
  
  
  # const #4:  Right inequality sum_i l_{j,i,t} <= I_{j,t} * sum (U.dt - L.bt)
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4*n.obs
  sum.bnd = sum(U.bt - L.bt)
  A4.a = diag(n.obs) %x% (diag(4)%x% t(rep(1,d.x))) 
  A4.b = (-sum.bnd)*diag(n.obs*4) 
  A.const[nrow.c:nrow.f ,] = cbind(A4.a, A4.b, matrix(0, nrow = nrow(A4.a), ncol=2*n.obs+ 2*d.z+4*d.x)) # 4*n.obs
  rm(A4.a, A4.b)#; gc()
  
  
  # const #4.5  Left inequality [- sum l_{j,i,t}] <= 0
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4*n.obs
  A4.c = diag(n.obs) %x% (diag(4)%x% t(rep(-1,d.x)))
  A.const[nrow.c:nrow.f ,] = cbind(A4.c, matrix(0,ncol = ncol.A - 4*d.x*n.obs,nrow = nrow(A4.c)))  # 4*n.obs
  rm(A4.c)#; gc()
  
  
  # const # 5 Right inequality sum (beta - l) <= (1-I)* sum (U.dt - L.bt)
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4*n.obs
  A5.a = diag(n.obs) %x% (diag(4)%x% t(rep(-1,d.x)))
  A5.b = sum.bnd*diag(n.obs*4)
  A5.c = rep(1,n.obs)%x%(diag(4)%x%t(rep(1,d.x))) 
  A.const[nrow.c:nrow.f ,] = cbind(A5.a, A5.b,matrix(0, nrow = nrow(A5.a), ncol = 2*n.obs+ 2*d.z), A5.c)  # 4*n.obs
  rm(A5.a, A5.b, A5.c)#; gc()
  
  
  # const # 5.5 Left inequality sum(l - beta) <=0
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4*n.obs
  A5.d = diag(n.obs) %x% (diag(4)%x% t(rep(1,d.x)))
  A5.e = rep(1,n.obs)%x%(diag(4)%x%t(rep(-1,d.x)))
  A.const[nrow.c:nrow.f ,] = cbind(A5.d,matrix(0, nrow = nrow(A5.d), ncol = 6*n.obs+ 2*d.z),A5.e )  # 4*n.obs
  rm(A5.d, A5.e)#; gc()
  
  # const #6.1 I < f(d_1)
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4*n.obs
  A6.a = diag(4*n.obs)
  A6.b = diag(n.obs)%x%c(-1,1,1,-1)
  A.const[nrow.c:nrow.f ,] = cbind(matrix(0,nrow = nrow(A6.a), ncol = 4*d.x*n.obs),
                                   A6.a, A6.b, matrix(0,nrow = nrow(A6.a), ncol = n.obs),
                                   matrix(0,nrow = nrow(A6.a), ncol = 2*d.z + 4*d.x)) # 4*n.obs
  rm(A6.a, A6.b)#; gc()
  
  # const #6.2 I < f(d_2)
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4*n.obs
  A6.c = diag(4*n.obs)
  A6.d = diag(n.obs)%x%c(-1,-1,1,1)
  A.const[nrow.c:nrow.f ,] = cbind(matrix(0,nrow = nrow(A6.c), ncol = 4*d.x*n.obs),
                                   A6.c,matrix(0,nrow = nrow(A6.c), ncol = n.obs),A6.d,
                                   matrix(0,nrow = nrow(A6.c), ncol = 2*d.z + 4*d.x)) # 4*n.obs
  rm(A6.c, A6.d)#; gc()
  
  # const #6.3 -I + f(d_1,d_2) <= 0
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4*n.obs
  A6.e = -diag(4*n.obs)
  A6.f = diag(n.obs)%x%matrix(c(1,-1,-1,1),ncol = 1,byrow = T)
  A6.g = diag(n.obs)%x%matrix(c(1,1,-1,-1),ncol = 1,byrow = T)
  A.const[nrow.c:nrow.f ,] = cbind(matrix(0,nrow = nrow(A6.e), ncol = 4*d.x*n.obs),
                                   A6.e, A6.f, A6.g,
                                   matrix(0,nrow = nrow(A6.e), ncol = 2*d.z + 4*d.x)) # 4*n.obs
  rm(A6.e, A6.f, A6.g); gc()
  # const #7 Left inequality of sum I_{j,t} 
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4
  A7.a =  t(rep(-1/n.obs,n.obs)) %x% diag(4) 
  A.const[nrow.c:nrow.f ,] =  cbind(matrix(0,ncol =4*d.x*n.obs, nrow = nrow(A7.a)),
                                    A7.a, matrix(0,ncol = 2*n.obs + 2*d.z + 4*d.x, nrow = nrow(A7.a))) # 4
  
  # const #7.5  Right inequality of sum I_{j,t} 
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4
  A7.b =  t(rep(1/n.obs,n.obs)) %x% diag(4)  
  A.const[nrow.c:nrow.f ,] = cbind(matrix(0,ncol =4*d.x*n.obs, nrow = nrow(A7.b)),
                                   A7.b, matrix(0,ncol = 2*n.obs + 2*d.z + 4*d.x, nrow = nrow(A7.a))) # 4
  
  # const # gamma upper and lower bound
  nrow.c = nrow.f + 1; nrow.f = nrow.f +  4*d.z
  A8.a = rbind(diag(2*d.z), -diag(2*d.z)) 
  A.const[nrow.c:nrow.f ,] = cbind(matrix(0,nrow = nrow(A8.a), ncol = 4*d.x*n.obs + 4*n.obs + 2*n.obs),
                                   A8.a, matrix(0,nrow = nrow(A8.a), ncol = 4*d.x)) # 4*d.z
  
  # A.const.old = rbind(A1, A2, A2.5, A3, A3.5, A4, A4.5,
  #                  A5, A5.5, A6.1, A6.2, A6.3, A7, A7.5,A8)
  
  # b.const # 1 beta tilde upper and bounds
  b.1 = c(rep(1,4)%x%(U.bt - L.bt),rep(0,4*d.x))
  b.const = b.1
  # b.const #2 l.tilde - beta.tilde <= 0 
  b.2 = rep(0,4*d.x*n.obs)
  b.const = c(b.const, b.2)
  # b.const #2.5 -l.tilde <= 0 
  b.2.5 = rep(0,4*d.x*n.obs)
  b.const = c(b.const, b.2.5) 
  # b.const #3  Left inequality of z_t'gamma
  b.3 = c(M+eta, M+eta)
  b.const = c(b.const, b.3)
  # b.const #3.5  Right inequality of z_t'gamma
  b.3.5 = rep(0,2*n.obs)
  b.const = c(b.const, b.3.5)
  # b.const #4:  Right inequality sum_i l_{j,i,t} <= I_{j,t} * sum (U.dt - L.bt)
  b.4 = rep(0,4*n.obs)
  b.const = c(b.const, b.4)
  # b.const #4.5  Left inequality [- sum l_{j,i,t}] <= 0
  b.4.5 = rep(0,4*n.obs)
  b.const = c(b.const, b.4.5)
  # b.const # 5 Right inequality sum (beta - l) <= (1-I)* sum (U.dt - L.bt)
  b.5 = rep(sum.bnd, 4*n.obs)
  b.const = c(b.const, b.5)
  # b.const # 5.5 Left inequality sum(l - beta) <=0
  b.5.5 = rep(0, 4*n.obs)
  b.const = c(b.const, b.5.5)
  # b.const #6.1 I < f(d_1)
  b.6.1 = rep(c(0,1,1,0), n.obs)
  # b.const #6.2 I < f(d_2)
  b.6.2 = rep(c(0,0,1,1), n.obs)
  # b.const #6.3 -I + f(d_1, d_2) <=0
  b.6.3 = rep(c(1,0,-1,0), n.obs)
  b.const = c(b.const, b.6.1, b.6.2, b.6.3)
  # b.const #7 Left inequality of sum I_{j,t} 
  b.7 = rep(-tau1, 4)
  # b.const #7.5 Right inequality of sum I_{j,t} 
  b.7.5 = rep(tau2, 4)
  b.const = c(b.const, b.7, b.7.5)
  # const # gamma upper and lower bound
  b.8 = c(rep(U.gm,2), rep(-L.gm,2))
  b.const  = c(b.const, b.8)
  rm(b.1, b.2, b.2.5, b.3, b.3.5, b.4, b.4.5, 
     b.5, b.5.5, b.6.1, b.6.2, b.6.3, b.7, b.7.5, b.8)
  gc()
  return(list(A.const=A.const,b.const=b.const))
  
}


estimate <- function (y, x, z, Q.obj, L.obj, objcon, A.const, b.const, L.bt, L.gm, 
                      params=list(OutputFlag=1, Cuts = 3, TimeLimit = 15*60)) {
  
  # Call Library
  library("gurobi")
  
  model <- list()
  n.obs = length(y)
  d.x = ncol(x)
  d.z = ncol(z)
  
  # Quadratic objective function
  model$Q       = Q.obj
  model$obj     = L.obj
  model$objcon  = objcon
  
  # Linear constraints
  model$A       = A.const
  model$rhs     = b.const
  model$sense   = rep('<=', length(model$rhs))
  model$lb      = c(rep(0, 4*d.x*n.obs), rep(0,4*n.obs), rep(0,2*n.obs), 
                    rep(L.gm,2),rep(0,4*d.x))            #Done
  #model$lb      = c(rep(-10^10,d.x), rep(-10^10,d.x*n.obs), rep(-10^10,n.obs),  rep(-10^10,d.x), rep(-10^10,d.z))         #Done
  
  model$vtype   = c(rep('C', 4*d.x*n.obs), rep('B', 4*n.obs), rep('B',2*n.obs), 
                    rep('C', 2*d.z), rep('C', 4*d.x))   
  model$modelsense = 'min'
  result <- gurobi(model, params=params)
  
  return(result)
}



MIQP_joint = function(x, y, z, U.gm, L.gm,U.bt, L.bt, const){
  n.obs = nrow(x)
  d.z = ncol(z); d.x = ncol(x)
  A.gm = c(1,-1) %x% diag(rep(1,d.z)) 
  b.gm = c(U.gm, -L.gm)
  M = rep(NA,n.obs)
  for (i in (1:n.obs)){
    M[i] = get_m(z=z[i,],A=A.gm,b=b.gm)$m
  }  
  
  constraint = joint_constraint(A = const$A.const, b = const$b.const, d.x = d.x, d.z = d.z, n.obs = n.obs, M = M, z = z)
  # 1. Construct the objective function  
  L.obj = get_L(x=x,y=y,L.bt=L.bt, d.z=d.z)
  Q.obj = get_Q(x=x, L.bt=L.bt, d.z=d.z)
  objcon = mean(y^2)
  gc()
  # const = build_constraint(L.bt, U.bt, L.gm, U.gm, M, d.x, d.z, n.obs, z, tau1=tau1, tau2=tau2,  eta=1e-4)
  A.const = constraint$A.const; b.const = constraint$b.const
  result = estimate(y=y, x=x,z=z, Q.obj=Q.obj, L.obj=L.obj, objcon=objcon, A.const=A.const, b.const=b.const, L.bt=L.bt, L.gm=L.gm)
  rm(A.const,b.const,L.obj,Q.obj,M);gc()
  beta = result$x[(4*d.x*n.obs + 4*n.obs + 2*n.obs +2*d.z+ 1) : (4*d.x*n.obs + 4*n.obs + 2*n.obs + 2*d.z+4*d.x)] + rep(L.bt, 4)  
  gamma = result$x[(4*d.x*n.obs + 4*n.obs + 2*n.obs +1) : (4*d.x*n.obs + 4*n.obs + 2*n.obs + 2*d.z)]
  time = result$runtime
  objval = result$objval
  #  a1 = -gamma[1]/gamma[3]; b1 = (a1 + b1*x, from = -20, to = 20, ylim = c(-20,20))-gamma[2]/gamma[3]; a2 = -gamma[4]/gamma[6]; b2 =- gamma[5]/gamma[6]
  # curve
  #  curve(a2 + b2*x, from = -20, to = 20, ylim = c(-20,20), add = T)
  gamma = cbind(gamma[1:d.z], gamma[(d.z+1):(2*d.z)])
  beta = cbind(beta[1:d.x], beta[(d.x+1):(2*d.x)], beta[(2*d.x+1):(3*d.x)], beta[(3*d.x+1):(4*d.x)] )
  return(list(gamma = gamma, beta = beta, time = time,objval=objval))
}



#-------------------------------------------
#          MIQP--Iterative optimization    #
#-------------------------------------------


Grid_gen= function(U, L, s=20, p){
  U = U.gm; L = L.gm; nG = s^(p+1)
  dim = length(U)
  #dim.x = dim(x)[2]
  Grid = matrix(0.0, nrow =  nG, ncol = 2*dim)
  for(n in 1:dim){
    for(g in 1:nG){
      Grid[g,n] = sample(seq(L[n], U[n], 0.5), 1)
      Grid[g,n+dim] = sample(seq(L[n], U[n], 0.5), 1)
    }}
  return(Grid)
}


step1_grid = function(y, x, z,  grid){
  nG = nrow(grid)
  dim.x = dim(x)[2]
  dim.z = dim(z)[2]
  result = matrix(NA, nrow=nG, ncol=(1+4*dim.x))
  for(i in (1:nG)){
    z.index.1 = as.numeric(z %*% grid[i, 1:dim.z])
    z.index.2 = as.numeric(z %*% grid[i, (dim.z+1):(2*dim.z)])
    x.reg = cbind(x*(z.index.1 > 0)*(z.index.2 >0),
                  x*(z.index.1 <= 0)*(z.index.2 >0),
                  x*(z.index.1 <= 0)*(z.index.2 <=0),
                  x*(z.index.1 > 0)*(z.index.2 <=0))
    m = lm(y~x.reg-1)
    result[i,] = c(sum(m$resid^2),coef(m))
  }
  # result = na.omit(result)
  result[is.na(result)] <- 0
  opt = which.min(result[,1])
  bt.hat = result[opt,-1]
  gm.hat = as.numeric(Grid[opt,])
  
  return(list(bt.hat=bt.hat, gm.hat=gm.hat))
}


A_build = function(n.obs=n.obs, d.x=d.x, d.z = d.z,  M, z, tau1=tau1, tau2=tau2,  eta=1e-6){
  # dim of parameters {I_k,t; d_j,t; gamma_j}
  ncol.A = 4*n.obs + 2*n.obs + 2*d.z 
  nrow.A =  4*n.obs + 12*n.obs + 8 + 4*d.z
  A.const = matrix(0, nrow = nrow.A, ncol = ncol.A)
  I.dx = diag(d.x); I.dx.all = diag(4*d.x)
  I.df = diag(d.z)
  
  
  
  # const #3:  Left inequality of z_t'gamma
  nrow.c = 1; nrow.f = 2*n.obs
  A3.a = diag(c(M + 2*eta, M+2*eta)) 
  A3.b = diag(2)%x% (-z)
  A.const[nrow.c:nrow.f ,]= cbind(matrix(0, nrow = nrow(A3.a), ncol =  4*n.obs), A3.a, A3.b) # 2*n.obs
  rm(A3.a, A3.b)
  
  
  # const #3.5:  Right inequality of z_t'gamma
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 2*n.obs
  A3.c = -diag(c(M,M))
  A3.d = diag(2)%x% (z)
  A.const[nrow.c:nrow.f ,] = cbind(matrix(0, nrow = nrow(A3.c), ncol = 4*n.obs), A3.c, A3.d) # 2*n.obs
  rm(A3.c, A3.d); gc()
  
  # const #6.1 I < f(d_1)
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4*n.obs
  A6.a = diag(4*n.obs)
  A6.b = diag(n.obs)%x%c(-1,1,1,-1)
  A.const[nrow.c:nrow.f ,] = cbind(A6.a, A6.b, matrix(0,nrow = nrow(A6.a), ncol = n.obs),
                                   matrix(0,nrow = nrow(A6.a), ncol = 2*d.z)) # 4*n.obs
  rm(A6.a, A6.b)
  
  # const #6.2 I < f(d_2)
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4*n.obs
  A6.c = diag(4*n.obs)
  A6.d = diag(n.obs)%x%c(-1,-1,1,1)
  A.const[nrow.c:nrow.f ,] = cbind(A6.c,matrix(0,nrow = nrow(A6.c), ncol = n.obs),A6.d,
                                   matrix(0,nrow = nrow(A6.c), ncol = 2*d.z)) # 4*n.obs
  rm(A6.c, A6.d); gc()
  
  # const #6.3 -I + f(d_1,d_2) <= 0
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4*n.obs
  A6.e = -diag(4*n.obs)
  A6.f = diag(n.obs)%x%matrix(c(1,-1,-1,1),ncol = 1,byrow = T)
  A6.g = diag(n.obs)%x%matrix(c(1,1,-1,-1),ncol = 1,byrow = T)
  A.const[nrow.c:nrow.f ,] = cbind(A6.e, A6.f, A6.g,
                                   matrix(0,nrow = nrow(A6.e), ncol = 2*d.z)) # 4*n.obs
  rm(A6.e, A6.f, A6.g)
  
  # const #7 Left inequality of sum I_{j,t} 
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4
  A7.a =  t(rep(-1/n.obs,n.obs)) %x% diag(4) 
  A.const[nrow.c:nrow.f ,] =  cbind(A7.a, matrix(0,ncol = 2*n.obs + 2*d.z, nrow = nrow(A7.a))) # 4
  rm(A7.a); gc()
  # const #7.5  Right inequality of sum I_{j,t} 
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 4
  A7.b =  t(rep(1/n.obs,n.obs)) %x% diag(4)  
  A.const[nrow.c:nrow.f ,] = cbind(A7.b, matrix(0,ncol = 2*n.obs + 2*d.z , nrow = nrow(A7.b))) # 4
  
  # const # gamma upper and lower bound
  nrow.c = nrow.f + 1; nrow.f = nrow.f +  4*d.z
  A8.a = rbind(diag(2*d.z), -diag(2*d.z)) 
  rm(A7.b)
  A.const[nrow.c:nrow.f ,] = cbind(matrix(0,nrow = nrow(A8.a), ncol = 4*n.obs + 2*n.obs),
                                   A8.a) # 4*d.z
  
  rm(A8.a); gc()
  #-------------------------------------------------------
  # b.const
  # b.const #3  Left inequality of z_t'gamma
  b.3 = c(M+eta, M+eta)
  b.const = b.3
  # b.const #3.5  Right inequality of z_t'gamma
  b.3.5 = rep(0,2*n.obs)
  b.const = c(b.const, b.3.5)
  # b.const #4:  Right inequality sum_i l_{j,i,t} <= I_{j,t} * sum (U.dt - L.bt)
  # b.const #6.1 I < f(d_1)
  b.6.1 = rep(c(0,1,1,0), n.obs)
  # b.const #6.2 I < f(d_2)
  b.6.2 = rep(c(0,0,1,1), n.obs)
  # b.const #6.3 -I + f(d_1, d_2) <=0
  b.6.3 = rep(c(1,0,-1,0), n.obs)
  b.const = c(b.const, b.6.1, b.6.2, b.6.3)
  # b.const #7 Left inequality of sum I_{j,t} 
  b.7 = rep(-tau1, 4)
  # b.const #7.5 Right inequality of sum I_{j,t} 
  b.7.5 = rep(tau2, 4)
  b.const = c(b.const, b.7, b.7.5)
  # const # gamma upper and lower bound
  b.8 = c(rep(U.gm,2), rep(-L.gm,2))
  b.const  = c(b.const, b.8)
  rm(b.3, b.3.5, b.6.1, b.6.2, b.6.3, b.7, b.7.5, b.8)
  gc()
  
  return(list(A.const = A.const, b.const = b.const))
}


iter_constraint = function(A, b, n.obs, d.z, d.x, M, z, eta=1e-6){
  # dim of parameters {I_k,t; d_j,t; gamma_j}
  ncol.A = 4*n.obs + 2*n.obs + 2*d.z 
  nrow.A =  4*n.obs 
  A.const = matrix(0, nrow = nrow.A, ncol = ncol.A)
  I.dx = diag(d.x); I.dx.all = diag(4*d.x)
  I.df = diag(d.z)
  
  
  
  # const #3:  Left inequality of z_t'gamma
  nrow.c = 1; nrow.f = 2*n.obs
  A3.a = diag(c(M + 2*eta, M+2*eta)) 
  A3.b = diag(2)%x% (-z)
  A.const[nrow.c:nrow.f ,]= cbind(matrix(0, nrow = nrow(A3.a), ncol =  4*n.obs), A3.a, A3.b) # 2*n.obs
  rm(A3.a, A3.b)
  
  
  # const #3.5:  Right inequality of z_t'gamma
  nrow.c = nrow.f + 1; nrow.f = nrow.f + 2*n.obs
  A3.c = -diag(c(M,M))
  A3.d = diag(2)%x% (z)
  A.const[nrow.c:nrow.f ,] = cbind(matrix(0, nrow = nrow(A3.c), ncol = 4*n.obs), A3.c, A3.d) # 2*n.obs
  rm(A3.c, A3.d); gc()
  
  A.const = rbind(A.const, A)
  
  b.3 = c(M+eta, M+eta)
  b.const = c(b.3, b)
  return(list(A.const=A.const,b.const=b.const))
}

estimate_gm <- function(y, x, z, bt , M, U.gm, L.gm, A.const, b.const,
                        params=list(OutputFlag=0, FeasibilityTol=1e-6, DualReductions = 0, 
                                    NodeFileStart = 0.5,  TimeLimit = 2*60)){
  
  # Call necessary libraries
  library('gurobi')   # Gurobi library for MIO
  
  # Data dictionary: declare parameter values.
  n.obs = length(y)         # Sample size
  dim.gm = ncol(z)         # The dimension of gamma
  gm.hat = rep(NA,dim.gm)  # Maximizer of the problem
  # do we need?  d.x = ncol(x)         # The dimension of delta
  
  
  # Declare a model 
  model <- list()
  
  A.1 = (t(x[1,]) %*% bt)^2 - 2*y[1]*t(x[1,]) %*% bt
  for (i in c(2:n.obs)){
    A.1 = cbind(A.1,  (t(x[i,]) %*% bt)^2 - 2*y[i]*t(x[i,]) %*% bt )
  }
  
  A = (2/n.obs) * A.1
  L = cbind(A, matrix(0,1,(2*d.z+2*n.obs))) 
  
  model$obj = L
  
  
  model$A          = A.const
  # rm(A.const); gc()
  model$rhs        = b.const
  model$sense      = c(rep('<=',nrow(model$A)))
  
  # Other model parameter setting
  model$vtype   = c( rep('B', 4*n.obs), rep('B',2*n.obs), 
                     rep('C', 2*d.z))   
  model$modelsense = 'min'
  model$lb      = c( rep(0,4*n.obs), rep(0,2*n.obs), 
                     rep(L.gm,2))       
  result = gurobi(model, params)
  
  # Return the estimate for gamma 
  return(list(obj=result$objval, sol=result$x, gm=result$x[-(1: (4*n.obs + 2*n.obs))], result=result))  
  
}

estimate_bt = function(y, x, z, gm){
  z.index.1 = as.numeric(z %*% gm[1: dim(z)[2]])
  z.index.2 = as.numeric(z %*% gm[(dim(z)[2]+1): (2*dim(z)[2])])
  x.reg = cbind(x*(z.index.1 > 0)*(z.index.2 >0),
                x*(z.index.1 <= 0)*(z.index.2 >0),
                x*(z.index.1 <= 0)*(z.index.2 <=0),
                x*(z.index.1 > 0)*(z.index.2 <=0))
  m = lm(y~x.reg-1)
  bt = coef(m)
  bt[is.na(bt)] = 0
  bt = matrix(bt, ncol = 4, byrow = F)
  return(list(bt = bt, all_result = m))
}


MIQP_iter = function(x, y, z,gs = F, bt, Grid, U.gm, L.gm,U.bt, L.bt,tau1=0, tau2=0.9, max.iter= 30,pres= 1e-4){
  start_time = Sys.time()
  n.obs = nrow(x)
  d.x = ncol(x); d.z = ncol(z)
  
  if(gs == F){
    step1.out = step1_grid(y=y,x=x,z = z,grid=Grid)
    rm(Grid); gc()
    bt.hat.step1 = step1.out$bt.hat  
  }else{bt.hat.step1 = bt}
  
  
  
  A.gm = c(1,-1) %x% diag(rep(1,d.z)) 
  b.gm = c(U.gm, -L.gm)
  M = rep(NA,n.obs)
  for (i in (1:n.obs)){
    M[i] = get_m(z=z[i,],A=A.gm,b=b.gm)$m
  }  
  # step1.out = step1_grid(y=y,x=x,z = z,grid=Grid)
  #  rm(Grid); gc()
  # bt.hat.step1 = step1.out$bt.hat
  #  gm.hat.step1 = step1.out$gm.hat
  # Bt = matrix(bt.hat.step1, ncol = 4, byrow = F)
  #  cat(' Grid search gamma =', gm.hat.step1, '\n')
  # Bt = matrix(c(1,1,1, 1,1,0, 2, 0, -1, -1, 3, 0), ncol = 4, byrow = F)
  const = A_build(n.obs=n.obs, d.x=d.x, d.z = d.z, z = z, M = M, tau1=tau1, tau2=tau2,  eta=1e-5)
  A.const = const$A.const; b.const = const$b.const
  rm(const)
  Bt = matrix(bt.hat.step1, ncol = 4, byrow = F)
  for (cnt.it in (1:max.iter)){
    # Estimate gm.hat by MIO and Update ap.hat
    step2_1.out = estimate_gm(y=y, x=x, z= z, bt=Bt, L.gm=L.gm, U.gm =U.gm, A.const = A.const, b.const = b.const)
    gm.hat = step2_1.out$gm
    
    # Update ap.hat
    step2_2.out = estimate_bt(y=y, x=x, z=z, gm=gm.hat)
    bt.hat = step2_2.out$bt
    
    bt.diff = bt.hat - Bt
    if(sum(bt.diff^2)< pres){break}
    # Update the initial values of bt.hat and dt.hat
    Bt = bt.hat
    # cat(' Iteration =', cnt.it, '\n')
  }
  gm.hat = cbind(gm.hat[1:d.z], gm.hat[(d.z+1):(2*d.z)])
  rm(A.const); gc()
  end_time   = Sys.time()
  run_time   = as.numeric(end_time - start_time, units="secs")
  return(list(gamma = gm.hat, beta = Bt,time = run_time))
}


MIQP_iter_fast = function(x, y, z,gs = F, bt, Grid, const, U.gm, L.gm,  max.iter= 30,pres= 1e-6){
  start_time = Sys.time()
  n.obs = nrow(x)
  d.x = ncol(x); d.z = ncol(z)
  
  if(gs == F){
    step1.out = step1_grid(y=y,x=x,z = z,grid=Grid)
    rm(Grid); gc()
    bt.hat.step1 = step1.out$bt.hat  
  }else{bt.hat.step1 = bt}
  
  
  
  A.gm = c(1,-1) %x% diag(rep(1,d.z)) 
  b.gm = c(U.gm, -L.gm)
  M = rep(NA,n.obs)
  for (i in (1:n.obs)){
    M[i] = get_m(z=z[i,],A=A.gm,b=b.gm)$m
  }  
  # step1.out = step1_grid(y=y,x=x,z = z,grid=Grid)
  #  rm(Grid); gc()
  # bt.hat.step1 = step1.out$bt.hat
  #  gm.hat.step1 = step1.out$gm.hat
  # Bt = matrix(bt.hat.step1, ncol = 4, byrow = F)
  #  cat(' Grid search gamma =', gm.hat.step1, '\n')
  # Bt = matrix(c(1,1,1, 1,1,0, 2, 0, -1, -1, 3, 0), ncol = 4, byrow = F)
  const = iter_constraint(A = const$A.const, b = const$b.const, n.obs = n.obs, d.z=d.z, d.x=d.x, M=M, z=z)
  A.const = const$A.const; b.const = const$b.const
  rm(const)
  Bt = matrix(bt.hat.step1, ncol = 4, byrow = F)
  for (cnt.it in (1:max.iter)){
    # Estimate gm.hat by MIO and Update ap.hat
    step2_1.out = estimate_gm(y=y, x=x, z= z, bt=Bt, L.gm=L.gm, U.gm =U.gm, A.const = A.const, b.const = b.const)
    gm.hat = step2_1.out$gm
    
    # Update ap.hat
    step2_2.out = estimate_bt(y=y, x=x, z=z, gm=gm.hat)
    bt.hat = step2_2.out$bt
    
    bt.diff = bt.hat - Bt
    if(sum(bt.diff^2)< pres){break}
    # Update the initial values of bt.hat and dt.hat
    Bt = bt.hat
    # cat(' Iteration =', cnt.it, '\n')
  }
  gm.hat = cbind(gm.hat[1:d.z], gm.hat[(d.z+1):(2*d.z)])
  rm(A.const); gc()
  end_time   = Sys.time()
  run_time   = as.numeric(end_time - start_time, units="secs")
  return(list(gamma = gm.hat, beta = Bt,time = run_time))
}


#-------------------------------------------
#          Accessment for results          #
#-------------------------------------------

dist_v_H = function(v, H){
  bias = H-v
  mse = sapply(1:ncol(H), function(x){sum(bias[,x]^2)})
  ii = which.min(mse)
  bias = bias[,ii]
  return(list(ii = ii, bias = bias, mse = min(mse)))
}

dist_T_H = function(V, H){
  e = 0
  for(c in 1:ncol(V)){
    d = min(colSums((V[,c] - H)^2))
    e = e+d
  }
  return(e)
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



SLM_est <- function(x, y, z, U.gm, L.gm,U.bt, L.bt, tau1, tau2, bt.hat.step1, max.iter= 30, method="iter"){
  if(method=="iter"){
    start_time = Sys.time()
    est        = MIQP_iter(x = x, y = y, z=  z, bt.hat.step1 = bt.hat.step1, U.gm = U.gm, L.gm = L.gm, L.bt = L.bt, U.bt = U.bt, 
                           tau1=tau1,tau2 =tau2,max.iter = max.iter)
    end_time   = Sys.time()
    run_time   = as.numeric(end_time - start_time, units="secs")
    
  }else{
    start_time = Sys.time()
    est        = MIQP_joint(x = x, y = y, z=  z, U.gm = U.gm, L.gm = L.gm, L.bt = L.bt, U.bt = U.bt, 
                            tau1=tau1,tau2 =tau2)
    end_time   = Sys.time()
    run_time   = as.numeric(end_time - start_time, units="secs")
  }
  output = list(beta = est$beta, gamma = est$gamma, time = run_time)
  return(output)
}


#-------------------------------------------
#          Data generation process         #
#-------------------------------------------




DGP <- function(ar, ma, level, n, sd = 1){
  if(ar != 0 & ma == 0){
    ts = arima.sim(n = n, list(ar = level))
  }else if(ma!=0 & ar == 0){
    ts = arima.sim(n = n, list(ma = level))
  }else{ts = rnorm(n = n, sd =sd)}
  return(ts)
}



Model_DGP <- function(model, n, p.x = 3, p.z = 2, 
                      ar = 0, ma = 0, level = 0,  sd = 1 ){
  
  noise = rnorm(n, mean = 0, sd = sd)
  p = p.x + p.z
  d.z = p.z+1; d.x = p.x+1; n.obs = n
  data = matrix(0.0, nrow = n, ncol = p)
  for(i in 1:p.x){data[,i] = DGP(n = n, ar = ar, ma = ma, level = level, sd = sd)}
  for(c in (p.x+1):p){data[,c] = DGP(ar = ar, ma = ma, level = level, n = n); 
  data[,c] = data[,c]/max(abs(data[,c]))}
  data = as.data.frame(data)
  names(data) = c(sprintf("x%s", 1:p.x),sprintf("z%s", 1:p.z))
  
  
  if(model == "3A"){
    
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
    
  }else if(model == "3B"){
    
    gamma.0 = cbind(c(-1/4,1,1),c(1/4,1,1))
    beta.0 = cbind(c(1,1,1,1), c(-3,2,-1,0),
                   c(0,1,3,-1))
    
    x = as.matrix(data[,1:p.x])
    x = cbind(rep(1, n), x)
    
    z = as.matrix(data[,(1+p.x):p])
    z = cbind(rep(1, n), z)
    regime = regime_class(Z=z, gamma = gamma.0)
    regime = regime[,1:3]
    
    
    y = rowSums((x%*% beta.0) * regime) + noise
    
    
  }else if(model == "2A"){
    
    gamma.0 = c(0,1,-1)
    beta.0 = cbind(c(1,1,1,1), c(-3,2,-1,0))
    
    x = as.matrix(data[,1:p.x])
    x = cbind(rep(1, n), x)
    
    z = as.matrix(data[,(1+p.x):p])
    z = cbind(rep(1, n), z)
    # regime.o = regime_class(Z=z, gamma = gamma.0)
    
    z.prod = z %*% gamma.0
    z.1 = as.numeric(z.prod>=0)
    z.2 = as.numeric(z.prod<0)
    regime = cbind(z.1, z.2)
    
    
    y = rowSums((x%*% beta.0) * regime) + noise
    
  }else if(model == "2B"){
    
    gamma.0 = cbind(c(0,1,-1),c(0,1,1))
    beta.0 = cbind(c(1,1,1,1), c(-3,2,-1,0))
    
    x = as.matrix(data[,1:p.x])
    x = cbind(rep(1, n), x)
    
    z = as.matrix(data[,(1+p.x):p])
    z = cbind(rep(1, n), z)
    regime.o = regime_class(Z=z, gamma = gamma.0)
    regime = matrix(0.0, nrow = nrow(regime.o), ncol = 2)
    regime[,1] = regime.o[,1]; regime[,2] = regime.o[,2]+ regime.o[,3]+ regime.o[,4]
    
    
    y = rowSums((x%*% beta.0) * regime) + noise
    
  }else if(model == "1"){
    
    
    beta.0 = matrix(c(1,1,1,1),ncol = 1)
    
    x = as.matrix(data[,1:p.x])
    x = cbind(rep(1, n), x)
    
    z = as.matrix(data[,(1+p.x):p])
    z = cbind(rep(1, n), z)
    # regime.o = regime_class(Z=z, gamma = gamma.0)
    
    y = (x%*% beta.0) + noise
    regime = matrix(1, nrow = n, ncol = 1)
    
  }else if(model == "4"){
    
    gamma.0 = cbind(c(0,1,-1),c(0,1,1))
    beta.0 = cbind(c(1,1,1,1), c(-3,2,-1,0),
                   c(0,1,3,-1), c(2,-1,0,2))
    
    x = as.matrix(data[,1:p.x])
    x = cbind(rep(1, n), x)
    
    z = as.matrix(data[,(1+p.x):p])
    z = cbind(rep(1, n), z)
    regime = regime_class(Z=z, gamma = gamma.0)
    
    y = rowSums((x%*% beta.0) * regime) + noise
  }
  return(list(x = x, z = z, y = y, regime = regime, beta = beta.0))
  
}


DGP_vec = function(n=100, ar = 0, ma = 1, level = 0.1){
  Sigma = matrix(0.1, nrow = 5, ncol = 5)
  for(i in 1:5){Sigma[i,i]= 1}
  
  VAR = diag(rep(level*ar, 5))
  VMA = diag(rep(level*ma, 5))
  
  library(multiwave)
  XZ = varma(N = n, k= 5, VAR =VAR, VMA = VMA )
  
  x = cbind(1, XZ[,1:3])
  z = cbind(1, XZ[,4:5])
  for(c in 2:3){
    z[,c] = z[,c] / max(abs(z[,c]))
  }
  
  gamma.0 = cbind(c(0,1,-1),c(0,1,1))
  beta.0 = cbind(c(1,1,1,1), c(-3,2,-1,0),
                 c(0,1,3,-1), c(2,-1,0,2))
  noise = c()
  for(i in 1:n){
    noise[i] = rnorm(n=1, mean = 0, sd = 1+0.1*x[n,2]^2 + 0.1*z[n,2]^2)
  }
  
  regime = regime_class(Z=z, gamma = gamma.0)
  
  y = rowSums((x%*% beta.0) * regime) + noise
  
  
  return(list(x = x, y = y, z= z))
  
}
