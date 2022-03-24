update_int = function(params, constants) {
  
  lambda = params$lambda; pi_int = params$pi_int
  gamma = params$gamma; gamma_int = params$gamma_int
  beta = params$beta; beta_z = params$beta_z
  CI_inv = params$CI_inv; Z = params$Z; X = params$X
  X_int = params$X_int; y = params$y
  heredity = constants$heredity
  p = constants$p; tau = constants$tau
  c = constants$c; lambda_prior = constants$lambda_prior
  
  # --- update lambda --- #
  G = matrix(0,p,p)
  if(heredity == "strong"){
    G = gamma%*%t(gamma)
  }else if(heredity == "weak"){
    
    for(h in 1:p){
      if(gamma[h] != 0){
        G[h,] = 1
        G[,h] = 1
      }
    }
  }else{
    stop("provide proper heredity condition ")
  }
  
  if(lambda_prior == "cont_spike"){
    # continuous spike
    H = G[lower.tri(G)]
    ind = which(H != 0)
    
    if (length(ind)>0){
      lambda[-ind] = 0 
      X_int_curr = X_int[,ind]    
      E_inv = diag(1/(gamma_int[ind]*(tau)^2 + (1-gamma_int[ind])*(c*tau)^2),
                   nrow = length(ind))
      V = solve(t(X_int_curr)%*%CI_inv%*%X_int_curr + E_inv)
      csi = y - (Z%*%beta_z + X%*%beta)
      m = V%*%(t(X_int_curr)%*%CI_inv%*%csi)
      lambda[ind] = as.vector(bayesSurv::rMVNorm(n = 1, mean = m, Sigma = V))
    }else{
      lambda = rep(0,length(gamma_int))
    }
    
    # --- update gamma_int --- #
    for (j in 1:length(lambda)){
      if (H[j] == 0){
        gamma_int[j] = 0
        lambda[j] = 0
      }else{
        b = dnorm(lambda[j],0,(c*tau))
        a = dnorm(lambda[j],0,tau)
        gamma_int[j] = rbinom(1,1,a/(a+b))
      }
    }
  }
  
  if(lambda_prior != "cont_spike"){
    # discrete spike
    # very slow especially with weak heredity
    H = G[lower.tri(G)]
    ind_int = which(H != 0)
    lambda[-ind_int] = 0
    gamma_int[-ind_int] = 0
    csi = y - (X%*%beta + Z%*%beta_z)
    
    if (length(ind_int)>0){
      for (j in ind_int){
        
        # j-th interaction variable not in the model
        gamma_int[j] = 0
        ind = which(gamma_int != 0)
        if(length(ind)==0){
          logpi_0 = 0
        }else{
          X_0 = X_int[,ind]
          A_delta_0_inv = t(X_0)%*%CI_inv%*%X_0 + diag(length(ind))
          m_0 = t(X_0)%*%CI_inv%*%csi
          S_0 = t(m_0)%*%solve(A_delta_0_inv)%*%m_0
          det0 = determinant(A_delta_0_inv,log = T)
          logpi_0 = -0.5*det0$modulus[1] +
            0.5*S_0-
            0.5*log(det(diag(length(ind))))
          logpi_0 = as.numeric(logpi_0)
        }
        
        # j-th variable in the model
        gamma_int[j] = 1
        ind = which(gamma_int != 0)
        X_1 = X_int[,ind]
        A_delta_1_inv = t(X_1)%*%CI_inv%*%X_1 + diag(length(ind))
        m_1 = t(X_1)%*%CI_inv%*%csi
        S_1 = t(m_1)%*%solve(A_delta_1_inv)%*%m_1
        det1 = determinant(A_delta_1_inv,log = T)
        logpi_1 = -0.5*det1$modulus[1] +
          0.5*S_1-
          0.5*log(det(diag(length(ind))))
        logpi_1 = as.numeric(logpi_1)
        
        # set gamma_j based on likelihood ratio
        R_j = exp(logpi_0 - logpi_1)
        prob_inv = 1+ R_j*(1-pi_int)/pi_int
        gamma_int[j] = rbinom(1,1,1/prob_inv)
        
        if(gamma_int[j] == 0){
          lambda[j] = 0
        }
      }
    }else{
      lambda = rep(0,length(gamma_int))
      gamma_int = rep(0,length(gamma_int))
    }
    
    # --- update lambda_j diff from zero ---
    
    ind = which(gamma_int != 0)
    if(length(ind) > 0){
      X_1 = X_int[,ind]
      D_inv = diag(length(ind))
      V = solve(t(X_1)%*%CI_inv%*%X_1 + D_inv)
      m = V%*%(t(X_1)%*%CI_inv%*%csi)
      lambda[ind] = as.vector(bayesSurv::rMVNorm(n = 1, mean = m, Sigma = V))
    }
  }
  
  # --- update pi_int ---
  d2 = sum(gamma_int)
  pi_int = rbeta(1,1+d2,5+ncol(X_int)-d2)
  

  params[["pi_int"]] = pi_int
  params[["lambda"]] = lambda
  params[["gamma_int"]] = gamma_int
  
  return(params)
}