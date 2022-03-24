# One Gibbs step for updating the parameters of the predictors in X
update_gamma_beta = function(params, constants) {
  
  y = params$y; Z = params$Z
  X = params$X; X_int = params$X_int
  gamma = params$gamma; lambda = params$lambda
  beta = params$beta; beta_z = params$beta_z
  pi_prob = params$pi_prob; CI_inv = params$CI_inv
  p = constants$p; heredity = constants$heredity
  
  # update gamma
  set = sample(1:p,p)
  for (j in set){
    csi = y - (X_int%*%lambda + Z%*%beta_z)
    
    # j-th variable NOT in the model
    gamma[j] = 0
    ind = which(gamma != 0)
    if(length(ind) == 0){
      logpi_0 = 0
    }else{
      X_0 = X[,ind]
      A_delta_0_inv = t(X_0)%*%CI_inv%*%X_0 + diag(length(ind))
      m_0 = t(X_0)%*%CI_inv%*%csi
      S_0 = t(m_0)%*%solve(A_delta_0_inv)%*%m_0
      det0 = determinant(A_delta_0_inv,log = T)
      logpi_0 = -0.5*det0$modulus[1] +
        0.5*S_0 -
        0.5*log(det(diag(length(ind))))
      
      logpi_0 = as.numeric(logpi_0)
    }
    
    # j-th variable in the model
    gamma[j] = 1
    ind = which(gamma != 0)
    X_1 = X[,ind]
    A_delta_1_inv = t(X_1)%*%CI_inv%*%X_1 + diag(length(ind))
    m_1 = t(X_1)%*%CI_inv%*%csi
    S_1 = t(m_1)%*%solve(A_delta_1_inv)%*%m_1
    det1 = determinant(A_delta_1_inv,log = T)
    
    logpi_1 = -0.5*det1$modulus[1] + 
      0.5*S_1 -
      0.5*log(det(diag(length(ind))))
    logpi_1 = as.numeric(logpi_1)
    
    # update gamma_j based on likelihood ratio
    R_j = exp(logpi_0 - logpi_1)
    prob_inv = 1 + R_j*(1-pi_prob)/pi_prob
    gamma[j] = rbinom(1,1,1/prob_inv)
    
    # update beta_j based on gamma_j
    if(gamma[j] == 0){
      beta[j] = 0
    }
    
    # set to zero the lambdas according to heredity condition
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
      stop("provide proper heredity condition")
    }
    H = G[lower.tri(G)]
    ind = which(H == 0)
    lambda[ind] = 0
  }
  
  # --- sample beta_j diff from zero ---
  csi = y - (X_int%*%lambda + Z%*%beta_z)
  ind = which(gamma != 0)
  if(length(ind) > 0){
    X_1 = X[,ind]
    D_inv = diag(length(ind))
    V = solve(t(X_1)%*%CI_inv%*%X_1 + D_inv)
    m = V%*%(t(X_1)%*%CI_inv%*%csi)
    beta[ind] = as.vector(bayesSurv::rMVNorm(n = 1, mean = m, Sigma = V))
    beta[-ind] = 0
  }
  
  # --- update pi --- #
  d1 = sum(gamma)
  pi_prob = rbeta(1,1+d1,20+p-d1)
  
  params[["gamma"]] = gamma
  params[["beta"]] = beta
  params[["lambda"]] = lambda
  params[["pi_prob"]] = pi_prob
  #params["G"] = G  # G is reset every time?

  return(params)
}

