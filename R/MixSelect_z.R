update_gamma_beta_z = function(params, constants) {
  
  y = params$y
  Z = params$Z
  X = params$X
  X_int = params$X_int
  beta = params$beta
  lambda = params$lambda
  gamma_z = params$gamma_z
  beta_z = params$beta_z
  CI_inv = params$CI_inv
  tau = constants$tau
  c = constants$c
  
  # --- update beta_z --- #
  D_inv = diag(1/(gamma_z*(tau)^2 + (1-gamma_z)*(tau*c)^2),nrow = length(beta_z))
  V = solve(t(Z)%*%CI_inv%*%Z + D_inv)
  csi = y - (X_int%*%lambda + X%*%beta)
  m = V%*%(t(Z)%*%CI_inv%*%csi)
  beta_z = as.vector(bayesSurv::rMVNorm(n = 1, mean = m, Sigma = V))
  
  # --- update gamma_z --- with continuos spike#
  for (j in 1:length(beta_z)){
    b = dnorm(beta_z[j],0,c*tau)
    a = dnorm(beta_z[j],0,tau)
    gamma_z[j] = rbinom(1,1,a/(a+b))
  }
  
  params[["beta_z"]] = beta_z
  params[["gamma_z"]] = gamma_z
  
  return(params)
}