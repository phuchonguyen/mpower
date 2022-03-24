# Imputes missing data due to NAs and LOD
# Augments the observed data with imputed data
update_na_lod = function(params, constants) {
  
  if(constants$na) {
    params = update_factor_na(params, constants)
  }
  
  if(constants$lod){
    params = update_lod(params, constants)
  }
  
  if(constants$na | constants$lod){
    y = params$y
    X = params$X
    X_int = params$X_int
    l = params$l
    rho = params$rho
    p = constants$p
    n = constants$n
    exp_C = constants$exp_C
    
    # update interactions matrix
    X_int = model.matrix(y~.^2 - 1,as.data.frame(X))
    X_int = X_int[,(p+1):ncol(X_int)]
    
    # update projection matrix
    P = diag(n) - X%*%solve(t(X)%*%X)%*%t(X)
    Pt = t(P)
    X_big_proj_train = rbind(X,X)
    P_big_train = diag(2*n) - X_big_proj_train%*%solve(t(X_big_proj_train)%*%X_big_proj_train)%*%t(X_big_proj_train)
    
    params[["P"]] = P
    params[["X_big_proj_train"]] = X_big_proj_train
    params[["P_big_train"]] = P_big_train
    params[["X_int"]] = X_int
  }
  
  return(params)
}

# Imputes the missing data in X, Z
# Updates the latent factor params
update_factor_na = function(params, constants) {
  
  eta_na = params$eta_na
  Lambda_na = params$Lambda_na
  ps_na = params$ps_na
  Plam_na = params$Plam_na
  W = params$W
  XZ = params$XZ
  
  list_na = constants$list_na
  W_na = constants$W_na
  p = constants$p
  p_na = constants$p_na
  n = constants$n
  k_na = constants$k_na
  
  # update missing data
  W_pred = eta_na %*% t(Lambda_na) + mvtnorm::rmvnorm(n, sigma = diag(1/ps_na))
  W[W_na] = W_pred[W_na]
  XZ[,list_na] = W
  
  # update Sigma 
  Wtil = W - eta_na %*% t(Lambda_na) 
  ps_na = rgamma(p_na, 1 + 0.5*n, 1+0.5*colSums(Wtil^2))
  Sigma_na = diag(1/ps_na)
  
  # update Lambda
  eta2 = t(eta_na) %*% eta_na
  for (j in 1:p_na) {
    Qlam = diag(Plam_na[j, ]) + ps_na[j] * eta2
    blam = ps_na[j] * (t(eta_na) %*% W[, j])
    Llam = t(chol(Qlam))
    zlam = rnorm(k_na)
    vlam = solve(Llam, blam)
    mlam = solve(t(Llam), vlam)
    ylam = solve(t(Llam), zlam)
    Lambda_na[j, ] = t(ylam + mlam)
  }
  
  
  # update Eta
  Lmsg = Lambda_na * ps_na
  Veta1 = diag(k_na) + t(Lmsg) %*% Lambda_na
  eigs = eigen(Veta1)
  if(all(eigs$values > 1e-6)) {
    Tmat = sqrt(eigs$values) * t(eigs$vectors)
  } else {
    Tmat = chol(Veta1)
  }
  R = qr.R(qr(Tmat))
  S = solve(R)
  Veta = S %*% t(S)                                              
  Meta = W %*% Lmsg %*% Veta                                    
  eta_na = Meta + matrix(rnorm(n*k_na), nrow = n, ncol = k_na) %*% t(S)
  
  params[["eta_na"]] = eta_na
  params[["Lambda_na"]] = Lambda_na
  params[["ps_na"]] = ps_na
  params[["Plam_na"]] = Plam_na
  params[["Sigma_na"]] = Sigma_na
  params[["W"]] = W
  params[["XZ"]] = XZ
  params[["X"]] = XZ[,1:p]
  params[["Z"]] = XZ[,(p+1):ncol(XZ)]
  
  return(params)  
}

# Use truncated Normal to impute data below LOD
update_lod = function(params, constants) {
  
  # --- update data under LOD --- #
  X_lod = constants$X_lod
  lod_vec = constants$lod_vec
  X = params$X
  Lambda_na = params$Lambda_na
  eta_na = params$eta_na
  ps_na = params$ps_na
  
  for (i in 1:nrow(X_lod)) {
    for(j in 1:ncol(X_lod)){
      if(X_lod[i,j]){
        X[i, j] = truncnorm::rtruncnorm(n = 1, a = -Inf, 
                                        b = lod_vec[j],
                                        mean = Lambda_na[j, ] %*% eta_na[i, ], 
                                        sd = sqrt(1/ps_na[j]) )
      }
    }
  }
  
  params[["X"]] = X
  
  return(params)
}