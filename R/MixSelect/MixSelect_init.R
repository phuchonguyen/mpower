# Initializes parameters
# Sets NAs in X, Z to zeros
# Sets LOD indicators of NAs to FALSE
# Sets measurements below LOD to zeros
# Creates interaction matrices

init_params = function(y, X , Z,
                       nrun, burn, thin,
                       tau, c, exp_C,
                       heredity, lambda_prior,
                       na, k_na, list_na, sd_lambda,
                       lod, X_lod, lod_vec) {
  
  n = nrow(X); p = ncol(X)
  
  # --- interactions --- #
  X_int = model.matrix(y~.^2 - 1,as.data.frame(X))
  X_int = X_int[,(p+1):ncol(X_int)]
  
  # --- projection matrix --- #
  P = diag(n) - X%*%solve(t(X)%*%X)%*%t(X)
  Pt = t(P)
  X_big_proj_train = rbind(X,X)
  P_big_train = diag(2*n) - X_big_proj_train%*%solve(t(X_big_proj_train)%*%X_big_proj_train)%*%t(X_big_proj_train)
  
  # --- initial value for rho --- #
  if(exp_C == 1){
    rho = 4
  }else if(exp_C == 2){
    rho = 1
  }else{
    rho = 1
  }
  
  # --- other initial values --- #
  pi_prob = 0.5
  pi_int = 0.1
  gamma = rbinom(p,1,pi_prob)
  gamma_int = rbinom(ncol(X_int),1,pi_int)
  beta = rep(0,p)
  lambda = rep(0,ncol(X_int))
  beta_z = rep(0,ncol(Z))
  gamma_z = rep(0,ncol(Z))
  sigmasq = 1
  l = rep(1,p)
  gamma_l = rep(1,p)
  
  # --- covariance function --- #
  X_sweep = sweep(X, 2, l, `*`)
  C = rho*fields::Exp.cov(X_sweep, x2=NULL, p = exp_C)
  PCPt = P%*%C%*%Pt
  CI = PCPt + sigmasq*diag(n)
  chol = chol(CI)
  CI_inv = chol2inv(chol)
  logdetCI = 2*as.numeric(sum(log((diag(chol)))))
  
  
  params = list(
    PCPt = PCPt,
    CI = CI,
    CI_inv = CI_inv,
    logdetCI = logdetCI,
    y = y,
    X = X,
    Z = Z,
    X_int = X_int,
    P = P,
    X_big_proj_train = X_big_proj_train,
    P_big_train = P_big_train,
    rho = rho,
    pi_prob = pi_prob,
    pi_int = pi_int,
    gamma = gamma,
    gamma_int = gamma_int,
    beta = beta,
    lambda = lambda,
    beta_z = beta_z,
    gamma_z = gamma_z,
    sigmasq = sigmasq,
    l = l,
    gamma_l = gamma_l
  )
  
  constants = list(
    n = n,
    p = p,
    tau = tau,
    c = c,
    exp_C = exp_C,
    heredity = heredity,
    lambda_prior = lambda_prior,
    na = na,
    lod = lod,
    k_na = k_na,
    list_na = list_na,
    sd_lambda = sd_lambda,
    lod_vec = lod_vec
  )
  
  # --- initial values for Factor model/NA --- #
  if(na){
    
    # W is the matrix of variables to be included in the Factor Model and NA-imputation
    XZ = cbind(X,Z)
    if(is.null(list_na)) list_na = rep(T,ncol(XZ))
    W = XZ[,list_na]
    p_na = ncol(W)
    if(is.null(k_na)) k_na = floor(2*log(p_na))
    W_na = W %>% is.na()
    
    # set na to zero just to compute some of the quantities later
    X_na = X %>% is.na()
    X[X_na] = 0
    Z_na = Z %>% is.na()
    Z[Z_na] = 0
    
    # parameters of factor model
    Lambda_na = matrix(0, nrow = p_na, ncol = k_na)
    Plam_na =  matrix(sd_lambda, nrow = p_na, ncol = k_na)   
    eta_na = matrix(0, nrow = n, ncol = k_na)
    ps_na = rep(1,p_na)
    
    params[["XZ"]] = XZ
    params[["W"]] = W
    params[["X"]] = X
    params[["Z"]] = Z
    params[["Lambda_na"]] = Lambda_na
    params[["Plam_na"]] = Plam_na
    params[["eta_na"]] = eta_na
    params[["ps_na"]] = ps_na
    
    constants[["W_na"]] = W_na
    constants[["p_na"]] = p_na
    constants[["k_na"]] = k_na
  }
  
  # --- LOD --- #
  if(lod){
    if (is.null(X_lod)) stop("argument X_lod missing")
    if (is.null(lod_vec)) stop("argument lod_vec missing")
    X[X_lod] = 0
    X_na = X %>% is.na()
    X_lod[X_na] = F
    
    params[["X"]] = X
    constants[["X_lod"]] = X_lod
  }
  
  return(list(constants = constants,
              params = params))
}

# Calculate mean vector
get_mu = function(params, constants) {
  return(params$y - (params$Z%*%params$beta_z + 
                       params$X%*%params$beta + 
                       params$X_int%*%params$lambda))
}