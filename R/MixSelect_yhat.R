predict_insample = function(params, constants) {
  # in sample prediction
  rho = params$rho
  P_big_train = params$P_big_train
  sigmasq = params$sigmasq
  beta_z = params$beta_z
  beta = params$beta
  y = params$y
  Z = params$Z
  X = params$X
  X_int = params$X_int
  lambda = params$lambda
  l = params$l
  p = constants$p
  n = constants$n
  exp_C = constants$exp_C
  
  X_sweep = sweep(X, 2, l, `*`)
  K_big_train = rho*fields::Exp.cov(rbind(X_sweep,X_sweep), p = exp_C)
  PKPt_big_train = P_big_train%*%K_big_train%*%t(P_big_train) + sigmasq*diag(2*n)
  PKPt_12_train = PKPt_big_train[1:n,(n+1):(2*n)]
  Sig_22_inv_train = solve(PKPt_big_train[(n+1):(2*n),(n+1):(2*n)])
  mu = Z%*%beta_z + X%*%beta + X_int%*%lambda
  
  y_hat = mu + PKPt_12_train%*%Sig_22_inv_train%*%(y-mu)
  
  return(y_hat)
}