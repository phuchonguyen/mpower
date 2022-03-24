update_sigmasq = function(params, constants) {
  
  sigmasq = params$sigmasq
  mu = params$mu
  CI_inv = params$CI_inv
  CI = params$CI
  logdetCI = params$logdetCI
  PCPt = params$PCPt
  n = constants$n
  
  csi = params$y - (params$X_int%*%params$lambda + params$X%*%params$beta + params$Z%*%params$beta_z)
  V = CI/sigmasq
  V_inv = CI_inv*sigmasq
  logdetV = logdetCI - n*log(sigmasq)
  wss = t(csi)%*%V_inv%*%csi
  # Update CI_inv
  sigmasq = 1/rgamma(1, (1+n)*0.5, (1+wss)*0.5)
  CI_inv = V_inv/sigmasq
  CI = V*sigmasq
  logdetCI = logdetV + n*log(sigmasq)
  
  params[["sigmasq"]] = sigmasq
  params[["CI"]] = CI
  params[["CI_inv"]] = CI_inv
  params[["logdetCI"]] = logdetCI
  
  return(params)
}