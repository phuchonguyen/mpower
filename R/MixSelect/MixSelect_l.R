update_l = function(params, constants) {
  
  rho = params$rho; sigmasq = params$sigmasq
  X = params$X; mu = params$mu; CI = params$CI
  CI_inv = params$CI_inv; P = params$P; Pt = t(P)
  l = params$l; logdetCI = params$logdetCI
  gamma_l = params$gamma_l; PCPt = params$PCPt
  p = constants$p; n = constants$n
  exp_C = constants$exp_C
  
  # ret = update_l_Rcpp(rho, sigmasq, X,
  #                     mu, CI, CI_inv,
  #                     P, l, logdetCI,
  #                     gamma_l, PCPt,
  #                     p, n, exp_C)
  
  # --- updateL's: one add or remove move ---
  if(rho == 0){

    l = rep(0,p)
    gamma_l = rep(0,p)
    X_sweep = sweep(X, 2, l, `*`)
    C = matrix(0,n,n)
    PCPt = matrix(0,n,n)
    CI = diag(n)*sigmasq
    chol = chol(CI)
    logdetCI = 2*as.numeric(sum(log((diag(chol)))))
    CI_inv = diag(n)*(1/sigmasq)

  }

  if (rho != 0){
    set = sample(1:p,1)
    for (j in set){

      if (gamma_l[j] == 0){
        #propose from the prior so that log prop cancels out with the
        l_starj = rgamma(1,1)
        gamma_star = gamma_l
        gamma_star[j] = 1
        l_star = l; l_star[j] = l_starj
        X_sweep_star = sweep(X, 2, l_star, `*`)

      }

      if (gamma_l[j] == 1){
        l_starj = 0
        gamma_star = gamma_l
        gamma_star[j] = 0
        l_star = l; l_star[j] = l_starj
        X_sweep_star = sweep(X, 2, l_star, `*`)
        logdet_prior = - log(9/10) + log(1/10)
      }

      C_star = rho*fields::Exp.cov(X_sweep_star, x2=NULL, p = exp_C)
      PCPt_star = P%*%C_star%*%t(P)
      CI_star = PCPt_star + diag(n)*sigmasq

      chol = chol(CI_star)
      logdetCI_star = 2*as.numeric(sum(log((diag(chol)))))
      CI_inv_star = chol2inv(chol)
      # logdetCI_star = det(CI_star)
      # CI_inv_star = solve(CI_star)

      logr = logdetCI_star - logdetCI +
        t(mu)%*%(CI_inv_star - CI_inv)%*%mu
      logr = -0.5*logr

      logu = log(runif(1))

      if (logr > logu){
        l = l_star
        X_sweep = X_sweep_star
        gamma_l = gamma_star
        C = (C_star + t(C_star))/2
        CI = (CI_star + t(CI_star))/2
        CI_inv = (CI_inv_star + t(CI_inv_star))/2
        logdetCI = logdetCI_star
        PCPt = PCPt_star
      }
    }
  }

  ind = which(l != 0)

  if (length(ind) > 0){

    if (length(ind) == 1){
      set = ind
    }else{
      set = sample(c(ind),1)
    }

    for (j in set){

      l_starj = rgamma(1,l[j])
      l_star = l; l_star[j] = l_starj
      X_sweep_star = sweep(X, 2, l_star, `*`)

      C_star = rho*fields::Exp.cov(X_sweep_star, x2=NULL, p = exp_C)
      PCPt_star = P%*%C_star%*%Pt
      CI_star = PCPt_star + diag(n)*sigmasq

      chol = chol(CI_star)
      logdetCI_star = 2*as.numeric(sum(log((diag(chol)))))
      CI_inv_star = chol2inv(chol)
      # logdetCI_star = det(CI_star)
      # CI_inv_star = solve(CI_star)

      loglik = logdetCI_star - logdetCI +
        t(mu)%*%(CI_inv_star - CI_inv)%*%mu
      loglik = -0.5*loglik

      logprop = dgamma(l[j],l_starj,log = T) - dgamma(l_starj,l[j],log = T)
      logprior =  0.5*(dgamma(l_starj,1,log = T) - dgamma(l[j],1,log = T))

      logr = loglik + logprop + logprior
      logu = log(runif(1))

      if (logr > logu){
        l = l_star
        C = (C_star + t(C_star))/2
        X_sweep = X_sweep_star
        CI = (CI_star + t(CI_star))/2
        CI_inv = (t(CI_inv_star) + CI_inv_star)/2
        logdetCI = logdetCI_star
        PCPt = PCPt_star
      }
    }
  }
  
  params[["l"]] = l
  params[["gamma_l"]] = gamma_l
  params[["CI"]] = CI
  params[["CI_inv"]] = CI_inv
  params[["logdetCI"]] = logdetCI
  params[["PCPt"]] = PCPt
  
  return(params)
}