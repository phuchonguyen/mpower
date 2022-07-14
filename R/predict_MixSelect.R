# Function to predict y on a test set given the output of the function MixSelect, X_test and Z_test
predict_MixSelect = function(gibbs, X_test = NULL, Z_test = NULL,
                             na = F, list_na = NULL, lod = F, X_lod = NULL){

  # X and y from gibbs object
  X = gibbs$X
  y = gibbs$y

  # interactions
  X_int = model.matrix(y~.^2 - 1,as.data.frame(X))
  X_int = X_int[,(p+1):ncol(X_int)]

  # checks
  if(is.null(X_test)) X_test = X
  p = ncol(X_test); n_test = nrow(X_test)
  if(ncol(gibbs$beta) != p) stop("Mismatching number of columns in X and X_test")

  y_test = numeric(n_test)
  X_int_test = model.matrix(y_test~.^2 - 1,as.data.frame(X_test))
  X_int_test = X_int_test[,(p+1):ncol(X_int)]

  if(is.null(Z_test)) Z_test = matrix(0, nrow = n_test, ncol = 1)
  if(gibbs$beta_z != ncol(Z_test)) stop("Mismatching number of columns in Z and Z_test")

  S = nrow(gibbs$beta)
  N = n_test + nrow(X)
  y_pred = matrix(0, nrow = S, ncol = n_test)

  X_big_proj = rbind(X_test,X)
  P_big = diag(nrow(X) + n_test) - X_big_proj%*%solve(t(X_big_proj)%*%X_big_proj)%*%t(X_big_proj)

  if(na){
    XZ = cbind(X,Z); # XZ_test = cbind(X_test,Z_test)
    list_na = gibbs$list_na
    W = XZ[,list_na]; # W_test = XZ_test[,list_na]
    p_na = ncol(W)
    k_na = gibbs$k_na
    W_na = W %>% is.na()
  }

  if(lod){

    if(is.null(X_lod)) stop("provide X_lod")
    X_na = X %>% is.na();
    X[X_lod] = 0
    X_lod[X_na] = F
    lod_vec = gibbs$lod

  }

  for(s in 1:S){

    if(na){

      # sample missing data
      Lambda_na = gibbs$Lambda[s,,]
      eta_na = gibbs$eta[s,,]
      ps_na = gibbs$ps[s,]
      W_pred = eta_na %*% t(Lambda_na) + mvtnorm::rmvnorm(nrow(X), sigma = diag(1/ps_na))
      W[W_na] = W_pred[W_na]
      XZ[,list_na] = W
      X = XZ[,1:p]
      Z = XZ[,(p+1):ncol(XZ)]

      # recompute quantities
      X_big_proj = rbind(X_test,X)
      P_big = diag(nrow(X) + n_test) - X_big_proj%*%solve(t(X_big_proj)%*%X_big_proj)%*%t(X_big_proj)
      X_int = model.matrix(y~.^2 - 1,as.data.frame(X))
      X_int = X_int[,(p+1):ncol(X_int)]

    }

    if(lod){

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
    }

    X_sweep_pred = sweep(X_test, 2, gibbs$l[s], `*`)
    X_sweep = sweep(X, 2, gibbs$l[s], `*`)

    K_big = gibbs$rho[s]*fields::Exp.cov(rbind(X_sweep_pred,X_sweep), p = gibbs$exp_C)
    PKPt_big = P_big%*%K_big%*%t(P_big) + gibbs$sigmasq[s]*diag(N)
    PKPt_12 = PKPt_big[1:n_test,(n_test+1):N]
    Sig_22_inv = solve(PKPt_big[(n_test+1):N,(n_test+1):N])

    mu_pred = Z_test%*%gibbs$beta_z + X_test%*%gibbs$beta[s,] + X_int_test%*%gibbs$lambda[s,]
    mu = Z%*%gibbs$beta_z + X%*%beta + X_int%*%gibbs$lambda

    y_pred[s,] = mu_pred + PKPt_12%*%Sig_22_inv%*%(y-mu)
  }

  return(y_pred)

}
