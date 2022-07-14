#' Gibbs sampler for MixSelect regression, based on a decomposition of the regression
#'           surface in:
#'           - main effects
#'           - interactions with heredity constraint
#'           - nonlinear deviation distributed as projected GP
#' ArXiv: https://arxiv.org/abs/1911.01910
#' Email: ff31@duke.edu for questions
#' @param y response vector;
#' @param X predictor matrix (n x p) to be included in projected GP;
#' @param Z covariate adjustments matrix (n x q);
#' @param nrun number of iterations;
#' @param burn burn-in period;
#' @param thin thinning interval; lambda_prior = "cont_spike",
#' @param tau variance hyperparameter for continuous spike and slab on main effects and interactions
#' @param c shrinkage parameter for continuous spike and slab on main effects and interactions
#' @param exp_C exponent of exponential covariance, default is  (2.2) in paper
#' @param heredity "strong" or "weak"
#' @param lambda_prior "cont_spike" or "discrete_spike", the  latter is extremely slow and  not  suggested
#' @param na T if analysis with missing data
#' @param list_na vector of dimension p+q indicating which variables in X and Z need to be included in  Factor Model,
#'                   this is because maybe we don't want to use a factor model for all the covariates, or maybe there are some variables
#'                   that are not missing
#'                   If null, use all columns in X and Z
#' @param k_na number of latent factors for analysis with NAs, refer to Section 5.3 in paper
#'                  Defaule is 2*log(p) where p is number of variables included in the Factor Model
#' @param sd_lambda hyperparamenter for element of Lambda matrix
#' @param lod T if analysis with data below the Limit of Detection (LOD)
#' @param X_lod Matrix of dimension nxp, 1 if observation i,j has been observed  below LOD and 1  otherwise
#' @param lod_vec vector or dimension  px1 containing  the LOD for each  variable  in  X
#' @param verbose logical
#' @importFrom stats dgamma dnorm model.matrix rbeta rbinom rgamma rnorm runif
#' @importFrom utils setTxtProgressBar txtProgressBar

MixSelect = function(y, X , Z = NULL,
                     nrun = 2000, burn = nrun/2, thin = 1,
                     tau = 1, c = 0.1, exp_C = 2,
                     heredity = "strong", lambda_prior = "cont_spike",
                     na = F, k_na = NULL, list_na = NULL, sd_lambda = 0.1,
                     lod = F, X_lod = NULL, lod_vec = NULL,
                     verbose = T){
  # --- default Z --- #
  if(is.null(Z)){
    Z = matrix(0,nrow = nrow(X), ncol = 1);
  }

  # --- verbose --- #
  at = ceiling(nrun/100)
  if(verbose) {
    pb = txtProgressBar(style = 3)
  }

  # --- initialize parameters --- #
  out_params = init_params(y = y, X = X, Z = Z,
                            nrun = nrun, burn = burn, thin = thin,
                            tau = tau, c = c, exp_C = exp_C,
                            heredity = heredity, lambda_prior = lambda_prior,
                            na = na, k_na = k_na, list_na = list_na, sd_lambda = sd_lambda,
                            lod = lod, X_lod = X_lod, lod_vec = lod_vec)
  const = out_params$constants
  curr_params = out_params$params

  # --- storage --- #
  p = const$p
  n = const$n
  S = floor((nrun - burn)/thin)
  beta_st = gamma_st = matrix(NA, nrow = S, ncol = p)
  beta_z_st = gamma_z_st = matrix(NA, nrow = S, ncol = ncol(Z))
  lambda_st = gamma_int_st = matrix(NA, nrow = S, ncol = ncol(curr_params$X_int))
  l_st = gamma_l_st = matrix(NA, nrow = S, ncol = p)
  rho_st = sigma_sq_st = numeric(S)
  y_hat = y_hat2 = y_hat3 =  matrix(NA, nrow = S, ncol = n)
  mu_st = matrix(NA, nrow = S, ncol = nrow(X))
  Omega_st = array(0,c(S, p, p))
  Omega_01 = array(0,c(S, p, p))
  Lambda_st = array(0,c(S, p, k_na))

  count = 1
  for (s in 1:nrun){

    curr_params = update_na_lod(params = curr_params, constants = const)
    curr_params = update_gamma_beta(params = curr_params, constants = const)
    curr_params = update_gamma_beta_z(params = curr_params, constants = const)
    curr_params = update_int(params = curr_params, constants = const)
    curr_params[["mu"]] = get_mu(params = curr_params, constants = const)
    curr_params = update_l(params = curr_params, constants = const)
    curr_params = update_rho(params = curr_params, constants = const)
    curr_params = update_sigmasq(params = curr_params, constants = const)

    if(verbose & (s %% at == 0)) setTxtProgressBar(pb, s / nrun)

    if(s > burn & s%%thin == 0){

      # --- storage of coeffs ---
      beta_st[count,] = curr_params$beta
      gamma_st[count,] = curr_params$gamma
      beta_z_st[count,] = curr_params$beta_z
      gamma_z_st[count,] = curr_params$gamma_z
      lambda_st[count,] = curr_params$lambda
      gamma_int_st[count,] = curr_params$gamma_int
      rho_st[count] = curr_params$rho
      sigma_sq_st[count] = curr_params$sigmasq
      l_st[count,] = curr_params$l
      gamma_l_st[count,] = curr_params$gamma_l

      Omega_curr = matrix(0,p,p)
      Omega_curr[lower.tri(Omega_curr)] = curr_params$lambda/2
      Omega_st[count,,] = Omega_curr + t(Omega_curr)

      Omega_curr = matrix(0,p,p)
      Omega_curr[lower.tri(Omega_curr)] = curr_params$gamma_int
      Omega_01[count,,] = Omega_curr + t(Omega_curr)

      y_hat[count,] = predict_insample(params = curr_params, constants = const)

      count = count + 1
    }
  }

  ret_list = list(beta = beta_st,
                  gamma_beta = gamma_st,
                  beta_z = beta_z_st,
                  gamma_z = gamma_z_st,
                  lambda = lambda_st,
                  gamma_int = gamma_int_st,
                  sigmasq = sigma_sq_st,
                  rho = rho_st,
                  exp_C = exp_C,
                  l = l_st,
                  gamma_l = gamma_l_st,
                  Omega = Omega_st,
                  Omega_01 = Omega_01,
                  y_hat = y_hat,
                  X = X)

  if(na){
    ret_list[["Lambda"]] = Lambda_st
  }

  return(ret_list)
}
