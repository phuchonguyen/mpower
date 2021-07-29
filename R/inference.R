#' @importFrom stats coef cov cov2cor model.matrix qnorm qt quantile rbeta rbinom rnorm summary.glm glm

#' @export
InferenceModel <- function(fun, ...) {
  mod <- list("fun_name" = fun, "args" = list(...))
  if (fun == "bma") {
    mod[["fun"]] <- bma_wrapper
  } else if (fun == "bkmr") {
    mod[["fun"]] <- bkmr_wrapper
  } else if (fun == "ms") {
    mod[["fun"]] <- mixselect_wrapper
  } else if (fun == "bws") {
    mod[["fun"]] <- bws_wrapper
  } else if (fun == "qgc") {
    mod[["fun"]] <- qgcomp_lin_wrapper
  } else if (fun == "fin") {
    mod[["fun"]] <- fin_wrapper
  } else if (fun == "glm") {
    mod[["fun"]] <- glm_wrapper
  } else if (is.function(fun)) {
    mod[["fun"]] <- fun
    mod[["fun_name"]] <- "custom"
  } else {
    stop("`fun` not implemented")
  }

  new_InferenceModel(mod)
}

new_InferenceModel <- function(l = list()) {
  stopifnot(is.list(l))
  structure(l, class = "mpower_InferenceModel")
}

#' @export
fit <- function(mod, X, y, alpha) {
  UseMethod("fit")
}

fit.mpower_InferenceModel <- function(mod, X, y, alpha) {
  mod$fun(y = y, X = X, alpha = alpha, args = mod$args)
}

mixselect_wrapper <- function(y, X, Z=NULL, alpha=0.05, args=list()) {
  s <- Sys.time()
  ms_out <- do.call(MixSelect, c(list(y = y, X = X, Z = Z), args))
  ms_time <- Sys.time() - s
  ms_beta <- cbind(ms_out$beta * ms_out$gamma_beta, ms_out$lambda * ms_out$gamma_int)
  ms_beta <- apply(ms_beta, 2, quantile, c(alpha/2, 0.5, 1-alpha/2), na.rm=T)
  ms_beta_pip <- cbind(ms_out$gamma_beta, ms_out$gamma_int)
  ms_beta_pip <- apply(ms_beta_pip, 2, mean, na.rm=T)
  ms_gp <- apply(ms_out$l,
                 2, quantile, c(alpha/2, 0.5, 1-alpha/2), na.rm=T)
  ms_gp_pip <- apply(ms_out$gamma_l, 2, mean, na.rm=T)
  ms_beta_gp_pip <- 1*((ms_out$gamma_beta + ms_out$gamma_l) > 0)
  ms_beta_gp_pip <- apply(ms_beta_gp_pip, 2, mean, na.rm=T)
  return(list(beta = ms_beta,               # CIs for linear terms
              beta_pip = ms_beta_pip,       # PIP for linear terms
              gp = ms_gp,                   # CIs for lengthscale terms
              gp_pip = ms_gp_pip,           # PIPs for lengthscale terms (detection of nonlinear terms)
              pip = ms_pip,                 # detection of effects
              time = ms_time))
}

bkmr_wrapper <- function(y, X, Z=NULL, alpha=0.05, args=list()) {
  s <- Sys.time()
  km_fit <- do.call(bkmr::kmbayes, c(list(y=y, Z=X, X=Z), args))
  km_time <- Sys.time() - s
  km_pip <- bkmr::ExtractPIPs(km_fit) # by default keeps the second half of all samples
  km_pip <- km_pip[, (ncol(km_pip)/2 + 1):ncol(km_pip)]
  pip <- km_pip$PIP
  if (is.null(pip)) pip <- km_pip$condPIP
  return(list(pip = pip,
              group_pip = km_pip$groupPIP,
              time = km_time))  # detect of effect
}

bma_wrapper <- function(y, X, Z=NULL, interact=FALSE, alpha=0.05, args=list()) {
  if (!interact) {
    X_int <- X
  } else {
    X_int <- model.matrix(~ .^2 - 1, as.data.frame(X))
  }
  if (is.null(Z)) {
    X_full <- X_int
  } else {
    X_full <- cbind(X_int, Z)
  }
  s <- Sys.time()
  bma_out <- do.call(BMA::bic.glm, c(list(x=X_full, y=y), args))
  bma_time <- Sys.time() - s
  # Model-averaged coef, the first one is the intercept term
  p <- ncol(X_int)
  bma_beta <- rbind(bma_out$postmean[2:(p+1)] - qnorm(1-alpha/2)*bma_out$postsd[2:(p+1)],
                   bma_out$postmean[2:(p+1)],
                   bma_out$postmean[2:(p+1)] + qnorm(1-alpha/2)*bma_out$postsd[2:(p+1)])
  # Significance of a variable is when either the CI of a main effect or interaction is significant
  bma_beta_sig <- 1*(bma_beta[1,]>0 | bma_beta[3,]<0)
  bma_sig <- rep(NA, ncol(X))
  for (i in 1:ncol(X)) {
    id <- which(grepl(i, colnames(X_int)))
    bma_sig[i] <- 1*(sum(bma_beta_sig[id]) > 0)
  }
  # Sum inclusion probability postprob of all models with a main effect of interaction of Xi
  bma_pip <- rep(NA, ncol(X))
  W <- bma_out$which
  p <- bma_out$postprob
  m <- nrow(W)
  for (i in 1:ncol(X)) {
    id <- which(grepl(i, colnames(X_full)))
    bma_pip[i] <- sum(sapply(1:m, function(j, W, p, id) {
      r <- W[j,]
      if (sum(r[id]) > 0) {p[j]}
      else {0}
    }, W, p, id))
  }

  return(list(beta = bma_beta,         # CIs for main effect and interactions
              ci = bma_sig,            # CIs significant for either main effect or interactions
              pip = bma_pip,           # PIPs of either main effect or interactions
              time = bma_time))
}

fin_wrapper <- function(y, X, Z=NULL, alpha=0.05, args=list()) {
  y <- matrix(y, length(y), 1)
  if (!is.null(Z)) warning("Inclusion of covariates not implemented in package infinitefactor for FIN. Ignoring X from regression")
  s <- Sys.time()
  fin_out <- do.call(infinitefactor::interactionDL, c(list(y=y, X=X), args))
  fin_time <- Sys.time() - s
  fin_beta <- apply(fin_out$mainEffectSamps, 1, quantile, c(alpha/2, 0.5, 1-alpha/2), na.rm=T)
  fin_int <-  apply(apply(fin_out$interactionSamps, 3, c),
                    1, quantile, c(alpha/2, 0.5, 1-alpha/2), na.rm=T)
  p <- ncol(fin_beta)
  colnames(fin_int) <- paste0(rep(1:p, each=p), rep(1:p, p))
  fin_int <- fin_int[, c(paste0(1:p, 1:p), combn(1:p, m=2, FUN=paste0, collapse=""))]
  # significance of a variable is when either the CI of a main effect or interaction is significant
  fin_beta_sig <- apply(fin_beta, 2, function(x) x[1] > 0 | x[3] < 0)
  fin_int_sig <- apply(fin_int, 2, function(x) x[1] > 0 | x[3] < 0)
  fin_sig <- rep(NA, p)
  for (i in 1:p) {
    id <- which(grepl(i, colnames(fin_int)))
    fin_sig[i] <- 1*(sum(fin_beta_sig[i], fin_int_sig[id]) > 0)
  }
  return(list(beta = fin_beta,        # CIs main effects
              beta_int = fin_int,     # CIs interactions
              ci = fin_sig,           # detection of effects
              time = fin_time))
}

qgcomp_lin_wrapper <- function(y, X, Z=NULL, alpha=0.05, args=list()) {
  if (is.null(Z)) {
    dat <- data.frame(y=y, x=X)
    colnames(dat) <- c("y", colnames(X))
  } else {
    dat <- data.frame(y=y, x=X, z=Z)
    colnames(dat) <- c("y", colnames(X), colnames(Z))
  }
  s <- Sys.time()
  qc.fit <- do.call(qgcomp::qgcomp.noboot, c(list(y ~ . , data = dat, expnms = colnames(X)), args))
  qgc_time <- Sys.time() - s
  p <- ncol(X)
  qgc_sig <- rep(qc.fit$pval[2] <= alpha, p)  # p-value of psi1
  # Get CI for each estimate, only available for linear qgcomp
  qc.fit$qx$y = qc.fit$fit$data$y
  newfit <- do.call(glm, c(list(y ~ ., data=dat), args))
  me <- coef(newfit)[2:(1+p)]
  sde <- summary.glm(newfit)$coef[,2][2:(1+p)]
  qgc_beta <- rbind(me - qt(1-alpha/2, nrow(X)-p)*sde,
                    me,
                    me + qt(1-alpha/2, nrow(X)-p)*sde)
  rownames(qgc_beta) <- c(paste0(alpha/2, "%"), "50.0%", paste0(1-alpha/2, "%"))
  return(list(ci = qgc_sig,    # detection  of effect
              pos_psi = qc.fit$pos.psi,
              neg_psi = qc.fit$neg.psi,
              pos_weights = qc.fit$pos.weights,
              neg_weights = qc.fit$neg.weights,
              beta = qgc_beta,  # glm CI
              time = qgc_time))
}

glm_wrapper <- function(y, X, Z=NULL, alpha=0.05, args=list()) {
  if (is.null(Z)) {
    dat <- data.frame(y=y, x=X) %>%
      set_colnames(c("y", colnames(X)))
  } else {
    dat <- data.frame(y=y, x=X, z=Z) %>%
      set_colnames(c("y", colnames(X), colnames(Z)))
  }
  p <- ncol(X)
  s <- Sys.time()
  fit <- do.call(glm, c(list(y ~ ., data=dat), args))
  time <- Sys.time() - s
  me <- coef(fit)[2:(1+p)]
  sde <- summary.lm(fit)$coef[,2][2:(1+p)]
  beta <- rbind(me - qt(1-alpha/2, nrow(X)-p)*sde, me,
                me + qt(1-alpha/2, nrow(X)-p)*sde)
  rownames(beta) <- c(paste0(alpha/2, "%"), "50.0%", paste0(1-alpha/2, "%"))
  return(list(ci = 1*(beta[1,] > 0 | beta[3,] < 0),
              beta = beta,
              time = time))
}

# TODO: suppress rstan log messages
bws_wrapper <- function(y, X, Z=NULL, alpha=0.05, args) {
  s <- Sys.time()
  fit <- do.call(bws::bws, c(list(y = y, X = X, Z = Z), args))
  bws_time <- Sys.time() - s
  samps <- rstan::extract(fit)
  quantile_beta <- quantile(samps$theta1, c(alpha/2, 0.5, 1-alpha/2))
  all_pos <- sum(quantile_beta > 0) == 3
  all_neg <- sum(quantile_beta < 0) == 3
  bws_sig <- (all_pos | all_neg) * 1
  bws_beta <- apply(sweep(samps$w, 1, samps$theta1, "*"), 2,
                    quantile, c(alpha/2, 0.5, 1-alpha/2))
  bws_weights <- apply(samps$w, 2, quantile, c(alpha/2, 0.5, 1-alpha/2))
  return(list(ci = rep(bws_sig, ncol(X)),
              beta = bws_beta,
              weights = bws_weights,
              time = bws_time))
}




