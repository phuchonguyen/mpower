#' Correlated predictors generator
#'
#' This function creates a generative model for the correlated, mixed-scale
#' predictors.
#'
#' @param method A string, one of the three options "resampling", "estimation",
#'   or "cvine". Default is "estimation". See Details.
#' @param data A dataframe or matrix, required for "resampling" and "estimation"
#'   method.
#' @param G A guesstimate pairwise correlation matrix for "cvine" method. See
#'   Details.
#' @param m A positive number indicating uncertainty in the guesstimate G,
#'   larger means more uncertainty. Default is 100.
#' @param sbg_args A list of named arguments for function [sbgcop.mcmc()].
#' @param cvine_marginals A list describing the univariate distribution of each predictor.
#' @param resamp_prob A vector of sampling probability for each observation in data. Must sums to 1.
#' @param nudge A number, default 10e-10 to add to the diagonal of the covariance
#'   matrix for numerical stability.
#'
#' @section Details:
#'
#'   There are three methods to generate data:
#'
#'   1. Resampling: if we have enough data of the predictors, we can
#'   resample to get realistic joint distributions and dependence amongst them.
#'
#'   2. Estimation: if we have a small sample from, for example, a pilot study,
#'   we can sample from a semiparametric copula model (Hoff 2007) after learning
#'   the dependence and univariate marginals of the predictors.
#'
#'   3. C-vine: if no pilot data exists, we can still set rough guesstimate of
#'   the dependence and univariate marginals. The C-vine algorithm (Joe 2006)
#'   generates positive semi-definite correlation matrix given the guesstimate G.
#'   The guesstimate G is a symmetric p x p matrix whose ij-th item is between
#'   -1 and 1 and is the guesstimate correlation between predictor ith and jth.
#'   G doesn't need to be a valid correlation matrix. The method works well when
#'   values in G are not extreme (i.e., 0.999, -0.999).
#'
#' @section References:
#'
#'   Hoff P (2007). "Extending the rank likelihood for
#'   semiparametric copula estimation." Ann. Appl. Stat, 1(1), 265-283.
#'
#'   Joe H (2006). “Generting random correlation matrices based on partial
#'   correlations.”Journal of Multivariate Analysis, 97, 2177-2189.
#' @export
MixtureModel <- function(method = "estimation", data = NULL, G = NULL, m = 100,
                         nudge = 10e-10, sbg_args = list(nsamp = 1000),
                         cvine_marginals = list(), resamp_prob = NULL) {
  mod <- NULL
  if (nudge < 0) stop("nudge must be positive")
  if (method == "estimation" | method == "resampling") {
    if (!is.matrix(data) & !is.data.frame(data)) {
      stop("Input must be a matrix or a data.frame",
           call. = FALSE)
    }
    var_name <- colnames(data)
    if (is.null(var_name)) var_name <- paste0("x", 1:ncol(data))
  }
  if (method == "estimation") {
    mod <- new_estimation_MixtureModel(data = data, nudge = nudge,
                                       var_name = var_name, args = sbg_args)
  } else if (method == "resampling") {
    mod <- new_resampling_MixtureModel(data = data, var_name = var_name,
                                       prob = resamp_prob)
  } else if (method == "cvine") {
    var_name <- colnames(G)
    if (is.null(var_name)) var_name <- paste0("x", 1:ncol(G))
    mod <- new_cvine_MixtureModel(G = G, m = m, var_name = var_name, p = ncol(G),
                                  marginals = cvine_marginals)
  } else {stop("Method not implemented", call. = FALSE)}

  validate_MixtureModel(mod)
}

validate_MixtureModel <- function(x) {
  stopifnot(is.character(x$var_name))
  stopifnot(is.numeric(x$p))
  x
}

new_resampling_MixtureModel <- function(data = data.frame(),
                                        var_name = character(), prob = numeric()) {
  if(!is.null(prob)) stopifnot(sum(prob)==1 & all(prob)>0)
  x <- list(data = data, var_name = var_name, p = ncol(data), prob = prob)
  structure(x, class = "mpower_resampling_MixtureModel")
}

new_estimation_MixtureModel <- function(data = numeric(), nudge = numeric(),
                                        var_name = character(), args = list()) {
  if(!is.numeric(as.matrix(data))) {
    stop("Data must be numeric only. Convert categories to integer but
         do not use One-hot-encoding.")
  }
  cat("\nEstimating the joint distribution using sbgcop\n")
  sbgcop_fit <- do.call(sbgcop::sbgcop.mcmc, c(list(Y = data), args))
  cat("\n*** MCMC Summary ***")
  cat("\nNumber of samples", sbgcop::summary.psgc(sbgcop_fit)$nsamp)
  cat("\nEffective sample sizes")
  print(sbgcop::summary.psgc(sbgcop_fit)$ESS)
  cat("\nSummary plots for univariate marginals, pair-wise correlation and regression parameters.")
  sbgcop::plot.psgc(sbgcop_fit)
  mtext("Univariate marignals (left), pairwise correlation (middle), pairwise regression parameter (right)",
        outer = T, cex = 0.7)
  R <- sbgcop_fit$C.psamp
  x <- list(data = data, R = R, var_name = var_name, p = ncol(data), nudge = nudge)
  structure(x, class = "mpower_estimation_MixtureModel")
}

new_cvine_MixtureModel <- function(G = numeric(), m = numeric(),
                                   var_name = character(),
                                   p = numeric(), nudge = numeric(),
                                   marginals = marginals) {
  if (!all(G <= 1 & G >= -1)) {
    stop("Values in G must be between -1 and 1", call. = FALSE)
  }
  if (ncol(G) != nrow(G)) stop("G is not symmetric", call. = FALSE)
  if (!is.numeric(m)) stop("Input `m` not numeric", call. = FALSE)
  if (m <= 0) {stop("Input `m` must be positive", call. = FALSE)}
  x <- list(G = G, m = m, var_name = var_name, nudge = nudge, p = p,
            marginals = marginals)
  structure(x, class = "mpower_cvine_MixtureModel")
}

#' Generates a matrix of n observations of p predictors
#' @param obj A MixtureModel object
#' @param n An integer, number of observations to generate
#' @return A (n x p) dataframe
#' @export
genx <- function(obj, n) {
  UseMethod("genx")
}

genx.mpower_resampling_MixtureModel <- function(obj, n) {
  stopifnot(n >= 0, n %% 1 == 0)
  if (n > nrow(obj$data)) {
    warning("Number of observations ", nrow(obj$data), " less than ", n, ". Samples are dependent.")
  }
  indices <- sample(1:nrow(obj$data), n, replace = TRUE, prob = obj$prob)
  obj$data[indices, ] %>%
    magrittr::set_colnames(obj$var_name)
}

genx.mpower_estimation_MixtureModel <- function(obj, n) {
  stopifnot(n >= 0, n %% 1 == 0)
  R <- obj$R[,,sample(dim(obj$R)[3], 1)]
  while (!is_positive_definite(R)) {
    R <- R + diag(obj$nudge, obj$p)
    warning("Adding a nudge to the diagonal because R is not positive definite.")
  }
  Z <- MASS::mvrnorm(n, mu = rep(0, obj$p), Sigma = R)
  X <- vapply(1:obj$p, function(j) {
    quantile(obj$data[,j],
             probs = pnorm(Z[,j], 0, sqrt(R[j,j])), na.rm=TRUE, type=1)},
    numeric(n))
  rownames(X) <- NULL
  colnames(X) <- obj$var_name
  return(as.data.frame(X))
}

# TODO: implement an qmultinom function that takes p=probabilities
genx.mpower_cvine_MixtureModel <- function(obj, n) {
  stopifnot(n >= 0, n %% 1 == 0)
  R <- cvine(d = obj$p, S = obj$G, m = obj$m)
  while (!is_positive_definite(R)) {
    R <- R + diag(obj$nudge, obj$p)
    warning("Adding a nudge to the diagonal because R is not positive definite.")
  }
  Z <- MASS::mvrnorm(n, mu = rep(0, obj$p), Sigma = R)
  X <- vapply(1:obj$p, function(j) {
    namej <- obj$var_name[j]
    p <- pnorm(Z[,j], 0, sqrt(R[j,j]))
    text <- gsub(")$", ", p=p)", obj$marginals[[namej]])
    eval(parse(text = text))
    },
    numeric(n)) %>%
    magrittr::set_colnames(obj$var_name)
  return(X)
}

