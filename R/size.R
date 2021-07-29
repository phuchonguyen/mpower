#' @export
estimate_snr <- function(f, xmod, sigma=NULL, m=5000, family="gaussian") {
  X <- genx(xmod, n = m)
  mu <- f(X)
  if (family == "gaussian" & !is.null(sigma)) {
    sigma_signal <- sum((mu - mean(mu))^2)/(m-1)
    return(sigma_signal / sigma)
  } else if (family == "bernoulli") {
    s <- n <- 1
    warning("Not implemented")
  } else if (family == "poisson") {
    s <- n <- 1
    warning("Not implemented")
  } else {
    stop("Family not implemented")
  }
}

#' @export
rsq2snr <- function(x) {
  (1/r - 1)^(-1)
}

#' @export
scale_f <- function(f, rho, sigma, xmod, m=5000, family="gaussian") {
  if (family != "gaussian") stop("Family not implemented")
  X <- genx(xmod, n = m)
  mu <- f(X)
  sigma_signal <- sum((mu - mean(mu))^2)/(m-1)
  s <- rho * (sigma / sigma_signal)
  f_new <- function(...) {sqrt(s) * f(...)}
  return(list(var_signal = sigma_signal, scale_s = sqrt(s)))
}

#' @export
scale_sigma <- function(f, rho, xmod, m=5000, family="gaussian") {
  if (family != "gaussian") stop("Family not implemented")
  X <- genx(xmod, n = m)
  mu <- f(X)
  sigma_signal <- sum((mu - mean(mu))^2)/(m-1)
  return(list(var_signal = sigma_signal, sigma = sqrt(sigma_signal/rho)))
}
