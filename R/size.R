#' @export
estimate_snr <- function(y, mu, family="gaussian") {
  n <- length(y)
  if (family == "gaussian") {
    s <- sum((mu - mean(mu, na.rm = TRUE))^2)/(n-1)
    n <- sum((y - mu)^2)/(n-1)
  } else if (family == "bernoulli") {
    s <- n <- 1
    warning("Not implemented")
  } else if (family == "poisson") {
    s <- n <- 1
    warning("Not implemented")
  } else {
    stop("Family not implemented")
  }
  return(s/n)
}

#' @export
scale_f <- function(f, rho, sigma, xmod, m=5000) {
  X <- genx(xmod, n = m)
  mu <- f(X)
  sigma_signal <- sum((mu - mean(mu))^2)/(m-1)
  s <- rho * (sigma / sigma_signal)
  f_new <- function(...) {sqrt(s) * f(...)}
  return(list(var_signal = sigma_signal, scale_s = sqrt(s)))
}

#' @export
scale_sigma <- function(f, rho, xmod, m=5000) {
  X <- genx(xmod, n = m)
  mu <- f(X)
  sigma_signal <- sum((mu - mean(mu))^2)/(m-1)
  return(list(var_signal = sigma_signal, sigma = sqrt(sigma_signal/rho)))
}
