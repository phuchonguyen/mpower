#' @export
estimate_snr <- function(f, xmod, sigma=NULL, m=5000, family="gaussian") {
  X <- genx(xmod, n = m)
  mu <- f(X)
  if (family == "gaussian" & !is.null(sigma)) {
    sigma_signal <- sum((mu - mean(mu))^2)/(m-1)
    return(sigma_signal / sigma)
  } else if (family == "binomial") {
    y <- rbinom(length(mu), size = 1, prob = mu)
    de_noise <- binomial_de(y, mu, 1)
    de_signal <- binomial_de(y, mean(y), 1) - de_noise
    return(de_signal / de_noise)
  } else if (family == "poisson") {
    s <- n <- 1
    warning("Not implemented")
  } else {
    stop("Family not implemented")
  }
}

binomial_de <- function(y, mu, m) {
  d1 <- y*(log(y) - log(mu))
  d1 <- ifelse(is.na(d1), 0, d1)
  d2 <- (m-y)*(log(m-y) - log(m-mu))
  d2 <- ifelse(is.na(d2), 0, d2)
  return(2*sum(d1+d2))
}

#' @export
rsq2snr <- function(r) {
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
