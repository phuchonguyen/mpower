#' Monte Carlo approximation of the SNR
#' @param ymod A OutcomeModel object
#' @param xmod A MixtureModel object
#' @param m Number of MC samples
#' @return An estimate SNR
#' @export
estimate_snr <- function(ymod, xmod, m = 5000) {
    X <- genx(xmod, n = m)
    mu <- ymod$f(X)
    if (ymod$family == "gaussian" & !is.null(ymod$sigma)) {
        sigma_signal <- get_sigma_signal(mu, m, ymod$sigma)
        return(sigma_signal/ymod$sigma^2)
    } else if (ymod$family == "binomial") {
        return(get_de_snr(mu))
    } else if (ymod$family == "poisson") {
        s <- n <- 1
        warning("Poisson not implemented. Returning NA.")
        return(NA)
    } else {
        stop("Family not implemented.")
    }
}

#' Convert R-squared value to the SNR
#' @param r R-squared value
#' @export
rsq2snr <- function(r) {
    (1/r - 1)^(-1)
}

#' Rescale the mean function of an OutcomeModel to meet a given SNR
#' @param snr A SNR
#' @param ymod A OutcomeModel object to modify
#' @param xmod A MixtureModel object
#' @param m Number of MC samples to estimate the SNR of a proposed noise variance
#' @return A new OutcomeModel object
#' @export
scale_f <- function(snr, ymod, xmod, m = 5000) {
    if (ymod$family != "gaussian")
        stop("Family not implemented")
    X <- genx(xmod, n = m)
    mu <- ymod$f(X)
    sigma_signal <- get_sigma_signal(mu, m, ymod$sigma)
    s <- snr * (ymod$sigma/sigma_signal)
    f_new <- function(...) {
        sqrt(s) * ymod$f(...)
    }
    return(OutcomeModel(f = f_new, family = ymod$family, sigma = ymod$sigma))
}

#' Rescale the noise variance of a Gaussian OutcomeModel to meet a given SNR
#' @param snr A SNR
#' @param ymod A OutcomeModel object to modify
#' @param xmod A MixtureModel object
#' @param m Number of MC samples to estimate the SNR of a proposed noise variance
#' @return A new OutcomeModel object
#' @export
scale_sigma <- function(snr, ymod, xmod, m = 5000) {
    if (ymod$family != "gaussian")
        stop("Family not implemented")
    X <- genx(xmod, n = m)
    mu <- ymod$f(X)
    sigma_signal <- get_sigma_signal(mu, m, ymod$sigma)
    return(OutcomeModel(f = ymod$f, family = ymod$family, sigma = sqrt(sigma_signal/snr)))
}

get_sigma_signal <- function(mu, m, sigma) {
    statistics <- function(mu, idx) {
        sum((mu[idx] - mean(mu[idx]))^2)/(m - 1)
    }
    out <- boot::boot(mu, statistics, R = 100)$t
    message("Estimated SNR is ", round(mean(out/sigma^2), 4), " with bootstrap s.e. ",
        round(stats::sd(out/sigma^2), 4))
    return(mean(out))
}

get_de_snr <- function(mu) {
    statistics <- function(mu, idx) {
        y <- stats::rbinom(length(mu[idx]), size = 1, prob = mu[idx])
        de_noise <- binomial_de(y, mu[idx], 1)
        de_signal <- binomial_de(y, mean(y), 1) - de_noise
        de_signal/de_noise
    }
    out <- boot::boot(mu, statistics, R = 100)$t
    message("Estimated SNR is ", round(mean(out), 4), " with bootstrap s.e. ", round(stats::sd(out),
        4))
    return(mean(out))
}

binomial_de <- function(y, mu, m) {
    d1 <- y * (log(y) - log(mu))
    d1 <- ifelse(is.na(d1), 0, d1)
    d2 <- (m - y) * (log(m - y) - log(m - mu))
    d2 <- ifelse(is.na(d2), 0, d2)
    return(2 * sum(d1 + d2))
}
