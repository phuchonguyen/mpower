#' Citation:
#' Daniel Lewandowski, Dorota Kurowicka, Harry Joe,
#' Generating random correlation matrices based on vines and extended onion method,
#' Journal of Multivariate Analysis, Volume 100, Issue 9, 2009, Pages 1989-2001,
#' ISSN 0047-259X, https://doi.org/10.1016/j.jmva.2009.04.008.
#' @param d Number of dimension
#' @param alpha Parameter for Beta distribution
#' @param beta Parameter for Beta distribution
#' @param S A 'guess' of the correlation matrix
#' @param m A number that indicates how much the random matrices vary from S
#' @return A random positive-definite correlation matrix
cvine <- function(d, alpha = 10, beta = 10, S = NULL, m = 100) {
    if (!is.null(S))
        P0 <- cor2partial(S)
    P <- matrix(0, d, d)
    R <- diag(1, d, d)
    for (k in 1:(d - 1)) {
        for (i in (k + 1):d) {
            if (!is.null(S)) {
                alpha <- m * (P0[k, i]/2 + 0.5)
                beta <- m - alpha
            }
            P[k, i] <- stats::rbeta(1, alpha, beta)  # sample partial correlation from Beta Distribution
            P[k, i] <- (P[k, i] - 0.5) * 2  # shift to [-1, 1]
            p <- P[k, i]
            if (k > 1) {
                for (l in seq(k - 1, 1, -1)) {
                  # converting partial correlation to raw correlation
                  p <- p * sqrt((1 - P[l, i]^2) * (1 - P[l, k]^2)) + P[l, i] * P[l,
                    k]
                }
            }
            R[k, i] <- p
            R[i, k] <- p
        }
    }

    return(R)
}

#' Convert a correlation matrix into a partial correlation matrix
#' @param r A correlation matrix
#' @return A partial correlation matrix
cor2partial <- function(r) {
    d <- nrow(r)
    if (d <= 2)
        return(r)
    pcor <- diag(1, d)
    for (k in 2:d) {
        pcor[1, k] <- r[1, k]
        pcor[k, 1] <- r[1, k]
    }
    for (k in 2:(d - 1)) {
        for (i in (k + 1):d) {
            pcor[k, i] <- partial(r, x = c(k, i), y = seq(1, k - 1))[1, 2]
            pcor[i, k] <- pcor[k, i]
        }
    }
    return(pcor)
}

#' Partial correlations between elements in x and elements in y
#' @param r A correlation matrix
#' @param x A vector of indices
#' @param y A vector of indices
#' @return A partial correlation matrix
partial <- function(r, x, y) {
    rr <- r[c(x, y), ][, c(x, y)]
    rx <- 1:length(x)
    ry <- (length(x) + 1):(ncol(rr))
    Cx <- rr[rx, rx] - rr[rx, ry] %*% MASS::ginv(rr[ry, ry]) %*% rr[ry, rx]

    return(stats::cov2cor(Cx))

}
