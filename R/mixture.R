#' Correlated predictors generator
#'
#' This function creates a generative model for the correlated, mixed-scale
#' predictors.
#'
#' @param method A string, one of the three options "resampling", "estimation",
#'   or "cvine". Default is "estimation". See Details.
#' @param data A dataframe or matrix, required for resampling" and "estimation"
#'   method.
#' @param G A guesstimate pairwise correlation matrix for "cvine" method. See
#'   Details.
#' @param m A positive number indicating uncertainty in the guesstimate G,
#'   larger means more uncertainty. Default is 100.
#' @param sbg_args A list of named arguments, except Y, for function `sbgcop.mcmc()`.
#' @param cvine_marginals A named list describing the univariate distribution
#'   of each predictor. See Details.
#' @param cvine_dtypes A named list describing the data type of each variable.
#' @param resamp_prob A vector of sampling probability for each observation in data. Must sums to 1.
#' @param nudge A number, default 10e-10 to add to the diagonal of the covariance
#'   matrix for numerical stability.
#'
#' @section Details:
#'
#'   There are three methods to generate data:
#'
#'   1. Resampling: if we have enough data of the predictors, we can
#'   resample to get realistic joint distributions and dependence among them.
#'
#'   2. Estimation: if we have a small sample from, for example, a pilot study,
#'   we can sample from a semi-parametric copula model (Hoff 2007) after learning
#'   the dependence and univariate marginals of the predictors.
#'
#'   3. C-vine: if no pilot data exists, we can still set rough guesstimate of
#'   the dependence and univariate marginals. The C-vine algorithm (Joe 2006)
#'   generates positive semi-definite correlation matrix given the guesstimate G.
#'   The guesstimate G is a symmetric p x p matrix whose ij-th item is between
#'   -1 and 1 and is the guesstimate correlation between predictor ith and jth.
#'   G doesn't need to be a valid correlation matrix. The method works well when
#'   values in G are not extreme (i.e., 0.999, -0.999). Built-in functions for
#'   univariate marginals include: `qbeta` , `qbinom`, `qcauchy`, `qchisq`, `qexp`, `qf`,
#'   `qgamma`, `qgeom`, `qhyper`, `qlogis`, `qlnorm`, `qmultinom`, `qnbinom`,
#'   `qnorm`, `qpois`, `qt`, `qunif`, `qweibull`.
#' @return A MixtureModel object.
#' @examples
#' data("nhanes1518")
#' xmod <- mpower::MixtureModel(nhanes1518, method = "resampling")
#'
#' @section References:
#'
#'   Hoff P (2007). 'Extending the rank likelihood for
#'   semiparametric copula estimation.' Ann. Appl. Stat, 1(1), 265-283.
#'
#'   Joe H (2006). “Generating random correlation matrices based on partial
#'   correlations.”Journal of Multivariate Analysis, 97, 2177-2189.
#' @export
MixtureModel <- function(method = "estimation", data = NULL, G = NULL, m = 100, nudge = 1e-09,
    sbg_args = list(nsamp = 1000), cvine_marginals = list(), cvine_dtypes = list(),
    resamp_prob = NULL) {
    mod <- NULL
    if (nudge < 0)
        stop("nudge must be positive")
    if (method == "estimation" | method == "resampling") {
        if (!is.matrix(data) & !is.data.frame(data)) {
            stop("Input must be a matrix or a data.frame", call. = FALSE)
        }
        var_name <- colnames(data)
        if (is.null(var_name))
            var_name <- paste0("x", 1:ncol(data))
    }
    if (method == "estimation") {
        sbg_args[["Y"]] <- NULL
        var_dtypes <- sapply(data, class)
        if ("character" %in% var_dtypes) {
            stop("Data must be numeric or factor with numeric categories", call. = FALSE)
        }
        if ("factor" %in% var_dtypes) {
            warning("Categorical data (factors) must have numeric categories.")
        }
        mod <- new_estimation_MixtureModel(data = data, nudge = nudge, var_name = var_name,
            var_dtypes = var_dtypes, args = sbg_args)
    } else if (method == "resampling") {
        mod <- new_resampling_MixtureModel(data = data, var_name = var_name, prob = resamp_prob)
    } else if (method == "cvine") {
        var_name <- names(cvine_marginals)
        if (is.null(var_name))
            var_name <- paste0("x", 1:ncol(G))
        mod <- new_cvine_MixtureModel(G = G, m = m, var_name = var_name, p = ncol(G),
            marginals = cvine_marginals, dtypes = cvine_dtypes)
    } else {
        stop("Method not implemented", call. = FALSE)
    }

    validate_MixtureModel(mod)
}

validate_MixtureModel <- function(x) {
    stopifnot(is.character(x$var_name))
    stopifnot(is.numeric(x$p))
    x
}

new_resampling_MixtureModel <- function(data = data.frame(), var_name = character(),
    prob = numeric()) {
    if (!is.null(prob))
        stopifnot(sum(prob) == 1 & all(prob) > 0)
    x <- list(data = data, var_name = var_name, p = ncol(data), prob = prob)
    structure(x, class = "mpower_resampling_MixtureModel")
}

new_estimation_MixtureModel <- function(data = numeric(), nudge = numeric(), var_name = character(),
    var_dtypes = character(), args = list()) {
    if (!is.numeric(as.matrix(data))) {
        stop("Data must be numeric only. Convert ordinal factors to integer or
    use One-hot-encoding.")
    }
    message("\nEstimating the joint distribution using sbgcop\n")
    sbgcop_fit <- do.call(sbgcop::sbgcop.mcmc, c(list(Y = data), args))
    R <- sbgcop_fit$C.psamp
    x <- list(data = data, R = R, var_name = var_name, var_dtypes = var_dtypes, p = ncol(data),
        nudge = nudge, sbgcop_summary = sbgcop::summary.psgc(sbgcop_fit))
    structure(x, class = "mpower_estimation_MixtureModel")
}

new_cvine_MixtureModel <- function(G = numeric(), m = numeric(), var_name = character(),
    p = numeric(), nudge = numeric(), marginals = list(), dtypes = list()) {
    if (!all(G <= 1 & G >= -1)) {
        stop("Values in G must be between -1 and 1", call. = FALSE)
    }
    if (ncol(G) != nrow(G))
        stop("G is not symmetric", call. = FALSE)
    if (!is.numeric(m))
        stop("Input `m` not numeric", call. = FALSE)
    if (m <= 0) {
        stop("Input `m` must be positive", call. = FALSE)
    }
    x <- list(G = G, m = m, var_name = var_name, nudge = nudge, p = p, marginals = marginals,
        dtypes = dtypes)
    structure(x, class = "mpower_cvine_MixtureModel")
}

#' Generates a matrix of n observations of p predictors
#' @param obj A MixtureModel object.
#' @param n An integer, number of observations to generate.
#' @return A (n x p) dataframe.
#' @export
genx <- function(obj, n) {
    UseMethod("genx")
}

#' @export
genx.mpower_resampling_MixtureModel <- function(obj, n) {
    stopifnot(n >= 0, n%%1 == 0)
    if (n > nrow(obj$data)) {
        warning("Number of observations ", nrow(obj$data), " less than ", n, ". Samples are dependent.")
    }
    indices <- sample(1:nrow(obj$data), n, replace = TRUE, prob = obj$prob)
    obj$data[indices, ] %>%
        magrittr::set_colnames(obj$var_name)
}

#' @export
genx.mpower_estimation_MixtureModel <- function(obj, n) {
    stopifnot(n >= 0, n%%1 == 0)
    R <- obj$R[, , sample(dim(obj$R)[3], 1)]
    while (!is_positive_definite(R)) {
        R <- R + diag(obj$nudge, obj$p)
        warning("Adding a nudge to the diagonal because R is not positive definite.")
    }
    Z <- MASS::mvrnorm(n, mu = rep(0, obj$p), Sigma = R)
    X <- vapply(1:obj$p, function(j) {
        stats::quantile(obj$data[, j], probs = stats::pnorm(Z[, j], 0, sqrt(R[j,
            j])), na.rm = TRUE, type = 1)
    }, numeric(n))
    rownames(X) <- NULL
    colnames(X) <- obj$var_name
    for (j in 1:obj$p) {
        if (obj$var_dtypes[j] == "factor")
            X[, j] <- as.factor(X[, j])
    }
    return(as.data.frame(X))
}

#' @export
genx.mpower_cvine_MixtureModel <- function(obj, n) {
    stopifnot(n >= 0, n%%1 == 0)
    R <- cvine(d = obj$p, S = obj$G, m = obj$m)
    while (!is_positive_definite(R)) {
        R <- R + diag(obj$nudge, obj$p)
        warning("Adding a nudge to the diagonal because R is not positive definite.")
    }
    Z <- MASS::mvrnorm(n, mu = rep(0, obj$p), Sigma = R)
    X <- vapply(1:obj$p, function(j) {
        namej <- obj$var_name[j]
        p <- stats::pnorm(Z[, j], 0, sqrt(R[j, j]))
        text <- gsub(")$", ", p=p)", obj$marginals[[namej]])
        eval(parse(text = text))
    }, numeric(n)) %>%
        as.data.frame() %>%
        magrittr::set_colnames(obj$var_name)
    for (namej in names(obj$dtypes)) {
        if (obj$dtypes[[namej]] == "factor") {
            X[namej] <- as.factor(X[[namej]])
        }
    }
    return(X)
}

#' Quantile function for the multinomial distribution, size = 1
#' @param p A quantile.
#' @param probs A vector of probabilities for each level.
#' @return Gives the quantile function
#' @export
qmultinom <- function(p, probs) {
    stopifnot(all(probs < 1) & all(probs > 0))
    stopifnot(abs(sum(probs) - 1) < 10^-8)
    cum_probs <- cumsum(probs)
    vapply(p, function(x) min(which(x <= cum_probs)), numeric(1))
}

#' Visualize marginals and Gaussian copula correlations of simulated data
#' @param obj A MixtureModel object.
#' @param split A logical, whether to display numbers on half of the covariance
#' matrix.
#' @return A 'ggplot2' graphics.
#' @export
mplot <- function(obj, split = TRUE) {
    UseMethod("mplot")
}

#' @export
mplot.mpower_resampling_MixtureModel <- function(obj, split = TRUE) {
    g1 <- plot_marginals(obj$data)
    if (!is.numeric(as.matrix(obj$data))) {
        warning("Data must be numeric only to plot correlation matrix. Convert ordinal factors to integer or
    use One-hot-encoding.")
    } else{
        g2 <- stats::cor(obj$data, method = "spearman") %>%
            reshape2::melt() %>%
            ggplot(aes(!!sym("Var2"), !!sym("Var1"), fill = !!sym("value"))) +
            geom_tile() +
            scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#d55E00",
                limit = c(-1, 1), name = "Spearman\nCorrelation") +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(angle = 90),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank()) +
            coord_fixed() + labs(x = "", y = "", title = "Spearman correlation matrix of original data")
        return(list(hist = g1, corr = g2))
    }
    return(list(hist = g1))
}

#' @export
mplot.mpower_estimation_MixtureModel <- function(obj, split = TRUE) {
    g1 <- plot_marginals(obj$data)
    R_mean <- apply(obj$R, c(1, 2), mean)
    g2 <- plot_corr(R_mean, title = "Gaussian copula correlation matrix", split = split)
    return(list(hist = g1, corr = g2))
}

#' @export
mplot.mpower_cvine_MixtureModel <- function(obj, split = TRUE) {
    n <- 100
    data <- genx(obj, n)
    g1 <- plot_marginals(data)
    R <- cvine(d = obj$p, S = obj$G, m = obj$m)
    g2 <- plot_corr(R, title = "An example correlation matrix by C-vine", split = split)
    return(list(hist = g1, corr = g2))
}

plot_marginals <- function(data) {
    g <- list()
    nums <- data %>%
        dplyr::select_if(is.numeric)
    if (ncol(nums) > 0) {
        nums <- nums %>%
            tidyr::pivot_longer(tidyr::everything(), values_to = "value", names_to = "name")
        temp <- list(ggplot(nums, aes(x = !!sym("value"))) + geom_histogram() + facet_wrap(stats::as.formula("~name")) +
            labs(title = "Univariate distributions"))
        g <- c(g, temp)
    }
    cat <- data %>%
        dplyr::select_if(purrr::negate(is.numeric))
    if (ncol(cat) > 0) {
        cat <- cat %>%
            tidyr::pivot_longer(tidyr::everything(), values_to = "value", names_to = "name")
        temp <- list(ggplot(cat, aes(x = !!sym("value"))) + geom_bar() + facet_wrap(stats::as.formula("~name")) +
            labs(title = "Univariate distributions") + theme(axis.text.x = element_text(angle = 90)))
        g <- c(g, temp)
    }
    return(g)
}

plot_corr <- function(C, title = "", split = TRUE) {
    if (!split) {
        g <- reshape2::melt(C) %>%
            ggplot(aes(!!sym("Var2"), !!sym("Var1"), fill = !!sym("value"))) + geom_tile() +
            scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#d55E00",
                limit = c(-1, 1), name = "Pearson\nCorrelation") + theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), panel.background = element_blank()) +
            coord_fixed() + labs(x = "", y = "", title = title)
        return(g)
    }
    upper_C <- lower_C <- C
    upper_C[lower.tri(C)] <- NA
    lower_C[upper.tri(C)] <- NA
    meltNum <- reshape2::melt(lower_C, na.rm = T)
    meltColor <- reshape2::melt(upper_C, na.rm = T)
    g <- ggplot() + labs(x = NULL, y = NULL, title = title) + geom_tile(data = meltColor,
        mapping = aes(!!sym("Var2"), !!sym("Var1"), fill = !!sym("value"))) + geom_text(data = meltNum,
        mapping = aes(!!sym("Var2"), !!sym("Var1"), label = round(!!sym("value"),
            digits = 2))) + scale_x_discrete(position = "top") + scale_fill_gradient2(low = "#0072B2",
        mid = "white", high = "#d55E00", limit = c(-1, 1), name = "Pearson\nCorrelation") +
        theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank()) + coord_fixed()
    return(g)
}
