#' Statistical model that returns significance criterion
#'
#' This function creates a wrapper function for a statistical model and its
#' applicable significance criterion. It finds relationships between a matrix
#' of predictors and a vector of outcomes using the statistical model, and
#' determines if the relationships are 'significant' according to pre-specified
#' criterion for that model.
#'
#' @param model A string of the name of a built-in statistical model or a function
#'   implements a statistical model and returns a list of significance
#'   criteria for each predictor. Built-in options include 'bma', 'bkmr', 'mixselect',
#'   'bws', 'qgc', 'fin', 'glm'.
#' @param name A string, name of the statistical model. Default to string input
#'   of model.
#' @param ... Additional keyword arguments for the statistical model
#' @return An InferenceModel object.
#' \item{model}{a function that takes matrices of predictors and outcomes and
#' returns a list of significance criteria.}
#' \item{model_name}{a string.}
#' @examples
#' imod <- mpower::InferenceModel(model = 'glm', family = 'gaussian',
#' formula = y ~ Poverty*(poly(Age, 2) + HHIncome + HomeOwn + Education))
#' @export
InferenceModel <- function(model, name = NULL, ...) {
    mod <- list(model_name = model, args = list(...))
    if (is.function(model)) {
        mod[["model"]] <- model
        mod[["model_name"]] <- name
    } else if (is.character(model)) {
        if (model == "bma") {
            if (!requireNamespace("BMA")) {
                stop("Package 'BMA' not installed.")
            }
            mod[["model"]] <- bma_wrapper
        } else if (model == "bkmr") {
            if (!requireNamespace("bkmr")) {
                stop("Package 'bkmr' not installed.")
            }
            mod[["model"]] <- bkmr_wrapper
        } else if (model == "bws") {
            if (!requireNamespace("bws")) {
                stop("Package 'bws' not installed.")
            }
            message("Estimating the power of detecting the overall effect based on posterior credible interval. Running sim_power() with cores > 1 is recommended.")
            mod[["model"]] <- bws_wrapper
        } else if (model == "qgc") {
            if (!requireNamespace("qgcomp")) {
                stop("Package 'qgcomp' not installed.")
            }
            mod[["model"]] <- qgcomp_lin_wrapper
        } else if (model == "fin") {
            if (!requireNamespace("infinitefactor")) {
                stop("Package 'infinitefactor' not installed.")
            }
            mod[["model"]] <- fin_wrapper
        } else if (model == "glm") {
            message("Estimating the power of conditional t-test on each regression coefficient.")
            mod[["model"]] <- glm_wrapper
        } else {
            stop(paste(model, "not implemented"))
        }
    } else {
        stop("Invalid input for `model`")
    }

    new_InferenceModel(mod)
}

new_InferenceModel <- function(l = list()) {
    stopifnot(is.list(l))
    structure(l, class = "mpower_InferenceModel")
}

#' Fits the model to given data and gets a list of significance criteria
#' @param mod An InferenceModel object
#' @param x A matrix of predictors
#' @param y A vector of outcome
#' @return A list of some of the following significance criteria:
#' \item{beta}{The smaller
#' posterior probability of being to one side of zero for linear term, given
#' either the main effect or interaction is non-zero. Applicable to 'bma',
#' 'bws', 'fin', and 'mixselect' model.}
#' \item{interact_beta}{Same as linear_beta but
#' for pair-wise interactions. Applicable to 'fin' model.}
#' \item{pip}{Posterior inclusion probability (PIP) of either a linear or
#' non-linear effect. Applicable to 'bma', 'bkmr', and 'mixselect' model.}
#' \item{group_pip}{PIP of either a linear or non-linear effect. Applicable to
#' 'bkmr' model.}
#' \item{linear_pip}{PIP of a linear effect. Applicable to 'mixselect' model.}
#' \item{gp_pip}{PIP of a non-linear effect. Applicable to 'mixselect' model.}
#' \item{pval}{The p-value of the combined effect, the same for all
#' predictors. Applicable to 'glm', and 'qgc' model.}
#' \item{time}{elapsed time to fit the model.}
#' @export
fit <- function(mod, x, y) {
    UseMethod("fit")
}

#' @export
fit.mpower_InferenceModel <- function(mod, x, y) {
    if (length(mod$args) == 0) {
        return(mod$model(y = y, x = x))
    }
    return(mod$model(y = y, x = x, args = mod$args))
}

#' Fits a BKMR model with significance criteria: PIP and group-wise PIP
#' @param x A matrix of predictors
#' @param y A vector of outcome
#' @param args A list of arguments, see R function `bkmr::kmbayes()`
#' @return A list of vectors whose values are between 0 and 1
#' \item{pip}{PIP for component-wise selection or conditional (with-in group)
#' PIP for hierarchical variable selection.}
#' \item{group_pip}{PIP for group-specific selection.}
#' \item{time}{elapsed time to fit the model.}
#' @section Reference:
#'
#'   Bobb JF, Henn BC, Valeri L, Coull BA (2018). “Statistical software for
#'   analyzing the health effects of multiple concurrent exposures via Bayesian
#'   kernel machine regression.”Environ-mental
#'   Health,17(67).doi:10.1186/s12940-018-0413-y.
#' @export
bkmr_wrapper <- function(y, x, args = list()) {
    args[["varsel"]] <- TRUE
    x <- stats::model.matrix(~. - 1, x)
    s <- Sys.time()
    km_fit <- do.call(bkmr::kmbayes, c(list(y = y, Z = x), args))
    km_time <- Sys.time() - s
    km_pip <- bkmr::ExtractPIPs(km_fit)  # by default keeps the second half of all samples
    # km_pip <- km_pip[, (ncol(km_pip)/2 + 1):ncol(km_pip)]
    pip <- km_pip$PIP
    if (is.null(pip))
        pip <- km_pip$condPIP
    names(pip) <- colnames(x)
    group_pip <- km_pip$groupPIP
    if (!is.null(group_pip))
        names(group_pip) <- colnames(x)
    return(list(pip = pip, group_pip = group_pip, time = km_time))
}

#' Fits a linear model with Bayesian model selection with significance criteria:
#' PIP and posterior probability of nonzero coefficients being on one side of
#' zero.
#' @param x A matrix of predictors
#' @param y A vector of outcome
#' @param args A list of arguments see R function `BMA::bic.glm()`.
#' @return A list of vectors whose values are between 0 and 1
#' \item{beta}{The smaller posterior probability of the coefficients being to
#' one side of zero: min(Pr(beta >0), Pr(beta<0)).}
#' \item{pip}{PIP of the
#' effect not being zero.}
#' \item{time}{elapsed time to fit the model.}
#' @section Reference:
#'
#'   Raftery A, Hoeting J, Volinsky C, Painter I, Yeung KY (2021).BMA: Bayesian
#'   model averaging. R package version 3.18.15.
#' @export
bma_wrapper <- function(y, x, args = list()) {
    x <- stats::model.matrix(~., x)[, -1]  # full rank matrix with dummy variables
    s <- Sys.time()
    bma_out <- do.call(BMA::bic.glm, c(list(x = x, y = y), args))
    bma_time <- Sys.time() - s
    p <- ncol(x)  # Model-averaged coef, the first one is the intercept term
    # Posterior probability of being less than zero.
    bma_zero <- stats::pnorm(0, bma_out$postmean[2:(p + 1)], bma_out$postsd[2:(p +
        1)])
    bma_zero <- vapply(bma_zero, function(x) min(x, 1 - x), numeric(1))
    # Sum inclusion probability postprob of all models with a main effec Xi
    bma_pip <- bma_out$probne0
    names(bma_pip) <- names(bma_zero) <- colnames(x)
    return(list(beta = bma_zero, pip = bma_pip, time = bma_time))
}

#' Fits a Bayesian factor model with interactions
#' @param x A matrix of predictors
#' @param y A vector of outcome
#' @param args A list of arguments see R function
#'   `infinitefactor::interactionDL()` in 'infinitefactor' package.
#' @return A list of vectors whose values are between 0 and 1
#' \item{beta}{The smallest posterior probability of the coefficients being to
#' one side of zero for either main effect or interaction: min(Pr(beta >0),
#' Pr(beta<0)).}
#' \item{linear_beta}{The smaller of posterior probability of
#' the main effects being to one side of zero.}
#' \item{interact_beta}{Same as
#' linear_beta but for pair-wise interactions.}
#' \item{time}{elapsed time to fit the model.}
#' @section Reference:
#'
#'   Ferrari F, Dunson DB (2020). “Bayesian factor analysis for inference on
#'   interactions.”Journal of the American Statistical Association, 0(0), 1–12.
#' @export
fin_wrapper <- function(y, x, args = list(nrun = 2000)) {
    y <- matrix(y, length(y), 1)
    x <- stats::model.matrix(~. - 1, x)
    if (!is.null(args$z))
        warning("Inclusion of linear effects for covariates not implemented in package infinitefactor for FIN")
    s <- Sys.time()
    fin_out <- do.call(infinitefactor::interactionDL, c(list(y = y, X = x), args))
    fin_time <- Sys.time() - s
    # Posterior probability of effects being to one side of zero
    fin_beta <- apply(fin_out$mainEffectSamps > 0, 1, mean, na.rm = T)
    fin_beta <- vapply(fin_beta, function(x) min(x, 1 - x), numeric(1))
    fin_int <- apply(apply(fin_out$interactionSamps, 3, c) > 0, 1, mean, na.rm = T)
    fin_int <- vapply(fin_int, function(x) min(x, 1 - x), numeric(1))
    p <- ncol(x)
    names(fin_int) <- paste0(rep(1:p, each = p), rep(1:p, p))
    fin_int <- fin_int[c(paste0(1:p, 1:p), utils::combn(1:p, m = 2, FUN = paste0,
        collapse = ""))]
    # Significance of a variable: either main effect or interaction is
    # significant
    fin_zero <- rep(NA, p)
    for (i in 1:p) {
        id <- which(grepl(i, colnames(fin_int)))
        fin_zero[i] <- min(fin_beta[i], fin_int[id])
    }
    names(fin_zero) <- names(fin_beta) <- colnames(x)
    return(list(beta = fin_zero, linear_beta = fin_beta, interact_beta = fin_int,
        time = fin_time))
}

#' Fits a linear Quantile G-Computation model with no interactions
#' @param x A matrix of predictors
#' @param y A vector of outcome
#' @param args A list of arguments see R function `qgcomp::qgcomp.noboot()`.
#' @return A list
#' \item{pval}{The p-value of the combined effect, the
#' same for all predictors.}
#' \item{pos_weights}{Positive weights. See 'qgcomp' package.}
#' \item{neg_weights}{Negative weights. See 'qgcomp' package.}
#' \item{time}{elapsed time to fit the model.}
#' @section Reference:
#'
#'   Keil AP, Buckley JP, O’Brien KM, Ferguson KK, Zhao S, White AJ (2020). “A
#'   Quantile-based g-computation  approach  to  addressing  the  effects of
#'   exposure  mixtures.”Environmental Health Perspectives, 128(4).
#' @export
qgcomp_lin_wrapper <- function(y, x, args = list()) {
    dat <- data.frame(y = y, x = x)
    colnames(dat) <- c("y", colnames(x))
    s <- Sys.time()
    qc.fit <- do.call(qgcomp::qgcomp.noboot, c(list(y ~ ., data = dat), args))
    qgc_time <- Sys.time() - s
    p <- ncol(x)
    qgc_pval <- qc.fit$pval[2]  #rep(qc.fit$pval[2], p)
    names(qgc_pval) <- "overall"
    return(list(pval = qgc_pval, pos_weights = qc.fit$pos.weights, neg_weights = qc.fit$neg.weights,
        time = qgc_time))
}

#' Fits a generalized linear model
#' @param x A matrix of predictors
#' @param y A vector of outcome
#' @param args A list of arguments see R `glm` function.
#' @return A list
#' \item{pval}{The p-value of the linear main effect.}
#' \item{time}{elapsed time to fit the model.}
#' @export
glm_wrapper <- function(y, x, args = list()) {
    if (is.null(args[["formula"]])) {
        args[["formula"]] <- y ~ .
    }
    dat <- data.frame(y = y, x = x) %>%
        set_colnames(c("y", colnames(x)))
    s <- Sys.time()
    fit <- do.call(stats::glm, c(list(data = dat), args))
    time <- Sys.time() - s
    var_names <- colnames(fit$model)[-1]
    p <- length(var_names)
    all_pval <- stats::summary.glm(fit)$coef[, 4]
    all_int_pval <- all_pval[grepl("\\:", names(all_pval))]
    all_main_pval <- all_pval[!grepl("\\:", names(all_pval))]
    pval <- int_pval <- main_pval <- rep(NA, p)
    names(pval) <- names(int_pval) <- names(main_pval) <- var_names
    # Save the smallest p-value of all terms with the variable name
    for (var in var_names) {
        int_pval[var] <- min(all_int_pval[grepl(var, names(all_int_pval))])
        main_pval[var] <- min(all_main_pval[grepl(var, names(all_main_pval))])
        pval[var] <- min(int_pval[var], main_pval[var])
    }
    return(list(pval = pval, main_pval = main_pval, int_pval = int_pval, time = time))
}

#' Fits a Bayesian weighted sums
#' @param x A matrix of predictors
#' @param y A vector of outcome
#' @param args A list of arguments see R `bws::bws()`` function.
#' @return A list
#' \item{beta}{The smaller posterior probability of
#' the combined overall effect being to one side of zero: min(Pr(beta >0),
#' Pr(beta<0)). The same for all predictor.}
#' \item{weights}{The 95\% CI of the
#' contribution of each predictor to the overall effect. Between 0 and 1.}
#' \item{time}{elapsed time to fit the model.}
#' @section Reference:
#'
#'   Hamra GB, MacLehose RF, Croen L, Kauffman EM, Newschaffer C (2021).
#'   “Bayesian weighted sums: a flexible approach to estimate summed mixture
#'   effects.” International Journal of Environmental Research and Public
#'   Health, 18(4), 1373.
#' @export
bws_wrapper <- function(y, x, args = list(iter = 2000)) {
    message("Sampling bws")
    s <- Sys.time()
    fit <- do.call(bws::bws, c(list(y = y, X = x), args))
    bws_time <- Sys.time() - s
    samps <- rstan::extract(fit)
    beta_zero <- mean(samps$theta1 > 0)
    beta_zero <- min(beta_zero, 1 - beta_zero)
    # beta_zero <- rep(beta_zero, ncol(x))
    weights <- apply(samps$w, 2, mean)
    names(beta_zero) <- "overall"
    names(weights) <- colnames(x)
    return(list(beta = beta_zero, weights = weights, time = bws_time))
}
