#' Outcome generator
#'
#' This function creates a generative model of the outcome given a matrix of
#' predictors.
#'
#' @param f A string that describes the relationships between the predictors and
#'   outcome or a function that takes an input matrix and returns a vector of
#'   outcome: \eqn{E(y|x) = g(f(x))} where g is a link function that depends on
#'   the family argument.
#' @param family A string, 'gaussian', 'binomial', or 'poisson' for continuous,
#'   binary, or count outcomes.
#' @param sigma A number, Gaussian noise standard deviation if applicable.
#' @param f_args A named list of additional arguments to f
#' @return An OutcomeModel object. Attributes:
#' \item{f}{mean function.}
#' \item{sigma}{a number for the Gaussian observation noise.}
#' \item{family}{a string 'gaussian' or 'binomial'.}
#' @examples
#' # Define BMI as a ratio of weight and height plus random Gaussian error with standard deviation 1.
#' bmi_model <- mpower::OutcomeModel(f = 'weight/(height^2)', sigma = 1, family = 'gaussian')
#' @export
OutcomeModel <- function(f, family = "gaussian", sigma = 1, f_args = list()) {
    mu <- NULL
    if (is.character(f)) {
        if (family == "gaussian") {
            mu <- function(x) {
                parse_text_to_f(x, f)
            }
        } else if (family == "binomial") {
            mu <- function(x) {
                fx <- parse_text_to_f(x, f)

                1/(1 + exp(-fx))
            }
        } else if (family == "poisson") {
            mu <- function(x) {
                fx <- parse_text_to_f(x, f)

                exp(fx)
            }
        } else {
            stop("Family not implemented")
        }
    } else if (is.function(f)) {
        mu <- function(x) {
            do.call(f, args = c(list(x), f_args))
        }
    } else {
        stop("Argument `f` must be a string or function")
    }
    new_OutcomeModel(list(f = mu, sigma = sigma, family = family))
}

parse_text_to_f <- function(x, f) {
    as.data.frame(x) %>%
        dplyr::mutate(`:=`(!!"fx", !!parse_expr(f))) %>%
        dplyr::select("fx") %>%
        as.matrix() %>%
        as.vector()
}

#' This function updates values in an OutcomeModel object
#' @param obj An OutcomeModel object
#' @param name A string for name of the attribute to be changed
#' @param value An appropriate data type
set_value <- function(obj, name, value) {
    UseMethod("set_value")
}

set_value.mpower_OutcomeModel <- function(obj, name, value) {
    for (i in seq_along(name)) {
        if (name[i] == "sigma" & value[i] < 0)
            stop("Noise sd must be nonnegative")
        if (name[i] == "snr" & value[i] < 0)
            stop("SNR must be nonnegative")
        obj[[name[i]]] <- value[i]
    }

    obj
}

set_value.default <- function(obj, name, value) {
    warning("Not implemented")
    obj
}

new_OutcomeModel <- function(y = list()) {
    stopifnot(is.list(y))
    stopifnot("f" %in% names(y))
    structure(y, class = "mpower_OutcomeModel")
}

#' Generates a vector of outcomes
#' @param obj An OutcomeModel object
#' @param x An (n x p) matrix of predictors
#' @return An n-vector of outcomes
#' @export
geny <- function(obj, x) {
    UseMethod("geny")
}

#' @export
geny.mpower_OutcomeModel <- function(obj, x) {
    if (obj$family == "gaussian") {
        obj$f(x) + stats::rnorm(nrow(x), 0, obj$sigma)
    } else if (obj$family == "binomial") {
        stats::rbinom(nrow(x), size = 1, prob = obj$f(x))
    } else if (obj$family == "poisson") {
        stats::rpois(nrow(x), obj$f(x))
    } else {
        stop("Family not implemented")
    }

}
