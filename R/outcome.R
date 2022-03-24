#' Outcome generator
#'
#' This function creates a generative model of the outcome given a matrix of
#' predictors.
#'
#' @param f A string that describes the relationships between the
#'   predictors and outcome or a function that takes an input matrix and returns
#'   a vector of outcome: \eqn{E(y|X) = g(f(X))} where g is a link function
#'   that depends on the family argument.
#' @param family A string, "gaussian" or "binomial" for continuous or binary
#'   outcomes.
#' @param sigma A number, Gaussian noise standard deviation if applicable.
#' @return An OutcomeModel object. Attributes: f: mean function, sigma: a number Gaussian
#'   observation noise, family: a string "gaussian" or "binomial".
#' @export
OutcomeModel <- function(f, family = "gaussian", sigma = 1) {
  mu <- NULL
  if (is.character(f)) {
    if (family == "gaussian") {
      mu <- function(X) {
        as.data.frame(X) %>%
          dplyr::mutate(!! "fx" := !! parse_expr(f)) %>%
          dplyr::select("fx") %>% as.matrix() %>% as.vector()
      }
    } else if (family == "binomial") {
      mu <- function(X) {
        fx <- as.data.frame(X) %>%
          dplyr::mutate(!! "fx" := !! parse_expr(f)) %>%
          dplyr::select("fx") %>% as.matrix() %>% as.vector()

        1/(1 + exp(-fx))
      }
    } else {
      stop("Family not implemented")
    }
  } else if (is.function(f)) {
    mu <- function(X) {
      f(X)
    }
  } else {
    stop("Argument `f` must be a string or function")
  }
  new_OutcomeModel(list(f = mu, sigma = sigma, family = family))
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
    if (name[i] == "sigma" & value[i] < 0) stop("Noise sd must be nonnegative")
    if (name[i] == "snr" & value[i] < 0) stop("SNR must be nonnegative")
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
#' @param X A matrix of predictors
#' @export
geny <- function(obj, X) {
  UseMethod("geny")
}

geny.mpower_OutcomeModel <- function(obj, X) {
  if (obj$family == "gaussian") {
    obj$f(X) + rnorm(nrow(X), 0, obj$sigma)
  } else if (obj$family == "binomial") {
    rbinom(nrow(X), size = 1, prob = obj$f(X))
  } else {
    stop("Family not implemented")
  }

}
