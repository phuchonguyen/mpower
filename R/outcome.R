#' @importFrom rlang parse_expr :=
#' @importFrom dplyr mutate select

#' @export
OutcomeModel <- function(f, family="gaussian") {
  mf <- NULL
  if (is.character(f)) {
    if (family == "gaussian") {
      mf <- function(X) {
        as.data.frame(X) %>%
          mutate(!! "y" := !! parse_expr(f)) %>%
          select("y") %>% as.matrix() %>% as.vector()
      }
    } else if (family == "binomial") {
      mf <- function(X) {
        mu <- as.data.frame(X) %>%
          mutate(!! "mu" := !! parse_expr(f)) %>%
          select("mu") %>% as.matrix() %>% as.vector()

        1/(1 + exp(-mu))
      }
    } else {
      stop("Family not implemented")
    }
  } else if (is.function(f)) {
    mf <- function(X) {
      f(X)
    }
  } else {
    stop("Argument `f` must be a string or function that returns a vector of outcome")
  }

  new_OutcomeModel(list(f=mf, s=1, sigma=1, family=family))

}

set_value <- function(object, name, value) {
  UseMethod("set_value")
}

set_value.mpower_OutcomeModel <- function(object, name, value) {
  for (i in seq_along(name)) {
    if (name[i] == "sigma" & value[i] < 0) stop("Noise variance `sigma` must be nonnegative")
    if (name[i] == "rho" & value[i] < 0) stop("SNR `rho` must be nonnegative")
    if (name[i] == "s" & value[i] < 0) stop("Mean function scale `s` must be nonnegative")
    object[[name[i]]] <- value[i]
  }

  object
}

set_value.default <- function(object, name, value) {
  warning("Not implemented")
  object
}

new_OutcomeModel <- function(y = list()) {
  stopifnot(is.list(y))
  stopifnot("f" %in% names(y))
  structure(y, class = "mpower_OutcomeModel")
}

geny <- function(y, X) {
  UseMethod("geny")
}

geny.mpower_OutcomeModel <- function(y, X) {
  if (y$family == "gaussian") {
    y$s * y$f(X) + rnorm(nrow(X), 0, y$sigma)
  } else if (y$family == "binomial") {
    rbinom(nrow(X), size = 1, prob = y$f(X))
  } else {
    stop("Family not implemented")
  }

}
