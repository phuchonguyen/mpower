#' @export
MixtureModel <- function(..., method = "estimation") {
  args <- list(...)
  mod <- NULL

  if (method == "estimation") {
    if (!is.matrix(args[[1]]) & !is.data.frame(args[[1]])) {
      stop("Input must be a matrix or a data.frame",
      call. = FALSE)
    }
    name <- args$var_name
    if (is.null(name)) name <- colnames(args[[1]])
    if (is.null(name)) name <- paste0("x", 1:ncol(args[[1]]))
    if (is.null(args$per)) {
      mod <- new_estimation_MixtureModel(data = args[[1]], name = name)
    } else {
      mod <- new_estimation_MixtureModel(data = args[[1]], per = args$per, name = name)
    }

  }

  if (method == "resampling") {
    if (!is.matrix(args[[1]]) & !is.data.frame(args[[1]])) {
      stop("Input must be a matrix or a data.frame",
      call. = FALSE)
    }
    name <- args$var_name
    if (is.null(name)) name <- colnames(args[[1]])
    if (is.null(name)) name <- paste0("x", 1:ncol(args[[1]]))

    mod <- new_resampling_MixtureModel(data = args[[1]], name = name)
  }

  if (method == "cvine") {
    if (is.matrix(args$S)) {
      S <- args$S
    } else if (length(args$S)==1 & is.numeric(args$p)) {
      S <- matrix(args$S, nrow = args$p, ncol = args$p)
      diag(S) <- 1
    } else {
      stop("Input `S` is missing or not a numeric matrix", call. = FALSE)
    }

    if (!all(S <= 1 & S >= -1)) stop("Values in S must be between -1 and 1", call. = FALSE)
    name <- args$var_name
    if (is.null(name)) name <- colnames(S)
    if (is.null(name)) name <- paste0("x", 1:ncol(S))

    if (is.numeric(args$m) & length(args$m)==1) {
      m <- args$m
    } else {
      stop("Input `m` is missing or not a scalar", call. = FALSE)
    }
    if (m <= 0) {stop("Input `m` must be positive", call. = FALSE)}

    mod <- new_cvine_MixtureModel(S = S, m = m, name = name)
  }

  validate_MixtureModel(mod)
}

validate_MixtureModel <- function(x) {
  if (!is.null(x$data) & !is.matrix(x$data) & !is.data.frame(x$data)) {
    stop("`data` must be a matrix or a data.frame",
    call. = FALSE)
  }

  if (!is.null(x$var_name) & !is.character(x$var_name)) {
    stop("Predictors' names `var_name` must be character",
    call. = FALSE)
  }

  if (!is.null(x$R) & !is_positive_definite(x$R)) {
    stop("Predictors' correlation matrix `R` is not positive definite",
    call. = FALSE)
  }

  if (!is.null(x$per)) {
    if (x$per < 0) stop("Perpturbation `per` is not positive", call. = FALSE)
  }

  if (!is.null(x$m)) {
    if (x$m <= 0) stop("Argument `m` is not positive", call. = FALSE)
  }

  x
}

new_resampling_MixtureModel <- function(data = numeric(), name = character()) {
  stopifnot(is.numeric(as.matrix(data)) & is.character(name))
  x <- list(data = data, var_name = name, p = ncol(data))
  structure(x, class = "mpower_resampling_MixtureModel")
}

new_estimation_MixtureModel <- function(data = numeric(), per = 10e-10, name = character()) {
  stopifnot(is.numeric(as.matrix(data)) & is.numeric(per) & is.character(name))
  C <- cov(data)
  while (!is_positive_definite(C)) {
    C <- C + diag(per, ncol(data))
    per <- 10*per
  }
  R <- cov2cor(C)
  x <- list(data = data, R = R, per = per, var_name = name, p = ncol(R))
  structure(x, class = "mpower_estimation_MixtureModel")
}

new_cvine_MixtureModel <- function(S = numeric(), m = 100, name = character(), per = 10e-10) {
  stopifnot(is.numeric(S) & is.numeric(m) & is.character(name) & is.numeric(per))
  R <- cvine(d = ncol(S), S = S, m = m)
  while (!is_positive_definite(R)) {
    R <- R + diag(per, ncol(S))
    per <- 10*per
  }
  x <- list(R = R, S = S, m = m, var_name = name, per = per, p = ncol(R))
  structure(x, class = "mpower_cvine_MixtureModel")
}

genx <- function(x, n) {
  UseMethod("genx")
}

genx.mpower_resampling_MixtureModel <- function(x, n) {
  stopifnot(n >= 0, n %% 1 == 0)
  x$data[sample(1:nrow(x$data), n, replace = TRUE), ] %>%
    magrittr::set_colnames(x$var_name)
}

genx.default <- function(x, n) {
  stopifnot(n >= 0, n %% 1 == 0)
  MASS::mvrnorm(n, mu = rep(0, ncol(x$R)), Sigma = x$R) %>%
    magrittr::set_colnames(x$var_name)
}
