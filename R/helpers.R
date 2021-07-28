# Using cholesky decomp to check if symmetric matrix S is positive.definite
is_positive_definite <- function(S) {
  tryCatch({
    L <- chol(S)
    TRUE
    },
    error = function(e) FALSE)
}

# Source: https://github.com/christophM/iml/blob/master/R/utils.R
has_predict <- function(o) {
  classes <- class(o)
  any(unlist(lapply(classes, function(x) {
    "predict" %in% attr(methods(class = x), "info")$generic
  })))
}

log_msg <- function(text, dopar=FALSE, socket=NULL, ...) {
  msg <- sprintf(paste0(as.character(Sys.time()), ": ", text, "\n"), ...)
  cat(msg)
  if (dopar & !is.null(socket)) {
    write.socket(socket, msg)
  }
}
