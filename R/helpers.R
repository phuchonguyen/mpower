# Using cholesky decomp to check if symmetric matrix S is positive.definite
is_positive_definite <- function(S) {
    tryCatch({
        L <- chol(S)
        TRUE
    }, error = function(e) FALSE)
}
