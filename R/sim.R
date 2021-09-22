#' @importFrom magrittr %>% set_colnames
#' @importFrom utils combn make.socket methods write.socket
#' @importFrom ggplot2 ggplot aes aes_string scale_colour_discrete geom_path geom_point labs
#' @importFrom purrr map
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom abind abind


#' @export
sim_power <- function(xmod, ymod, imod, s=1000, n=100,
                sigma=NULL, rho=NULL, alpha=0.05,
                cores=1, m=5000, file=NULL) {
  # scale effect size
  if (is.null(sigma) & !is.null(rho)) {
    sigma <- scale_sigma(f = ymod$f, rho = rho, xmod = xmod, m = m)$sigma
    ymod <- set_value(ymod, c("sigma", "rho"), c(sigma, rho))
  } else if (!is.null(sigma) & !is.null(rho)) {
    scale_s <- scale_f(f = ymod$f, rho = rho, sigma = sigma, xmod = xmod, m = m)$scale_s
    ymod <- set_value(ymod, c("sigma", "rho", "s"), c(sigma, rho, scale_s))
  } else if (!is.null(sigma) & is.null(rho)) {
    ymod <- set_value(ymod, c("sigma", "s"), c(sigma, 1))
  } else {
    stop("sigma or rho missing")
  }

  # run MC sims
  prog_int <- max(1, s %/% 5)
  if (cores < 2) {
    out <- as.list(rep(NA, s))
    for (i in seq_along(out)) {
      if (i %% prog_int == 0) log_msg(paste0("Processing ", i," simulations of ", s))
      X <- genx(xmod, n)
      y <- geny(ymod, X)
      out[[i]] <- fit(imod, X, y, alpha)
    }
  } else {
    # Set up parallel backend to use many processors
    cores <- detectCores()
    cl <- makeCluster(cores[1]-1) #not to overload your computer
    registerDoParallel(cl)
    out <- foreach(i = 1:s) %dopar% {
      if (i %% prog_int == 0) log_msg(paste0("Processing ", i," simulations of ", s))
      X <- genx(xmod, n)
      y <- geny(ymod, X)
      fit(imod, X, y, alpha)
    }
    stopCluster(cl)
  }

  if (!is.null(file)) saveRDS(out, file = file)

  new_Sim(list(s = s, n = n, alpha = alpha,
               xmod = xmod, ymod = ymod, imod = imod,
               sims = out))
}

new_Sim <- function(x = list()) {
  stopifnot(is.list(x))
  structure(x, class = "mpower_Sim")
}

#' @export
summarize_power <- function(object, crit, thres=NULL, digits=3) {
  if (crit == "ci") message(paste("For a", 1-object$alpha, "% CI", sep = " "))
  if (!(crit %in% names(object$sims[[1]]))) stop(paste("Criterion", crit, "not implemented for", object$imod$fun_name))
  temp <- object$sims %>%
    map(crit) %>%
    map(function(e) if(!is.null(thres)){e >= thres} else {e}) %>%
    abind(along = -1)
  if (length(dim(temp))==3) {
    temp %>%
      colMeans() %>% round(digits) %>%
      set_colnames(object$xmod$var_name)
  } else if (is.matrix(temp)) {
    temp %>%
      set_colnames(object$xmod$var_name) %>%
      colMeans() %>% round(digits)
  } else {
    temp %>% colMeans() %>% round(digits)
  }
}

#' @export
plot_pip_thresholds <- function(object, thres=seq(0, 1, 0.1), digits=3) {
  p <- object$xmod$p
  res <- matrix(NA, length(thres), p)
  if (!("pip" %in% objects$sims[[1]])) stop(paste("Criterion PIP not implemented for", object$imod$fun_name, sep=" "))
  for (i in seq_along(thres)) {
    res[i, ] <- object$sims %>%
      map("pip") %>%
      map(function(x) x >= thres[i]) %>%
      abind(along = -1) %>%
      colMeans() %>% round(digits)
  }
  tibble::as_tibble(cbind(thres, res)) %>%
    set_colnames(c("thres", object$xmod$var_name)) %>%
    tidyr::pivot_longer(cols = tidyselect::starts_with("x")) %>%
    ggplot(aes_string(x = "thres", y = "value", group = "name", color = "name")) +
    geom_point() + geom_path() + labs(y = "Type I error rate/ Power", x = "PIP threshold")
}

#' @export
sim_curve <- function(xmod, ymod, imod, s=1000, n=100, sigma=1, rho=NULL, alpha=0.05,
                     cores=1, m=5000, files=NULL) {
  powr <- list()
  k <- 1
  for (i in seq_along(n)) {
    nsize <- n[i]
    for (j in seq_along(sigma)) {
      vnoise <- sigma[j]
      vrho <- rho[j]
      powr[[k]] <- sim_power(xmod = xmod, ymod = ymod, imod = imod,
                           s = s, n = nsize,
                           sigma = vnoise, rho = vrho,
                           alpha = alpha, cores = cores, m = m, file = files[k])
      k <- k+1
    }
  }

  powr
}

#' @export
plot_power_curve <- function(x, crit, var, thres=NULL, digits=3) {
  n <- rep(NA, length(x))
  s <- rep(NA, length(x))
  powr <- matrix(NA, nrow = length(x), ncol = x[[1]]$xmod$p)
  k <- 1
  for (l in x) {
    n[k] <- l$n
    s[k] <- l$ymod$rho
    temp <- summarize_power(object = l, crit = crit, thres = thres, digits = digits)
    if (!is.null(nrow(temp))) temp <- temp[2,]
    powr[k,] <- temp
    k <- k+1
  }
  data.frame(cbind(n, s, powr)) %>%
    set_colnames(c("n", "s", x[[1]]$xmod$var_name)) %>%
    mutate(n = as.factor(n)) %>%
    ggplot(aes_string(x = "s", y = var, colour = "n")) + geom_path() + geom_point() +
    labs(x="SNR", y="Type I error rate/ Power") +
    scale_colour_discrete(name = "Sample size")
}






