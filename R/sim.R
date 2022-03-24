#' Power analysis for multiple settings using Monte Carlo simulation
#'
#' This function can be used to create power curves by calling
#'   sim_power() on combinations of its inputs.
#'
#' @param xmod A MixtureModel object.
#' @param ymod One or a list of OutcomeModel object(s).
#' @param imod An InferenceModel object.
#' @param s An integer for number of Monte Carlo simulation.
#' @param n An interger or a vector of sample sizes.
#' @param cores An integer for number of processing cores. When cores > 1,
#'   parallelism is automatically applied.
#' @param file A string, a file name with no extension to write samples to
#'   periodically. By default write to an RDS file.
#' @param errorhandling A string "remove", "stop", or "pass". If an error occurs
#'   in any iteration, remove that iteration (remove), return the error message
#'   verbatim in the output (pass), or terminate the loop (stop). Default is
#'   "remove". See R package `foreach` for more details.
#' @return A SimCurve object. Attributes: s: a number, snr: a number or list of
#'   numbers, n: a number or list of numbers, xmod: a MixtureModel, ymod: one or
#'   a list of OutcomeModels, imod: an InferenceModel, sims: a list of simulation
#'   output matrices.
#' @export
sim_curve <- function(xmod, ymod, imod, s=100, n=100,
                      cores = 1, file = NULL, errorhandling = "remove",
                      snr_iter = 10000) {
  stopifnot(class(xmod) %in% c("mpower_estimation_MixtureModel",
                               "mpower_cvine_MixtureModel",
                               "mpower_resampling_MixtureModel"))
  stopifnot(class(imod) == "mpower_InferenceModel")
  if (class(ymod) == "mpower_OutcomeModel") {
    ymod <- list(ymod)
  }
  sim_list <- list()
  ymod_list <- list()
  n_vec <- c()
  snr_vec <- rep(0, length(ymod))
  cur_file <- NULL
  k <- 1
  for (i in seq_along(n)) {
    cur_n <- n[i]
    for (j in seq_along(ymod)) {
      cur_ymod <- ymod[[j]]
      cur_file <- if(!is.null(file)) {paste(file, i, j, sep = "-")}
      print(paste("Simulation for n =", cur_n, "and the", j, "th outcome model"))
      out <- sim_power(xmod = xmod, ymod = cur_ymod, imod = imod,
                       s = s, n = cur_n, snr_iter = snr_iter,
                       cores = cores, file = cur_file)
      sim_list[[k]] <- out$sims
      ymod_list[[k]] <- out$ymod
      snr_vec[j] <- snr_vec[j] + out$snr
      n_vec <- c(n_vec, out$n)
      k <- k+1
    }
  }
  new_SimCurve(list(s = s, n = n_vec,
                    snr = rep(snr_vec/length(n), length(n)),
                    xmod = xmod, ymod = ymod_list, imod = imod,
                    sims = sim_list))
}

new_SimCurve <- function(x = list()) {
  stopifnot(is.list(x))
  structure(x, class = "mpower_SimCurve")
}


#' Power analysis using Monte Carlo simulation
#'
#' @param xmod A MixtureModel object.
#' @param ymod An OutcomeModel object.
#' @param imod An InferenceModel object.
#' @param s An integer for number of Monte Carlo simulation.
#' @param n An interger for sample size in each simulation.
#' @param cores An integer for number of processing cores. When cores > 1,
#'   parallelism is automatically applied.
#' @param file A string, a file name with no extension to write samples to
#'   periodically. By default write to an RDS file.
#' @param errorhandling A string "remove", "stop", or "pass". If an error occurs
#'   in any iteration, remove that iteration (remove), return the error message
#'   verbatim in the output (pass), or terminate the loop (stop). Default is
#'   "remove". See R package `foreach` for more details.
#' @param snr_iter An integer for number of Monte Carlo samples to estimate SNR
#' @return A PowerSim object. Attributes: xmod, ymod, imod, s, n, simulation
#'   samples.
#' @export
sim_power <- function(xmod, ymod, imod, s = 100, n = 100,
                      cores = 1, file = NULL, errorhandling = "remove",
                      snr_iter = 10000) {
  stopifnot(class(xmod) %in% c("mpower_estimation_MixtureModel",
                               "mpower_cvine_MixtureModel",
                               "mpower_resampling_MixtureModel"))
  stopifnot(class(ymod) == "mpower_OutcomeModel")
  stopifnot(class(imod) == "mpower_InferenceModel")
  pb <- utils::txtProgressBar(max = s, style = 3)
  progress <- function(i) utils::setTxtProgressBar(pb, i)
  if (cores < 2) {
    out <- as.list(rep(NA, s))
    for (i in seq_along(out)) {
      out[[i]] <- tryCatch({
        X <- genx(xmod, n)
        y <- geny(ymod, X)
        fit(imod, X, y)
        },
        error=function(e) {
          message(e)
          return(NA)
        })
      progress(i)
      if (!is.null(file) & (i%%50==0)) saveRDS(out, file = paste0(file, ".rds"))
    }
  } else {
    # Set up parallel backend to use many processors
    cores_max <- parallel::detectCores()
    cl <- parallel::makeCluster(min(cores_max[1]-1, cores))
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress=progress)
    out <- foreach::foreach(i = 1:s, .errorhandling=errorhandling, .options.snow=opts,
                            .export = c("genx.mpower_estimation_MixtureModel",
                                        "genx.mpower_resampling_MixtureModel",
                                        "genx.mpower_cvine_MixtureModel",
                                        "geny", "fit.mpower_InferenceModel")) %dopar% {
      X <- genx(xmod, n)
      y <- geny(ymod, X)
      r <- fit(imod, X, y)
      ## write results to disk
      if (!is.null(file)) saveRDS(r, file = paste0(file, i, ".rds"))
      return(r)
    }
    parallel::stopCluster(cl)
  }
  close(pb)
  snr <- estimate_snr(ymod, xmod, m = snr_iter)
  new_Sim(list(s = s, n = n, snr = snr,
               xmod = xmod, ymod = ymod, imod = imod,
               sims = out))
}

#TODO: Add SNR attribute to this!
new_Sim <- function(x = list()) {
  stopifnot(is.list(x))
  structure(x, class = "mpower_Sim")
}


#' Tabular and graphical summaries of power simulation
#'
#' @param sim A Sim or a list of Sim objects, output from [sim_power()] or
#'   [sim_curve()].
#' @param crit A string specifying the significance criteria.
#' @param thres A number or vector of numbers spefifying the thresholds of
#'   "significance".
#' @param digits An integer for the number of decimal points to display.
#' @return A data.frame summary of power for each predictor for each
#'   combination of thresholds, sample size, signal-to-noise ratios.
#' @export
summary <- function(sim, crit, thres, digits = 3) {
  UseMethod("summary")
}

summary.mpower_Sim <- function(sim, crit, thres, digits = 3, how = "greater") {
  if(seq_along(crit) > 1)  stop("Summary limited to one criteria", call. = FALSE)
  cat("\n\t*** POWER ANALYSIS SUMMARY ***")
  cat("\nNumber of Monte Carlo simulations:", sim$s)
  cat("\nNumber of observations in each simulation:", sim$n)
  cat("\nData generating process estimated SNR:", round(sim$snr,2))
  cat("\nInference model:", sim$imod$model_name)
  cat("\nSignificance criterion:", crit)
  if (length(thres) == 1) {
    # For one threshold and one Sim object
    res <- summary_one_sim_one_thres(sim$sims, crit, thres, digits, how, pivot = T)
  } else if (length(thres) > 1) {
    # For many threshold and one Sim object
    res <- summary_one_sim_many_thres(sim$sims, crit, thres, digits, how)
  }
  colnames(res) <- c("thres", " ", "power")
  for (i in 1:length(thres)) {
    cat("\n\nSignificance threshold: ", thres[i])
    cat("\n", res %>% dplyr::filter(thres == thres[i]) %>%
          dplyr::select(-thres) %>% knitr::kable(), sep="\n")
  }
  return(res)
}

summary.mpower_SimCurve <- function(sim, crit, thres, digits = 3, how = "greater") {
  cat("\n\t*** POWER CURVE ANALYSIS SUMMARY ***")
  cat("\nNumber of Monte Carlo simulations:", sim$s)
  cat("\nNumber of observations in each simulation:", unique(sim$n))
  cat("\nData generating process estimated SNR (for each outcome model):",
      round(unique(sim$snr),2))
  cat("\nInference model:", sim$imod$model_name)
  cat("\nSignificance criterion:", crit)
  cat("\nSignificance threshold: ", thres[1])
  if (length(thres) == 1) {
    # For one threshold and many Sim objects
    res <- summary_many_sim_one_thres(sim$sims, sim$n, sim$snr, crit, thres, digits, how)
  } else {
    warning("Using the first threshold only")
    res <- summary_many_sim_one_thres(sim$sims, sim$n, sim$snr, crit, thres[1], digits, how)
  }
  cat("\n", res %>%
        dplyr::select(test, power, n, snr, -thres) %>%
        magrittr::set_colnames(c(" ", "power", "n", "snr")) %>%
        knitr::kable(), sep="\n")
  return(res)
}

summary_many_sim_one_thres <- function(sim, n, snr, crit, thres, digits, how) {
  stopifnot(is.list(sim))
  num_sims <- length(sim)
  res <- tibble::tibble()
  for (i in 1:num_sims) {
    res <- dplyr::bind_rows(res,
                            summary_one_sim_one_thres(sim[[i]], crit, thres, digits, how, pivot = F))
  }
  res <- res %>%
    dplyr::mutate(n = n, snr = round(snr, 3)) %>%
    tidyr::pivot_longer(-c("n", "snr", "thres"), names_to = "test", values_to = "power")
  g <- res %>%
    dplyr::mutate(snr = as.factor(snr)) %>%
    ggplot(aes_string(x = "n", y = "power", colour = "snr")) + geom_path() + geom_point() +
    labs(x="Sample size", y="Type I error rate/ Power") +
    scale_colour_discrete(name = "SNR") +
    facet_wrap(as.formula("~test"))
  print(g)
  return(res)
}

summary_one_sim_many_thres <- function(sim, crit, thres, digits, how) {
  res <- data.frame()
  for (i in seq_along(thres)) {
    res <- dplyr::bind_rows(res,
                            summary_one_sim_one_thres(sim, crit, thres[i], digits, how, pivot = T))
  }
  g <- res %>%
    ggplot(aes_string(x = "thres", y = "power", group = "test", color = "test")) +
    geom_point() + geom_path() +
    labs(y = "Type I error rate/ Power", x = "Significance threshold") +
    scale_colour_discrete(name = "")
  print(g)
  return(res)
}


#' See summary()
#' Returns a data.frame. Each column is power of a predictor. A column for threshold.
summary_one_sim_one_thres <- function(sim, crit, thres, digits, how, pivot=FALSE) {
  res <- sim %>%
    purrr::map(crit) %>%
    purrr::map(function(e) {
      if (how == "greater") {e >= thres}
      else {e <= thres}}) %>%
    abind::abind(along = -1) %>%
    tibble::as_tibble() %>%
    dplyr::summarise(across(everything(), mean)) %>% round(digits) %>%  #TODO: correct colnames???
    dplyr::mutate(thres = thres)
  if (pivot) {
    res <- res %>%
      tidyr::pivot_longer(-thres, names_to = "test", values_to = "power")
  }
  return(res)
}

