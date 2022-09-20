#' Power curve using Monte Carlo simulation
#'
#' This function can be used to create power curves by calling sim_power() on
#' combinations of many sample sizes and signal-to-noise ratio (SNR).
#'
#' @param xmod A MixtureModel object.
#' @param ymod One or a list of OutcomeModel object(s).
#' @param imod An InferenceModel object.
#' @param s An integer for the number of Monte Carlo simulations.
#' @param n An integer or a vector of sample sizes.
#' @param cores An integer for the number of processing cores. When cores > 1,
#'   parallelism is automatically applied.
#' @param file A string, a file name with no extension to write samples to
#'   periodically. By default write to an RDS file.
#' @param errorhandling A string "remove", "stop", or "pass". If an error occurs
#'   in any iteration, remove that iteration (remove), return the error message
#'   verbatim in the output (pass), or terminate the loop (stop). Default is
#'   "remove". See R package 'foreach' for more details.
#' @param snr_iter An integer for number of Monte Carlo samples to estimate SNR.
#' @param cluster_export A vector of functions to pass to the
#'   parallel-processing clusters.
#' @return A SimCurve object with the following attributes:
#'   \item{s}{a number of simulations.}
#'   \item{snr}{a real number or array of real numbers for SNR of each OutcomeModel.}
#'   \item{n}{a number or vector of sample sizes.}
#'   \item{xmod}{the MixtureModel used.}
#'   \item{ymod}{the OutcomeModel used.}
#'   \item{imod}{the InferenceModel used.}
#'   \item{sims}{a list of simulation output matrices.}
#' @examples
#' data("nhanes1518")
#' chems <- c("URXCNP", "URXCOP", "URXECP", "URXHIBP", "URXMBP", "URXMC1",
#' "URXMCOH", "URXMEP","URXMHBP", "URXMHH", "URXMHNC", "URXMHP", "URXMIB",
#' "URXMNP", "URXMOH", "URXMZP")
#' chems_mod <- mpower::MixtureModel(nhanes1518[, chems], method = "resampling")
#' bmi_mod <- mpower::OutcomeModel(f = "0.2*URXCNP + 0.15*URXECP +
#' 0.1*URXCOP*URXECP", family = "binomial")
#' logit_mod <- mpower::InferenceModel(model = "glm", family = "binomial")
#' logit_out <- mpower::sim_curve(xmod=chems_mod, ymod=bmi_mod, imod=logit_mod,
#' s=20, n=c(500, 1000), cores=2, snr_iter=1000)
#' logit_df <- summary(logit_out, crit="pval", thres=0.05, how="lesser")
#' plot_summary(logit_out, crit="pval", thres=0.05, how="lesser")
#' @export
sim_curve <- function(xmod,
                      ymod,
                      imod,
                      s = 100,
                      n = 100,
                      cores = 1,
                      file = NULL,
                      errorhandling = "stop",
                      snr_iter = 10000,
                      cluster_export = c()
                      ) {
  stopifnot(inherits(xmod, "mpower_estimation_MixtureModel") |
              inherits(xmod, "mpower_cvine_MixtureModel") |
              inherits(xmod, "mpower_resampling_MixtureModel"))
  stopifnot(inherits(imod, "mpower_InferenceModel"))
  if (inherits(ymod, "mpower_OutcomeModel")) {
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
      message(paste("Simulation for n =", cur_n, "and the", j, "th outcome model"))
      out <- sim_power(xmod = xmod,
                       ymod = cur_ymod,
                       imod = imod,
                       s = s,
                       n = cur_n,
                       snr_iter = snr_iter,
                       cores = cores,
                       file = cur_file,
                       errorhandling = errorhandling,
                       cluster_export = cluster_export)
      sim_list[[k]] <- out$sims
      ymod_list[[k]] <- out$ymod
      snr_vec[j] <- snr_vec[j] + out$snr
      n_vec <- c(n_vec, out$n)
      k <- k+1
    }
  }
  new_SimCurve(list(s = s, n = n_vec,
                    snr = rep(snr_vec/length(n), length(n)),
                    xmod = xmod,
                    ymod = ymod_list,
                    imod = imod,
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
#' @param s An integer for number of Monte Carlo simulations.
#' @param n An integer for sample size in each simulation.
#' @param cores An integer for number of processing cores. When cores > 1,
#'   parallelism is automatically applied.
#' @param file A string, a file name with no extension to write samples to
#'   periodically. By default write to an RDS file.
#' @param errorhandling A string "remove", "stop", or "pass". If an error occurs
#'   in any iteration, remove that iteration (remove), return the error message
#'   verbatim in the output (pass), or terminate the loop (stop). Default is
#'   "remove". See R package `foreach` for more details.
#' @param snr_iter An integer for number of Monte Carlo samples to estimate SNR
#' @param cluster_export A vector of functions to pass to the parallel-processing clusters
#' @return A PowerSim object. Attributes:
#'   \item{s}{a number of simulations.}
#'   \item{snr}{a real number for SNR of the OutcomeModel.}
#'   \item{n}{a number of sample sizes.}
#'   \item{xmod}{the MixtureModel used.}
#'   \item{ymod}{the OutcomeModel used.}
#'   \item{imod}{the InferenceModel used.}
#'   \item{sims}{an output matrices.}
#' @examples
#' data("nhanes1518")
#' chems <- c("URXCNP", "URXCOP", "URXECP", "URXHIBP", "URXMBP", "URXMC1",
#' "URXMCOH", "URXMEP","URXMHBP", "URXMHH", "URXMHNC", "URXMHP", "URXMIB",
#' "URXMNP", "URXMOH", "URXMZP")
#' chems_mod <- mpower::MixtureModel(nhanes1518[, chems], method = "resampling")
#' bmi_mod <- mpower::OutcomeModel(f = "0.2*URXCNP + 0.15*URXECP +
#' 0.1*URXCOP*URXECP", family = "binomial")
#' logit_mod <- mpower::InferenceModel(model = "glm", family = "binomial")
#' logit_out <- mpower::sim_power(xmod=chems_mod, ymod=bmi_mod, imod=logit_mod,
#' s=100, n=2000, cores=2, snr_iter=2000)
#' logit_df <- summary(logit_out, crit="pval", thres=0.05, how="lesser")
#' plot_summary(logit_out, crit="pval", thres=0.05, how="lesser")
#' @export
sim_power <- function(xmod,
                      ymod,
                      imod,
                      s = 100,
                      n = 100,
                      cores = 1,
                      file = NULL,
                      errorhandling = "stop",
                      snr_iter = 10000,
                      cluster_export = c()
                      ) {
  stopifnot(inherits(xmod, "mpower_estimation_MixtureModel") |
              inherits(xmod, "mpower_cvine_MixtureModel") |
              inherits(xmod, "mpower_resampling_MixtureModel"))
  stopifnot(inherits(ymod, "mpower_OutcomeModel"))
  stopifnot(inherits(imod, "mpower_InferenceModel"))
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
    cl <- snow::makeCluster(min(cores_max-1, cores), type="SOCK")
    snow::clusterExport(cl, c(c("genx", "geny", "fit"), cluster_export))
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress=progress)
    out <- foreach::foreach(i = 1:s,
                            .errorhandling = errorhandling,
                            .options.snow = opts) %dopar% {
      X <- genx(xmod, n)
      y <- geny(ymod, X)
      r <- fit(imod, X, y)
      ## write results to disk
      if (!is.null(file)) saveRDS(r, file = paste0(file, i, ".rds"))
      return(r)
    }
    snow::stopCluster(cl)
  }
  close(pb)
  snr <- estimate_snr(ymod, xmod, m = snr_iter)$est
  new_Sim(list(s = s,
               n = n,
               snr = snr,
               xmod = xmod,
               ymod = ymod,
               imod = imod,
               sims = out))
}

new_Sim <- function(x = list()) {
  stopifnot(is.list(x))
  structure(x, class = "mpower_Sim")
}

#' Tabular summaries of power simulation
#' @param sim A Sim or a SimCurve object, output from `sim_power()` or
#'   `sim_curve()`.
#' @param crit A string specifying the significance criteria.
#' @param thres A number or vector of numbers specifying the thresholds of
#'   "significance".
#' @param digits An integer for the number of decimal points to display.
#' @param how A string, whether to compare the criterion 'greater' or 'lesser'
#'   than the threshold.
#' @return A data.frame summary of power for each predictor for each combination
#'   of thresholds, sample size, signal-to-noise ratios.
#' @export
summary <- function(sim, crit, thres, digits = 3, how = "greater") {
  UseMethod("summary", sim)
}

#' Plot summaries of power simulation
#' @param sim A Sim or a SimCurve object, output from `sim_power()` or
#'   `sim_curve()`.
#' @param crit A string specifying the significance criteria.
#' @param thres A number or vector of numbers specifying the thresholds of
#'   "significance".
#' @param digits An integer for the number of decimal points to display.
#' @param how A string, whether to compare the criterion 'greater' or 'lesser'
#'   than the threshold.
#' @return A 'ggplot2' graphics.
#' @export
plot_summary <- function(sim, crit, thres, digits = 3, how = "greater") {
  UseMethod("plot_summary", sim)
}

#' @export
plot_summary.mpower_Sim <- function(sim, crit, thres, digits = 3, how = "greater") {
  if(seq_along(crit) > 1)  stop("Summary limited to one criteria", call. = FALSE)
  if (length(thres) == 1) {
    res <- summary_one_sim_one_thres(sim$sims, crit, thres, digits, how, pivot = T)
  } else if (length(thres) > 1) {
    # For many threshold and one Sim object
    res <- summary_one_sim_many_thres(sim$sims, crit, thres, digits, how)
  }
  g <- res %>%
    ggplot(aes(x = !! sym("thres"), y = !! sym("power"),
               group = !! sym("test"), color = !! sym("test"))) +
    geom_hline(yintercept = 0.8, linetype = 2, colour = "gray") +
    geom_point() + geom_path() +
    labs(y = "Type I error rate/ Power", x = "Significance threshold") +
    scale_colour_discrete(name = "")
  return(g)
}

#' @export
plot_summary.mpower_SimCurve <- function(sim, crit, thres, digits = 3, how = "greater") {
  if (length(thres) == 1) {
    # For one threshold and many Sim objects
    res <- summary_many_sim_one_thres(sim$sims, sim$n, sim$snr, crit, thres, digits, how)
  } else {
    warning("Using the first threshold only")
    res <- summary_many_sim_one_thres(sim$sims, sim$n, sim$snr, crit, thres[1], digits, how)
  }
  if (length(unique(res$snr)) > 1) {
    g <- res %>%
      dplyr::mutate(snr = as.factor(!! sym("snr"))) %>%
      ggplot(aes(x = !! sym("n"), y = !! sym("power"), colour = !! sym("snr"))) +
      geom_hline(yintercept = 0.8, linetype = 2, colour = "grey") +
      geom_path() + geom_point() +
      labs(x="Sample size", y="Type I error rate/ Power") +
      scale_colour_discrete(name = "SNR") +
      facet_wrap(stats::as.formula("~test"))
  } else {
    g <- res %>%
      ggplot(aes(x = !! sym("n"), y = !! sym("power"), colour = !! sym("test"))) +
      geom_hline(yintercept = 0.8, linetype = 2, colour = "grey") +
      geom_path() + geom_point() +
      labs(x="Sample size", y="Type I error rate/ Power") +
      scale_colour_discrete(name = " ")
  }
  return(g)
}

#' @export
summary.mpower_Sim <- function(sim, crit, thres, digits = 3, how = "greater") {
  if(seq_along(crit) > 1)  stop("Summary limited to one criteria", call. = FALSE)
  if (length(thres) == 1) {
    # For one threshold and one Sim object
    res <- summary_one_sim_one_thres(sim$sims, crit, thres, digits, how, pivot = T)
  } else if (length(thres) > 1) {
    # For many threshold and one Sim object
    res <- summary_one_sim_many_thres(sim$sims, crit, thres, digits, how)
  }
  colnames(res) <- c("thres", " ", "power")
  cat("\n\t*** POWER ANALYSIS SUMMARY ***")
  cat("\nNumber of Monte Carlo simulations:", sim$s)
  cat("\nNumber of observations in each simulation:", sim$n)
  cat("\nData generating process estimated SNR:", round(sim$snr,2))
  cat("\nInference model:", sim$imod$model_name)
  cat("\nSignificance criterion:", crit)
  for (tt in thres) {
    cat("\n\nSignificance threshold: ", tt)
    cat("\n", res %>% dplyr::filter(thres == tt) %>%
          dplyr::select(-thres) %>% knitr::kable(), sep="\n")
  }
  return(res)
}

#' @export
summary.mpower_SimCurve <- function(sim, crit, thres, digits = 3, how = "greater") {
  if (length(thres) == 1) {
    # For one threshold and many Sim objects
    res <- summary_many_sim_one_thres(sim$sims, sim$n, sim$snr, crit, thres, digits, how)
  } else {
    warning("Using the first threshold only")
    res <- summary_many_sim_one_thres(sim$sims, sim$n, sim$snr, crit, thres[1], digits, how)
  }
  cat("\n\t*** POWER CURVE ANALYSIS SUMMARY ***")
  cat("\nNumber of Monte Carlo simulations:", sim$s)
  cat("\nNumber of observations in each simulation:", unique(sim$n))
  cat("\nData generating process estimated SNR (for each outcome model):",
      round(unique(sim$snr),2))
  cat("\nInference model:", sim$imod$model_name)
  cat("\nSignificance criterion:", crit)
  cat("\nSignificance threshold:", thres[1])
  cat("\n", res %>%
        dplyr::select(tidyselect::all_of(c("test", "power", "n", "snr"))) %>%
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
    dplyr::mutate(n = n, snr = round(snr, 2)) %>%
    tidyr::pivot_longer(cols = -c("n", "snr", "thres"), names_to = "test", values_to = "power")
  return(res)
}

summary_one_sim_many_thres <- function(sim, crit, thres, digits, how) {
  res <- data.frame()
  for (i in seq_along(thres)) {
    res <- dplyr::bind_rows(res,
                            summary_one_sim_one_thres(sim, crit, thres[i], digits, how, pivot = T))
  }
  return(res)
}

# Returns a data.frame. Each column is power of a predictor. A column for threshold.
summary_one_sim_one_thres <- function(sim, crit, thres, digits, how, pivot=FALSE) {
  res <- sim %>%
    purrr::map(crit) %>%
    purrr::map(function(e) {
      if (how == "greater") {e >= thres}
      else {e <= thres}}) %>%
    abind::abind(along = -1) %>%
    tibble::as_tibble() %>%
    dplyr::summarise(dplyr::across(tidyselect::everything(), mean)) %>%
    round(digits) %>%
    dplyr::mutate(thres = thres)
  if (pivot) {
    res <- res %>%
      tidyr::pivot_longer(-thres, names_to = "test", values_to = "power")
  }
  return(res)
}
