#' Inferring the model parameters
#'
#' @param data The data is fitted to infer the model parameters.
#' @param ngener The number of the generations of the algorithm.
#' @param nsample The number of the samples for each generation.
#' @param acc_rate The acceptance rate that determines the best portion of the samples to be accepted.

#' @return The abundance matrix for all plots and species through time.
#' @examples
#' plantsim(nplot = 5, nspe = 2, t = 10, ini_abundance = matrix(1:10, 5, 2),
#'  growth_rate = 1, interaction_matrix = matrix(0.001, 2, 2), st_portion = 0.7, surv_rate = 1)
#' @export
#'

ParaInfer <- function(data, ngener, nsample, acc_rate = 0.25) {
  load(data)
  dim_data = dim(plot_abundance)
  nplot = dim_data[1]
  nspe = dim_data[2]
  nt = dim_data[3]

  nparas <- nspe ** 2
  paras <- array(0, dim = c(ngener, nsample, nparas))
  fitness <- array(1, dim = c(ngener, nsample))

  for (para_ind in c(1:nparas)) {
    paras[1, , para_ind] <- runif(nsample, 0, 1)
  }


  for (gener_ind in c(2:ngener)) {
    print(paste0("Iteration ", gener_ind, "..."))
    # progress bar
    pb <-
      progress::progress_bar$new(format = "[:bar] :current/:total (:percent)", total = nsample)
    pb$tick(0)
    Sys.sleep(3)
    for (sample_ind in c(1:nsample)) {
      sim_result <- plantsim(
        nplot = nplot,
        nspe = nspe,
        t = nt,
        ini_abundance = runif(nspe, 1, 100),
        growth_rate = 1,
        interaction_matrix = matrix(paras[gener_ind -
                                            1, sample_ind, ], nspe, nspe),
        st_portion = 0.7,
        surv_rate = 1
      )
      fitness[gener_ind, sample_ind] <-
        1 / sqrt(sum((sim_result - plot_abundance) ** 2))
      pb$tick(1)
      Sys.sleep(1 / 100)
    }

    # pick up the best fit simulations
    best_fit_index <-
      sort(fitness[gener_ind, ], index.return = TRUE)$ix
    best_fit_index_rate <-
      best_fit_index[floor(nsample * (1 - acc_rate) + 1):nsample]

    best_fit_paras_rate <- paras[gener_ind - 1, best_fit_index_25, ,drop = FALSE]
    best_fit_rate <- fitness[gener_ind, best_fit_index_25]
    weight_fit <- best_fit_rate / sum(best_fit_rate)

    # means and vars of the best fit paras
    #mean_best_fit_paras <- colMeans(best_fit_paras_rate)
    var_best_fit_paras <- apply(best_fit_paras_rate, 1, var)

    # sample all paras from the best fit paras according to their weights
    sample_paras <- NULL
    for (para_ind in c(1:nparas)) {
      sample_paras <-
        cbind(
          sample_paras,
          sample(
            best_fit_paras_rate[,, para_ind],
            nsample,
            prob = weight_fit,
            replace = TRUE
          )
        )
    }

    # disturb the samples a bit
    proposed_paras <- NULL
    for (para_ind in c(1:nparas)) {
      proposed_paras <- cbind(proposed_paras,
                              unlist(
                                lapply(sample_paras[, para_ind], rnorm, n = 1, sd = var_best_fit_paras[para_ind])
                              ))
    }

    # use the new proposed samples
    paras[gener_ind, , ] <- proposed_paras
  }
}
