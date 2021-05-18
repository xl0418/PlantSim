#' Simulate a plant community
#'
#' @param nplot The number of the plots in the community
#' @param nspe The number of the species
#' @param t The total time steps
#' @param ini_abundance The initial abundance of the species
#' @param growth_rate The growth rate for each plot and each species. It could be universally same.
#' @param interaction_matrix The interaction matrix among species
#' @param st_portion The portion of the seedlings that stay at the parental plot.
#' @param surv_rate The survival rate for each plot and each species. It could be universally same.
#' @param filesave The file name to be saved.
#' @param vary_k Set vary_k be a vector of two entries providing the mean and the sd. Otherwise, it is FALSE giving the same k = 1 to all plots.
#' @param kill_rate Kill all plants in a percentage of plots.
#' @return The abundance matrix for all plots and species through time.
#' @importFrom stats rnorm rpois runif var
#' @examples
#' plantsim(nplot = 5, nspe = 2, t = 10, ini_abundance = matrix(1:10, 5, 2),
#'  growth_rate = 1, interaction_matrix = matrix(0.001, 2, 2), st_portion = 0.7, surv_rate = 1)
#' @export
#'
plantsim <-
  function(nplot,
           nspe,
           t,
           ini_abundance,
           growth_rate,
           interaction_matrix,
           st_portion,
           surv_rate,
           filesave = NULL,
           vary_k = FALSE,
           kill_rate = 0.05) {
    # Check if the initial abundance is compatible with the dim (nplot, nspe) or a constant that applies to every plot and species
    if (is.matrix(ini_abundance) &&
        all(dim(ini_abundance) == c(nplot, nspe))) {

    } else {
      return (print(
        "Provide the initial abundances that is at the dim of c(nplot, nspe)."
      ))
    }

    # Check the growth rate dimension. Either universally same or customized by users
    if (is.null(dim(growth_rate))) {
      growth_rate <- matrix(growth_rate, nrow = nplot, ncol = nspe)
    } else if (any(dim(growth_rate) != c(nplot, nspe))) {
      return (
        print(
          "The dim of the growth rate is not equal to the dim of the number of plots X the number of species."
        )
      )
    }
    # Check the survival rate dimension.
    if (is.null(dim(surv_rate))) {
      surv_rate <- matrix(surv_rate, nrow = nplot, ncol = nspe)
    } else if (any(dim(surv_rate) != c(nplot, nspe))) {
      return (
        print(
          "The dim of the survival rate is not equal to the dim of the number of plots X the number of species."
        )
      )
    }
    # Check the dimension of the interaction matrix
    if (dim(interaction_matrix)[1] != nspe) {
      return(print("Wrong dim for the interaction matrix."))
    }

    # Set a variable K for each plot
    if (length(vary_k) == 2) {
      k = round(rnorm(nplot, mean = vary_k[1], sd = vary_k[2]))
    } else {
      k = rep(1, nplot)
    }

    # initialize the community matrix
    plot_abundance <- array(0, dim = c(nplot, nspe, t+1))
    stay_seeds <- array(0, dim = c(nplot, nspe, t))
    dispersal_seeds <- array(0, dim = c(nplot, nspe, t))
    seeds_before_disp <- array(0, dim = c(nplot, nspe, t))
    plot_abundance[, , 1] <- round(ini_abundance)

    # Ricker model
    for (ts in c(1:t )) {
      # initialize the new seeds gain matrix
      new_seeds <- matrix(0, nrow = nplot, ncol = nspe)
      # Producing the seeds for each spec and in each plot following the Ricker model
      # and use Poisson draw to calculate out how many seed are produced.
      for (plo in c(1:nplot)) {
        for (spe in c(1:nspe)) {
          temp_seeds <-
              plot_abundance[plo, spe, ts] * growth_rate[plo, spe] * exp(1 - plot_abundance[plo, , ts] %*% interaction_matrix[spe,] / k[nplot])
          new_seeds[plo, spe] <- temp_seeds
          if (is.nan(new_seeds[plo, spe])) {
            print("Overflow numbers generated!")
            break
          }
        }
      }

      seeds_before_disp[,,ts] <- new_seeds
      # initialize the update_seeds matrix
      update_seeds <- matrix(0, nrow = nplot, ncol = nspe)

      # stay seeds and dispersal seeds
      stay_seeds[,,ts] <-  rpois(length(new_seeds), st_portion * new_seeds)

      dis_seeds <- (st_portion != 1) * pmax(new_seeds - stay_seeds[,,ts], 0)

      # the seeds rain for each species by Poisson draws
      seeds_rain <- colSums(dis_seeds)

      # apply survival rate to the seeds rain
      actual_seeds_rain <- matrix(0, nrow = nplot, ncol = nspe)
      for (col_ind in c(1:nspe)) {
        actual_seeds_rain[ , col_ind] <- surv_rate[, col_ind] * seeds_rain[col_ind] / nplot
      }

      # kill plants in selected plots
      kill_plots <-
        sort(sample(x = c(1:nplot), size = round(kill_rate * nplot)))
      stay_seeds[kill_plots,,ts] <- 0
      # dispersal seeds
      for (spe_ind in c(1:nspe)) {
        dispersal_seeds[, spe_ind, ts] <- rpois(nplot, actual_seeds_rain[spe_ind])
      }

      # seeds rain joins the local seeds
      update_seeds <-
        stay_seeds[,,ts] + dispersal_seeds[,,ts]
      # apply survival rate to seeds
      plot_abundance[, , ts + 1] <- update_seeds
    }

    if (is.null(filesave)) {

    }
    else {
      save(plot_abundance, file = filesave)
    }
    return(list(all = plot_abundance[,,1:t, drop = FALSE],
                stay = stay_seeds,
                dispersal = dispersal_seeds,
                bef_dispersal = seeds_before_disp))
  }
