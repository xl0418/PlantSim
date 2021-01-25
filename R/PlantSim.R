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
#' @return The abundance matrix for all plots and species through time.
#' @examples
#' PlantSim(nplot = 5, nspe = 2, t = 10, ini_abundance = matrix(1:10, 5, 2),
#'  growth_rate = 1, interaction_matrix = matrix(0.001, 2, 2), st_portion = 0.7, surv_rate = 1)
#' @export
#'
PlantSim <-
  function(nplot,
           nspe,
           t,
           ini_abundance,
           growth_rate,
           interaction_matrix,
           st_portion,
           surv_rate,
           filesave=NULL) {
    if (is.matrix(ini_abundance) && all(dim(ini_abundance) == c(nplot, nspe))) {
    } else {
      return (print(
        "Provide the initial abundances that is at the dim of c(nplot, nspe)."
      ))
    }


    if (is.null(dim(growth_rate))) {
      growth_rate <- matrix(growth_rate, nrow = nplot, ncol = nspe)
    } else if (any(dim(growth_rate) != c(nplot, nspe))) {
      return (
        print(
          "The dim of the growth rate is not equal to the dim of the number of plots X the number of species."
        )
      )
    }

    if (is.null(dim(surv_rate))) {
      surv_rate <- matrix(surv_rate, nrow = nplot, ncol = nspe)
    } else if (any(dim(surv_rate) != c(nplot, nspe))) {
      return (
        print(
          "The dim of the survival rate is not equal to the dim of the number of plots X the number of species."
        )
      )
    }

    if (dim(interaction_matrix)[1] != nspe) {
      return(
        print("Wrong dim for the interaction matrix.")
      )
    }
    # initialize the community matrix
    plot_abundance <- array(0, dim = c(nplot, nspe, t))
    plot_abundance[, , 1] <- ini_abundance

    for (ts in c(1:(t - 1))) {
      new_seedlings <- matrix(0, nrow = nplot, ncol = nspe)
      for (plo in c(1:nplot)) {
        for (spe in c(1:nspe)) {
          new_seedlings[plo, spe] <-
            plot_abundance[plo, spe, ts] * growth_rate[plo, spe] * exp(1 - plot_abundance[plo, , ts] %*% interaction_matrix[spe,])
          if (is.nan(new_seedlings[plo, spe])) {
            print("Overflow numbers generated!")
            break
          }
          }
      }
      update_seedlings <- matrix(0, nrow = nplot, ncol = nspe)
      for (plo in c(1:nplot)) {
        update_seedlings[plo,] <-
          st_portion * new_seedlings[plo, ] + (1 - st_portion) * (colSums(new_seedlings) - new_seedlings[plo,]) / (nplot - 1)
      }

      plot_abundance[, , ts + 1] <- surv_rate * update_seedlings
    }

    if (is.null(filesave)){}
    else {save(plot_abundance, file = filesave)}
    return(plot_abundance)
  }
