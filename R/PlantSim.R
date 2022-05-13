#' Simulate a plant community
#'
#' @param nplot The number of the plots in the community
#' @param nspe The number of the species
#' @param t The total time steps
#' @param ini_abundance The initial abundance of the species
#' @param growth_rate The growth rate for each plot and each species. It could be universally same.
#' @param interaction_matrix The interaction matrix among species
#' @param st_portion The portion of the seedlings that stay at the parental plot.
#' @param filesave The file name to be saved.
#' @param vary_k Set vary_k be a vector of two entries providing the mean and the sd. Otherwise, it is FALSE giving the same k = 1 to all plots.
#' @param distribution Distribution mode. "uniform" stands for the uniform dispersal; "random" stands for the random dispersal; "Gaussian" stands for the Gaussian dispersal.
#' @param sig_disp The variance of the dispersal for the Gaussian dispersal kernel.
#' @param model Specify local population models, i.e. Ricker, BH, PL, for the Ricker model, the Beverton-Holt model and the power law model.
#' @param boundary If TRUE, the boundary effect is considered. Otherwise, the grid is a torus.
#' @param cell_kill_rate The percentage of cells eliminates plant individuals.
#' @return The abundance matrix for all plots and species through time.
#' @importFrom stats rnorm rpois runif var
#' @examples
#' plantsim(nplot = 5, nspe = 2, t = 10, ini_abundance = matrix(1:10, 5, 2),
#'  growth_rate = 1, interaction_matrix = matrix(0.001, 2, 2), st_portion = 0.7)
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
           filesave = NULL,
           vary_k = FALSE,
           distribution = "uniform",
           model = "Ricker",
           sig_disp = 1,
           boundary = FALSE,
           cell_kill_rate = 0) {
    # Check if the initial abundance is compatible with the dim (nplot, nspe) or a constant that applies to every plot and species
    if (is.matrix(ini_abundance) &&
        all(dim(ini_abundance) == c(nplot, nspe))) {

    } else {
      return (print(
        "Provide the initial abundances that is at the dim of c(nplot, nspe)."
      ))
    }

    if (model == "PL") {
      if (nspe > 1) {
        return (print(
          "Only one species is simulated under the power law model."
        ))
      }
    } else {

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

    # sample the cells that will eliminate individuals
    num_cells_kill <- ceiling(nplot * cell_kill_rate)
    sample_cells <- sample(c(1:nplot), size = num_cells_kill)

    # initialize the community matrix
    plot_abundance <- array(0, dim = c(nplot, nspe, t + 1))
    stay_seeds <- array(0, dim = c(nplot, nspe, t))
    dis_seeds <- array(0, dim = c(nplot, nspe, t))

    dispersal_seeds <- array(0, dim = c(nplot, nspe, t))
    seeds_before_disp <- array(0, dim = c(nplot, nspe, t + 1))
    plot_abundance[, , 1] <- round(ini_abundance)
    seeds_before_disp[, , 1] <- round(ini_abundance)
    if (distribution == "Gaussian") {
      x <- sqrt(nplot)
      y <- sqrt(nplot)
      if (x - floor(x) != 0) {
        return(print("Pls make nplot = x^2..."))
      }
      #generate a grid of coordinates
      xy <- expand.grid(seq(1, x, by =1), seq(1, y, by =1))



      if (boundary){
        #calculating the distance between all populations.

        dmat <- sqrt(outer(xy[,1], xy[,1], "-")^2+outer(xy[,2], xy[,2], "-")^2)

        Dmat <- exp(-(dmat/sig_disp)^2) #gaussian seed dispersal kernel with variance sig_disp
        Dmat <- sweep(Dmat, 1, apply(Dmat, 2, sum), "/")
      } else {
        dmat <- as.matrix(som.nn::dist.torus(xy))
        Dmat <- exp(-(dmat/sig_disp)^2)
        Dmat <- sweep(Dmat, 1, apply(Dmat, 2, sum), "/")
      }

      st_portion_gaussian <- Dmat[1,1]

    }

    # Ricker model
    for (tt in c(1:t)) {
      # initialize the new seeds gain matrix
      new_seeds <- matrix(0, nrow = nplot, ncol = nspe)
      # Producing the seeds for each spec and in each plot following the Ricker model
      # and use Poisson draw to calculate out how many seeds stay.
      for (plo in c(1:nplot)) {
        for (spe in c(1:nspe)) {
          if (model == "Ricker"){
            temp_seeds <-
              plot_abundance[plo, spe, tt] * exp(growth_rate[plo, spe] - plot_abundance[plo, , tt] %*% interaction_matrix[spe, ] / k[nplot])
          } else if (model == "PL") {
            temp_seeds <-
               growth_rate[plo, spe] * plot_abundance[plo, spe, tt] ** interaction_matrix[1, 1]
          } else if (model == "BH") {
            temp_seeds <-
              plot_abundance[plo, spe, tt] * growth_rate[plo, spe] / (1 + plot_abundance[plo, , tt] %*% interaction_matrix[spe, ] / k[nplot])
          } else {
            return(print("Pls specify a local population model..."))
          }

          new_seeds[plo, spe] <- temp_seeds
          if (is.nan(new_seeds[plo, spe])) {
            print("Overflow numbers generated!")
            break
          }
        }
      }


      if (distribution == "Gaussian") {
        for (spe.id in c(1:nspe)) {
          total_pop <-as.vector(new_seeds[, spe.id] %*% Dmat)
          stay_seeds[, spe.id, tt] <- rpois(nplot, new_seeds[, spe.id] * st_portion_gaussian)
          dispersal_seeds[, spe.id, tt] <- rpois(nplot, pmax(total_pop - stay_seeds[, spe.id, tt], 0))
          seeds_before_disp[, spe.id, tt] <- new_seeds[, spe.id]
          plot_abundance[, spe.id, tt + 1] <- stay_seeds[, spe.id, tt] + dispersal_seeds[, spe.id, tt]
        }
      } else {
        # initialize the update_seeds matrix
        update_seeds <- matrix(0, nrow = nplot, ncol = nspe)

        # stay seeds and dispersal seeds
        stay_seeds[, , tt] <-
          rpois(length(new_seeds), st_portion * new_seeds)


        dis_seeds[, , tt] <-
          rpois(length(new_seeds), (1 - st_portion) * new_seeds)

        seeds_before_disp[, , tt] <-
          stay_seeds[, , tt] + dis_seeds[, , tt]

        # the seeds rain for each species by Poisson draws
        seeds_rain <- colSums(dis_seeds[, , tt, drop = FALSE])

        # dispersal seeds
        if (distribution == "uniform") {
          prob_dis = rep(1 / nplot, nplot)
        } else {
          probs <- runif(nplot, 0, 1)
          prob_dis = probs / sum(probs)
        }

        for (spe_ind in c(1:nspe)) {
          dispersal_seeds[, spe_ind, tt] <-
            if (seeds_rain[spe_ind] < .Machine$integer.max) {
              stats::rmultinom(1, seeds_rain[spe_ind], prob = prob_dis)
            }
          else {
            rep(seeds_rain[spe_ind] / nplot, nplot) # large integer check
          }
        }
        # seeds rain joins the local seeds
        update_seeds <-
          stay_seeds[, , tt] + dispersal_seeds[, , tt]

        # kill individuals in bad cells
        update_seeds[sample_cells, ] <- 0

        # apply survival rate to seeds
        plot_abundance[, , tt + 1] <- update_seeds
      }

    }

    if (is.null(filesave)) {

    }
    else {
      save(plot_abundance, file = filesave)
    }
    return(
      list(
        all = plot_abundance[, , 1:t, drop = FALSE],
        stay = stay_seeds,
        dispersal = dispersal_seeds,
        bef_dispersal = seeds_before_disp[, , 1:t, drop = FALSE],
        st_rate = ifelse(distribution == "Gaussian", st_portion_gaussian, st_portion)
      )
    )
  }
