#' @title Derive the beta priors on all of the primary event probabilities
#' @description Elicit the beta prior distributions for each primary event probability by the
#' pairwise comparison approach.
#' @param cornerstone_event_index: the index of the primary event that is the cornerstone event.
#' This had to be defined when you defined the tree structure in the file.
#' @param cornerstone_event_interval: the range over which the cornerstone event probability is elicited from the expert.
#' It's mapped to the central 95\% probability interval of a beta distribution.
#' @param pairwise_matrix: the matrix of pairwise comparisons. A square matrix whose dimension is the number of primary events.
#' If a pairwise comparison is not done then the entry is NA.
#' @return A matrix of dimension number of primary events by 2; each row is the
#' @examples beta_parameters <- derive_primary_event_priors(cornerstone_event_index = tree$corner,
#'                                             cornerstone_event_interval  = c(0.01,0.05),
#'                                             pairwise_matrix             = matrix(c(1,   1/4,  2, 1/2,
#'                                                                                   4, 1,  2, 4,
#'                                                                                   1/2, 1/2, 1, 1/5,
#'                                                                                  (2,   1/4, 5, 1), nrow=4,byrow=TRUE))
#' Or here is an example of an incomplete comparison matrix that you could also use as the pairwise_matrix
#' argument in the derive_primary_event_priors function:
#' pairwise_matrix = matrix(c(1,   1/4,  2,  1/2,
#'                            4,   1,    NA, NA,
#'                            1/2, NA,   1,  NA,
#'                            2,   NA,   NA, 1), nrow=4,byrow=TRUE)
#' @export derive_primary_event_prior

# #########################################################################################################
# These are functions that derive the beta priors on all of the primary event probabilities
# #########################################################################################################

# #########################################################################################################
# derive_primary_event_priors.R
# #########################################################################################################
# Takes in a matrix of pairwise comparisons and the cornerstone event upper and lower bounds as an input
# and computes the beta parameters for the other primary events.
# Inputs:
# cornerstone_event_index: integer, the index, from 1 to the number of primary events, of the cornerstone event
# cornerstone_event_interval: a 2-vector, the upper and lower probability bounds for the cornerstone event, to be mapped to the
#                   central 95% probability of a beta distribution
# pairwise_matrix: a square matrix of pairwise comparisons.  Assumes either NA or -1 for comparisons that are not made
# Output:
# beta_params: a matrix of dimension number of primary events by 2.

derive_primary_event_priors <- function(cornerstone_event_index,cornerstone_event_interval, pairwise_matrix) {
  # Initializes the number of primary events
  if (length(dim(pairwise_matrix)[1])==0) {
    n_primary_events <- 1
  } else {
    n_primary_events <- dim(pairwise_matrix)[1]
  }
  # define beta_shape matrix
  beta_shape <- matrix(data=NA,nrow=n_primary_events,ncol=2)
  # Derive the parameters for the cornerstone event
  beta_shape[cornerstone_event_index,] <- find_beta_parameters(cornerstone_event_interval)
  # Derive the event weights from the pairwise comparison matrix
  weights <- derive_weights(pairwise_matrix)
  # cat(" Weights are ",weights)
  # Derive the parameters for the other events
  non_cornerstone_events <- 1:n_primary_events
  non_cornerstone_events <- non_cornerstone_events[-cornerstone_event_index]
  for (i in non_cornerstone_events) {
    beta_shape[i,] <- find_beta_parameters(c(weights[i]*cornerstone_event_interval[1]/weights[cornerstone_event_index],
                                             min(c(0.9999999,weights[i]*cornerstone_event_interval[2]/weights[cornerstone_event_index]))))
  }
  return(beta_shape)
}
# #########################################################################################################

# #########################################################################################################
# find_beta_parameters.R
# #########################################################################################################
# Takes an interval in (0,1) and identifies the beta distribution with those parameters.
# Restricted to both beta parameters > 1
find_beta_parameters <- function(interval) {
  if ((min(interval) <= 0.0) | (max(interval) >= 1.0)) {
    cat("WARNING: find_beta_parameters: interval is outside range (0,1).\n")
  }
  if (interval[1] >= interval[2]) {
    cat("WARNING: find_beta_parameters: lower bound not less than upper bound.\n ")
  }
  shape <- log(log(c(2,2)))
  objective_function <- function(shape,bounds) {
    return(log(1+(bounds[1]-qbeta(0.025,exp(exp(shape[1])),exp(exp(shape[2]))))^2 +
                 (bounds[2]-qbeta(0.975,exp(exp(shape[1])),exp(exp(shape[2]))))^2))
  }
  beta_shape <- optim(shape,objective_function, bounds=interval, method="BFGS")
  return(exp(exp(beta_shape$par)))
}
# #########################################################################################################

# #########################################################################################################
# derive_weights.R
# #########################################################################################################
# Takes in a matrix of pairwise comparisons and produces weights that are the normalised geometric mean of valid
# entries of each row in the matrix.  Valid entries are any strictly positive numbers.  It is assumed that negative
# numbers or NA indicates that that comparison was not made.
# Input:
# pairwise_matrix: a square matrix of pairwise comparisons.  Assumes either NA or -1 for comparisons that are not made
# Output:
# weights: weights that are the normalised geometric mean of valid entries of each row in the matrix
derive_weights <- function(pairwise_matrix) {
  log_gm_mean = function(a) { # Computes geometric mean for vector a of strictly positive elements
    if (length(a)==1) {
      return(log(a))
    } else {
      return(sum(log(a))/(length(a)-1))
    }
  }
  n_primary_events <- dim(pairwise_matrix)[1]
  log_weights <- vector(mode="numeric",length=n_primary_events)
  for (i in 1:n_primary_events) {
    valid_comparisons <- (pairwise_matrix[i,] > 0) & (is.na(pairwise_matrix[i,])==FALSE)
    log_weights[i] = log_gm_mean(pairwise_matrix[i,valid_comparisons])
  }
  weights_norm <- exp(log_weights - max(log_weights))
  return(weights_norm/sum(weights_norm))
}
# #########################################################################################################
