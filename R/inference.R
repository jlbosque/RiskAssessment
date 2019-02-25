#' @title Inference
#' @description This function computes the Posterior distribution of event probabilities
#' by MCMC and store the results.
#' It takes the prior beta_parameters (output from derive_primary_event_priors),
#' the data (might be simulated or real), and
#' a list of 4 parameters for the MCMC: n_MCMC,  logit_p_proposal_sd, burn_in and thinning.
#' The output is matrix p_samples.
#' posterior is a matrix that contains samples from the posterior distribution of the event probabilities.
#' Each row of the matrix is one sample.
#' Each column is samples for one of the event probabilities, ordered as specified in tree.
#' @param beta_parameters: A matrix of dimension number of primary events by 2.
#' It has the beta parameters computed by the derive_primary_event_priors function.
#' @param data: A Matrix with number of rows equal to the number of simulated data and number of columns equal to the number of nodes in the tree.
#' It has the values computed by the Simulate_data function.
#' @param data_type: is one of "complete" (all nodes are observed), "top_only" (only the top node is
#'        observed, all others are NA) or "incomplete" (each node is randomly observed with
#'        probability p_obs, otherwise it is NA)
#' @param MCMC_parameters: It is a list with the following variables:
#'
#' n_MCMC:              number of MCMC samples to generate
#'
#' burn_in:             initial proportion of samples to discard as burn_in
#'
#' thinning:            save 1 in every thinning samples to use in construction of posterior distributions
#'
#' logit_p_proposal_sd: the std. dev. of the proposal on the logit p_i
#'
#' @param tree: Tree object.
#' @return A matrix that contains samples from the posterior distribution of the event probabilities.
#' Each row of the matrix is one sample.
#' Each column is samples for one of the event probabilities, ordered as specified in tree.
#' @examples MCMC_parameters <- list(n_MCMC=10000, burn_in=0.2, thinning=10, logit_p_proposal_sd=0.2)
#' posterior <- inference(beta_parameters, data, data_type="complete", MCMC_parameters, tree)
#' @export inference

# ##########################################################################################################
# COMPUTE POSTERIOR DISTEIBUTION OF EVENT PROBABILITUES BY MCMC AND STORE RESULTS
# ##########################################################################################################
# posterior is a matrix that contains samples from the posterior distribution of the event probabilities.
# Each row of # the matrix is one sample.  Each column is samples for one of the event probabilities,
# ordered as specified in tree.
# ##########################################################################################################
inference <- function (beta_parameters, data, data_type, MCMC_parameters, tree) {

  if (length(dim(data)[1])==0) {
    data <- matrix(data,nrow=1)
  }
  # ###############################################################################################################
  # Compute prior means and modes
  # ###############################################################################################################
  prior_means <- vector(mode="numeric",length=tree$n_primary)
  prior_modes <- vector(mode="numeric",length=tree$n_primary)
  for (i in 1:tree$n_primary) {
    prior_means[i] =  beta_parameters[i,1]/sum(beta_parameters[i,])
    prior_modes[i] = (beta_parameters[i,1]-1)/(sum(beta_parameters[i,])-2)
  }

  # ###############################################################################################################
  # Set up vectors and matrices to store the MCMC output
  # ###############################################################################################################
  logit_p_samples = matrix(NA,nrow=MCMC_parameters$n_MCMC, ncol=tree$n_primary)
  logit_p_current <- log(prior_means) - log(1-prior_means)
  logit_p_samples[1,] <- logit_p_current # Start the MCMC at the prior mean values
  p_current <- exp(logit_p_current)/(1+exp(logit_p_current))
  # If we're analysing incomplete data then compute the likelihood_encoding otherwise leave it null
  if ((data_type=="incomplete") | (data_type=="intermediate")) {
    likelihood_encoding <- compute_likelihood_encoding(data, tree)
  } else {
    likelihood_encoding <- matrix()
  }
  log_likelihood_current <- log_likelihood(data, p_current, data_type, likelihood_encoding, tree)
  accepted <- 0
  # ###############################################################################################################
  # Do the MCMC
  # Computed via a random walk Metropolis on logit p_i
  # ###############################################################################################################
  for (i in 2:MCMC_parameters$n_MCMC) {
    logit_p_proposal        <- propose_logit_p(logit_p_current, tree$n_primary, MCMC_parameters$logit_p_proposal_sd,limits=c(-500,16))
    p_proposal              <- exp(logit_p_proposal)/(1+exp(logit_p_proposal))
    log_prior_ratio         <- log_beta_prior_ratio(p_proposal,p_current,beta_parameters)
    log_likelihood_proposal <- log_likelihood(data,p_proposal,data_type,likelihood_encoding,tree)
    log_accept_ratio = log_prior_ratio + log_likelihood_proposal - log_likelihood_current
    if (log(runif(1,0,1)) < log_accept_ratio) { # Accept the proposal
      logit_p_samples[i,] <- logit_p_proposal
      logit_p_current <- logit_p_proposal
      p_current <- p_proposal
      log_likelihood_current <- log_likelihood_proposal
      accepted <- accepted + 1
    } else {
      logit_p_samples[i,] <- logit_p_samples[i-1,]
    }
#    if (i%%100000==0) { cat("\n Iteration",i,"of",n_MCMC)}
  }
#  cat('\n',accepted,' proposals out of ',n_MCMC,' accepted.\n')
  # Remove burn in and thin samples
  exp_logit_p_samples <- exp(logit_p_samples)
  p_samples <- exp_logit_p_samples/(1+exp_logit_p_samples)
  p_samples <- p_samples[round(1+MCMC_parameters$burn_in*MCMC_parameters$n_MCMC):MCMC_parameters$n_MCMC,]
  thinning_indices = seq(from=1,to=dim(p_samples)[1],by=MCMC_parameters$thinning)
  p_samples <- p_samples[thinning_indices,]

  # Derive the probability for all of the non-primary events
  p_samples <- derive_nonprimary_event_sample_values(p_samples,tree)

 return (p_samples)
}

# ####################################################################################################################################
# compute_likelihood_encoding
# ####################################################################################################################################
# Derives a list that encodes the likelihood for a set of data that includes incomplete observations
# It does this for the 2 examples in the paper but obviously we'd want a completely
# general version of the function that worked on any given fault tree.
# The likelihood is a product of terms, one for each observation. Each of these terms is a sum of terms that
# are products of p_i and (1-p_i).
# The list encodes this as follows:
#   Each element of the list encodes the likelihood for 1 observation
#   Each of these elements consists of a matrix with number of columns equal to the number of primary events.
#   Each row of the matrix encodes one of the product terms with 1 (for p_i) and 0 (for 1-p_i) and NA (for does not appear)
# ####################################################################################################################################
compute_likelihood_encoding <- function(data,tree) {
  if (length(dim(data)[1])==0) {
    n_obs <- 1  # data consists of a single observation
    likelihood_encoding <- matrix(ncol=1+length(data))
  } else {
    n_obs <- dim(data)[1]
    likelihood_encoding <- matrix(ncol=1+dim(data)[2])
  }
  for (obs_index in 1:n_obs) {
    # Infer data that can be logically inferred from the observed values.
    observation <- logically_complete (data[obs_index,],tree)
    # Are any primary events unobserved? If so then likelihood is sum of terms over possible combinations of
    # primary events. Otherwise it's just the likelihood of all primary events.
    if (any(is.na(observation[1:tree$n_primary]))) {
      # Derive all possible combinations of primary events (with observed ones fixed at their observed value)
      unobserved_combinations <- compute_all_combinations(observation,tree)
      obs_encoder_matrix <- c()
      # For each combination see if it is consistent with the tree logic and if so then add it to the matrix of combinations
      for (j in 1:dim(unobserved_combinations)[1]) {
        if (feasible_complete_observation(unobserved_combinations[j,],tree)) {
          obs_encoder_matrix <- rbind(obs_encoder_matrix,unobserved_combinations[j,])
        }
      }
      likelihood_encoding <- rbind(likelihood_encoding,cbind(matrix(rep(obs_index,dim(obs_encoder_matrix)[1]),ncol=1),
                                                             obs_encoder_matrix))
    } else {
      likelihood_encoding <- rbind(likelihood_encoding,c(obs_index,observation))
    }
  }
  # Remove the first row of NAs in likelihood_encoding that was just there to set it up with the right number of columns
  likelihood_encoding <- likelihood_encoding[-1,]
  return(likelihood_encoding)
}
# ####################################################################################################################################

# ####################################################################################################################################
# compute_all_combinations
# ####################################################################################################################################
# This function takes an observation and returns all combinations of the observation with values for the unobserved primary events. The unobserved
# non-primary events are then logically deduced.
# Note that that is no guarantee that the observation is consistent with the logic of the fault tree.  Observed non-primary event values may not be
# consistent with a particular combination of primary event values.
# ####################################################################################################################################
compute_all_combinations <- function(observation,tree) {
  primary_events                  <- observation[1:tree$n_primary]
  unobserved_primary_events_index <- which(is.na(primary_events))
  number_unobserved_primary_events <- length(unobserved_primary_events_index)
  # The function expand.grid gets us all the unobserved combinations
  unobserved_combinations <- expand.grid(rep(list(0:1), number_unobserved_primary_events))
  # Combinations will contain all the combinations. Each column is a component of the observation
  n_combinations <- dim(unobserved_combinations)[1]
  combinations <- matrix(NA,nrow=n_combinations,ncol=length(observation))
  unobserved_combination_index <- 1
  # Go through the primary events. Add in a column from the unobserved combinations
  # if it's unobserved. Add in the observed value otherwise
  for (i in 1:length(primary_events)) {
    if (is.na(observation[i])) {
      combinations[,i] <- unobserved_combinations[,unobserved_combination_index]
      unobserved_combination_index <- unobserved_combination_index + 1
    } else {
      combinations[,i] <- rep(observation[i],n_combinations)
    }
  }
  # Go through the non-primary events and add in the logically implied value if unobserved.
  start <- tree$n_primary+1
  end <- tree$n_nodes-1
  for (i in start:tree$n_nodes){
    if (is.na(observation[i])) { # if the event is not observed
      combinations[,i] <- 1
      child <- tree$children[i-tree$n_primary,]
      if (tree$nodes[[i]]$Logic  == "and") {
        for (j in 1:end) {
          if (child[j] == 1) {
            combinations[,i] <- combinations[,i] * combinations[,j]
          }
        }
      }
      else {
        for (j in 1:end) {
          if (child[j] == 1) {
            combinations[,i] <- combinations[,i] * (1-combinations[,j])
          }
        }
        combinations[,i] <- 1 - combinations[,i]
      }
    } else {
      combinations[,i] <- rep(observation[i],n_combinations)
    }
  }

  return(combinations)
}
# ####################################################################################################################################

# ####################################################################################################################################
# logically_complete
# ####################################################################################################################################
# This function takes a vector of observations from the tree (1, 0 or NA if unobserved) and derives the values of any unobserved
# nodes that are logically implied from the observed ones.
# ####################################################################################################################################
logically_complete <- function(observation,tree) {
  # Cycle up the tree. Do this twice to catch all possible logical implications.
  # For each non-primary node: if it's observed then does it imply values for its parents?  If it's not observed then do its
  # parents imply its value?
  # For an unobserved OR node:  it's 1 if at least 1 parent is observed to be 1; it's 0 if all parents are observed and are 0.
  # For an observed OR node:    if it's 0 then all parents must be 0; if it's 1 and all but one parent are observed to be 0 then the
  #                             1 unobserved parent must be 1.
  # For an unobserved AND node: it's 0 if at least 1 parent is observed to be 0; it's 1 if all parents are observed and are 1.
  # For an observed AND node:   if it's 1 then all parents must be 1; if it's 0 and all but one parent are observed to be 1 then the
  #                             1 unobserved parent must be 0.

  for (jj in 1:2) {  # Cycle twice
    # For every Intermediate event
    for (ii in 1:tree$n_inter) {
        # ii<-3
      parents <- observation[tree$children[ii,]==1]
      if (length(tree$children[ii,tree$children[ii,]==1]) < length(parents)) parents <- parents[1:length(tree$children[ii,tree$children[ii,]==1])]
      if (is.na(observation[ii+tree$n_primary]))  { # If it's not observed then check to see if it can be inferred from its parent
          if (tree$nodes[[ii+tree$n_primary]]$Logic  == "or")  { # OR gate
             if (any(parents==1,na.rm=TRUE))   { observation[ii+tree$n_primary] <- 1  }  # There is an observed parent with value 1 => = 1
             if ((all(parents==0,na.rm=TRUE)) &
                (any(is.na(parents))==FALSE)) { observation[ii+tree$n_primary] <- 0 } # All parents observed with value 0 => = 0
          } else { # AND gate
             if (any(parents==0,na,rm=TRUE))  {observation[ii+tree$n_primary] <- 0 }  # There is an observed parent with value 0  => = 0
             if ((all(parents==1,na.rm=TRUE)) &
                 (any(is.na(parents))==FALSE)) { observation[ii+tree$n_primary] <- 1 } # All parents observed with value 1 => = 1
          }
      } else { # if it is observed then can it infer any unobserved parents?
        if (tree$nodes[[ii+tree$n_primary]]$Logic  == "or")  { # OR gate
            if (observation[ii+tree$n_primary]==0) {
              parents <- rep(0,length(parents)) # All parents must be 0
            } else {
              if (all(parents==0,na.rm=TRUE) & (sum(is.na(parents))==1)) {
                 # 1 unobserved parent and it has to be 1 because all other parents are observed to be 0.
                 parents[which(is.na(parents))] <- 1
              }
            }
        } else {  # AND gate
            if (observation[ii+tree$n_primary]==1) {
                parents <- rep(1,length(parents)) # All parents must be 0
            } else {
                if (all(parents==1,na.rm=TRUE) & (sum(is.na(parents))==1)) {
                    # 1 unobserved parent and it has to be 1 because all other parents are observed to be 0.
                    parents[which(is.na(parents))] <- 0
                }
            }
        }
        observation[tree$children[ii,]==1]<-parents
      }
    } # end for ii
  } # end for jj
  return(observation)
}
# ####################################################################################################################################

# ####################################################################################################################################
# feasible_complete_observation
# ####################################################################################################################################
# This function takes a list of complete observed events (as a binary vector)
# and returns whether they are consistent with the simple example used in the paper
# ####################################################################################################################################
feasible_complete_observation <- function(list_of_events,tree) {
  feasible <- TRUE
  logic <- c(0)
  end <- tree$n_nodes-1

  for (i in 1:tree$n_inter){
    logic[i] <- 1
    child <- tree$children[(i),]
    if (tree$nodes[[i+tree$n_primary]]$Logic  == "and") {
      for (j in 1:end) {
        if (child[j] == 1) {
          logic[i] <- logic[i] * list_of_events[j]
        }
      }
      feasible[i] <- (list_of_events[i+tree$n_primary] == logic[i])
    }
    else {
      for (j in 1:end) {
        if (child[j] == 1) {
          logic[i] <- logic[i] * (1-list_of_events[j])
        }
      }
      feasible[i] <- (list_of_events[i+tree$n_primary] == (1 -logic[i]))
    }
  }

  return(all(feasible))
}
# ####################################################################################################################################
# Now the functions that evaluate the log likelihood for a set of data with a given set of values for p.
# It also contains a function that evaluates the log prior ratio for two values of p.
# ####################################################################################################################################
# log_likelihood
# ####################################################################################################################################
# ####################################################################################################################################
log_likelihood  <- function(data, p, data_type, likelihood_encoding,tree) {
  if (data_type=="complete") {
    ll <- log_complete_likelihood(data,p)
  } else if (data_type=="top_only") {
    ll <- log_top_event_likelihood(data,p,tree)
  } else {
    ll <- log_incomplete_likelihood(data,p,likelihood_encoding)
  }
  return(ll)
}
# ####################################################################################################################################

# ####################################################################################################################################
# log_complete_likelihood
# ####################################################################################################################################
# Evaluates the likelihood in the case where we observe all events
# ####################################################################################################################################
log_complete_likelihood <- function(data,p) {
  loglike <- 0
  n_primary_events <- length(p)
  for (k in 1:dim(data)[1]) {
    loglike <- loglike + sum(data[k,1:n_primary_events]*log(p) + (1-data[k,1:n_primary_events])*log(1-p))
  }
  return(loglike)
}
# ####################################################################################################################################

# ####################################################################################################################################
# log_top_event_likelihood
# ####################################################################################################################################
# Evaluates the likelihood in the case where we just observe the top event (again just for the case of the simple example in the paper)
# ASSUMES THAT THE TOP EVENT OBSERVATIONS ARE THE LAST COMPONENT IN AN OBSERVATION, AND HENCE THE LAST COLUMN OF THE MATRIX data
# ####################################################################################################################################
# JLB. The function is generalized to any tree structure.
# ####################################################################################################################################
log_top_event_likelihood <- function(data,p,tree) {

    if (length(dim(data)[2])==0) {
    n_events <- length(data)  # data consists of a single observation
  } else {
    n_events <- dim(data)[2]
  }
  # If values of p_top are very close to 1 or 0 then p_top can be evaluated as 1 or 0, leading to a NaN in the log-likelihood.  So threshold
  # values of p_top to be strictly between 0 and 1.
  start <- tree$n_primary+1
  for (i in start:tree$n_nodes){
    child <- tree$children[(i-tree$n_primary),]
    if (tree$nodes[[i]]$Logic  == "and") {
      p[i] <- calculate_and_scalar(child,p)
    }
    else {
      p[i] <- calculate_or_scalar(child,p)
    }
  }
  p_top <- p[tree$n_nodes]

  if (p_top==1) {
    p_top <- 0.999999
  }
  if (p_top==0) {
    p_top <- 0.000000000000001
  }
  loglike <- sum(data[,n_events]*log(p_top) + (1-data[,n_events])*log(1-p_top))
  return(loglike)
}
# ####################################################################################################################################

# ####################################################################################################################################
# log_incomplete_likelihood
# ####################################################################################################################################
# Evaluates the likelihood in the case where data are incomplete (arbitrarily unobserved events)
# The matrix likelihood_encoding has been computed by the compute_likelihood_encoding function and it is explained there.
# The log likelihood is the sum of terms, one for each observation (being that observation's log likelihood).
# The log likelihood for each observation is the log of a sum of terms; each of those terms is a product whose log we calculate.
# Hence we want the log of a sum of terms that we know the log of!  We use log-sum-exp type rule to do this in a computationally
# stable way
log_incomplete_likelihood <- function(data,p,likelihood_encoding) {
  which_term <- function(e,obs_index) {
    if (is.na(e)) {
      return(1)
    } else {
      return((p[obs_index]^e)*((1-p[obs_index])^(1-e)))
    }
  }
  loglike <- 0
  if (length(dim(data)[1])==0) {
    n_obs <- 1  # data consists of a single observation
  } else {
    n_obs <- dim(data)[1]
  }
  for (i in 1:n_obs) {
    # X is the log likelihoods for each combination of primary events that are consistent with the data
    # for observation i, evaluated at p
    # The first column of likelihood_encoding contains the observation number, the next length(p) columns
    # contain values for all of the primary events
    X <-    likelihood_encoding[likelihood_encoding[,1]==i,2:(1+length(p))]%*%log(p) +
      (1-likelihood_encoding[likelihood_encoding[,1]==i,2:(1+length(p))])%*%log(1-p)
    # The log likelihood for observation i is the log_sum_exp of these terms, which we do in a numerically stable way
    # by using the maximum value
    maxX <- max(X)
    loglike <- loglike + maxX + log(sum(exp(X-maxX)))
  }
  return(loglike)
}
# ####################################################################################################################################
# log_beta_prior_ratio
# ####################################################################################################################################
# Evaluates the log ratio of 2 beta priors
# log_prior_ratio returns the log ratio of beta distributions at values p1 and p2 with parameters in beta_params
# ####################################################################################################################################
log_beta_prior_ratio <- function(p1,p2,beta_params) {
  ratio <- 0
  for (j in 1:length(p1)) {
    ratio <- ratio + dbeta(p1[j], shape1=beta_parameters[j,1],shape2=beta_params[j,2], log = TRUE)
    - dbeta(p2[j], shape1=beta_parameters[j,1],shape2=beta_params[j,2], log = TRUE)
  }
  return(ratio)
}
# ####################################################################################################################################
# propose_logit_p
# ####################################################################################################################################
# Generate random proposal for vector of logit primary event probabilities
# current is the vector of current logit probabilities, n is the number of them (=length of current) and sd is the bandwidth of the
# proposal move.
# propose_logit_p returns a vector or proposed values of logit of the primary event probabilities
# ####################################################################################################################################
propose_logit_p <- function(current,n,sd_p,limits) {
  proposal <- current + rnorm(n,mean=0,sd=sd_p)
  # If any value in proposal is less than -500 or more than 16 then there's potentially a problem when it's converted
  # to p in that it gets mapped to 0 or 1, which then leads to an error when calculating the log likelihood (which has
  # terms like log(p) and log(1-p).  So threshold proposals forced to be between these 2 values (as specified in limits)
  proposal[proposal < limits[1]] <- limits[1]
  proposal[proposal > limits[2]] <- limits[2]
  return(proposal)
}
# ####################################################################################################################################

# ####################################################################################################################################
# derive_nonprimary_event_sample_values.R
# ####################################################################################################################################
# This function takes the matrix of sampled primary event probabilities and adds in columns for the implied value of the
# probabilities for all of the non-primary events.
# Inputs:
#    p_samples: a matrix of all sampled values of the primary events; one column for each primary event
#    tree: The structure of the tree
# Output is p_samples with a column for each of the non-primary probabilities added on
# ################################################################################################################################
derive_nonprimary_event_sample_values <- function(p_samples,tree)  {

  paux <- matrix (1, nrow = dim(p_samples)[1], ncol = tree$n_inter)
  end <- tree$n_nodes-1

  for (i in 1:tree$n_inter){
    child <- tree$children[i,]
    if (tree$nodes[[i+tree$n_primary]]$Logic  == "and") {
      for (j in 1:end) {
        if (child[j] == 1) {
          if ((tree$nodes[[j]]$Type == "corner") || (tree$nodes[[j]]$Type == "leaf")) {
            paux[,i] <- paux[,i] * p_samples[,j]
          }
          else {
            paux[,i] <- paux[,i] * paux[,(j-tree$n_primary)]
          }
        }
      }
    }
    else {
      for (j in 1:end) {
        if (child[j] == 1) {
          if ((tree$nodes[[j]]$Type == "corner") || (tree$nodes[[j]]$Type == "leaf")) {
            paux[,i] <- paux[,i] * (1-p_samples[,j])
          }
          else {
            paux[,i] <- paux[,i] * (1-paux[,j-tree$n_primary])
          }
        }
      }
      paux[,i] <- (1 - paux[,i])
    }
  }
  p_samples <- cbind(p_samples, paux)

  return(p_samples)
}

calculate_and_scalar <- function(child, p){
  aux <- 1
  for (j in 1:length(child)) {
    if (child[j]==1) {
      aux <- aux * p[j]
    }
  }
  return (aux)
}

calculate_or_scalar <- function(child, p){
  aux <- 1
  for (j in 1:length(child)) {
    if (child[j]==1) {
      aux <- aux * (1-p[j])
    }
  }
  aux <- (1-aux)
  return (aux)
}
