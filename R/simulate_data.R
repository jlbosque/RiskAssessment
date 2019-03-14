#' @title Simulate Data
#' @description Function that simulate complete data, incomplete data and also approximate the top event probability distribution by simulation.
#' The data is a matrix type in R with the number of rows equal to the number of observation and the number
#'        of columns equal to the number of nodes in the tree.  The primary nodes must be listed as
#'        the first columns and the top event node must be the last column.
#'
#'     Each element of the matrix is either a 1 (that event was observed to occur), 0 (that event
#'        was observed not to occur) or NA (that event was not observed)
#' @param n_simulated_data: number of observations to be generated.
#' @param tree_definition is the fault tree structure to be used, as defined above in the create_fault_tree function.
#' @param data_type: is one of "complete" (all nodes are observed), "top_only" (only the top node is
#'        observed, all others are NA) or "incomplete" (each node is randomly observed with
#'        probability p_obs, otherwise it is NA)
#' @param true_primary_p is an array of the true probabilities of each primary event occurring.
#'        Obviously must be of the same dimension as the number of primary events in tree.
#' @param p_obs: in the case that data_type="incomplete", it is the probability that a node is observed.
#' @return A Matrix with number of rows equal to the number of simulated data and number of columns equal to the number of nodes in the tree.
#' @examples data <- simulate_data(n_simulated_data = 5,
#'                       tree_definition  = tree,
#'                       data_type        = "complete",
#'                       true_primary_p   = c(0.02,0.05,0.05,0.1),
#'                       p_obs            = 0.5)
#' @export simulate_data

# ####################################################################################################################################
# Monte Carlo simulation functions
# Contains routines that simulate complete data, incomplete data and also approximate the top event probability distribution by simulation
# ####################################################################################################################################

# ###############################################################################################################
# Data
# Either it's in the matrix data or, if we're simulating data, simulate it now
# ###############################################################################################################
simulate_data <- function (n_simulated_data, tree_definition, data_type, true_primary_p, p_obs) {
  if (data_type == "top_only") {        # Simulate data that consists of the top event only
    data <- simulate_top_event_data(n=n_simulated_data,p=true_primary_p,tree_definition)
  } else if (data_type=="incomplete") { # Simulate data where each event is observed independently with probability p_obs
    data <- simulate_incomplete_data(n=n_simulated_data,p=true_primary_p,p_obs=p_obs,tree_definition)
  } else if (data_type=="intermediate") { # Simulate data that consists of the intermediate and top events only
    data <- simulate_intermediate_data(n=n_simulated_data,p=true_primary_p,tree_definition)
  } else {                              # Simulate data where all events observed
    data <- simulate_complete_data(n=n_simulated_data,p=true_primary_p,tree_definition)
  }
return (data)
}

# ####################################################################################################################################
# simulate_complete_data
# ####################################################################################################################################
# Simulate complete data
# Simulate n values. Each value is a set of indept. Bernoulli values with sucess probabilities given in the vector p
# ####################################################################################################################################
simulate_complete_data <- function (n, p, tree) {

  data <- matrix(NA,nrow=n,ncol=tree$n_nodes)
  # Simulate the primary events
  for (i in 1:tree$n_primary) {
    data[,i] <- rbinom(n,size=1,prob=p[i])
  }

  # Infer the values of the intermediate and top events
  start <- tree$n_primary+1
  for (j in start:tree$n_nodes) {
    if (tree$nodes[[j]]$Logic == "and") {
      data[,j] <- and_logic (d=data, child=tree$children[(j-tree$n_primary),])
    }
    else {
      data[,j] <-  or_logic (d=data, child=tree$children[(j-tree$n_primary),])
    }
  }
  return(data)
}
# ####################################################################################################################################

# ####################################################################################################################################
# simulate_intermediate_data
# ####################################################################################################################################
# Simulate data from all non-primary e.g. intermediate and top events
# Simulate n values. Each value is a set of indept. Bernoulli values with sucess probabilities given in the vector p
# ####################################################################################################################################
simulate_intermediate_data <- function(n, p, tree) {

  data <- matrix(NA,nrow=n,ncol=tree$n_nodes)
  # Simulate the primary events
  for (i in 1:tree$n_primary) {
    data[,i] <- rbinom(n,size=1,prob=p[i])
  }

  # Infer the values of the intermediate and top events
  start <- tree$n_primary+1
  for (j in start:tree$n_nodes) {
    if (tree$nodes[[j]]$Logic == "and") {
      data[,j] <- and_logic (d=data, child = tree$children[(j-tree$n_primary),])
    }
    else {
      data[,j] <- or_logic (d=data, child = tree$children[(j-tree$n_primary),])
    }
  }
  data[,1:tree$n_primary] <- NA

  return(data)
}

# ####################################################################################################################################
# simulate_incomplete_data
# ####################################################################################################################################
# Simulate incomplete data with each observation observed with probability p_obs
# Simulate n values. Each value is a set of indept. Bernoulli values with sucess
# probabilities given in the vector p
# ####################################################################################################################################
simulate_incomplete_data <- function(n,p,p_obs,tree) {
  # Generate complete data
  data <- matrix(NA,nrow=n,ncol=tree$n_nodes)
  for (i in 1:length(p)) {
    data[,i] <- rbinom(n,size=1,prob=p[i])
  }

  # Each element of the matrix is observed with probability p_obs, otherwise it is replaced by NA
  start <- tree$n_primary+1
  for (j in start:tree$n_nodes) {
    if (tree$nodes[[j]]$Logic == "and") {
      data[,j] <- and_logic (d=data, child = tree$children[(j-tree$n_primary),])
    }
    else {
      data[,j] <- or_logic (d=data, child = tree$children[(j-tree$n_primary),])
    }
  }

  # Each element of the matrix is observed with probability p_obs, otherwise it is replaced by NA
  observed_flag <- matrix(rbinom(n*tree$n_nodes,size=1,prob=p_obs),nrow=n)

  # Replace unobserved elements with NA
  data[observed_flag==0] <- NA
  return(data)
}
# ####################################################################################################################################

# ####################################################################################################################################
# simulate_top_event_data
# ####################################################################################################################################
# Simulate data where the top event only is observed
# Simulate n values. Each value is a set of indept. Bernoulli values with sucess probabilities given in the vector p
# ####################################################################################################################################
simulate_top_event_data <- function(n, p, tree) {

  # Generate complete data
  data <- matrix(NA,nrow=n,ncol=tree$n_nodes)
  # Simulate the primary events
  for (i in 1:length(p)) {
    data[,i] <- rbinom(n,size=1,prob=p[i])
  }

  # Infer the values of the intermediate and top events
  start <- tree$n_primary+1
  for (j in start:tree$n_nodes) {
    if (tree$nodes[[j]]$Logic == "and") {
      data[,j] <- and_logic (d=data, child=tree$children[(j-tree$n_primary),])
    }
    else {
      data[,j] <- or_logic (d=data, child = tree$children[(j-tree$n_primary),])
    }
  }

  if (length(dim(data)[2])==0) { # Is the tree just the top event?
    data_top_event <- data
  } else {
    # Return n observations of the top event only with all other events recorded as NA
    n_events <- dim(data)[2]
    data_top_event <- cbind(matrix(NA,nrow=n,ncol=n_events-1),data[,n_events])
  }
  return(data_top_event)
}
# ####################################################################################################################################

# ####################################################################################################################################
# simulate_top_event_probability
# ####################################################################################################################################
# Simulate the top event probability prior in a fault tree
# Assume independent beta(2,2) priors on primary event probabilities
# Function inputs are: the number of simulations n_simulations to use, and a n x 2 matrix consisting of n beta parameter pairs,
# one for each of the n primary event probabilities.  Output is a KDE of the prior density of the top event probability
# ####################################################################################################################################
simulate_top_event_probability <- function(n_simulations, beta_params, tree, x_limits, lty) {
  x_min <- 0
  x_max <- 1
  simulated_p <- matrix(nrow=tree$n_nodes,ncol=n_simulations)
  # Simulate the primary event probabilities from their beta prior
  for (i in 1:tree$n_primary) {
    simulated_p[i,] <- rbeta(n_simulations,shape1=beta_params[i,1],shape2=beta_params[i,2])
  }

  start <- tree$n_primary+1
  for (j in start:tree$n_nodes){
    if (tree$nodes[[j]]$Logic  == "and") {
      simulated_p[j,] <- calculate_and_scalar(child=tree$children[(j-tree$n_primary),], mysim=simulated_p)
    }
    else {
      simulated_p[j,] <- calculate_or_scalar(child=tree$children[(j-tree$n_primary),], mysim=simulated_p)
    }
  }

  # Create an estimate of the density of the top event probability (which we assume is the last-indexed event)
  p_top <- density(simulated_p[tree$n_nodes,],from=x_min,to=x_max)
  plot_y_max <- max(p_top$y)
  par(mfrow=c(1,1))
  plot(p_top$x,p_top$y,type="l",xlab="TOP EVENT PROBABILITY",ylab="DENSITY",main="",
       xlim=x_limits,ylim=c(0,1.05*plot_y_max),lwd=2,lty=lty)
  grid()
  
  cat("Prior probability of top event:","\n")
  cat("   Mean is", mean(simulated_p[tree$n_nodes,]), "\n") # modified for ATV example
  cat("   Standard deviation is", sd(simulated_p[tree$n_nodes,]), "\n")
  cat("   Central 95% probability interval is (",quantile(simulated_p[tree$n_nodes,],0.025),
      ", ",quantile(simulated_p[tree$n_nodes,],0.975),")", "\n",sep="")
}
# ###########################################################################################################

# ###########################################################################################################
# and_logic. Performs the AND function of all the children of a node.
# ###########################################################################################################
and_logic <- function(d, child) {

  first<-TRUE
  for (k in 1:(ncol(d)-1)) {
    if (child[k]==1) {
      if (first==TRUE) {
        aux<-d[,k]
        first<-FALSE
      }
      else {
        aux <- aux & d[,k]
      }
    }
  }
  return (as.integer(aux))
}

# ################################################################################
# or_logic. Performs the OR function of all the children of a node.
# ################################################################################
or_logic <- function(d, child) {
  first<-TRUE
  for (k in 1:(ncol(d)-1)) {
    if (child[k]==1) {
      if (first==TRUE) {
        aux<-d[,k]
        first<-FALSE
      }
      else {
        aux <- aux | d[,k]
      }
    }
  }
  return (as.integer(aux))
}

# #####################################################################################################
# and_scalar. Calculates the probability of a "and" tree node based its children's probabilities
# #####################################################################################################
calculate_and_scalar <- function(child, mysim){
  if (length(dim(mysim))==0) {  # just passing a vector of probabilities
    just_a_vector <- TRUE 
  } else {                      # else we're passing a matrix 
    just_a_vector <- FALSE
  }
  aux <- 1
  for (j in 1:length(child)) {
    if (child[j]==1) {
      if (just_a_vector) {
         aux <- aux * mysim[j]
      } else {
         aux <- aux * mysim[j,]
      }     
    }
  }
  return (aux)
}

# #####################################################################################################
# and_scalar. Calculates the probability of a "or" tree node based its children's probabilities.
# #####################################################################################################
calculate_or_scalar <- function(child, mysim){
  if (length(dim(mysim))==0) {  # just passing a vector of probabilities
    just_a_vector <- TRUE 
  } else {                      # else we're passing a matrix 
    just_a_vector <- FALSE
  }  
  aux <- 1
  for (j in 1:length(child)) {
    if (child[j]==1) {
      if (just_a_vector) {
        aux <- aux * (1-mysim[j])
      } else {
        aux <- aux * (1-mysim[j,])
      }    
    }
  }
  aux <- (1-aux)
  return (aux)
}
