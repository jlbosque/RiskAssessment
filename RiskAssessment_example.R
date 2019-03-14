# #######################################################################################################################
# This script takes you through an analysis of the simple example that appears in "Quantitative system risk assessment 
# from incomplete  data with belief networks and pairwise comparison elicitation" by De Persis et al.
# Any queries on this code can be sent to Simon Wilson at simon.wilson@tcd.ie
# #######################################################################################################################

# #######################################################################################################################
# LOAD IN ALL FUNCTIONS #################################################################################################
# #######################################################################################################################
source("./R/inference.R")
source("./R/derive_primary_event_priors.R")
source("./R/simulate_data.R")
source("./R/tree_definition.R")
source("./R/report_posterior.R")
# #######################################################################################################################

# #######################################################################################################################
# LOAD IN THE FAULT TREE STRUCTURE FROM THE SPECIFIED FILE ##############################################################
# #######################################################################################################################
# The tree is specified in the file tree.conf and a tree object is created.
tree <- create_fault_tree("tree.conf")
# #######################################################################################################################

# #######################################################################################################################
# DO THE ELICITATION ON THE PRIMARY EVENT PROBABILITIES #################################################################
# #######################################################################################################################
# Elicit the beta prior distributions for each primary event probability by the 
# pairwise comparison approach.  For this we call the derive_primary_event_priors function
# with arguments:
#   cornerstone_event_index: the index of the primary event that is the cornerstone event.  This
#                            had to be defined when you defined the tree structure in the tree.conf
#                            file.
#   cornerstone_event_interval: the range over which the cornerstone event probability is elicited
#                            from the expert.  Technically, it's mapped to the central 95% prob.
#                            interval of a beta distribution.
#   pairwise_matrix: the matrix of pairwise comparisons, as described in the paper. A square matrix
#                            whose dimension is the number of primary events.  If a particular
#                            pairwise comparison is not done then the entry is NA
# Here's a complete pairwise matrix (as used in the paper)
pairwise_matrix            = matrix(c(1.00, 0.21, 1.04, 0.53,
                                      1.52, 1.00, 1.04, 1.52,
                                      0.53, 0.53, 1.00, 0.17,
                                      1.04, 0.21, 2.55, 1.00), nrow=4,byrow=TRUE)
# And here's an incomplete one (also used in the paper):
# pairwise_matrix = matrix(c(1.00, 0.21, 1.04, 0.53,
#                            1.52, 1.00, NA,   NA,
#                            0.53, NA,   1.00, NA,
#                            1.04, NA,   NA,   1.00), nrow=tree$n_primary,byrow=TRUE)
beta_parameters <- derive_primary_event_priors(cornerstone_event_index    = tree$corner, 
                                               cornerstone_event_interval = c(0.01,0.05),
                                               pairwise_matrix)
# Use these beta parameters to simulate the top event probability and plot a KDE of its prior.  The other arguments to
# this function are:
#   n_simulations: the number of simulations to generate
#   x_limits: the x-axis limits to use on the plot of the top event prior probability
simulate_top_event_probability(n_simulations=100000, beta_parameters, tree, x_limits=c(0,0.2))
# #######################################################################################################################

# #######################################################################################################################
# LOAD IN THE DATA ######################################################################################################
# The data has the following form:
#     It's a matrix type in R with the number of rows equal to the number of observation and the number 
#        of columns equal to the number of nodes in the tree.  The primary nodes must be listed as
#        the first columns and the top event node must be the last column.
#     Each element of the matrix is either a 1 (that event was observed to occur), 0 (that event
#        was observed not to occur) or NA (that event was not observed)
# You can specify data in 3 ways.
# The first is to simulate the data as follows: 
data <- simulate_data(n_simulated_data = 5, 
                      tree_definition  = tree, 
                      data_type        = "complete",
                      true_primary_p   = c(0.02,0.05,0.05,0.1), 
                      p_obs            = 0.5)
#   where:
#     n_simulated_data is the number of observations to be generated
#     tree_definition is the fault tree structure to be used, as defined above in the create_fault_tree function.
#     data_type is one of "complete" (all nodes are observed), "top_only" (only the top node is
#        observed, all others are NA) or "incomplete" (each node is randomly observed with 
#        probability p_obs, otherwise it is NA)
#     true_primary_p is an array of the true probabilities of each primary event occurring.  Obviously
#        must be of the same dimension as the number of primary events in tree.
#     p_obs is, in the case that data_type="incomplete", the probability that a node is observed
# The second is to load the data from a file (here a csv file).  The file format can be anything really, as long
#    as it is formatted to an R matrix object as decsribed above.  Here we input the data from a CSV file
data <- as.matrix(read.csv(file="sample_data.csv",sep=","))
colnames(data) <- NULL
# Thirdly, you can just directly enter the matrix of data
data <- matrix(c(1,0,0,0,1,0,1,
                 0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,
                 1,0,0,1,1,0,1), 
               ncol=tree$n_nodes, byrow=TRUE)
# ########################################################################################################################

# ########################################################################################################################
# SET UP THE MCMC ########################################################################################################
# n_MCMC:              number of MCMC samples to generate
# burn_in:             initial proportion of samples to discard as burn_in
# thinning:            save 1 in every thinning samples to use in construction of posterior distributions
# logit_p_proposal_sd: the std. dev. of the proposal on the logit p_i
MCMC_parameters <- list(n_MCMC=100000, burn_in=0.2, thinning=10, logit_p_proposal_sd=0.2)
# ########################################################################################################################

# ########################################################################################################################
# COMPUTE POSTERIOR DISTEIBUTION OF EVENT PROBABILITUES BY MCMC AND STORE RESULTS ########################################
# posterior is a matrix that contains samples from the posterior distribution of the event probabilities.  Each row of
# the matrix is one sample.  Each column is samples for one of the event probabilities, ordered as specified in tree.
# If using simulated data, the data_type value must match that used in the simulate_data function
posterior <- inference(beta_parameters, data, data_type="complete", MCMC_parameters, tree)
# ########################################################################################################################

# ########################################################################################################################
# REPORT THE POSTERIOR DISTRIBUTION IN VARIOUS WAYS ######################################################################
# ########################################################################################################################
# This includes a text summary for each event probability and 2 plots: one for the top event posterior and one for the 
# primary event posteriors.  There is an option to include the prior distribution in the top event plot
# This is still under development and you'll have to go into the code to alter plot parameters such as axes limits!
report_posterior(posterior, tree,
                 text_summary = TRUE,                   # Print summary posterior statistics to the screen?
                 plot_top_event_posterior = TRUE,       # Plot the top event probability posterior?
                 plot_primary_event_posterior = FALSE,  # Plot the primary event probabilities posteriors?
                 include_priors = TRUE,                 # Include the prior distribution in the top event plot?
                 save_figures = FALSE)                  # Save the 2 figures to file?
# ########################################################################################################################


