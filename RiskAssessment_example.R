# library("RiskAssessment", lib.loc="~/my_stuff/research/esa/code/Jose_Luis_code")

# #######################################################################################################################
# LOAD IN ALL FUNCTIONS #################################################################################################
source("RiskAssessment/R/inference.R")
source("RiskAssessment/R/derive_primary_event_priors.R")
source("RiskAssessment/R/simulate_data.R")
source("RiskAssessment/R/tree_definition.R")
# #######################################################################################################################

# #######################################################################################################################
# LOAD IN THE FAULT TREE STRUCTURE FROM THE SPECIFIED FILE ##############################################################
# The tree is specified in the file tree.conf and a tree object is created.
tree <- create_fault_tree("~/my_stuff/research/esa/code/Jose_Luis_code/tree.conf")
# #######################################################################################################################

# #######################################################################################################################
# DO THE ELICITATION ON THE PRIMARY EVENT PROBABILITIES #################################################################
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
beta_parameters <- derive_primary_event_priors(cornerstone_event_index    = tree$corner, 
                                               cornerstone_event_interval = c(0.01,0.05), 
                                               pairwise_matrix            = matrix(c(1,   1/4,  2, 1/2,
                                                                                     4, 1,  2, 4,
                                                                                     1/2, 1/2, 1, 1/5,
                                                                                     2,   1/4, 5, 1), nrow=4,byrow=TRUE))
# Or here's an example of an incomplete comparison matrix that you could also use as the pairwise_matrix 
# argument in the derive_primary_event_priors function:
# pairwise_matrix = matrix(c(1,   1/4,  2,  1/2,
#                             4,   1,    NA, NA,
#                             1/2, NA,   1,  NA,
#                             2,   NA,   NA, 1), nrow=4,byrow=TRUE)
# #######################################################################################################################

# #######################################################################################################################
# LOAD IN THE DATA ######################################################################################################
# The data has the following form:
#     It's a matrix type in R with the number of rows equal to the number of observation and the number 
#        of columns equal to the number of nodes in the tree.  The primary nodes must be listed as
#        the first columns and the top event node must be the last column.
#     Each element of the matrix is either a 1 (that event was observed to occur), 0 (that event
#        was observed not to occur) or NA (that event was not observed)
# You can simulate the data as follows: 
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
# Or you can load the data from a file (here a csv file).  The file format can be anything really, as long
#    as it is formatted to an R matrix object as decsribed above.  Here we input the data from a CSV file
data <- as.matrix(read.csv(file="~/my_stuff/research/esa/code/Jose_Luis_code/sample_data.csv",sep=","))
colnames(data) <- NULL
# Or you can just directly enter the matrix of data
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
MCMC_parameters <- list(n_MCMC=10000, burn_in=0.2, thinning=10, logit_p_proposal_sd=0.2)
# ########################################################################################################################

# ########################################################################################################################
# COMPUTE POSTERIOR DISTEIBUTION OF EVENT PROBABILITUES BY MCMC AND STORE RESULTS ########################################
# posterior is a matrix that contains samples from the posterior distribution of the event probabilities.  Each row of
# the matrix is one sample.  Each column is samples for one of the event probabilities, ordered as specified in tree.
posterior <- inference(beta_parameters, data, data_type="complete", MCMC_parameters, tree)
# ########################################################################################################################

# ########################################################################################################################
# REPORT THE POSTERIOR DISTRIBUTION IN VARIOUS WAYS ######################################################################
report_posterior(posterior, tree,
                 text_summary = TRUE,                   # Print summary posterior statistics to the screen?
                 plot_top_event_posterior = TRUE,       # Plot the top event probability posterior?
                 plot_primary_event_posterior = TRUE,  # Plot the primary event probabilities posteriors?
                 include_priors = TRUE,                # Include the prior distribution in these plots?
                 save_figures = TRUE)                  # Save the figures to file?


