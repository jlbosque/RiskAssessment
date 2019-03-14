#' @title Report_Posterior
#' @description This function reports the posterior distribution in varios ways. 
#' This includes a text summary for each event probability and 2 plots: one for the top event posterior and one for the 
#' primary event posteriors.  There is an option to include the prior distribution in the top event plot
#' You'll have to go into the code to alter plot parameters such as axes limits!
#' @param posterior 
#' @param tree: 
#' @param text_summary: Print summary posterior statistics to the screen?
#' @param plot_top_event_posterior: Plot the top event probability posterior?
#' @param plot_primary_event_posterior: Plot the primary event probabilities posteriors?
#' @param include_priors: Include the prior distribution in the top event plot?
#' @param save_figures: Save the 2 figures to file?
#' @return 
#' @examples report_posterior(posterior, tree,
#'                              text_summary = TRUE,                   
#'                              plot_top_event_posterior = TRUE,       
#'                              plot_primary_event_posterior = FALSE,  
#'                              include_priors = TRUE,                 
#'                              save_figures = FALSE)                 
#' @export report_posterior

report_posterior <- function(posterior, tree, text_summary, plot_top_event_posterior, 
                             plot_primary_event_posterior, include_top_event_prior,
                             save_figures) 
{
  # ########################################################################################################################
  # Print to screen a set summary statistics
  # ########################################################################################################################
  if (text_summary) {
    cat("Summary statistics of each p_i = P(E_i = 1):\n\n")
    colnames(posterior) <- sapply(as.character(1:dim(posterior)[2]), function(x) paste("p_",x,sep=""),simplify="array")
    summary.matrix(posterior)
    cat("\n")
  }
  # ########################################################################################################################

  # ########################################################################################################################
  # Plot the probability of the top event
  # ########################################################################################################################
  if (plot_top_event_posterior) {
    p_top <- density(posterior[,tree$n_nodes],from=0,to=1)
    plot_y_max <- max(p_top$y)
    if (save_figures) pdf("top_event_inference.pdf")
    par(mfrow=c(1,1))
    if (include_top_event_prior) {
      simulate_top_event_probability(n_simulations=100000, beta_parameters, tree, x_limits=c(0,0.3), lty=2)
      lines(p_top$x,p_top$y,lwd=2,lty=1)
      legend("topright",c("PRIOR","POSTERIOR"),lty=c(2,1),lwd=c(2,2))
    }  else {
       plot(p_top$x,p_top$y,type="l",xlab="TOP EVENT PROBABILITY",ylab="POSTERIOR DENSITY",main="",
          xlim=c(0,0.3),ylim=c(0,13),lwd=2,lty=1)   
    }
    grid()
    if (save_figures) dev.off()
  }
  # ########################################################################################################################
  
  # ########################################################################################################################
  # Plot the probabilities for the primary events
  # ########################################################################################################################
  if (plot_primary_event_posterior) {
    n_primary_events = tree$n_primary
    cornerstone_event_index = tree$corner 
    non_cornerstone_events = 1:tree$n_primary
    non_cornerstone_events = non_cornerstone_events[-tree$corner]
    p_corner <- density(posterior[,cornerstone_event_index],from=0,to=1)
    if (save_figures) pdf("primary_event_inference.pdf")
    par(mfrow=c(1,1))
    plot(p_corner$x,p_corner$y,
         type="l",lwd=1.5,xlab="EVENT PROBABILITY",ylab="POSTERIOR DENSITY",lty=cornerstone_event_index,
         xlim=c(0,0.25),ylim=c(0,40))
    for (i in non_cornerstone_events) {
      p_i <- density(posterior[,i],from=0,to=1)
      lines(p_i$x,p_i$y,lty=i,lwd=1.5)
    }
    grid()
    legend_lty <- 1:n_primary_events
    legend_lwd <- rep(1.5,n_primary_events)
    legend_text <- rep("Event ",n_primary_events)
    for (i in 1:n_primary_events) {
      legend_text[i] <- paste(legend_text[i],as.character(i))
    }
    legend("topright",legend=legend_text,lty=legend_lty,lwd=legend_lwd,cex=0.8)
    if (save_figures) dev.off()
  }
  # ########################################################################################################################
}  
  