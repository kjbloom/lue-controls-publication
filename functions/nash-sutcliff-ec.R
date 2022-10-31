###############################################################################
### A function for calculating the Nash Sutcliff model efficiency criterion ###
###############################################################################

# see also hydroGOF::NSE gives identical results


NSec <- function(obs, sims) {
  
  numerator <- sum((obs - sims)^2)
  denominator <- sum((obs - mean(obs))^2)
  EC <- 1 - (numerator / denominator)
  
  return(EC)
  }


metRics <- function(obs, sims) {
  
  linmod <- lm( obs ~ sims )
  rsq    <- summary( linmod )$adj.r.squared
  rmse   <- Metrics::rmse( obs, sims )
  bias   <- mean( sims - obs )
  EC <- NSec(obs, sims)
  
  out <- list(Rsquared = rsq, RMSE = rmse, Bias = bias, nsEC = EC)
  
  return(out)
  }




