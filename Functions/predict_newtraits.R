require(tidyverse)
require(gllvm)
library(mvtnorm)

# function to predict distribution of a hypothetical species using 
# counterfactual traits
predict_newtraits <-
  function(mod,
           env_name,
           trait_name,
           env_new,
           trait_new,
           quadratic = FALSE,
           intervals = c("confidence", "prediction"),
           nsim = 1000) {

    # extract point estimates
    mu <- c(beta0 = unique(mod$params$beta0),
            mod$params$B)
    
    # variance covariance
    Sigma <- vcov(mod)
    Sigma <- (Sigma + t(Sigma)) / 2 # make sure covariance matrix is perfectly symetric
    Sigma <- Sigma[str_detect(rownames(Sigma), "^b$|^B"), 
                   str_detect(colnames(Sigma), "^b$|^B")]
    
    # simulate parameter "posterior"
    param_draws <- rmvnorm(nsim, mu, Sigma)
    
    # calculate trait-moderated environmental slopes
    mm_kappa    <- c(1, trait_new)
    
    if (quadratic) {
      kappa_terms <- 
        c("beta0",
          paste0(env_name, "1"), 
          paste0(env_name, "2"), 
          trait_name,
          paste0(env_name, "1:", trait_name), 
          paste0(env_name, "2:", trait_name))
      kappa_draws <- param_draws[, kappa_terms]
      beta_draws  <- 
        cbind(kappa_draws[, c(1, 4)] %*% mm_kappa,
              kappa_draws[, c(2, 5)] %*% mm_kappa,
              kappa_draws[, c(3, 6)] %*% mm_kappa)
    } else {
      
      kappa_terms <- 
        c("beta0",env_name, 
          trait_name,
          paste0(env_name, ":", trait_name))
      
      kappa_draws <- param_draws[, kappa_terms]
      
      beta_draws  <- 
        cbind(kappa_draws[, c(1, 3)] %*% mm_kappa,
              kappa_draws[, c(2, 4)] %*% mm_kappa)
      
    }

    # calculate the linear predictor
    # final model matrix
    if (quadratic) {
      mm_beta <- cbind(1, env_new, env_new^2)
    } else {
      mm_beta <- cbind(1, env_new)
    }

    # add residuals to betas?
    if (intervals == "confidence") {
      beta_draws <- t(beta_draws)
    } else if (intervals == "prediction") {
      sigmaB_terms <- 
        c('int', env_name)
      sigmaB <- mod$params$sigmaB[sigmaB_terms, sigmaB_terms]
      beta_draws <- apply(beta_draws, 1, function(x) {
        rmvnorm(1, x, sigmaB)
      })
    }
    
    # calculate the linear predictor
    eta_draws <- t(mm_beta %*% beta_draws)
    
    # inverse link
    mu_draws  <- exp(eta_draws)
    
    # summarise
    mu_median <- apply(mu_draws, 2, median)
    mu_lci_95    <- apply(mu_draws, 2, quantile, probs = 0.025)
    mu_uci_95    <- apply(mu_draws, 2, quantile, probs = 0.975)
    mu_lci_50    <- apply(mu_draws, 2, quantile, probs = 0.25)
    mu_uci_50    <- apply(mu_draws, 2, quantile, probs = 0.75)
    
    return(
      data.frame(trait = trait_new, 
                 env = env_new, 
                 fit = mu_median, 
                 lower_95 = mu_lci_95, 
                 upper_95 = mu_uci_95,
                 lower_50 = mu_lci_50, 
                 upper_50 = mu_uci_50))
}

