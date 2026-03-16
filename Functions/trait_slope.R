require(tidyverse)
require(gllvm)
library(mvtnorm)


# function to predict the slope for a specific environmental variable and trait combination
# or return baseline effect when trait_new is NULL
# function to predict the slope for a specific environmental variable and trait combination
# function to predict the slope for a specific environmental variable and trait combination
trait_effects <- function(mod,
                          env_name,
                          trait_prefix = NULL,  # Prefix to identify trait columns (e.g., "body_size" for "body_size_large", "body_size_small")
                          trait_type = c("categorical", "numeric"),  # Type of trait
                          trait_value = NULL,  # For numeric: "min", "mean", "max", or specific value; For categorical: specific category
                          traits,
                          nsim = 1000) {
  
  # Default trait_type if not specified
  trait_type <- match.arg(trait_type)
  
  # extract point estimates
  mu <- c(beta0 = unique(mod$params$beta0),
          mod$params$B)
  
  # variance covariance
  Sigma <- vcov(mod)
  Sigma <- (Sigma + t(Sigma)) / 2
  Sigma <- Sigma[str_detect(rownames(Sigma), "^b$|^B"), 
                 str_detect(colnames(Sigma), "^b$|^B")]
  
  # simulate parameter "posterior"
  param_draws <- rmvnorm(nsim, mu, Sigma)
  
  # Mean / SD trait data
  traits_mean <- colMeans(traits)
  traits_sd <- apply(traits, 2, sd)
  
  # Find all columns that match the trait prefix
  trait_cols <- colnames(traits)[grepl(paste0("^", trait_prefix), colnames(traits))]
  
  if (length(trait_cols) == 0) {
    stop(paste("No traits found with prefix", trait_prefix))
  }
  
  # Process based on trait type
  if (trait_type == "categorical") {
    # Get all categories of this trait
    categories <- gsub(paste0(trait_prefix, "_"), "", trait_cols)
    
    # If trait_value is specified, use only that category
    if (!is.null(trait_value)) {
      if (!trait_value %in% categories) {
        stop(paste("Category", trait_value, "not found for trait", trait_prefix))
      }
      
      trait_col <- paste0(trait_prefix, "_", trait_value)
      category_param <- paste0(env_name, ":", trait_cols)
      
      # main effect
      B_e <- param_draws[, env_name]
      # Fourth corner coef
      B_I <- param_draws[,grepl(paste0(env_name, ':'), colnames(param_draws))] 
      
      # Set all trait categories for ttrait to 0
      t_j_base <- traits_mean
      t_j_base[trait_cols] <- 0 
      
      # set trait category to 1
      t_j <- t_j_base
      t_j[trait_col] <- 1 
      
      # scale traits 
      t_j <- (t_j - traits_mean) / traits_sd
      
      # Check if the interaction term exists
      if (!all(category_param %in% colnames(param_draws))) {
        stop(paste("Parameter", category_param, "not found in model"))
      }
      
      # Calculate effect
      B_j <- B_e +  B_I %*% t_j
      
      result <- data.frame(
        env_var = env_name,
        trait = trait_prefix,
        trait_category = trait_value,
        trait_value = 1,
        median = quantile(B_j, 0.5),
        upper_95 = quantile(B_j, 0.975),
        lower_95 = quantile(B_j, 0.025),
        upper_50 = quantile(B_j, 0.75),
        lower_50 = quantile(B_j, 0.25)
      )
      
      return(rbind(result))
    } 
    
    else {
      
      # Return results for each category
      results <- data.frame()
      
      # Calculate the baseline effect 
      trait_col <- paste0(trait_prefix, "_", trait_value)
      category_param <- paste0(env_name, ":", trait_cols)
      
      # main effect
      B_e <- param_draws[, env_name]
      # Fourth corner coef
      B_I <- param_draws[,grepl(paste0(env_name, ':'), colnames(param_draws))] 
      
      # Set all trait categories for trait to 0
      t_j_base <- traits_mean
      t_j_base[trait_cols] <- 0 
      
      # Calculate effect
      B_j <- B_e +  B_I %*% ((t_j_base - traits_mean)/traits_sd)
      
      baseline_result <- data.frame(
        env_var = env_name,
        trait = trait_prefix,
        trait_category = "baseline",
        trait_value = "all_zero",
        median = quantile(B_j, 0.5),
        upper_95 = quantile(B_j, 0.975),
        lower_95 = quantile(B_j, 0.025),
        upper_50 = quantile(B_j, 0.75),
        lower_50 = quantile(B_j, 0.25)
      )
      
      results <- rbind(results, baseline_result)
      
      # Calculate the effect for each category
      for (category in categories) {
        
        trait_col <- paste0(trait_prefix, "_", category)
        
        t_j <- t_j_base
        t_j[trait_col] <- 1 
        t_j <- (t_j - traits_mean) / traits_sd
        B_j <- B_e +  B_I %*% t_j
        
        category_result <- data.frame(
          env_var = env_name,
          trait = trait_prefix,
          trait_category = category,
          trait_value = 1,
          median = quantile(B_j, 0.5),
          upper_95 = quantile(B_j, 0.975),
          lower_95 = quantile(B_j, 0.025),
          upper_50 = quantile(B_j, 0.75),
          lower_50 = quantile(B_j, 0.25)
        )
        
        results <- rbind(results, category_result)
      }
      
      return(results)
    }
  } 
  
  else if (trait_type == "numeric") {
    # For numeric traits, there should be only one column
    if (length(trait_cols) > 1) {
      warning("Multiple columns found for numeric trait. Using the first one:", trait_cols[1])
      trait_col <- trait_cols[1]
    } else {
      trait_col <- trait_cols
    }
    
    # Calculate min, max, mean for the trait
    trait_min <- min(traits[, trait_col])
    trait_max <- max(traits[, trait_col])
    trait_mean <- mean(traits[, trait_col])
    # Get the main effect 
    B_e <- param_draws[, env_name]
    # Fourth corner coef
    B_I <- param_draws[,grepl(paste0(env_name, ':'), colnames(param_draws))] 

    # Set all trait to mean
    t_j_base <- traits_mean
    
    # If trait_value is specified, use only that value
    if (!is.null(trait_value)) {
      if (is.numeric(trait_value)) {
        # Use the specific value provided
        raw_values <- trait_value
        value_labels <- as.character(trait_value)
      } else if (trait_value %in% c("min", "mean", "max")) {
        raw_values <- trait_value
        value_labels <- trait_value
      } else {
        stop("For numeric traits, trait_value must be NULL, a numeric value, or one of 'min', 'mean', 'max'")
      }
    } else {
      # If no value specified, return min, mean, and max
      raw_values <- c("min", "mean", "max")
      value_labels <- c("min", "mean", "max")
    }
    
    # Initialize results data frame
    results <- data.frame()
    
    
    # Calculate effect for each requested value
    for (i in seq_along(raw_values)) {
      value_label <- value_labels[i]
      
      # Determine the actual value based on the label
      if (value_label == "min") {
        raw_value <- trait_min
      } else if (value_label == "mean") {
        raw_value <- trait_mean
      } else if (value_label == "max") {
        raw_value <- trait_max
      } else {
        raw_value <- as.numeric(value_label)
      }
      
      t_j_base[trait_cols] <- raw_value
      t_j <- (t_j_base - traits_mean)/traits_sd
      B_j <- B_e +  B_I %*% t_j
      
      trait_result <- data.frame(
        env_var = env_name,
        trait = trait_prefix,
        trait_category = value_label,
        trait_value = raw_value,
        median = quantile(B_j, 0.5),
        upper_95 = quantile(B_j, 0.975),
        lower_95 = quantile(B_j, 0.025),
        upper_50 = quantile(B_j, 0.75),
        lower_50 = quantile(B_j, 0.25)
      )
      
      results <- rbind(results, trait_result)
    }
    
    return(results)
  }
}

