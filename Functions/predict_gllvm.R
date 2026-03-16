pred.gllvm <- function(object, seed = NULL, Xnew = NULL, LV = NULL,
                       ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!is.null(Xnew)) {
    Xnew = Xnew
    nRows = dim(Xnew)[1]
    nCols = dim(object$y)[2]
  }else{
    nRows = 1
    nCols = dim(Xnew)[2]
    Xnew <- object$X[rep(1:nRows, nsim), ]
  }
  # PREDICTIONS
  if (is.null(object$X)) {
    prs = predict.gllvm(object, level = 0, type = "response")
    print('NULL X')
  }
  else if (is.null(object$TR)) {
    prs = predict.gllvm(object, newX = Xnew, level = 1,newLV = LV,
                        type = "response")
  }
  else {
    prs = predict.gllvm(object, newX = Xnew, level = 1, newLV = LV,
                        type = "response")
  }

    return(prs)
}



sim.gllvm <- function(object, nsim = 1, seed = NULL, conditional = FALSE, Xnew = NULL, LV = NULL, sim_dist = FALSE,
                      ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!is.null(Xnew)) {
    Xnew = Xnew
    nRows = dim(Xnew)[1]
    nCols = dim(object$y)[2]
  }else{
    nRows = 1
    nCols = dim(Xnew)[2]
    Xnew <- object$X[rep(1:nRows, nsim), ]
  }
  # PREDICTIONS
  if (is.null(object$X)) {
    prs = predict.gllvm(object, level = 0, type = "response")
    print('NULL X')
  }
  else if (is.null(object$TR)) {
    
    prs = predict.gllvm(object, newX = Xnew, level = 1, newLV = LV,
                        type = "response")
  }
  else {
    prs = predict.gllvm(object, newX = Xnew, level = 1, newLV = LV,
                        type = "response")
  }

  nTot = nsim * nRows * nCols
  
  if (object$family == "negative.binomial") {
    invPhis = matrix(rep(object$params$inv.phi, each = nsim *
                           nRows), ncol = nCols)
  }else if (object$family == "ZINB") {
    invPhis = matrix(rep(object$params$ZINB.inv.phi, each = nsim *
                           nRows), ncol = nCols)
    phis = matrix(rep(object$params$phi, each = nsim * nRows),
                  ncol = nCols)
  }else if (object$family == "tweedie") {
    phis = matrix(rep(object$params$phi, each = nsim * nRows),
                  ncol = nCols)
  }else if (object$family %in% c("gaussian", "gamma", "beta", "ZIP",
                                 "ZINB")) {
    phis = matrix(rep(object$params$phi, each = nsim * nRows),
                  ncol = nCols)
  }
  
  newDat = switch(
    object$family,
    binomial = rbinom(nTot, size = rep(object$Ntrials,each = nsim * nRows), prob = prs),
    poisson = rpois(nTot,prs),
    negative.binomial = rnbinom(nTot, size = invPhis, mu = prs),
    gaussian = rnorm(nTot, mean = prs, sd = phis),
    gamma = rgamma(nTot, shape = phis, scale = prs/phis),
    exponential = rexp(nTot, rate = 1/prs),
    tweedie = fishMod::rTweedie(nTot, mu = c(prs), phi = c(phis), p = object$Power),
    ordinal = sims,
    beta = rbeta(nTot, shape1 = phis * prs, shape2 = phis * (1 - prs)),
    ZIP = ifelse(rbinom(nTot, size = 1, prob = phis) > 0, 0, rpois(nTot, lambda = prs)),
    ZINB = ifelse(rbinom(nTot,size = 1, prob = phis) > 0, 0, rnbinom(nTot, size = invPhis, mu = prs)),
    stop(gettextf("family '%s' not implemented ", object$family), domain = NA))
  newDat = as.data.frame(matrix(newDat, ncol = nCols))
  try(colnames(newDat) <- colnames(prs), silent = TRUE)
  try(rownames(newDat) <- rownames(prs), silent = TRUE)
  try(dimnames(newDat) <- dimnames(prs), silent = TRUE)
  
  if (sim_dist) {
    return(newDat)
  } else {
    return(prs)
  }
}



pred_01.gllvm <- function(object, nsim = 1, seed = NULL, Xnew = NULL, LV = NULL,
                      ...) {

  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!is.null(Xnew)) {
    Xnew = Xnew
    nRows = dim(Xnew)[1]
    nCols = dim(object$y)[2]
  }else{
    nRows = 1
    nCols = dim(Xnew)[2]
    Xnew <- object$X[rep(1:nRows, nsim), ]
  }
  # PREDICTIONS
  if (is.null(object$X)) {
    prs = predict.gllvm(object, level = 0, type = "response")
    print('NULL X')
  } else if (is.null(object$TR)) {
    
    prs = predict.gllvm(object, newX = Xnew, level = 1, newLV = LV,
                        type = "response")
  } else {
    prs = predict.gllvm(object, newX = Xnew, level = 0, newLV = LV,
                        type = "response")
  }
  
  nTot = nsim * nRows * nCols

  if (object$family == "negative.binomial") {
    invPhis = matrix(rep(object$params$inv.phi, each = nsim *
                           nRows), ncol = nCols)
  }else if (object$family == "ZINB") {
    invPhis = matrix(rep(object$params$ZINB.inv.phi, each = nsim *
                           nRows), ncol = nCols)
    phis = matrix(rep(object$params$phi, each = nsim * nRows),
                  ncol = nCols)
  }else if (object$family == "tweedie") {
    phis = matrix(rep(object$params$phi, each = nsim * nRows),
                  ncol = nCols)
  }else if (object$family %in% c("gaussian", "gamma", "beta", "ZIP")) {
    phis = matrix(rep(object$params$phi, each = nsim * nRows),
                  ncol = nCols)
  }
  
  prs_01 <- switch(
    object$family,
    binomial = {
      1 - dbinom(0, size = object$Ntrials, prob = prs)
    },
    poisson = {
      1 - dpois(0, lambda = prs)
    },
    negative.binomial = {
      1 - dnbinom(0, size = invPhis, mu = prs)
    },
    ZIP = {
      1 - (phis + (1 - phis) * dpois(0, lambda = prs))
    },
    ZINB = {
      1 - (phis + (1 - phis) * dnbinom(0, size = invPhis, mu = prs))
    },
    stop(gettextf("family '%s' not implemented", object$family), domain = NA)
  )
  
  prs_01 = as.data.frame(matrix(prs_01, ncol = nCols))
  
  try(colnames(prs_01) <- colnames(prs), silent = TRUE)
  try(rownames(prs_01) <- rownames(prs), silent = TRUE)
  try(dimnames(prs_01) <- dimnames(prs), silent = TRUE)
  return(prs_01)
}
