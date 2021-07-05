
######################################################
########## MCMC diagnostics of kstar #################
######################################################

## Name: kstar_MCMCdiag()
## Purpose: Perform MCMC diagnostics for kstar
## Input:
# kstar: the vector output of kstar from running the DPMPM model
# nrun: number of mcmc iterations used in running the DPMPM model
# burn: number of burn-in iterations used in running the DPMPM model
# thin: number of thinning used in running the DPMPM model
## Output:
# Traceplot: the traceplot of kstar post burn-in and thinning
# Autocorrplot: the autocorrelation plot of kstar post burn-in and thinning

kstar_MCMCdiag <- function(kstar, nrun, burn, thin){
  iters <- (nrun - burn) / thin
  kstar_thinned <- rep(NA, iters)
  for (s in 1:iters){
    kstar_thinned[s] <- kstar[(s - 1) * thin + 1]
  }
  kstar_thinned <- data.frame(kstar = kstar_thinned)
  
  p1 <- mcmc_trace(kstar_thinned)
  p2 <- mcmc_acf(kstar_thinned)
  
  return(list(Traceplot = p1, Autocorrplot = p2))
}


######################################################
############# DPMPM_nozeros_syn() ####################
######################################################

## Name: DPMPM_nozeros_syn()
## Purpose: Use DPMPM models to synthesize data where there are no structural zeros
## Input:
# X: data frame for the original data
# dj: a vector recording the number of categories of the variables
# nrun: number of mcmc iterations
# burn: number of burn-in iterations
# K: number of latent classes
# alpha and balpha: the hyperparameters in stick-breaking prior distribution for alpha
# m: number of imputations
# vars: the names of variables to be synthesized
# seed: choice of seed
## Output:
# impdata: m imputed datasets
# origdata: original data containing missing values
# alpha: save posterior draws of alpha, which can be used to check MCMC convergence
# kstar: saved number of occupied mixture components, which can be used to track whether K is large enough

DPMPM_nozeros_syn <- function(X, dj, nrun, burn, thin, K, aalpha, balpha, m, vars, seed, silent = TRUE){
  vars_all <- names(X)
  num_obs <- nrow(X)
  p <- dim(X)[2]
  origdata <- X
  
  ## start DP model
  model <- CreateModel(origdata, NULL, K, 0, aalpha, balpha, seed = seed)
  model$EnableTracer <- TRUE   ## this will allow the tracer to be enabled for the model project
  model$Run(0,burn,1, silent)
  
  ## initalize matrices to save parameter draws:
  ## alpha, and the number of occupied latent classes
  eff.sam <- (nrun-burn)/thin
  alpha <- matrix(rep(0,eff.sam),nrow=eff.sam,ncol=1)
  uniquezout <- matrix(rep(0,eff.sam),nrow=eff.sam,ncol=1)
  
  #model$traceable: gives a list of tracable objects     
  model$SetTrace(c("k_star"), nrun - burn) ## will reset tracer every time it is called
  ## alpha, and the number of occupied latent classes
  for (i in 1:(eff.sam-m)){
    ## run model for "thin" more iterations
    model$Run(0,thin,1,silent)
    
    ## store results
    output <- model$snapshot
    alpha[i]<-output$alpha
    uniquezout[i]<-length(unique(output$z))
  }
  
  #####################################################
  ###      generate synthetic data with DPMPM       ###
  #####################################################
  
  ## initiate matrix to store results
  synX <- data.frame(matrix(NA, nrow = dim(origdata)[1], ncol = dim(origdata)[2]))
  names(synX) <- vars_all
  
  synX_unsyn <- origdata %>% select(-vars)
  
  ## initiate list to store results
  syndata <- vector("list",m)
  
  set.seed(seed)
  ## loop from 1 to m, each iteration generates a synthetic dataset
  for (i in 1:m){
    ## run model for "thin" more iterations
    model$Run(0,thin,1,silent)
    
    ## store results
    output <- model$snapshot
    alpha[eff.sam-(m-i)] <- output$alpha
    uniquezout[eff.sam-(m-i)] <- length(unique(output$z))
    
    ## save multinomial probabilities and latent class assignments
    psiout_i <- output$psi
    zout_i <- output$z
    
    ## generate synthetic data record by record with DPMPM
    for (j in 1:num_obs){
      zj <- zout_i[j] + 1
      
      ## synthesize variables in vars
      for (k in vars){
        k_index <- which(vars_all == k)
        Xprob_k <- psiout_i[1:dj[k_index], zj, k_index]
        synX[j, k_index] <- which(rmultinom(1, 1, Xprob_k) == 1)
      }
    }
    
    ## transform synX to data frame syndata0
    syndata0 <- GetDataFrame(synX[, vars] - 1, origdata[, vars], cols = 1:length(vars))
    
    ## add unsynthesized columns
    syndata0[, names(synX_unsyn)] <- synX_unsyn
    
    ## store i-th synthetic dataset to the list: syndata
    syndata[[i]] = syndata0[, names(X)]
  }
  
  traced <- model$GetTrace() # return the lists of tracked values
  kstar <- traced[[1]]
  
  ## return results: syndata, origdata, alpha, uniquezout
  res <- list(syndata=syndata,
              origdata = origdata,
              alpha = alpha,
              #uniquezout = uniquezout,
              kstar = kstar
  )
  return(res)
}