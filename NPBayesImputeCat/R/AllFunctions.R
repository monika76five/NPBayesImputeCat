
######################################################
############# Function 1 #############################
######################################################

## Name: DPMPM_nozeros_imp()
## Purpose: Use DPMPM models to impute missing data where there are no structural zeros
## Input:
# X: data frame for the data containing missing values
# nrun: number of mcmc iterations
# burn: number of burn-in iterations
# K: number of latent classes
# alpha and balpha: the hyperparameters in stick-breaking prior distribution for alpha
# m: number of imputations
# seed: choice of seed
## Output:
# impdata: m imputed datasets
# origdata: original data containing missing values
# alpha: save posterior draws of alpha, which can be used to check MCMC convergence
# kstar: saved number of occupied mixture components, which can be used to track whether K is large enough


DPMPM_nozeros_imp <- function(X, nrun, burn, thin, K, aalpha, balpha, m, seed, silent = TRUE){
    origdata <- X
    
    ## start DP model
    model <- CreateModel(origdata, NULL, K, 0, aalpha, balpha, seed = seed)
    model$EnableTracer <- TRUE   ## this will allow the tracer to be enabled for the model project
    model$Run(0,burn,1,silent)
    
    ## initalize matrices to save parameter draws:
    ## alpha, and the number of occupied latent classes
    eff.sam <- (nrun-burn)/thin
    alpha <- matrix(rep(0,eff.sam),nrow=eff.sam,ncol=1)
    uniquezout <- matrix(rep(0,eff.sam),nrow=eff.sam,ncol=1)
    
    ## run the MCMC and save parameter draws:
    #model$traceable: gives a list of tracable objects     
    model$SetTrace(c("k_star"), nrun - burn) ## will reset tracer every time it is called
    ## alpha, and the number of occupied latent classes
    for (i in 1:(eff.sam-m)){
        ## run model for "thin" more iterations
        model$Run(0,thin,1,silent)
        
        ## extract output from model$snapshot
        output <- model$snapshot
        
        ## save alpha and uniquezout
        alpha[i] <- output$alpha
        uniquezout[i] <- length(unique(output$z))
    }
    
    #####################################################
    #####      impute missing data with DPMPM       #####
    #####################################################
    
    ## initiate matrix to store results
    impX <- matrix(NA, dim(origdata)[1], dim(origdata)[2])
    
    ## initiate list to store results
    impdata<-vector("list",m)
    
    ## loop from 1 to m, each iteration generates a synthetic dataset
    for (i in 1:m){
        ## run model for "thin" more iterations
        model$Run(0,thin,1,silent)
        
        ## store results
        output <- model$snapshot
        alpha[eff.sam-(m-i)] <- output$alpha
        uniquezout[eff.sam-(m-i)] <- length(unique(output$z))
        
        #retrieve parameters from the final iteration
        result <- model$snapshot
        
        #convert ImputedX matrix to dataframe, using proper factors/names etc.
        impX <- GetDataFrame(result$ImputedX,X)
        
        ## store i-th synthetic dataset to the list: impdata
        impdata[[i]] <- impX
    }
    
    traced <- model$GetTrace() # return the lists of tracked values
    kstar <- traced[[1]]
    
    ## return results: impdata, origdata, alpha, uniquezout
    res <- list(impdata=impdata,
    origdata = origdata,
    alpha = alpha,
    #uniquezout = uniquezout,
    kstar = kstar
    )
    return(res)
}

######################################################
############# Function 2 #############################
######################################################

## Name: DPMPM_zeros_imp()
## Purpose: Use DPMPM models to impute missing data where there are structural zeros
## Input:
# X: data frame for the data containing missing values
# MCZ: data frame containing the structural zeros definition
# Nmax: an upper truncation limit for the augmented sample size;
# nrun: number of mcmc iterations
# burn: burn-in
# K: number of latent classes
# m: number of imputations
# seed: choice of seed
## Output:
# impdata: m imputed datasets
# origdata: original data containing missing values
# alpha: saved posterior draws of alpha, which can be used to check MCMC convergence
# kstar: saved number of occupied mixture components, which can be used to track whether K is large enough
# Nmis: saved posterior draws of the augmented sample size, which can be used to check MCMC convergence


DPMPM_zeros_imp <- function(X, MCZ, Nmax, nrun, burn, thin, K, aalpha, balpha, m, seed, silent = TRUE){
    origdata <- X
    
    ## start DP model
    model <- CreateModel(origdata, MCZ, K, Nmax, aalpha, balpha, seed = seed)
    model$EnableTracer <- TRUE   ## this will allow the tracer to be enabled for the model project
    model$Run(0,burn,1,silent)
    
    ## initalize matrices to save parameter draws:
    ## alpha, and the number of occupied latent classes
    eff.sam <- (nrun-burn)/thin
    alpha <- matrix(rep(0,eff.sam),nrow=eff.sam,ncol=1)
    uniquezout <- matrix(rep(0,eff.sam),nrow=eff.sam,ncol=1)
    
    ## run the MCMC and save parameter draws:
    #model$traceable: gives a list of tracable objects     
    model$SetTrace(c("k_star","Nmis"), nrun - burn) ## will reset tracer every time it is called
    ## alpha, and the number of occupied latent classes
    for (i in 1:(eff.sam-m)){
        ## run model for "thin" more iterations
        model$Run(0,thin,1,silent)
        
        ## extract output from model$snapshot
        output <- model$snapshot
        
        ## save alpha and uniquezout
        alpha[i] <- output$alpha
        uniquezout[i] <- length(unique(output$z))
    }
    
    #####################################################
    #####      impute missing data with DPMPM       #####
    #####################################################
    
    ## initiate matrix to store results
    impX <- matrix(NA, dim(origdata)[1], dim(origdata)[2])
    
    ## initiate list to store results
    impdata <- vector("list",m)
    
    ## loop from 1 to m, each iteration generates a synthetic dataset
    for (i in 1:m){
        ## run model for "thin" more iterations
        model$Run(0,thin,1,silent)
        
        ## store results
        output <- model$snapshot
        alpha[eff.sam-(m-i)] <- output$alpha
        uniquezout[eff.sam-(m-i)] <- length(unique(output$z))
        
        #retrieve parameters from the final iteration
        result <- model$snapshot
        
        #convert ImputedX matrix to dataframe, using proper factors/names etc.
        impX <- GetDataFrame(result$ImputedX,X)
        
        ## store i-th synthetic dataset to the list: syndata
        impdata[[i]] <- impX
    }
    
    traced <- model$GetTrace() # return the lists of tracked values
    kstar <- traced[[1]]
    Nmis <- traced[[2]]
    
    ## return results: syndata, origdata, alpha, uniquezout
    res <- list(impdata=impdata,
    origdata = origdata,
    alpha = alpha,
    #uniquezout = uniquezout,
    kstar = kstar,
    Nmis = Nmis
    )
    return(res)
}

######################################################
############# Function 3 #############################
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
        syndata[[i]] = syndata0
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

######################################################
############# Function 4 #############################
######################################################

## Name: compute_probs()
## Purpose: Estimating marginal and joint probabilities in imputed or synthetic datasets
## Input:
# InputData: a list of imputed or synthetic datasets
# varlist: a list of variable names (or combination of names) to evaluate (marginal or joint) probabilities for
## Output:
# Results: a list of marginal and joint probabilities in imputed or synthetic datasets

compute_probs <- function(InputData, varlist){
    m <- length(InputData)
    Results <- vector("list", m)
    
    for (l in 1:m){
        Data_l <- InputData[[l]]
        prob_table <- vector("list",length(varlist))
        for(pp in 1:length(varlist)){
            prob_table[[pp]] <- as.data.frame(table(Data_l[,varlist[[pp]]]))
            if(length(varlist[[pp]])==1){
                colnames(prob_table[[pp]]) <- c(varlist[[pp]],"Freq")
            }
            prob_table[[pp]]$Qm <- prob_table[[pp]]$Freq/sum(prob_table[[pp]]$Freq)
            prob_table[[pp]]$Um <- prob_table[[pp]]$Qm*(1-prob_table[[pp]]$Qm)/sum(prob_table[[pp]]$Freq)
            prob_table[[pp]] <- prob_table[[pp]][,!colnames(prob_table[[pp]])=="Freq"]
        }
        Results[[l]]  <- prob_table
    }
    return(Results)
}

######################################################
############# Function 5 #############################
######################################################

## Name: fit_GLMs()
## Purpose: Fit GLM models for imputed or synthetic datasets
## Input:
# InputData: a list of imputed or synthetic datasets
# exp: GLM expression (for polr and nnet, those libraries should be loaded first)
## Output:
# Results: a list of GLM results

fit_GLMs <- function(InputData, exp){
    m <- length(InputData)
    Results <- as.list(seq_len(m))
    for (l in 1:m){
        Data_l <- InputData[[l]]
        Model_l <- eval_tidy(enquo(exp), Data_l)
        
        if(is.matrix(coef(Model_l))){ #check class instead?
            Results[[l]] <- cbind(expand.grid(Levels=rownames(coef(Model_l)),Parameter=colnames(coef(Model_l))),
            Qm=c(coef(Model_l)),
            Um=c(summary(Model_l)$standard.errors)^2)
            Results[[l]] <- (Results[[l]])[order(Results[[l]][,1]),]
            rownames(Results[[l]]) <- NULL
        } else {
            Results[[l]] <- cbind(Qm=c(coef(Model_l)),
            Um=c(coef(summary(Model_l))[,"Std. Error"])^2)
            rownames(Results[[l]]) <- names(coef(Model_l))
        }
    }
    return(Results)
}

######################################################
############# Function 6 #############################
######################################################

## Name: pool_estimated_probs()
## Purpose: Pool probability estimates from imputed or synthetic datasets
## Input:
# ComputeProbsResults: output from the compute_probs function
# method: choose between "imputation", "synthesis_full", "synthesis_partial"
## Output:
# Results: a list of marginal and joint probability results after combining rules

pool_estimated_probs <- function(ComputeProbsResults,
method = c("imputation", "synthesis_full", "synthesis_partial")){
    m <- length(ComputeProbsResults)
    Results <- lapply(1:length(ComputeProbsResults[[m]]),
    function(x) (ComputeProbsResults[[m]][[x]])[,!is.element(colnames(ComputeProbsResults[[m]][[x]]),c("Um"))])
    
    if (method == "imputation"){
        for(pp in 1:length(Results)){
            Results[[pp]]$Qm <- 0
            colnames(Results[[pp]])[colnames(Results[[pp]]) == "Qm"] <- "Estimate"
            
            Qm <- NULL
            Um <- NULL
            for (l in 1:m){
                Qm <- cbind(Qm,(ComputeProbsResults[[l]][[pp]])[,"Qm"])
                Um <- cbind(Um,(ComputeProbsResults[[l]][[pp]])[,"Um"])
            }
            
            Qbarm <- rowMeans(Qm)
            Bm <- apply(Qm,1,var)
            Ubarm <- rowMeans(Um)
            Tm <- (1+1/m)*Bm + Ubarm
            v <- (m-1)*(1+Ubarm/(1+1/m)*Bm)^2
            CI_Lower <- Qbarm - qt(0.975, v)*sqrt(Tm)
            CI_Upper <- Qbarm + qt(0.975, v)*sqrt(Tm)
            
            Results[[pp]]$Estimate <- Qbarm
            Results[[pp]] <- cbind(Results[[pp]],
            Std.Error=sqrt(Tm),
            Df=v,
            Statistic=(Qbarm/sqrt(Tm)),
            CI_Lower=CI_Lower, CI_Upper=CI_Upper)
        }
    }
    
    else if (method == "synthesis_full"){
        for(pp in 1:length(Results)){
            Results[[pp]]$Qm <- 0
            colnames(Results[[pp]])[colnames(Results[[pp]]) == "Qm"] <- "Estimate"
            
            Qm <- NULL
            Um <- NULL
            for (l in 1:m){
                Qm <- cbind(Qm,(ComputeProbsResults[[l]][[pp]])[,"Qm"])
                Um <- cbind(Um,(ComputeProbsResults[[l]][[pp]])[,"Um"])
            }
            
            Qbarm <- rowMeans(Qm)
            Bm <- apply(Qm,1,var)
            Ubarm <- rowMeans(Um)
            Tf <- (1+1/m)*Bm - Ubarm
            Tf[Tf < 0] <- Ubarm[Tf < 0]
            v <- (m-1)*(1 - (m*Ubarm/((m+1)*Bm)) )^2
            CI_Lower <- Qbarm - qt(0.975, v)*sqrt(Tf)
            CI_Upper <- Qbarm + qt(0.975, v)*sqrt(Tf)
            
            Results[[pp]]$Estimate <- Qbarm
            Results[[pp]] <- cbind(Results[[pp]],
            Std.Error=sqrt(Tf),
            Df=v,
            Statistic=(Qbarm/sqrt(Tf)),
            CI_Lower=CI_Lower, CI_Upper=CI_Upper)
        }
    }
    
    else {
        for(pp in 1:length(Results)){
            Results[[pp]]$Qm <- 0
            colnames(Results[[pp]])[colnames(Results[[pp]]) == "Qm"] <- "Estimate"
            
            Qm <- NULL
            Um <- NULL
            for (l in 1:m){
                Qm <- cbind(Qm,(ComputeProbsResults[[l]][[pp]])[,"Qm"])
                Um <- cbind(Um,(ComputeProbsResults[[l]][[pp]])[,"Um"])
            }
            
            Qbarm <- rowMeans(Qm)
            Bm <- apply(Qm,1,var)
            Ubarm <- rowMeans(Um)
            Tf <- Bm/m + Ubarm
            v <- (m-1)*(1+Ubarm/(Bm/m))^2
            CI_Lower <- Qbarm - qt(0.975, v)*sqrt(Tf)
            CI_Upper <- Qbarm + qt(0.975, v)*sqrt(Tf)
            
            Results[[pp]]$Estimate <- Qbarm
            Results[[pp]] <- cbind(Results[[pp]],
            Std.Error=sqrt(Tf),
            Df=v,
            Statistic=(Qbarm/sqrt(Tf)),
            CI_Lower=CI_Lower, CI_Upper=CI_Upper)
        }
    }
    return(Results)
}

######################################################
############# Function 7 #############################
######################################################

## Name: pool_fitted_GLMs()
## Purpose: Pool estimates of fitted GLM models in imputed or synthetic datasets
## Input:
# GLMResults: output from the fit_GLMs function
# method: choose between "imputation", "synthesis_full", "synthesis_partial"
## Output:
# Results: a list of GLM results after combining rules

pool_fitted_GLMs <- function(GLMResults,
method = c("imputation", "synthesis_full", "synthesis_partial")){
    m <- length(GLMResults)
    Results <- (GLMResults[[m]])[,!is.element(colnames(GLMResults[[m]]),c("Qm","Um"))]
    
    if (method == "imputation"){
        Qm <- NULL
        Um <- NULL
        for (l in 1:m){
            Qm <- cbind(Qm,(GLMResults[[l]])[,"Qm"])
            Um <- cbind(Um,(GLMResults[[l]])[,"Um"])
        }
        
        Qbarm <- rowMeans(Qm)
        Bm <- apply(Qm,1,var)
        Ubarm <- rowMeans(Um)
        Tm <- (1+1/m)*Bm + Ubarm
        v <- (m-1)*(1+Ubarm/(1+1/m)*Bm)^2
        CI_Lower <- Qbarm - qt(0.975, v)*sqrt(Tm)
        CI_Upper <- Qbarm + qt(0.975, v)*sqrt(Tm)
        
        Results <- cbind(Results,
        Estimate=Qbarm,
        Std.Error=sqrt(Tm),
        Df=v,
        Statistic=(Qbarm/sqrt(Tm)),
        CI_Lower=CI_Lower, CI_Upper=CI_Upper)
    }
    
    else if (method == "synthesis_full"){
        Qm <- NULL
        Um <- NULL
        for (l in 1:m){
            Qm <- cbind(Qm,(GLMResults[[l]])[,"Qm"])
            Um <- cbind(Um,(GLMResults[[l]])[,"Um"])
        }
        
        Qbarm <- rowMeans(Qm)
        Bm <- apply(Qm,1,var)
        Ubarm <- rowMeans(Um)
        Tf <- (1+1/m)*Bm - Ubarm
        Tf[Tf < 0] <- Ubarm[Tf < 0]
        v <- (m-1)*(1-m*Ubarm/((m+1)*Bm))^2
        CI_Lower <- Qbarm - qt(0.975, v)*sqrt(Tf)
        CI_Upper <- Qbarm + qt(0.975, v)*sqrt(Tf)
        
        Results <- cbind(Results,
        Estimate=Qbarm,
        Std.Error=sqrt(Tf),
        Df=v,
        Statistic=(Qbarm/sqrt(Tf)),
        CI_Lower=CI_Lower, CI_Upper=CI_Upper)
    }
    
    else {
        Qm <- NULL
        Um <- NULL
        for (l in 1:m){
            Qm <- cbind(Qm,(GLMResults[[l]])[,"Qm"])
            Um <- cbind(Um,(GLMResults[[l]])[,"Um"])
        }
        
        Qbarm <- rowMeans(Qm)
        Bm <- apply(Qm,1,var)
        Ubarm <- rowMeans(Um)
        Tf <- Bm/m + Ubarm
        v <- (m-1)*(1+Ubarm/(Bm/m))^2
        CI_Lower <- Qbarm - qt(0.975, v)*sqrt(Tf)
        CI_Upper <- Qbarm + qt(0.975, v)*sqrt(Tf)
        
        Results <- cbind(Results,
        Estimate=Qbarm,
        Std.Error=sqrt(Tf),
        Df=v,
        Statistic=(Qbarm/sqrt(Tf)),
        CI_Lower=CI_Lower, CI_Upper=CI_Upper)
    }
    return(Results)
}

######################################################
############# Function 8 #############################
######################################################

## Name: marginal_compare_all_imp()
## Purpose: Plot estimated marginal probabilities from observed data vs imputed datasets
## Input:
# obsdata: the observed data
# impdata: the list of m imputed datasets
# vars: the variable of interest
## Output:
# Plot: the barplot
# Comparison: a table of marginal probabilies from observed data vs imputed datasets

marginal_compare_all_imp <- function(obsdata, impdata, vars){
    n <- dim(obsdata)[1]
    m <- length(impdata)
    per.tables <- data.frame(table(obsdata[, vars]) / sum(table(obsdata[, vars])) * 100)
    names(per.tables) <- c("level", "observed")
    level <-NULL
    value <- NULL
    variable <- NULL
    
    for (l in 1:m){
        impdata_l <- impdata[[l]]
        old_names <- names(per.tables)
        per.tables <- data.frame(cbind(per.tables, as.numeric(table(impdata_l[, vars])) / n * 100))
        names(per.tables) <- c(old_names, paste0("imputed", l))
    }
    
    per.tables_long <- melt(per.tables)
    per.tables_long$name <- rep(vars, dim(per.tables_long)[1])
    p <- ggplot(data = per.tables_long, aes(x = level, y = value, fill = variable)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme_bw() + theme(legend.position="top") + theme(legend.title = element_blank()) +
    theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1, size = 10)) +
    scale_fill_manual(values=c("#E69F00", rep("#999999", m))) +
    facet_wrap(~ name) +
    xlab("") + ylab("Percent")
    
    return(list(Plot = p, Comparison = per.tables))
}

######################################################
############# Function 9 #############################
######################################################

## Name: marginal_compare_all_syn()
## Purpose: Plot estimated marginal probabilities from observed data vs synthetic datasets
## Input:
# obsdata: the observed data
# syndata: the list of m imputed datasets
# vars: the variable of interest
## Output:
# Plot: the barplot
# Comparison: a table of marginal probabilies from observed data vs imputed datasets

marginal_compare_all_syn <- function(obsdata, syndata, vars){
    
    n <- dim(obsdata)[1]
    m <- length(syndata)
    per.tables <- data.frame(table(obsdata[, vars]) / n * 100)
    names(per.tables) <- c("level", "original")
    level <-NULL
    value <- NULL
    variable <- NULL
    
    for (l in 1:m){
        syndata_l <- syndata[[l]]
        old_names <- names(per.tables)
        per.tables <- data.frame(cbind(per.tables, as.numeric(table(syndata_l[, vars])) / n * 100))
        names(per.tables) <- c(old_names, paste0("synthetic", l))
    }
    
    per.tables_long <- melt(per.tables)
    per.tables_long$name <- rep(vars, dim(per.tables_long)[1])
    p <- ggplot(data = per.tables_long, aes(x = level, y = value, fill = variable)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme_bw() + theme(legend.position="top") + theme(legend.title = element_blank()) +
    theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1, size = 10)) +
    scale_fill_manual(values=c("#E69F00", rep("#999999", m))) +
    facet_wrap(~ name) +
    xlab("") + ylab("Percent")
    
    return(list(Plot = p, Comparison = per.tables))
}

