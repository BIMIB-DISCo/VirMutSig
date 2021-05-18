# perform signatures discovery and rank estimation for a range of K signatures given a set of observations x
"signatures.decomposition" <- function( x, K, nmf_runs = 1000, num_processes = Inf, seed = NULL ) {
    
    # set the seed
    set.seed(seed)

    # setting up parallel execution
    pbackend <- NULL
    close_parallel <- FALSE
    if(is.na(num_processes) || is.null(num_processes)) {
        num_processes <- 1
    }
    else if(num_processes==Inf) {
        cores <- as.integer((detectCores()-1))
        num_processes <- min(cores,nmf_runs)
    }
    else {
        num_processes <- min(num_processes,nmf_runs)
    }
    if(num_processes==1) {
        pbackend <- "seq"
    }
    else {
        pbackend <- makeCluster(num_processes)
        res_clusterEvalQ <- clusterEvalQ(pbackend,library("NMF",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE))
        res_clusterEvalQ <- clusterEvalQ(pbackend,library("nnls",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE))
        clusterExport(pbackend,varlist=c("nmf.seed","nmf.nnls"),envir=environment())
        clusterSetRNGStream(pbackend,iseed=round(runif(1)*100000))
        close_parallel <- TRUE
        rm(res_clusterEvalQ)
        gc(verbose=FALSE)
    }

    cat("Performing signatures discovery and rank estimation...","\n")
    cat("Executing",num_processes,"processes in parallel...","\n")

    # perform signatures discovery and rank estimation
    alpha <- list()
    beta <- list()
    measures <- NULL
    for(i in 1:length(K)) {

        cat(paste0("Performing inference for K=",K[i],"..."),"\n")

        # perform the inference for current K
        results <- NULL
        while(is.null(results)) {
            results <- tryCatch({
                results <- nmf(x=x,rank=K[i],method=nmf.nnls,seed=nmf.seed,rng=round(runif(1)*10000),nrun=nmf_runs,.pbackend=pbackend)
                gc(verbose=FALSE)
                results
            }, error = function(e) {
                cat(paste0("An error has occurred: ",e$message),"\n")
                gc(verbose=FALSE)
                NULL
            }, finally = {
                gc(verbose=FALSE)
            })
        }

        alpha[[paste0(K[i],"_signatures")]] <- basis(results)
        beta[[paste0(K[i],"_signatures")]] <- coef(results)

        # compute and save quality measures
        curr_rss <- rss(results,x)
        curr_evar <- evar(results,x)
        curr_silhouette_alpha <- tryCatch({
            mean(silhouette(results,what="features")[,"sil_width"])
        }, error = function(e) {
            NA
        })
        curr_silhouette_beta <- tryCatch({
            mean(silhouette(results,what="samples")[,"sil_width"])
        }, error = function(e) {
            NA
        })
        curr_sparseness_alpha <- sparseness(as.vector(alpha[[paste0(K[i],"_signatures")]]))
        curr_sparseness_beta <- sparseness(as.vector(beta[[paste0(K[i],"_signatures")]]))
        if(nmf_runs>1) {
            curr_silhouette_consensus <- tryCatch({
                mean(silhouette(results,what="chc")[,"sil_width"])
            }, error = function(e) {
                NA
            })
            curr_measures <- matrix(c(K[i],curr_rss,curr_evar,curr_silhouette_alpha,curr_silhouette_beta,curr_sparseness_alpha,curr_silhouette_beta,cophcor(results),dispersion(results),curr_silhouette_consensus),nrow=1)
            colnames(curr_measures) <- c("Rank","Residual_Sum_Squares","Explained_Variance","Silhouette_Alpha","Silhouette_Beta","Sparseness_Alpha","Sparseness_Beta","Cophenetic_Coefficient","Dispersion_Coefficient","Silhouette_Consensus")
        }
        else {
            curr_measures <- matrix(c(K[i],curr_rss,curr_evar,curr_silhouette_alpha,curr_silhouette_beta,curr_sparseness_alpha,curr_silhouette_beta),nrow=1)
            colnames(curr_measures) <- c("Rank","Residual_Sum_Squares","Explained_Variance","Silhouette_Alpha","Silhouette_Beta","Sparseness_Alpha","Sparseness_Beta")
        }
        measures <- rbind(measures,curr_measures)

    }
    rownames(measures) <- 1:nrow(measures)

    # save results
    results <- list(alpha=alpha,beta=beta,measures=measures)
    
    # close pbackend
    if(close_parallel) {
        stopCluster(pbackend)
    }

    # return the discovered signatures
    return(results)
    
}

# initialize alpha and beta for nmf.nnls function
"nmf.seed" <- function( model, target ) {

    # initialize alpha with an empty matrix
    alpha <- array(NA,c(nrow(target),nbasis(model)))
    rownames(alpha) <- rownames(target)
    colnames(alpha) <- paste0("S",1:ncol(alpha))

    # randomly initialize beta
    beta <- matrix(runif(nbasis(model)*ncol(target)),nrow=nbasis(model),ncol=ncol(target))
    beta <- (beta/rowSums(beta)) # beta rows (i.e., signatures) must sum to 1
    if(any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))),] <- 0
    }
    rownames(beta) <- colnames(alpha)
    colnames(beta) <- colnames(target)

    # update results
    basis(model) <- alpha
    coef(model) <- beta

    # return updated model
    return(model)

}

# perform NMF by Non-negative least squares
"nmf.nnls" <- function( x, seed ) {

    # initialization
    alpha <- basis(seed) # exposures matrix
    beta <- coef(seed) # signatures matrix
    n <- nrow(x) # n is the number of observations in x, i.e., the samples
    J <- ncol(x) # J is the number of categories, i.e., the contexts
    K <- nrow(beta) # K is the number of signatures to be fitted

    # iteratively fit alpha and beta by Non-negative least squares (nnls)
    for(i in 1:20) {

        # update alpha, beta is kept fixed
        for(j in 1:n) {
            alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
        }

        # update beta, alpha is kept fixed
        for(k in 1:J) {
            beta[,k] <- nnls(alpha,as.vector(x[,k]))$x
        }

    }

    # normalize beta and perform final fit of alpha
    beta <- (beta/rowSums(beta))
    if(any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))),] <- 0
    }
    for(j in 1:n) {
        alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
    }
    colnames(alpha) <- paste0("S",1:ncol(alpha))
    rownames(beta) <- colnames(alpha)

    # update results
    basis(seed) <- alpha
    coef(seed) <- beta

    # return updated seed
    return(seed)

}
