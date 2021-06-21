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


# function to plot 6-contexts signatures
"plot.6.signatures" = function(beta) {
    x = data.table(
        signature = rep(rownames(beta), ncol(beta)), 
        cat = rep(colnames(beta), each = nrow(beta)), 
        value = c(beta)
    )
    
    glist <- list()
    for (i in 1:nrow(beta)) {
        plt <- ggplot(x[signature == rownames(beta)[i]]) + geom_bar(aes(x = cat, 
                                                                        y = value, fill = cat), stat = "identity", position = "identity") + 
            facet_wrap(~cat, nrow = 1, scales = "free_x") + theme(axis.text.x = element_text(angle = 90, 
                                                                                             hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
            ggtitle(rownames(beta)[i]) + theme(text = element_text(size=12),legend.position = "none") + 
            ylab("Frequency of mutations") + xlab("") + ylim(0,1)
        plt <- plt + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        glist[[i]] <- plt
    }
    p = grid.arrange(grobs = glist, ncol = ceiling(nrow(beta)/3))
    return(p)
}

# function to plot 96-contexts signatures
"plot.96.signatures" = function(beta) {
    
    x = data.table(
        signature = rep(rownames(beta), ncol(beta)), 
        cat = rep(colnames(beta), each = nrow(beta)), 
        value = c(beta)
    )
    
    
    x[,Context:=paste0(substr(cat,1,1),".",substr(cat,7,7))]
    x[,alt:=paste0(substr(cat,3,3),">",substr(cat,5,5))]
    glist <- list()
    for (i in 1:nrow(beta)) {
        plt <- ggplot(x[signature==rownames(beta)[i]]) + 
            geom_bar(aes(x=Context,y=value,fill=alt),stat="identity",position="identity") + 
            facet_wrap(~alt,nrow=1,scales="free_x") + 
            theme(axis.text.x=element_text(angle=90,hjust=1),panel.background=element_blank(),axis.line=element_line(colour="black")) + 
            ggtitle(rownames(beta)[i]) + theme(legend.position="none") + ylab("Frequency of mutations")
        plt <- plt + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        glist[[i]] <- plt
    }
    p = grid.arrange(grobs = glist, ncol = ceiling(nrow(beta)/3))
    return(p)
}

trinucleotides <- c('A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'A[C>G]A', 'A[C>G]C', 
                    'A[C>G]G', 'A[C>G]T', 'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T', 
                    'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T', 'A[T>C]A', 'A[T>C]C', 
                    'A[T>C]G', 'A[T>C]T', 'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T', 
                    'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T', 'C[C>G]A', 'C[C>G]C', 
                    'C[C>G]G', 'C[C>G]T', 'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T', 
                    'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T', 'C[T>C]A', 'C[T>C]C', 
                    'C[T>C]G', 'C[T>C]T', 'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T', 
                    'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T', 'G[C>G]A', 'G[C>G]C', 
                    'G[C>G]G', 'G[C>G]T', 'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T', 
                    'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T', 'G[T>C]A', 'G[T>C]C', 
                    'G[T>C]G', 'G[T>C]T', 'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T', 
                    'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T', 'T[C>G]A', 'T[C>G]C', 
                    'T[C>G]G', 'T[C>G]T', 'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T', 
                    'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T', 'T[T>C]A', 'T[T>C]C', 
                    'T[T>C]G', 'T[T>C]T', 'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T')