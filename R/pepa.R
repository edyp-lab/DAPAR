##' fudge2LRT: heuristic to choose the value of the hyperparameter
##' (fudge factor) used to regularize the variance estimator in the
##' likelihood ratio statistic (as implemented in samLRT). We follow
##' the heuristic described in [1] and adapt the code of the fudge2
##' function in the siggene R package.

##' [1] Tusher, Tibshirani and Chu, Significance analysis of
##' microarrays applied to the ionizing radiation response, PNAS 2001
##' 98: 5116-5121, (Apr 24).
##' 
##' @title Heuristic to choose the value of the hyperparameter
##' (fudge factor) used to regularize the variance estimator in the
##' likelihood ratio statistic
##' @param lmm.res.h0 a vector of object containing the estimates (used to
##' compute the statistic) under H0 for each connected component. If
##' the fast version of the estimator was used (as implemented in this
##' package), lmm.res.h0 is a vector containing averages of squared
##' residuals. If a fixed effect model was used, it is a vector of lm
##' objects and if a mixed effect model was used it is a vector or lmer
##' object.
##' @param lmm.res.h1 similar to lmm.res.h0, a vector of object containing
##' the estimates (used to compute the statistic) under H1 for each
##' protein.
##' @param cc a list containing the indices of peptides and proteins
##' belonging to each connected component.
##' @param n the number of samples used in the test
##' @param p the number of proteins in the experiment
##' @param s a vector containing the maximum likelihood estimate of the
##' variance for the chosen model. When using the fast version of the
##' estimator implemented in this package, this is the same thing as
##' the input lmm.res.h1. For other models (e.g. mixed models) it can
##' be obtained from samLRT.
##' @param alpha A vector of proportions used to build candidate values for
##' the regularizer. We use quantiles of s with these
##' proportions. Default to seq(0, 1, 0.05)
##' @param include.zero logical value indicating if 0 should be included in
##' the list of candidates. Default to TRUE.
##' @return (same as the fudge2 function of siggene):
##' s.zero: the value of the fudge factor s0.
##' alpha.hat: the optimal quantile of the 's' values. If s0=0, 'alpha.hat'
##' will not be returned.
##' vec.cv: the vector of the coefficients of variations.  Following
##'           Tusher et al. (2001), the optimal 'alpha' quantile is given
##'           by the quantile that leads to the smallest CV of the modified
##'           test statistics.
##' msg: a character string summarizing the most important information
##'           about the fudge factor.
##' @author Thomas Burger, Laurent Jacob
fudge2LRT <- function(lmm.res.h0, lmm.res.h1, cc, n, p, s, alpha = seq(0, 1, 0.05), include.zero = TRUE)
{
    
    if (max(alpha) > 1 || min(alpha) < 0)
        stop("alpha has to be between 0 and 1")
    
    if (any(round(100 * alpha, 10) != round(100 * alpha, 0))) {
        warning("At least one alpha is not a percentile. Only the first two decimal digits",
                " are retained.")
        alpha <- signif(alpha, 2)
    }
    
    if (length(alpha) == 1) {
        s.zero <- quantile(s, alpha)
        msg <- paste("s0 =", round(s.zero, 4), " (The", 100 *
                         alpha, "% quantile of the s values.) \n \n")
        
        invisible(return(list(s.zero = s.zero, alpha.hat = alpha,
                              vec.cv = NULL, msg = msg)))
        
    }
    
    fudge.quan <- quantile(s, alpha)
    fudge.quan <- sort(c(fudge.quan, 2*fudge.quan), decreasing=FALSE)
    if (include.zero)
        fudge.quan <- c(0, fudge.quan)
    
    n.alpha <- length(fudge.quan)
    ## n.pep <- rep(NA, p)
    ## for(ee in cc){
    ##     ee <- as.numeric(ee)
    ##     pp <- ee[ee <= p]   
    ##     n.pep[pp] <- sum(ee > p)
    ## }
    
    d.mat <- sapply(fudge.quan, FUN=function(s1){
        sam.res <- samLRT(lmm.res.h0, lmm.res.h1, cc, n, p, s1)
        ## (n*n.pep) * exp(llr)^(-(2/(n*n.pep))) # Hotelling T2
        exp(sam.res$llr.sam)^(-(2/(sam.res$sample.sizes))) # Hotelling T2
    })
    
    # r/outer(s, fudge.quan, "+")
    ## print(str(d.mat))

    n.uni.s <- length(unique(s))

    if (n.uni.s < 25)
        stop("For the computation of the fugde factor,", "\n",
             "there should be at least 25 genes with differing standard deviations.")

    n.int <- ifelse(n.uni.s > 500, 101, floor(n.uni.s/5))
    quan <- quantile(s, seq(0, 1, le = n.int))
    quan <- unique(round(quan, 8))
    n.int <- length(quan)
    int.s <- as.numeric(cut(s, quan, include.lowest = TRUE, right = FALSE))
    ## print(cbind(int.s, s))
    
    mad.mat <- matrix(0, n.int - 1, ncol(d.mat))
    for (i in 1:(n.int - 1)) {
        mad.mat[i, ] <- apply(d.mat[which(int.s == i), , drop = FALSE],2, mad)
    }
    
    cv <- function(x) {
        sd(x)/mean(x)
    }
    
    vec.cv <- apply(mad.mat, 2, cv)
    ## print(fudge.quan)
    ## x11()
    ## plot(fudge.quan, vec.cv)
    which.min <- which(vec.cv == min(vec.cv))
    if (include.zero & which.min == 1) {
        msg <- "s0 = 0 \n \n"
        s.zero <- 0
        invisible(return(list(s.zero = s.zero, vec.cv = vec.cv,msg = msg)))
    }
    
    s.zero <- fudge.quan[which.min]
    print(s.zero)
    if (include.zero)
        which.min <- which.min - 1
    
    alpha.hat <- alpha[which.min]
    msg <- paste("s0 =", round(s.zero, 4), " (The", 100 * alpha.hat,
                 "% quantile of the s values.)", "\n", "\n")
    invisible(return(list(alpha.hat = alpha.hat, s.zero = s.zero,
                          vec.cv = vec.cv, msg = msg)))
}



##' This function is xxxxxxx
##' 
##' @title xxxxxx
##' @param X an n.pep*n.prot indicator matrix.
##' @param y1 n.pep*n.samples matrice giving the observed counts for
## each peptide in each sample from the condition 1
##' @param y2  n.pep*n.samples matrice giving the observed counts for
## each peptide in each sample from the condition 2 
##' @return xxxxxxxxxx..
##' @author Thomas Burger, Laurent Jacob
LH0 <- function(X, y1, y2){
    
    n <- ncol(y1)+ncol(y2)
    
    ## Ytilde <- matrix(c(as.vector(y1), as.vector(y2)), ncol=1)
    ## ss <- norm(Ytilde, 'F')^2 - n*norm(matrix(rowMeans(cbind(y1, y2)), ncol=1), 'F')^2
    
    ss <- norm(y1, 'F')^2 + norm(y2, 'F')^2 - n*norm(matrix(rowMeans(cbind(y1, y2)), ncol=1), 'F')^2
    
    return(list(ss=ss))
}




##' This function is xxxxxxx
##' 
##' @title xxxxxx
##' @param X an n.pep*n.prot indicator matrix.
##' @param y1 n.pep*n.samples matrice giving the observed counts for
## each peptide in each sample from the condition 1
##' @param y2  n.pep*n.samples matrice giving the observed counts for
## each peptide in each sample from the condition 2 
##' @param j the index of the protein being tested, ie which has different
## expression in the two conditions under H1
##' @return xxxxxxxxxx..
##' @author Thomas Burger, Laurent Jacob
LH1 <- function(X, y1, y2, j){
    n1 <- ncol(y1)
    n2 <- ncol(y2)
    n <- n1 + n2
    xj <- X[, j, drop=FALSE]
    
    ## Ytilde <- matrix(c(as.vector(y1), as.vector(y2)), ncol=1)    
    ## ss <- norm(Ytilde, 'F')^2 - n*sum(rowMeans(cbind(y1, y2))^2)- (sum((xj/norm(xj, 'F')) * (rowMeans(y1) - rowMeans(y2)))^2)*(n1*n2/n)
    
    ss <- norm(y1, 'F')^2 + norm(y2, 'F')^2 - n*sum(rowMeans(cbind(y1, y2))^2)- (sum((xj/norm(xj, 'F')) * (rowMeans(y1) - rowMeans(y2)))^2)*(n1*n2/n)
    
    return(list(ss=ss))
}




##' This function computes a regularized version of the likelihood ratio
##' statistic. The regularization adds a user-input fudge factor s1 to
##' the variance estimator. This is straightforward when using a fixed
##' effect model (cases 'numeric' and 'lm') but requires some more care
##' when using a mixed model.

##' 
##' @title xxxxxx
##' @param lmm.res.h0 a vector of object containing the estimates (used to
##' compute the statistic) under H0 for each connected component. If
##' the fast version of the estimator was used (as implemented in this
##' package), lmm.res.h0 is a vector containing averages of squared
##' residuals. If a fixed effect model was used, it is a vector of lm
##' objects and if a mixed effect model was used it is a vector or lmer
##' object.
##' @param lmm.res.h1 similar to lmm.res.h0, a vector of object containing
##' the estimates (used to compute the statistic) under H1 for each
##' protein.
##' @param cc a list containing the indices of peptides and proteins
##' belonging to each connected component. 
##' @param n the number of samples used in the test
##' @param p  the number of proteins in the experiment 
##' @param s1 the fudge factor to be added to the variance estimate
##' @return llr.sam: a vector of numeric containing the regularized log
##' likelihood ratio statistic for each protein.
##' s: a vector containing the maximum likelihood estimate of the
##' variance for the chosen model. When using the fast version of the
##' estimator implemented in this package, this is the same thing as
##' the input lmm.res.h1.
##' lh1.sam: a vector of numeric containing the regularized log
##' likelihood under H1 for each protein.
##' lh0.sam: a vector of numeric containing the regularized log
##' likelihood under H0 for each connected component.
##' sample.sizes: a vector of numeric containing the sample size
##' (number of biological samples times number of peptides) for each
##' protein. This number is the same for all proteins within each
##' connected component.
##' @author Thomas Burger, Laurent Jacob
samLRT <- function(lmm.res.h0, lmm.res.h1, cc, n, p, s1){
    s <- lh1.sam <- llr.sam <- rep(NA, p)
    lh0.sam <- rep(NA, length(cc))
    sample.sizes <- rep(NA, p)
    ii <- 0
    for(ee in cc){
        ii <- ii + 1
        ee <- as.numeric(ee)
        local.prot <- ee[ee <= p]    
        local.pep <- ee[ee > p] - p
        if(class(lmm.res.h0) == 'numeric'){
            lh0.sam[ii] <- (lmm.res.h0[ii] + s1) * (length(local.pep)*n)
        }else{
            if(class(lmm.res.h0[[ii]]) == 'lm'){
                lh0.sam[ii] <- (mean(residuals(lmm.res.h0[[ii]])^2) + s1) * (length(local.pep)*n)
            }else{
                devcomp <- lme4::getME(lmm.res.h0[[ii]], 'devcomp')
                sigmaML2 <- devcomp$cmp['pwrss']/devcomp$dims['n']
                lh0.sam[ii] <- -(devcomp$cmp['ldL2'] + devcomp$dims['n'] * (1 + log(2 * pi * (sigmaML2 + s1))))/2
            }
        }
        for(jj in local.prot){
            sample.sizes[jj] <- (length(local.pep)*n)
            if(class(lmm.res.h1) == 'numeric'){
                lh1.sam[jj] <- (lmm.res.h1[jj] + s1) * (length(local.pep)*n)
                llr.sam[jj] <- (length(local.pep)*n) * (log(lh0.sam[ii]) - log(lh1.sam[jj]))
                s[jj] <- lmm.res.h1[jj]
            }else{
                if(class(lmm.res.h1[[jj]]) == 'lm'){
                    s[jj] <- mean(residuals(lmm.res.h1[[jj]])^2)
                    lh1.sam[jj] <- (mean(residuals(lmm.res.h1[[jj]])^2) + s1) * (length(local.pep)*n)
                    llr.sam[jj] <- (length(local.pep)*n) * (log(lh0.sam[ii]) - log(lh1.sam[jj]))
                }else{
                    devcomp <- lme4::getME(lmm.res.h1[[jj]], 'devcomp')
                    sigmaML2 <- devcomp$cmp['pwrss']/devcomp$dims['n']
                    s[jj] <- sigmaML2
                    lh1.sam[jj] <- -(devcomp$cmp['ldL2'] + devcomp$dims['n'] * (1 + log(2 * pi * (sigmaML2 + s1))))/2
                    llr.sam[jj] <- 2 * (lh1.sam[jj] - lh0.sam[ii])
                }
            }
        }
    }
    return(list(llr.sam=llr.sam, s=s, lh1.sam=lh1.sam, lh0.sam=lh0.sam, sample.sizes=sample.sizes))
}






##' This function is PEptide based Protein differential Abundance test
##' 
##' @title PEptide based Protein differential Abundance test
##' @param X Binary q x p design matrix for q peptides and p
## proteins. X_(ij)=1 if peptide i belongs to protein j, 0 otherwise.
##' @param y q x n matrix representing the log intensities of q peptides
## among n MS samples.
##' @param n1  number of samples under condition 1. It is assumed that the first n1 columns of y
## correspond to observations under condition 1. 
##' @param n2 number of samples under condition 2.
##' @return A list of the following elements:
## llr: log likelihood ratio statistic (maximum likelihood version).
## llr.map: log likelihood ratio statistic (maximum a posteriori version).
## llr.pv: p-value for llr.
## llr.map.pv: p-value for llr.map.
## mse.h0: Mean squared error under H0
## mse.h1: Mean squared error under H1
## s: selected regularization hyperparameter for llr.map.
## wchi2: weight used to make llr.map chi2-distributed under H0.
##' @author Thomas Burger, Laurent Jacob


pepa.test <- function(X, y, n1, n2){
    
    n <- n1+n2
    
    q <- nrow(X) # Number of peptides
    p <- ncol(X) # Number of proteins
    
    if(nrow(y) != q){
        stop('[pepa.test] y must have the same number of rows as X (one per peptide)')
    }
    
    y1 <- y[, 1:n1]
    y2 <- y[, -(1:n1)]
    
    ## Connected components: which peptides are connected to which proteins?
    print('[pepa.test] Identifying connected components in the peptide-protein graph')
    A <- matrix(0, nrow=p+q, ncol=p+q)
    A[1:p, (p+1):(p+q)] <- as.matrix(t(X))
    A[(p+1):(p+q), 1:p] <- as.matrix(X)
    g <- graphAM(A, edgemode='undirected', values=NA)
    nodes(g) <- as.character(1:(p+q))
    cc <- connComp(as(g, 'graphNEL'))
    
    ## Test proteins in each connected component
    protein.cc <- sapply(cc, FUN=function(ee){min(as.numeric(ee)) <= p})
    cc <- cc[protein.cc]
    
    ## Compute 'local' likelihoods, using only peptides belonging to
    ## protein of interest, or connected proteins.    
    llr <-  llr.scalable <- lh1 <- lh0 <- n.pep <- rep(NA, p)
    mse.h0 <- mse.h1 <- c()
    ii <- 0
    for(ee in cc){
        ii <- ii + 1
        print(sprintf('[pepa.test] Computing H0 likelihood for connected component %d/%d', ii, length(cc)))                
        ee <- as.numeric(ee)
        local.prot <- ee[ee <= p]
        if(length(local.prot) == 0){next()}
        local.pep <- ee[ee > p] - p
        local.X <- X[local.pep, local.prot, drop=FALSE]
        prot.barcode <- apply(local.X, 2, FUN=function(v) paste(as.character(v), collapse='-'))
        bc.dup <- duplicated(prot.barcode)
        equiv.X <- local.X[, !bc.dup, drop=FALSE]
        
        lik0 <- LH0(equiv.X, y1[local.pep, , drop=FALSE], y2[local.pep, , drop=FALSE])
        mse.h0[ii] <- lik0$ss/(length(local.pep)*n)
        
        for(jj in local.prot){
            print(sprintf('[pepa.test] Computing H1 likelihood for protein %d', jj))
            prot.equiv.idx <- match(names(prot.barcode[!bc.dup])[prot.barcode[!bc.dup] == prot.barcode[colnames(X)[jj]]],
                                    colnames(equiv.X))
            lik1 <- LH1(equiv.X, y1[local.pep, , drop=FALSE], y2[local.pep, , drop=FALSE], prot.equiv.idx)            
            mse.h1[jj] <- lik1$ss/(length(local.pep)*n)    
            llr[jj] <- (length(local.pep)*n) * (log(lik0$ss) - log(lik1$ss))
        }
    }
    s <- mse.h1   
    if(p < 500){# Following SAM paper recommendations, we don't attempt to estimate s0 if p<500
        s1 <- quantile(s, 0.05)
    }else{
        s1 <- fudge2LRT(mse.h0, mse.h1, cc, n, p, s)$s.zero
    }
    samLRT.res <- samLRT(mse.h0, mse.h1, cc, n, p, s1)
    llr.map <- samLRT.res$llr.sam
    
    ## Compute p-values
    llr.pv <- 1 - pchisq(llr, 1)
    ## Reweight LLR-sam for calibration
    s2.est <- mean(s)
    wchi2 <- rep((s2.est + s1)/s2.est, length(de))
    llr.map.pv <- 1 - pchisq(wchi2*llr.map, 1)
    
    return(list(llr=llr, llr.map=llr.map,
                llr.pv=llr.pv, llr.map.pv=llr.map.pv,
                mse.h0=mse.h0, mse.h1=mse.h1,
                s=s, s1=s1, wchi2=wchi2))
}
