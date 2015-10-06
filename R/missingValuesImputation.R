##' This method is a wrapper to the \code{imputeLCMD} package adapted to 
##' objects of class \code{\link{MSnSet}}.
##' 
##' @title Missing values imputation
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param method The imputation method to be used.
##' Choices are QRILC, KNN, BPCA and MLE. 
##' @return The object \code{obj} which has been imputed
##' @author Samuel Wieczorek
##' @examples data(UPSprotx2)
##' mvImputation(UPSprotx2, "QRILC")
mvImputation <- function(obj, method){
    #Check parameters
    param<-c("BPCA", "KNN", "MLE", "QRILC")
    if (sum(is.na(match(method, param)==TRUE))>0){
        warning("Param method is not correct.")
        return (NULL)
    }

    ## BPCA impute imputation at peptide level

    if (method == "BPCA"){
        ## replace all 0's by NAs
        tmp.data <- exprs(obj)
        nSamples <- dim(exprs(obj))[2]
        resultBPCA <- pca(exprs(obj), method="bpca", nPcs=(nSamples-1))
        exprs(obj) <- completeObs(resultBPCA)
        msg <- "Missing values imputation using bpca algorithm"
        obj@processingData@processing <- c(obj@processingData@processing,msg)
    } else if (method == "KNN"){
        exprs(obj) <- impute.wrapper.KNN(exprs(obj), 3)
        msg <- "Missing values imputation using KNN algorithm"
    obj@processingData@processing <- c(obj@processingData@processing,msg)
    } else if (method == "MLE"){
        exprs(obj) <- impute.wrapper.MLE(exprs(obj))
        msg <- "Missing values imputation using MLE algorithm"
        obj@processingData@processing <- c(obj@processingData@processing,msg)
    } else if (method == "QRILC"){
        exprs(obj) <- impute.QRILC(exprs(obj))[[1]]
        msg <- "Missing values imputation using QRILC algorithm"
        obj@processingData@processing <- c(obj@processingData@processing,msg)
    }
    return (obj)
}


