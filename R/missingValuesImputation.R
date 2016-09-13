##' This method is a wrapper to the \code{imputeLCMD} package adapted to 
##' objects of class \code{\link{MSnSet}}.
##' 
##' @title Missing values imputation from a \code{\link{MSnSet}} object
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param method The imputation method to be used.
##' Choices are QRILC, KNN, BPCA and MLE. 
##' @return The object \code{obj} which has been imputed
##' @author Alexia Dorffer
##' @examples
##' data(UPSpep25)
##' wrapper.mvImputation(UPSpep25, "QRILC")
wrapper.mvImputation <- function(obj, method){
qData <- Biobase::exprs(obj)
Biobase::exprs(obj) <- mvImputation(qData, method)
msg <- paste("Missing values imputation using ", method,  sep="")
obj@processingData@processing <- c(obj@processingData@processing,msg)

obj@experimentData@other$imputation.method <- method
return(obj)
}




##' This method is a wrapper to the \code{imputeLCMD} package adapted to 
##' a matrix.
##' 
##' @title Missing values imputation from a matrix
##' @param qData A dataframe that contains quantitative data.
##' @param method The imputation method to be used.
##' Choices are QRILC, KNN, BPCA and MLE. 
##' @return The matrix imputed
##' @author Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' qData <- Biobase::exprs(UPSpep25)
##' mvImputation(qData, "QRILC")
mvImputation <- function(qData, method){
#Check parameters
param<-c("BPCA", "KNN", "MLE", "QRILC")
if (sum(is.na(match(method, param)==TRUE))>0){
    warning("Param method is not correct.")
    return (NULL)
}

qData <- as.matrix(qData)
## BPCA impute imputation at peptide level

if (method == "BPCA"){
    if (getNumberOfEmptyLines(qData) > 0) {
        stop("Data contains rows in which 
                all elements are 'NA'. Remove them first")}
    ## replace all 0's by NAs
    tmp.data <- qData
    nSamples <- dim(qData)[2]
    resultBPCA <- pca(qData, method="bpca", nPcs=(nSamples-1))
    exprs <- completeObs(resultBPCA)
    #     msg <- "Missing values imputation using bpca algorithm"
    #     obj@processingData@processing <- c(obj@processingData@processing,msg)
} else if (method == "KNN"){
    if (getNumberOfEmptyLines(qData) > 0) {
        stop("Data contains rows in which all elements are 'NA'. 
            Remove them first")}
    exprs <- impute.wrapper.KNN(qData, 3)
    #     msg <- "Missing values imputation using KNN algorithm"
    #     obj@processingData@processing <- c(obj@processingData@processing,msg)
} else if (method == "MLE"){
    exprs <- impute.wrapper.MLE(qData)
    #     msg <- "Missing values imputation using MLE algorithm"
    #     obj@processingData@processing <- c(obj@processingData@processing,msg)
} else if (method == "QRILC"){
    exprs <- impute.QRILC(qData)[[1]]
    #     msg <- "Missing values imputation using QRILC algorithm"
    #     obj@processingData@processing <- c(obj@processingData@processing,msg)
}


return (exprs)
}




# Functions to use with the imp4p package

wrapper.identifyMCAR_MNAR <- function(obj, dat.slsa){
    #Estimation of the mixture model
    res <- estim.mix(tab=Biobase::exprs(obj), tab.imp=dat.slsa, conditions=as.factor(Biobase::pData(obj)$Label))
    #Computing probabilities to be MCAR
    born <- estim.bound(tab=Biobase::exprs(obj),conditions=as.factor(Biobase::pData(obj)$Label))
    proba <- prob.mcar.tab(born$tab.lower,born$tab.upper,res)
    return(proba)
}



wrapper.imputeImp4p <- function(obj, dat.slsa, proba, nb.iter=1, selec=100,ind.comp=0){
    
    #Multiple imputation strategy with 5 iterations (can be time consuming!)
    data.mi <- mi.mix(tab = Biobase::exprs(obj), 
                      tab.imp = dat.slsa, 
                      prob.MCAR = proba, 
                      conditions = as.factor(Biobase::pData(obj)$Label), 
                      repbio = as.factor(Biobase::pData(obj)$Bio.Rep), 
                      reptech = as.factor(Biobase::pData(obj)$Tech.Rep),
                      nb.iter = nb.iter,
                      selec = selec,
                      ind.comp = ind.comp)
    colnames(data.mi) <- colnames(Biobase::exprs(obj))
    Biobase::exprs(obj) <- data.mi
    
    msg <- paste("Missing values imputation using imp4p")
    obj@processingData@processing <- c(obj@processingData@processing,msg)
    
    obj@experimentData@other$imputation.method <- "imp4p"
    
    return(obj)
}


wrapper.impute.slsa <- function(obj, selec = 100){
    #Imputation of missing values with the slsa algorithm
    dat.slsa <- impute.slsa(Biobase::exprs(obj),
                            conditions=as.factor(Biobase::pData(obj)$Label),
                            repbio=as.factor(Biobase::pData(obj)$Bio.Rep),
                            reptech=as.factor(Biobase::pData(obj)$Tech.Rep),
                            selec = selec)
    return(dat.slsa)
}


