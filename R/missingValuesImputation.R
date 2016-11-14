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




##' This method is a wrapper to the function \code{impute.pa} of the package
##' \code{imp4p} adapted to an object of class \code{\link{MSnSet}}.
##'
##' @title Imputation of peptides having no values in a biological condition.
##' @param obj An object of class \code{\link{MSnSet}}.
##' @return The \code{exprs(obj)} matrix with imputed values instead of missing values.
##' @author Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' dat <- mvFilter(UPSpep25, type="allCond", th = 1)
##' dat <- wrapper.impute.pa(dat)
wrapper.impute.pa <- function(obj){
    cond <- as.factor(Biobase::pData(obj)$Label)
    Biobase::exprs(obj) <- impute.pa(Biobase::exprs(obj), conditions=cond)

    return (obj)
}


##' This method is a wrapper to the function \code{impute.mi} of the package \code{imp4p} adapted to
##' an object of class \code{\link{MSnSet}}.
##'
##' @title Missing values imputation using the LSimpute algorithm.
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param selec An integer which corresponds to the \code{selec} parameter of the
##' function \code{impute.slsa}.
##' @return The \code{exprs(obj)} matrix with imputed values instead of missing values.
##' @author Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' dat <- mvFilter(UPSpep25, type="allCond", th = 1)
##' dat <- wrapper.dapar.impute.mi(dat)
wrapper.dapar.impute.mi <- function (obj, nb.iter = 3, 
                               nknn = 15, selec = 600, siz = 500, weight = 1, ind.comp = 1, 
                               progress.bar = TRUE, x.min = 10, x.max = 35, x.step.mod = 300, 
                               x.step.pi = 300, nb.rei = 100, method = 4, gridsize = 300, 
                               q = 0.95, q.min = 0, q.norm = 3, eps = 2, methodi = "slsa",
                               lapala = TRUE) 
{
    
    conditions <- as.factor(Biobase::pData(obj)$Label)
    repbio <- as.factor(Biobase::pData(obj)$Bio.Rep)
    reptech <-as.factor(Biobase::pData(obj)$Tech.Rep)
    
    tab <- Biobase::exprs(obj)

    if (progress.bar == TRUE) {
        cat(paste("\n 1/ Initial imputation under the MCAR assumption with impute.rand ... \n  "))
    }
    dat.slsa = impute.rand(tab = tab, conditions = conditions)
    if (progress.bar == TRUE) {
        cat(paste("\n 2/ Estimation of the mixture model in each sample... \n  "))
    }
    res = estim.mix(tab = tab, tab.imp = dat.slsa, conditions = conditions, 
                    x.min = x.min, x.max = x.max, x.step.mod = x.step.mod, 
                    x.step.pi = x.step.pi, nb.rei = nb.rei, method = method, 
                    gridsize = gridsize)
    if (progress.bar == TRUE) {
        cat(paste("\n 3/ Estimation of the probabilities each missing value is MCAR... \n  "))
    }
    born = estim.bound(tab = tab, conditions = conditions, q = q)
    proba = prob.mcar.tab(born$tab.lower, born$tab.upper, res)
    if (progress.bar == TRUE) {
        cat(paste("\n 4/ Multiple imputation strategy with mi.mix ... \n  "))
    }
    data.mi = mi.mix(tab = tab, tab.imp = dat.slsa, prob.MCAR = proba, 
                     conditions = conditions, repbio = repbio, reptech = reptech, 
                     nb.iter = nb.iter, nknn = nknn, weight = weight, selec = selec, 
                     siz = siz, ind.comp = ind.comp, methodi = methodi, q = q, 
                     progress.bar = progress.bar)
    
    if (lapala == TRUE){
        if (progress.bar == TRUE) {
            cat(paste("\n\n 5/ Imputation of rows with only missing values in a condition with impute.pa ... \n  "))
        }
        data.final = impute.pa(tab = data.mi, conditions = conditions, 
                           q.min = q.min, q.norm = q.norm, eps = eps)
    }
    else {
        data.final <- data.mi
    }

    colnames(data.final) <- colnames(Biobase::exprs(obj))
    Biobase::exprs(obj) <- data.final
    
    msg <- paste("Missing values imputation using imp4p")
    obj@processingData@processing <- c(obj@processingData@processing,msg)
    
    obj@experimentData@other$imputation.method <- "imp4p"
    
    return(obj)
    
}
