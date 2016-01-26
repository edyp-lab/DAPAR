##' This method is a wrapper to the \code{imputeLCMD} package adapted to 
##' objects of class \code{\link{MSnSet}}.
##' 
##' @title Missing values imputation from a \code{\link{MSnSet}} object
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param method The imputation method to be used.
##' Choices are QRILC, KNN, BPCA and MLE. 
##' @return The object \code{obj} which has been imputed
##' @author Alexia Dorffer
##' @examples data(UPSprotx2)
##' wrapper.mvImputation(UPSprotx2, "QRILC")
wrapper.mvImputation <- function(obj, method){
  
  qData <- exprs(obj)
  exprs(obj) <- mvImputation(qData, method)
  msg <- paste("Missing values imputation using ", method,  sep="")
  obj@processingData@processing <- c(obj@processingData@processing,msg)
  return(obj)
  
}




##' This method is a wrapper to the \code{imputeLCMD} package adapted to a matrix.
##' 
##' @title Missing values imputation from a matrix
##' @param qData A dataframe that contains quantitative data.
##' @param method The imputation method to be used.
##' Choices are QRILC, KNN, BPCA and MLE. 
##' @return The matrix imputed
##' @author Samuel Wieczorek
##' @examples data(UPSprotx2)
##' qData <- exprs(UPSprotx2)
##' mvImputation(qData, "QRILC")
mvImputation <- function(qData, method){
  #Check parameters
  param<-c("BPCA", "KNN", "MLE", "QRILC")
  if (sum(is.na(match(method, param)==TRUE))>0){
    warning("Param method is not correct.")
    return (NULL)
  }
  
  ## BPCA impute imputation at peptide level
  
  if (method == "BPCA"){
    ## replace all 0's by NAs
    tmp.data <- qData
    nSamples <- dim(qData)[2]
    resultBPCA <- pca(qData, method="bpca", nPcs=(nSamples-1))
    exprs <- completeObs(resultBPCA)
    #     msg <- "Missing values imputation using bpca algorithm"
    #     obj@processingData@processing <- c(obj@processingData@processing,msg)
  } else if (method == "KNN"){
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
