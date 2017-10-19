##' This method is a wrapper to the \code{imputeLCMD} package adapted to 
##' objects of class \code{MSnSet}.
##' 
##' @title Missing values imputation from a \code{MSnSet} object
##' @param obj An object of class \code{MSnSet}.
##' @param method The imputation method to be used.
##' Choices are QRILC, KNN, BPCA and MLE. 
##' @return The object \code{obj} which has been imputed
##' @author Alexia Dorffer
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.mvImputation(Exp1_R25_pept[1:1000], "QRILC")
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)[1:1000]
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



##' This method is xxxxxxxx
##' 
##' @title xxxxxxxxx
##' @param obj xxxxxxxxx
##' @param qval xxxxxxxxx
##' @param factor xxxxxxxxx. 
##' @return xxxxxxxxx
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.impute.detQuant(Exp1_R25_pept)
wrapper.impute.detQuant <- function(obj, qval=0.025, factor=1){
  if (is.null(obj)){return(NULL)}
  
  qData <- Biobase::exprs(obj)
  values <- getQuantile4Imp(qData, qval, factor)
  
  Biobase::exprs(obj) <- impute.detQuant(qData, values$shiftedImpVal)
  msg <- "Missing values imputation using deterministic quantile"
  obj@processingData@processing <- c(obj@processingData@processing,msg)
  
  obj@experimentData@other$imputation.method <- "detQuantile"
  return(obj)
}



  
##' This method returns the q-th quantile of each colum of an expression set, up to a scaling factor
##'
##' @title Quantile imputation value definition
##' @param qData An expression set containing quantitative values of various replicates
##' @param qval The quantile used to define the imputation value
##' @param factor A scaling factor to multiply the imputation value with
##' @return A list of two vectors, respectively containing the imputation values and the rescaled imputation values
##' @author Thomas Burger
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' getQuantile4Imp(qData) 
getQuantile4Imp <- function(qData, qval=0.025, factor=1){
  r1 <- apply(qData, 2, quantile, qval, na.rm=TRUE)
  r2 <- r1*factor
  return(list(ImpVal = r1, shiftedImpVal = r2))
}


##' This method replaces each missing value by a given value
##'
##' @title Deterministic imputation
##' @param qData An expression set containing quantitative or missing values
##' @param values A vector with as many elements as the number of colums of qData
##' @return An imputed dataset
##' @author Thomas Burger
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' values <- getQuantile4Imp(qData)$shiftedImpVal
##' impute.detQuant(qData, values) 
impute.detQuant <- function(qData, values){
  for(i in 1:dim(qData)[2]){
    col <- qData[,i]
    col[which(is.na(col))] <- values[i]
    qData[,i] <- col
  }
  return(qData)
}









##' This method is a wrapper to the function \code{impute.pa} of the package
##' \code{imp4p} adapted to an object of class \code{MSnSet}.
##'
##' @title Imputation of peptides having no values in a biological condition.
##' @param obj An object of class \code{MSnSet}.
##' @param q.min Same as the function \code{impute.pa} in the package \code{imp4p}
##' @return The \code{exprs(obj)} matrix with imputed values instead of missing values.
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' dat <- mvFilter(Exp1_R25_pept[1:1000], type="allCond", th = 1)
##' dat <- wrapper.impute.pa(dat)
wrapper.impute.pa <- function(obj, q.min = 0.025){
    cond <- as.factor(Biobase::pData(obj)$Label)
    Biobase::exprs(obj) <- impute.pa(Biobase::exprs(obj), conditions=cond, q.min = q.min, q.norm=3,  eps=0)

    return (obj)
}


##' This method is a wrapper to the function \code{impute.mi} of the package \code{imp4p} adapted to
##' an object of class \code{MSnSet}.
##'
##' @title Missing values imputation using the LSimpute algorithm.
##' @param obj An object of class \code{MSnSet}.
##' @param nb.iter Same as the function \code{mi.mix} in the package \code{imp4p}
##' @param nknn Same as the function \code{mi.mix} in the package \code{imp4p}
##' @param selec Same as the function \code{mi.mix} in the package \code{imp4p}
##' @param siz Same as the function \code{mi.mix} in the package \code{imp4p}
##' @param weight Same as the function \code{mi.mix} in the package \code{imp4p}
##' @param ind.comp Same as the function \code{mi.mix} in the package \code{imp4p}
##' @param progress.bar Same as the function \code{mi.mix} in the package \code{imp4p}
##' @param x.step.mod Same as the function \code{estim.mix} in the package \code{imp4p}
##' @param x.step.pi Same as the function \code{estim.mix} in the package \code{imp4p}
##' @param nb.rei Same as the function \code{estim.mix} in the package \code{imp4p}
##' @param method Same as the function \code{estim.mix} in the package \code{imp4p}
##' @param gridsize Same as the function \code{estim.mix} in the package \code{imp4p}
##' @param q Same as the function \code{mi.mix} in the package \code{imp4p}
##' @param q.min Same as the function \code{impute.pa} in the package \code{imp4p}
##' @param q.norm Same as the function \code{impute.pa} in the package \code{imp4p}
##' @param eps Same as the function \code{impute.pa} in the package \code{imp4p}
##' @param methodi Same as the function \code{mi.mix} in the package \code{imp4p}
##' @param lapala xxxxxxxxxxx
##' @param distribution xxxxxxxxxxx
##' @return The \code{exprs(obj)} matrix with imputed values instead of missing values.
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' dat <- mvFilter(Exp1_R25_pept[1:1000], type="allCond", th = 1)
##' dat <- wrapper.dapar.impute.mi(dat, nb.iter=1)
wrapper.dapar.impute.mi <- function (obj, nb.iter = 3, 
                               nknn = 15, selec = 600, siz = 500, weight = 1, ind.comp = 1, 
                               progress.bar = TRUE, x.step.mod = 300, 
                               x.step.pi = 300, nb.rei = 100, method = 4, gridsize = 300, 
                               q = 0.95, q.min = 0, q.norm = 3, eps = 0, methodi = "slsa",
                               lapala = TRUE,
                               distribution="unif") 
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
                    x.step.mod = x.step.mod, 
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
        data.final = impute.pa2(tab = data.mi, conditions = conditions, 
                           q.min = q.min, q.norm = q.norm, eps = eps, distribution = distribution)
        } else {
        data.final <- data.mi
    }

    colnames(data.final) <- colnames(Biobase::exprs(obj))
    Biobase::exprs(obj) <- data.final
    
    msg <- paste("Missing values imputation using imp4p")
    obj@processingData@processing <- c(obj@processingData@processing,msg)
    
    obj@experimentData@other$imputation.method <- "imp4p"
    
    return(obj)
    
}


################################################
##' This method xxxxxxxxxxx
##' 
##' @title xxxxxxx
##' @param n An integer which is the number of simulation (same as in rbeta)
##' @param min An integer that corresponds to the lower bound of the interval
##' @param max An integer that corresponds to the upper bound of the interval
##' @param param1 An integer that is the first parameter of rbeta function.
##' @param param2 An integer that is second parameter of rbeta function.
##' @return A vector of n simulated values
##' @author Thomas Burger
##' @examples
##' translatedRandomBeta(1000, 5, 10, 1, 1)
translatedRandomBeta <- function(n, min, max, param1=3, param2=1){
    scale <- max-min
    simu <- rbeta(n,param1,param2)
    res <- (simu*scale) + min
    return(res)
}


################################################
##' This method is a wrapper to the function \code{impute.pa} from the package 
##' \code{imp4p} adapted to objects of class \code{MSnSet}.
##' 
##' @title Missing values imputation from a \code{MSnSet} object
##' @param obj An object of class \code{MSnSet}.
##' @param q.min A quantile value of the observed values allowing defining the 
##' maximal value which can be generated. This maximal value is defined by the
##' quantile q.min of the observed values distribution minus eps. 
##' Default is 0 (the maximal value is the minimum of observed values minus eps).
##' @param q.norm A quantile value of a normal distribution allowing defining 
##' the minimal value which can be generated. Default is 3 (the minimal value 
##' is the maximal value minus qn*median(sd(observed values)) where sd is the 
##' standard deviation of a row in a condition).
##' @param eps A value allowing defining the maximal value which can be 
##' generated. This maximal value is defined by the quantile q.min of the 
##' observed values distribution minus eps. Default is 0.
##' @param distribution The type of distribution used. Values are unif or beta.
##' @return The object \code{obj} which has been imputed
##' @author Thomas Burger, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.impute.pa2(Exp1_R25_pept[1:1000], distribution="beta")
wrapper.impute.pa2 <- function (obj, q.min = 0, q.norm = 3, eps = 0, distribution = "unif"){
    tab <- Biobase::exprs(obj)
    conditions <- as.factor(Biobase::pData(obj)$Label)

    tab_imp <- impute.pa2(tab, conditions, q.min, q.norm, eps, distribution)
    Biobase::exprs(obj) <- tab_imp
    
     
     return(obj)
    }





################################################
##' This method is a variation to the function \code{impute.pa} from the package 
##' \code{imp4p}.
##' 
##' @title Missing values imputation from a \code{MSnSet} object
##' @param tab An object of class \code{MSnSet}.
##' @param conditions xxxxxxxxxxxx
##' @param q.min A quantile value of the observed values allowing defining the 
##' maximal value which can be generated. This maximal value is defined by the
##' quantile q.min of the observed values distribution minus eps. 
##' Default is 0 (the maximal value is the minimum of observed values minus eps).
##' @param q.norm A quantile value of a normal distribution allowing defining 
##' the minimal value which can be generated. Default is 3 (the minimal value 
##' is the maximal value minus qn*median(sd(observed values)) where sd is the 
##' standard deviation of a row in a condition).
##' @param eps A value allowing defining the maximal value which can be 
##' generated. This maximal value is defined by the quantile q.min of the 
##' observed values distribution minus eps. Default is 0.
##' @param distribution The type of distribution used. Values are unif or beta.
##' @return The object \code{obj} which has been imputed
##' @author Thomas Burger, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.impute.pa2(Exp1_R25_pept[1:1000], distribution="beta")
impute.pa2 <- function (tab, conditions, q.min = 0, q.norm = 3, eps = 0, distribution = "unif"){
    tab_imp = tab
    qu = apply(tab_imp, 2, quantile, na.rm = TRUE, q.min)
    nb_cond = length(levels(conditions))
    nb_rep = rep(0, nb_cond)
    k = 1
    j = 1
    for (i in 1:nb_cond) {
        nb_rep[i] = sum((conditions == levels(conditions)[i]))
        sde = apply(tab_imp[, (k:(k + nb_rep[i] - 1))], 1, sd, 
                    na.rm = TRUE)
        while (j < (k + nb_rep[i])) {
            if(distribution == "unif")
            {
                tab_imp[which(is.na(tab_imp[, j])), j] = runif(n = sum(is.na(tab_imp[,j])), 
                                                               min = qu[j] - eps - q.norm * median(sde, na.rm = TRUE), 
                                                               max = qu[j] - eps)
            } else if (distribution == "beta"){
                tab_imp[which(is.na(tab_imp[, j])), j] = translatedRandomBeta(n = sum(is.na(tab_imp[,j])), 
                                                                              min = qu[j] - eps - q.norm * median(sde, na.rm = TRUE), 
                                                                              max = qu[j] - eps,
                                                                              param1 = 3,
                                                                              param2 = 1)
            }
            j = j + 1
        }
        k = k + nb_rep[i]
    }
    
    return(tab_imp)
}










