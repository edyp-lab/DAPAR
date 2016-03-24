##' Provides several methods to normalize quantitative data from
##' a \code{\link{MSnSet}} object.
##' They are organized in four main families : Strong Rescaling, 
##' Median Centering, Mean Centering, Mean CenteringScaling.
##' For the first family, two sub-categories are available : the sum by columns
##' and the quantiles method.
##' For the three other families, two categories are available : 
##' "Overall" which means that the value for each protein 
##' (ie line in the expression data tab) is computed over all the samples ;
##' "within conditions" which means that the value for each protein 
##' (ie line in the \code{exprs()} data tab) is computed condition
##' by condition.
##' 
##' @title Normalisation
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param family One of the following : Global Rescaling, 
##' Median Centering, Mean Centering, Mean Centering Scaling.
##' @param method "Overall" or "within conditions".
##' @return An instance of class \code{\link{MSnSet}} 
##' where the quantitative data in the \code{exprs()} tab
##' has been normalized.
##' @author Alexia Dorffer
##' @examples
##' data(UPSpep25)
##' wrapper.normalizeD(UPSpep25, "Median Centering", "within conditions")
wrapper.normalizeD <- function(obj, family, method){

qData <- Biobase::exprs(obj)
labels <- Biobase::pData(obj)[,"Label"]
Biobase::exprs(obj) <- normalizeD(qData, labels, family, method)
msg <- paste("Normalisation using family ", family,  sep="")
msg2 <- paste("With method ", method,  sep="")
obj@processingData@processing <- c(obj@processingData@processing,
                                    msg, msg2)
return(obj)
}


##' Provides several methods to normalize data from
##' a matrix.
##' They are organized in four main families : Strong Rescaling, 
##' Median Centering, Mean Centering, Mean CenteringScaling.
##' For the first family, two sub-categories are available : the sum by columns
##' and the quantiles method.
##' For the three other families, two categories are available : 
##' "Overall" which means that the value for each protein 
##' (ie line in the expression data tab) is computed over all the samples ;
##' "within conditions" which means that the value for each protein 
##' (ie line in the matrix) is computed condition
##' by condition.
##' 
##' @title Normalisation
##' @param qData A dataframe that contains quantitative data.
##' @param labels A vector of strings containing the column "Label" of 
##' the \code{pData()}.
##' @param family One of the following : Global Rescaling, 
##' Median Centering, Mean Centering, Mean Centering Scaling.
##' @param method "Overall" or "within conditions".
##' @return A matrix normalized
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' qData <- exprs(UPSpep25)
##' labels <- pData(UPSpep25)[,"Label"]
##' normalizeD(qData, labels, "Median Centering", "within conditions")
normalizeD <- function(qData, labels, family, method){
#Verification des parametres
paramfamily<-c("Global Rescaling", "Median Centering", "Mean Centering", 
                "Mean Centering Scaling")
if (sum(is.na(match(family, paramfamily)==TRUE))>0){
    warning("Parameter family is not correct")
    return (NULL)
}

parammethod<-c("sum by columns", "quantiles", "overall", 
                "within conditions")
if (sum(is.na(match(method, parammethod)==TRUE))>0){
    warning("Parameter method is not correct")
    return (NULL)
}


.temp <- qData
if (!is.null(.temp)){
    data <- .temp
    if (family == "Global Rescaling"){
    if (method == "sum by columns"){
        
        t <- 2^(.temp)
        s <- colSums(t, na.rm=TRUE)
        data <- t
        for ( i in 1:nrow(data)) { data[i,] <- t[i,] / s}
        
        .temp <- log2(data)
    }
    else if (method == "quantiles"){
        .temp <- normalize.quantiles(.temp)
        dimnames(.temp) <- list(rownames(qData),
                                colnames(qData))
    }
    }else if (family =="Median Centering"){
    medianOverSamples <- apply(.temp, 2, median, na.rm = TRUE)
    
    if (method == "overall"){
        cOverall <- median(medianOverSamples)
        .temp <- sweep(.temp, 2, medianOverSamples)
        .temp <- .temp + cOverall
        
    }
    else if (method == "within conditions"){
        .temp <- sweep(.temp, 2, medianOverSamples)
        
        cCond <- NULL
        for (l in unique(labels))
        {
        indices <- which(labels== l)
        cCond[l] <- median(medianOverSamples[indices])
        .temp[,indices] <- .temp[,indices] + cCond[l]
        }
        
    }
    } else if (family =="Mean Centering"){
    meanOverSamples <- apply(.temp, 2, mean, na.rm = TRUE)
    
    if (method == "overall"){
        cOverall <- mean(meanOverSamples)
        .temp <- sweep(.temp, 2, meanOverSamples)
        .temp <- .temp + cOverall
        
    } else if (method == "within conditions"){
        .temp <- sweep(.temp, 2, meanOverSamples)
        cCond <- NULL
        for (l in unique(labels))
        {
        indices <- which(labels== l)
        cCond[l] <- mean(meanOverSamples[indices])
        .temp[,indices] <- .temp[,indices] + cCond[l]
        }
    }
    } else if (family == "Mean Centering Scaling"){
    meanOverSamples <- apply(.temp, 2, mean, na.rm = TRUE)

    if (method == "overall"){
        cOverall <- mean(meanOverSamples)
        .temp <- sweep(.temp, 2, meanOverSamples)
        .temp <- scale(.temp,center=FALSE,
                        scale=TRUE)
        attr(.temp,"scaled:scale")<-NULL
        .temp <- .temp + cOverall
        
        
    } else if (method == "within conditions"){
        .temp <- sweep(.temp, 2, meanOverSamples)
        .temp <- scale(.temp,center=FALSE,
                        scale=TRUE)
        attr(.temp,"scaled:scale")<-NULL
        
        cCond <- NULL
        for (l in unique(labels))
        {
        indices <- which(labels== l)
        cCond[l] <- mean(meanOverSamples[indices])
        .temp[,indices] <- .temp[,indices] + cCond[l]
        }
        
    }
    }
    
    #     msg <- paste("Normalisation using ", method,  sep="")
    #     .temp@processingData@processing <- c(.temp@processingData@processing,
    #                                          msg)
}
return(.temp)
}


