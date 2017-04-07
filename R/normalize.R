##' Provides several methods to normalize quantitative data from
##' a \code{\link{MSnSet}} object.
##' They are organized in four main families : Global Alignment, 
##' Median Centering, Mean Centering, Mean Centering Scaling.
##' For the first family, two sub-categories are available : the sum by columns
##' and the quantiles method.
##' For the three other families, two categories are available : 
##' "Overall" which means that the value for each protein 
##' (ie line in the expression data tab) is computed over all the samples ;
##' "within conditions" which means that the value for each protein 
##' (ie line in the \code{exprs()} data tab) is computed condition
##' by condition.
##' 
##' @title Normalization
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param family One of the following : Global Alignment, 
##' Median Centering, Mean Centering, Mean Centering Scaling.
##' @param method "overall" or "within conditions".
##' @return An instance of class \code{\link{MSnSet}} 
##' where the quantitative data in the \code{exprs()} tab
##' has been normalized.
##' @author Alexia Dorffer
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.normalizeD(Exp1_R25_pept[1:1000], "Median Centering", "within conditions")
wrapper.normalizeD <- function(obj, family, method){

qData <- Biobase::exprs(obj)
labels <- Biobase::pData(obj)[,"Label"]
Biobase::exprs(obj) <- normalizeD(qData, labels, family, method)
msg <- paste("Normalisation using family ", family,  sep="")
msg2 <- paste("With method ", method,  sep="")
obj@processingData@processing <- c(obj@processingData@processing,
                                    msg, msg2)

obj@experimentData@other$normalizationFamily <- family
obj@experimentData@other$normalizationMethod <- method


return(obj)
}


##' Provides several methods to normalize data from
##' a matrix.
##' They are organized in four main families : Strong Rescaling, 
##' Median Centering, Mean Centering, Mean CenteringScaling.
##' For the first family, two sub-categories are available : the sum by columns
##' and the quantiles method.
##' For the three other families, two categories are available : 
##' "overall" which means that the value for each protein 
##' (ie line in the expression data tab) is computed over all the samples ;
##' "within conditions" which means that the value for each protein 
##' (ie line in the matrix) is computed condition
##' by condition.
##' 
##' @title Normalisation
##' @param qData A dataframe that contains quantitative data.
##' @param labels A vector of strings containing the column "Label" of 
##' the \code{pData()}.
##' @param family One of the following : Global Alignment, 
##' Median Centering, Mean Centering, Mean Centering Scaling.
##' @param method "Overall" or "within conditions".
##' @return A matrix normalized
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept[1:1000])
##' labels <- Biobase::pData(Exp1_R25_pept[1:1000])[,"Label"]
##' normalizeD(qData, labels, "Median Centering", "within conditions")
normalizeD <- function(qData, labels, family, method){
#Verification des parametres
paramfamily<-c("Global Alignment", "Median Centering", "Mean Centering", 
                "Mean Centering Scaling")
if (sum(is.na(match(family, paramfamily)==TRUE))>0){
    warning("Parameter family is not correct")
    return (NULL)
}

parammethod<-c("sum by columns", "quantile alignment", "overall", 
                "within conditions")
if (sum(is.na(match(method, parammethod)==TRUE))>0){
    warning("Parameter method is not correct")
    return (NULL)
}


.temp <- qData
if (!is.null(.temp)){
    data <- .temp
    
    
    
    ###############
    if (family == "Global Alignment"){
        if (method == "sum by columns"){
            t <- 2^(.temp)
            s <- colSums(t, na.rm=TRUE)
            data <- t
            for ( i in 1:nrow(data)) { data[i,] <- t[i,] / s}
            .temp <- log2(data)
        }
    else if (method == "quantile alignment"){
        .temp <- normalize.quantiles(.temp)
        dimnames(.temp) <- list(rownames(qData),colnames(qData))
        }
    }
    
    
    
    ####################
    else if (family =="Median Centering"){
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
    }
    
    
    
    ######################
    else if (family =="Mean Centering"){
        meanOverSamples <- apply(.temp, 2, mean, na.rm = TRUE)
    
        if (method == "overall"){
            cOverall <- mean(meanOverSamples)
            .temp <- sweep(.temp, 2, meanOverSamples)
            .temp <- .temp + cOverall
            } 
        else if (method == "within conditions"){
            .temp <- sweep(.temp, 2, meanOverSamples)
            cCond <- NULL
            for (l in unique(labels))
                {
                indices <- which(labels== l)
                cCond[l] <- mean(meanOverSamples[indices])
                .temp[,indices] <- .temp[,indices] + cCond[l]
                }
         }
    } 
    
    
    ###############
    else if (family == "Mean Centering Scaling"){
        meanOverSamples <- apply(.temp, 2, mean, na.rm = TRUE)

        if (method == "overall"){
            cOverall <- mean(meanOverSamples)
            .temp <- sweep(.temp, 2, meanOverSamples)
            .temp <- scale(.temp,center=FALSE,scale=TRUE)
            attr(.temp,"scaled:scale")<-NULL
            .temp <- .temp + cOverall
            } 
        else if (method == "within conditions"){
            .temp <- sweep(.temp, 2, meanOverSamples)
            .temp <- scale(.temp,center=FALSE, scale=TRUE)
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
##' @param method One of the following : Global Alignment (for normalizations 
##' of important magnitude), Quantile Centering, Mean Centering.  
##' @param type For the method "Global Alignment", the parameters are:
##' "sum by columns": operates on the original scale (not the log2 one) and propose 
##' to normalize each abundance by the total abundance of the sample (so as to focus
##' on the analyte proportions among each sample).
##' "Alignment on all quantiles": proposes to align the quantiles of all the 
##' replicates; practically it amounts to replace abundances by order statistics.
##' For the two other methods, the parameters are "overall" (shift all the 
##' sample distributions at once) or "within conditions" (shift the sample 
##' distributions within each condition at a time).
##' @param scaling A boolean that indicates if the variance of the data have to 
##' be forced to unit (variance reduction) or not.
##' @param quantile A float that corresponds to the quantile used to align the
##' data.
##' @return An instance of class \code{\link{MSnSet}} where the quantitative 
##' data in the \code{exprs()} tab has been normalized.
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.normalizeD2(Exp1_R25_pept[1:1000], "Quantile Centering", "within conditions")
wrapper.normalizeD2 <- function(obj, method, type, scaling=FALSE, quantile=0.15){
    
    qData <- Biobase::exprs(obj)
    labels <- Biobase::pData(obj)[,"Label"]
    Biobase::exprs(obj) <- normalizeD2(qData, labels, method, type, scaling, quantile)
    
    msg_method <- paste("Normalisation using method =", method,  sep="")
    msg_type <- paste("With type =", type,  sep="")
    
    if (method == "Global Alignment"){
        obj@processingData@processing <- c(obj@processingData@processing, msg_method, msg_type)
        obj@experimentData@other$normalizationMethod <- method
        obj@experimentData@other$normalizationType <- type
    }
    else if (method =="Quantile Centering"){
        msg_quantile <- paste("With quantile =", quantile,  sep="")
        obj@processingData@processing <- c(obj@processingData@processing, msg_method, msg_type, msg_quantile)
        obj@experimentData@other$normalizationMethod <- method
        obj@experimentData@other$normalizationType <- type
        obj@experimentData@other$normalizationQuantile <- quantile
    }   
    else if (method =="Mean Centering"){
        msg_scaling <- paste("With scaling =", scaling,  sep="")
        obj@processingData@processing <- c(obj@processingData@processing, msg_method, msg_type, msg_scaling)
        obj@experimentData@other$normalizationMethod <- msg_method
        obj@experimentData@other$normalizationType <- type
        obj@experimentData@other$normalizationScaling <- scaling
    } 

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
##' @param method One of the following : Global Alignment, 
##' Quantile Centering, Mean Centering.
##' @param type For the method "Global Alignment", the parameters are:
##' "sum by columns": operates on the original scale (not the log2 one) and propose 
##' to normalize each abundance by the total abundance of the sample (so as to focus
##' on the analyte proportions among each sample).
##' "Alignment on all quantiles": proposes to align the quantiles of all the 
##' replicates; practically it amounts to replace abundances by order statistics.
##' For the two other methods, the parameters are "overall" (shift all the 
##' sample distributions at once) or "within conditions" (shift the sample 
##' distributions within each condition at a time).
##' @param scaling A boolean that indicates if the variance of the data have to 
##' be forced to unit (variance reduction) or not.
##' @param quantile A float that corresponds to the quantile used to align the
##' data.
##' @return A matrix normalized
##' @author Samuel Wieczorek, Thomas Burger
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept[1:1000])
##' labels <- Biobase::pData(Exp1_R25_pept[1:1000])[,"Label"]
##' normalizeD2(qData, labels, "Quantile Centering", "within conditions", quantile = 0.15)
normalizeD2 <- function(qData, labels, method, type, scaling=FALSE, quantile=0.15){
    #Check of parameters
    parammethod<-c("Global Alignment", "Quantile Centering", "Mean Centering")
    paramtype<-c("sum by columns", "Alignment on all quantiles", "overall","within conditions")
    
    if (!(method %in% parammethod)){
        warning("Parameter method is not correct")
        return (NULL)
    }
    
    if (!(type %in% paramtype)){
        warning("Parameter type is not correct")
        return (NULL)
    }
    
        .temp <- qData
    if (!is.null(.temp)){
        data <- .temp
        
        
        
        ###############
        if (method == "Global Alignment"){
            if (type == "sum by columns"){
                t <- 2^(.temp)
                s <- colSums(t, na.rm=TRUE)
                data <- t
                for ( i in 1:nrow(data)) { data[i,] <- t[i,] / s}
                .temp <- log2(data)
            }
            else if (type == "Alignment on all quantiles"){
                .temp <- preprocessCore::normalize.quantiles(.temp)
                dimnames(.temp) <- list(rownames(qData),colnames(qData))
            }
        }
        
        
        
        ####################
        else if (method =="Quantile Centering"){
            q <- function(x) { quantile(x, probs=quantile, na.rm=TRUE) }
            medianOverSamples <- apply(.temp, 2, q)
            
            if (type == "overall"){
                cOverall <- q(medianOverSamples)
                .temp <- sweep(.temp, 2, medianOverSamples)
                .temp <- .temp + cOverall
                
            }
            else if (type == "within conditions"){
                .temp <- sweep(.temp, 2, medianOverSamples)
                
                cCond <- NULL
                for (l in unique(labels))
                {
                    indices <- which(labels== l)
                    cCond[l] <- q(medianOverSamples[indices])
                    .temp[,indices] <- .temp[,indices] + cCond[l]
                }
            }
        }
        
        
        
        ######################
        else if (method =="Mean Centering"){
            meanOverSamples <- apply(.temp, 2, mean, na.rm = TRUE)
            
            if (type == "overall"){
                cOverall <- mean(meanOverSamples)
                .temp <- sweep(.temp, 2, meanOverSamples)
                if (scaling){
                    .temp <- scale(.temp,center=FALSE,scale=TRUE)
                    attr(.temp,"scaled:scale")<-NULL 
                    }
                .temp <- .temp + cOverall
            } 
            else if (type == "within conditions"){
                .temp <- sweep(.temp, 2, meanOverSamples)
                if (scaling){
                    .temp <- scale(.temp,center=FALSE, scale=TRUE)
                    attr(.temp,"scaled:scale")<-NULL
                    }
                cCond <- NULL
                for (l in unique(labels))
                {
                    indices <- which(labels== l)
                    cCond[l] <- mean(meanOverSamples[indices])
                    .temp[,indices] <- .temp[,indices] + cCond[l]
                }
            }
        } 
        
# 
#              msg <- paste("Normalisation using ", method, " and ", type, sep="")
#              .temp@processingData@processing <- c(.temp@processingData@processing,
#                                                   msg)
    }
    return(.temp)
}

