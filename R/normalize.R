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
##' @param obj An object of class \code{\link{MSnSet}}
##' @param family One of the following : Global Rescaling, 
##' Median Centering, Mean Centering, Mean Centering Scaling
##' @param method "Overall" or "within conditions"
##' @return An instance of class \code{\link{MSnSet}} 
##' where the quantitative data in the \code{exprs()} tab
##' has been normalized.
##' @author Florence Combes, Samuel Wieczorek
##' @examples data(UPSprotx2)
##' normalizeD(UPSprotx2, "Median Centering", "within conditions")
normalizeD <- function(obj, family, method){
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


    
    .temp <- obj
    if (!is.null(.temp)){
        data <- exprs(.temp)
        if (family == "Global Rescaling"){
            if (method == "sum by columns"){
                
                t <- 2^(exprs(.temp))
                s <- colSums(t, na.rm=TRUE)
                data <- t
                for ( i in 1:nrow(data)) { data[i,] <- t[i,] / s}
                
                exprs(.temp) <- log2(data)
                }
            else if (method == "quantiles"){
                exprs(.temp) <- normalize.quantiles(exprs(.temp))
                dimnames(exprs(.temp)) <- list(rownames(exprs(obj)),
                                                colnames(exprs(obj)))
                }
        }else if (family =="Median Centering"){
            medianOverSamples <- apply(exprs(.temp),2,median,na.rm = TRUE)
            
            if (method == "overall"){
                cOverall <- median(medianOverSamples)
                exprs(.temp) <- sweep(exprs(.temp), 2, medianOverSamples)
                exprs(.temp) <- exprs(.temp) + cOverall
                
            }
            else if (method == "within conditions"){
                exprs(.temp) <- sweep(exprs(.temp), 2, medianOverSamples)
                
                labels <- unique(pData(.temp)$Label)
                cCond <- NULL
                for (l in labels)
                {
                    indices <- which(pData(obj)$Label== l)
                    cCond[l] <- median(medianOverSamples[indices])
                    exprs(.temp)[,indices] <- exprs(.temp)[,indices] + cCond[l]
                }
                
                }
        } else if (family =="Mean Centering"){
            meanOverSamples <- apply(exprs(.temp),2,mean,na.rm = TRUE)
            
            if (method == "overall"){
                cOverall <- mean(meanOverSamples)
                exprs(.temp) <- sweep(exprs(.temp), 2, meanOverSamples)
                exprs(.temp) <- exprs(.temp) + cOverall
                
            } else if (method == "within conditions"){
                exprs(.temp) <- sweep(exprs(.temp), 2, meanOverSamples)
                
                labels <- unique(pData(.temp)$Label)
                cCond <- NULL
                for (l in labels)
                {
                    indices <- which(pData(obj)$Label== l)
                    cCond[l] <- mean(meanOverSamples[indices])
                    exprs(.temp)[,indices] <- exprs(.temp)[,indices] + cCond[l]
                }
            }
        } else if (family == "Mean Centering Scaling"){
            meanOverSamples <- apply(exprs(.temp),2,mean,na.rm = TRUE)
            
            if (method == "overall"){
                cOverall <- mean(meanOverSamples)
                exprs(.temp) <- sweep(exprs(.temp), 2, meanOverSamples)
                exprs(.temp) <- scale(exprs(.temp),center=FALSE,
                                    scale=TRUE)
                attr(exprs(.temp),"scaled:scale")<-NULL
                exprs(.temp) <- exprs(.temp) + cOverall
                
            
            } else if (method == "within conditions"){
                exprs(.temp) <- sweep(exprs(.temp), 2, meanOverSamples)
                exprs(.temp) <- scale(exprs(.temp),center=FALSE,
                                    scale=TRUE)
                attr(exprs(.temp),"scaled:scale")<-NULL
                labels <- unique(pData(.temp)$Label)
                cCond <- NULL
                for (l in labels)
                {
                    indices <- which(pData(obj)$Label== l)
                    cCond[l] <- mean(meanOverSamples[indices])
                    exprs(.temp)[,indices] <- exprs(.temp)[,indices] + cCond[l]
                }
                
            }
        }

        msg <- paste("Normalisation using ", method,  sep="")
        .temp@processingData@processing <- c(.temp@processingData@processing,
                                            msg)
    }
    return(.temp)
}
