##' Provides several methods to normalize quantitative data from
##' a \code{MSnSet} object.
##' They are organized in six main families : GlobalQuantileAlignement, sumByColumns, QuantileCentering, MeanCentering, LOESS, vsn
##' For the first family, there is no type.
##' For the five other families, two type categories are available :
##' "Overall" which means that the value for each protein
##' (ie line in the expression data tab) is computed over all the samples ;
##' "within conditions" which means that the value for each protein
##' (ie line in the \code{exprs()} data tab) is computed condition
##' by condition.
##'
##' @title Normalisation
##' @param obj An object of class \code{MSnSet}.
##' @param method One of the following : "GlobalQuantileAlignment" (for
##' normalizations of important magnitude), "SumByColumns", "QuantileCentering",
##' "Mean Centering", "LOESS" and "vsn".
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
##' @param ... span parameter for LOESS method
##' @return An instance of class \code{MSnSet} where the quantitative
##' data in the \code{exprs()} tab has been normalized.
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' wrapper.normalizeD(Exp1_R25_pept[1:1000], "QuantileCentering", "within conditions")
wrapper.normalizeD <- function(obj, method, type=NULL, scaling=FALSE, quantile=0.15, ...){
  
  qData <- Biobase::exprs(obj)
  conds <- Biobase::pData(obj)[,"Condition"]
  Biobase::exprs(obj) <- normalizeD(qData, conds, method, type, scaling, quantile, ...)
  
  msg_method <- paste("Normalisation using method =", method,  sep="")
  msg_type <- paste("With type =", type,  sep="")
  
  switch(method,
  GlobalQuantileAlignment ={
    obj@processingData@processing <- c(obj@processingData@processing, msg_method, msg_type)
    obj@experimentData@other$normalizationMethod <- method
            },
  SumByColumns= {
    obj@processingData@processing <- c(obj@processingData@processing, msg_method, msg_type)
    obj@experimentData@other$normalizationMethod <- method
    obj@experimentData@other$normalizationType <- type
            },
  QuantileCentering = {
    msg_quantile <- paste("With quantile =", quantile,  sep="")
    obj@processingData@processing <- c(obj@processingData@processing, msg_method, msg_type, msg_quantile)
    obj@experimentData@other$normalizationMethod <- method
    obj@experimentData@other$normalizationType <- type
    obj@experimentData@other$normalizationQuantile <- quantile
            },
  MeanCentering = {
    msg_scaling <- paste("With scaling =", scaling,  sep="")
    obj@processingData@processing <- c(obj@processingData@processing, msg_method, msg_type, msg_scaling)
    obj@experimentData@other$normalizationMethod <- msg_method
    obj@experimentData@other$normalizationType <- type
    obj@experimentData@other$normalizationScaling <- scaling
          },
  LOESS = {
    msg_loess <- paste("With span =", ...,  sep="")
    obj@processingData@processing <- c(obj@processingData@processing, msg_method, msg_type, msg_loess)
    obj@experimentData@other$normalizationMethod <- msg_method
    obj@experimentData@other$normalizationType <- type
          },
  vsn ={
    obj@processingData@processing <- c(obj@processingData@processing, msg_method, msg_type)
    obj@experimentData@other$normalizationMethod <- method
    obj@experimentData@other$normalizationType <- type
  }
  )
  return(obj)
}




##' Provides several methods to normalize data from
##' a matrix.
##' They are organized in 6 main families : GlobalQuantileAlignment, sumByColumns, QuantileCentering, MeanCentering, LOESS and vsn.
##' For the first famlily there is no type.
##' For the five other families, two categories are available :
##' "Overall" which means that the value for each protein
##' (ie line in the expression data tab) is computed over all the samples ;
##' "within conditions" which means that the value for each protein
##' (ie line in the matrix) is computed condition
##' by condition.
##'
##' @title Normalisation
##' @param qData A dataframe that contains quantitative data.
##' @param conds A vector of strings containing the column "Condition" of
##' the \code{pData()}.
##' @param method One of the following : "GlobalQuantileAlignment", "SumByColumns", "QuantileCentering", "MeanCentering",
##' "LOESS", "vsn".
##' @param type For the method "GlobalQuantileAlignment" there is no associated type.
##' For the six other methods, the parameters are "overall" (shift all the
##' sample distributions at once) or "within conditions" (shift the sample
##' distributions within each condition at a time).
##' @param scaling A boolean that indicates if the variance of the data have to
##' be forced to unit (variance reduction) or not.
##' @param quantile A float that corresponds to the quantile used to align the
##' data.
##' @param ... span parameter for LOESS method
##' @return A matrix normalized
##' @author Samuel Wieczorek, Thomas Burger, Helene Borges
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept[1:1000])
##' conds <- Biobase::pData(Exp1_R25_pept[1:1000])[,"Condition"]
##' normalizeD(qData, conds, "QuantileCentering", "within conditions", quantile = 0.15)
normalizeD <- function(qData, conds, method, type=NULL, scaling=FALSE, quantile=0.15, ...){
  parammethod<-c("GlobalQuantileAlignment",
                 "SumByColumns",
                 "QuantileCentering",
                 "MeanCentering",
                 "LOESS",
                 "vsn")
  
  if (sum(is.na(match(method, parammethod)==TRUE))>0){
    warning("Parameter method is not correct")
    return (NULL)
  }
  
  paramtype<-c(NULL,
               "overall",
               "within conditions")
  if (sum(is.na(match(type, paramtype)==TRUE))>0){
    warning("Parameter type is not correct")
    return (NULL)
  }
  
  .temp <- qData
  if (!is.null(.temp)){
    data <- .temp
    
    ###############
    
    switch(method,
      GlobalQuantileAlignment = {
              .temp <- normalize.quantiles(.temp)
              dimnames(.temp) <- list(rownames(qData),colnames(qData))
              },
    SumByColumns = {
              t <- 2^(.temp)
      
              if (type == "overall"){
                  sums_cols <- colSums(t, na.rm=TRUE)
                  #normalisation
                  for ( i in 1:nrow(t)) {
                    t[i, ] <- (t[i, ] / sums_cols)*median(sums_cols)
                    }
              } else if (type == "within conditions"){
                  for (l in unique(conds)) {
                    indices <- which(conds== l)
                    sums_cols <- colSums(t[,indices], na.rm=TRUE)
                    for (i in 1:nrow(t)){
                      t[i,indices] <- (t[i,indices]/sums_cols) * median(sums_cols)
                    }
                  }
              }
      
          .temp <- log2(t)
        },
    
    QuantileCentering = {
      q <- function(x) { quantile(x, probs=quantile, na.rm=TRUE) }
      medianOverSamples <- apply(.temp, 2, q)
      
      if (type == "overall"){
        cOverall <- q(medianOverSamples)
        .temp <- sweep(.temp, 2, medianOverSamples)
        .temp <- .temp + cOverall
      } else if (type == "within conditions"){
        .temp <- sweep(.temp, 2, medianOverSamples)
         cCond <- NULL
        for (l in unique(conds))
          {
          indices <- which(conds== l)
          cCond[l] <- q(medianOverSamples[indices])
          .temp[,indices] <- .temp[,indices] + cCond[l]
          }
        }
    },
    MeanCentering = {
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
        for (l in unique(conds))
        {
          indices <- which(conds== l)
          cCond[l] <- mean(meanOverSamples[indices])
          .temp[,indices] <- .temp[,indices] + cCond[l]
        }
      }
    },
    vsn = {
      if(type == "overall"){
        vsn.fit <- vsn::vsnMatrix(2^(.temp))
        .temp <- vsn::predict(vsn.fit, 2^(.temp))
      }else if(type == "within conditions"){
        for (l in unique(labels)) {
          indices <- which(labels == l)
          vsn.fit <- vsn::vsnMatrix(2^(.temp[,indices]))
          .temp[,indices] <- vsn::predict(vsn.fit, 2^(.temp[,indices]))
          
        }
      }
      
    },
    
    ###############
    # data must be log-expressed.
    LOESS = {
      if(type == "overall"){
        # "..." to pass the span parameter to the loess function
        .temp <- normalizeCyclicLoess(x = .temp, method = "fast", ...)
        
      }else if(type == "within conditions"){
        for (l in unique(labels)) {
          indices <- which(labels == l)
          .temp[,indices] <- normalizeCyclicLoess(x = .temp[,indices],
                                                  method = "fast", ...)
        }
      }
      
    }
    )
    
    #
    #              msg <- paste("Normalisation using ", method, " and ", type, sep="")
    #              .temp@processingData@processing <- c(.temp@processingData@processing,
    #                                                   msg)
  }
  return(.temp)
}