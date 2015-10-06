##' Returns the percentage of missing values in the quantitative
##' data (\code{exprs()} table of the dataset).
##' 
##' @title Percentage of missing values
##' @param obj An object of class \code{\link{MSnSet}}.
##' @return A floating number
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' data(UPSprotx2)
##' getPourcentageOfMV(UPSprotx2)
getPourcentageOfMV <- function(obj){

NA.count<-apply(data.frame(exprs(obj)), 2, 
                function(x) length(which(is.na(data.frame(x))==TRUE)) )
pourcentage <- 100 * round(sum(NA.count)
                            /(dim(exprs(obj))[1]*dim(exprs(obj))[2]), 
                            digits=4)

return(pourcentage)
}

##' Returns the number of empty lines in the quantitative data
##' (i.e. \code{exprs()} table).
##' 
##' @title Returns the number of empty lines in the data
##' @param obj An object of class \code{\link{MSnSet}}.
##' @return An integer
##' @author Samuel Wieczorek
##' @examples
##' data(UPSprotx2)
##' getNumberOfEmptyLines(UPSprotx2)
getNumberOfEmptyLines <- function(obj){

n <- sum(apply(is.na(as.matrix(exprs(obj))), 1, all))
return (n)
}



##' Draws a pie chart of the missing values in the quantitative data
##' (i.e. \code{exprs()} table).
##' 
##' @title Pie chart of the missing values
##' @param obj An object of class \code{\link{MSnSet}}.
##' @return A pie chart
##' @author Samuel Wieczorek
##' @examples
##' data(UPSprotx2)
##' mvPieChart(UPSprotx2)
mvPieChart <- function(obj){
n <- getPourcentageOfMV(obj)

slices <- c(100-n, n) 
lbls <- c("Quant. data", "Missing values")
pct <- c(100-n, n)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls)
}



##' Filters the lines of \code{exprs()} table with conditions on the number
##' of missing values.
##' The user chooses the minimum amount of intensities that is acceptable and
##' the filter delete lines that do not respect this condition.
##' The condition may be on the whole line or condition by condition.
##' 
##' The different methods are :
##' "wholeMatrix": given a threshold \code{th}, only the lines that contain
##' at least \code{th} values are kept.
##' "allCond": given a threshold \code{th}, only the lines which contain
##' at least \code{th} values for each of the conditions are kept.
##' "atLeastOneCond": given a threshold \code{th}, only the lines that contain
##' at least \code{th} values, and for at least one condition, are kept.
##' 
##' @title Filter lines in the matrix of intensities w.r.t. some criteria
##' @param obj An object of class \code{\link{MSnSet}} containing
##' quantitative data.
##' @param type Method used to choose the lines to delete.
##' Values are : "none", "wholeMatrix", "allCond", "atLeastOneCond"
##' @param th An integer value of the threshold
##' @param processText A string to be included in the \code{\link{MSnSet}}
##' object for log. 
##' @return An instance of class \code{\link{MSnSet}} that have been filtered.
##' @author Florence Combes, Samuel Wieczorek
##' @examples data(UPSprotx2)
##' mvFilter(UPSprotx2, "wholeMatrix", 2)
mvFilter <- function(obj,type, th, processText=NULL )
{
    #Check parameters
    paramtype<-c("none", "wholeMatrix", "allCond", "atLeastOneCond") 
    if (sum(is.na(match(type, paramtype)==TRUE))>0){
        warning("Param type is not correct.")
        return (NULL)
    }

    paramth<-c(seq(0, nrow(pData(obj)), 1))
    if (sum(is.na(match(th, paramth)==TRUE))>0){
        warning("Param th is not correct.")
        return (NULL)
    }

    keepThat <- NULL

    if (type == "none"){
        keepThat <- seq(1:dim(exprs(obj))[1])
    } else if (type == "wholeMatrix"){
        for (i in 1:length(rownames(exprs(obj))))
            {
            if(sum(!is.na(exprs(obj)[i,])) >= th)  {
                keepThat <- c(keepThat,i) 
                }
            }
        } else if (type == "atLeastOneCond" || type == "allCond"){
        conditions <- unique(pData(obj)$Label)
        nbCond <- length(conditions)
        keepThat <- NULL
        for (i in 1:length(rownames(exprs(obj))))
            {
            t <- (type == "allCond")
            for (cond in 1:nbCond){
                ind <- which( pData(obj)$Label == conditions[cond])
                if (type == "atLeastOneCond") {
                    t <- t || (sum(!is.na(exprs(obj)[i,ind])) >= th)
                } else if (type == "allCond") {
            t <- t && (sum(!is.na(exprs(obj)[i,ind])) >= th)}
            }
        if(t)  {keepThat <- c(keepThat,i) }
        }
    }

obj <- obj[keepThat,]

    obj@processingData@processing <- 
        c(obj@processingData@processing, processText)
    return(obj)
}
