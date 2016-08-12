##' Returns the percentage of missing values in the quantitative
##' data (\code{exprs()} table of the dataset).
##' 
##' @title Percentage of missing values
##' @param obj An object of class \code{\link{MSnSet}}.
##' @return A floating number
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' getPourcentageOfMV(UPSpep25)
getPourcentageOfMV <- function(obj){

NA.count<-apply(data.frame(Biobase::exprs(obj)), 2, 
                function(x) length(which(is.na(data.frame(x))==TRUE)) )
pourcentage <- 100 * round(sum(NA.count)
                            /(dim(Biobase::exprs(obj))[1]*dim(Biobase::exprs(obj))[2]), 
                            digits=4)

return(pourcentage)
}

##' Returns the number of lines, in a given column, where content matches 
##' the prefix.
##' 
##' @title Number of lines with prefix
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param name The name of a column.
##' @param prefix A string
##' @return An integer
##' @author Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' getNumberOf(UPSpep25, "Potential.contaminant", "+")
getNumberOf <- function(obj, name=NULL, prefix=NULL){
if (is.null(name) || is.null(prefix) || (name=="") || (prefix=="")){
    return(0)}
if (!(is.null(name) || !is.null(name=="")) 
    && (is.null(prefix) || (prefix==""))){return(0)}

if(nchar(prefix) > 0){
    count <- length(which(substr(Biobase::fData(obj)[,name], 0, 1) == prefix))
} else { count <- 0}

return(count)
}


##' Plots a barplot of proportion of contaminants and reverse
##' 
##' @title Barplot of proportion of contaminants and reverse
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param idContaminants The name of a column of Contaminants
##' @param prefixContaminants The prefix to identify contaminants
##' @param idReverse The name of a column of Reverse
##' @param prefixReverse The prefix to identify Reverse
##' @return A barplot
##' @author Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' pref <- "+"
##' proportionConRev(UPSpep25, "Potential.contaminant", pref, "Reverse", pref)
proportionConRev <- function(obj, idContaminants=NULL, 
                            prefixContaminants=NULL, 
                            idReverse=NULL, prefixReverse=NULL){
#if (is.null(prefixContaminants) && is.null(prefixReverse) ){return(NULL)}
if (is.null(obj) ){return(NULL)}
nContaminants <- nReverse <- 0

nContaminants <- length(getIndicesOfLinesToRemove(obj, idContaminants, prefixContaminants))
nReverse <- length(getIndicesOfLinesToRemove(obj, idReverse, prefixReverse))

pctContaminants <- 100 * round(nContaminants/nrow(Biobase::fData(obj)),  digits=4)
pctReverse <- 100 * round(nReverse/nrow(Biobase::fData(obj)),  digits=4)

counts <- c(nrow(Biobase::fData(obj))-nContaminants-nReverse, nContaminants, 
            nReverse )
slices <- c(100-pctContaminants-pctReverse, pctContaminants, pctReverse ) 
lbls <- c("Quantitative data", "Contaminants", "Reverse")
pct <- c(100-pctContaminants-pctReverse, pctContaminants, pctReverse )
lbls <- paste(lbls, " : ", counts, " lines (", pct, "%)", sep="") 
#lbls <- paste(lbls,"%",sep="") # ad % to labels 

bp <- barplot(slices,
        #names.arg = lbls, 
        horiz = TRUE,
        col=c("lightgrey", "green", "blue"),
        axes = FALSE,
        cex.names = 1.5)

graphics::text(x= 20, y= bp, labels=as.character(lbls), xpd=TRUE, cex=1.5)

}





##' This function removes lines in the dataset based on a prefix string.
##' 
##' @title Removes lines in the dataset based on a prefix string.
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param idLine2Delete The name of the column that correspond to the 
##' data to filter
##' @param prefix A character string that is the prefix to find in the data
##' @return An object of class \code{\link{MSnSet}}.
##' @author Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' removeLines(UPSpep25, "Potential.contaminant")
##' removeLines(UPSpep25, "Reverse")
removeLines <- function(obj, idLine2Delete=NULL, prefix=NULL){
if ((prefix == "") || is.null(prefix)) {
    #warning ("No change was made")
    return (obj)}
    t <- (prefix == substring(Biobase::fData(obj)[,idLine2Delete],1,nchar(prefix)))
    ind <- which( t== TRUE)
    obj <- obj[-ind ]

return(obj)
}

##' This function returns the indice of the lines to delete, based on a 
##' prefix string
##' 
##' @title Get the indices of the lines to delete, based on a prefix string
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param idLine2Delete The name of the column that correspond to the data 
##' to filter
##' @param prefix A character string that is the prefix to find in the data
##' @return A vector of integers.
##' @author Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' getIndicesOfLinesToRemove(UPSpep25, "Potential.contaminant", prefix="+")
getIndicesOfLinesToRemove <- function(obj, idLine2Delete=NULL, prefix=NULL)
{
if ((prefix == "") || is.null(prefix)) {
   # warning ("No change was made")
    return (NULL)}
t <- (prefix == substring(Biobase::fData(obj)[,idLine2Delete],1,nchar(prefix)))
ind <- which( t== TRUE)
return(ind)
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
##' @examples
##' data(UPSpep25)
##' mvFilter(UPSpep25, "wholeMatrix", 2)
mvFilter <- function(obj,type, th, processText=NULL )
{
    #Check parameters
    paramtype<-c("none", "wholeMatrix", "allCond", "atLeastOneCond") 
    if (sum(is.na(match(type, paramtype)==TRUE))>0){
        warning("Param type is not correct.")
        return (NULL)
    }

    paramth<-c(seq(0, nrow(Biobase::pData(obj)), 1))
    if (sum(is.na(match(th, paramth)==TRUE))>0){
        warning("Param th is not correct.")
        return (NULL)
    }
    
    if(!is.integer(th)){th <- as.integer(th)}

    keepThat <- mvFilterGetIndices(obj,type, th)

obj <- obj[keepThat,]

    obj@processingData@processing <- 
        c(obj@processingData@processing, processText)
    return(obj)
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
##' @param keepThat A vector of integers which are the indices of lines to 
##' keep.
##' @param processText A string to be included in the \code{\link{MSnSet}}
##' object for log. 
##' @return An instance of class \code{\link{MSnSet}} that have been filtered.
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' mvFilter(UPSpep25, c(1:10))
mvFilterFromIndices <- function(obj,keepThat=NULL, processText=NULL )
{

if (is.null(keepThat)) {return(obj)}
obj <- obj[keepThat,]
obj@processingData@processing <- 
    c(obj@processingData@processing, processText)

txt <- unlist(strsplit(processText, split=" "))
obj@experimentData@other$mvFilter.method <- txt[3]
obj@experimentData@other$mvFilter.threshold <- txt[6]
return(obj)
}

##' Delete the lines of \code{exprs()} table identified by their indice.
##' 
##' @title Delete the lines in the matrix of intensities and the metadata table
##' given their indice.
##' @param obj An object of class \code{\link{MSnSet}} containing
##' quantitative data.
##' @param deleteThat A vector of integers which are the indices of lines to 
##' delete.
##' @param processText A string to be included in the \code{\link{MSnSet}}
##' object for log. 
##' @return An instance of class \code{\link{MSnSet}} that have been filtered.
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' mvFilter(UPSpep25, c(1:10))
deleteLinesFromIndices <- function(obj,deleteThat=NULL, processText=NULL )
{
    
    if (is.null(deleteThat)) {return(obj)}
    obj <- obj[-deleteThat,]
    obj@processingData@processing <- 
        c(obj@processingData@processing, processText)
    txt <- unlist(strsplit(processText, split=" "))
    if (txt[2] == "contaminants"){
    obj@experimentData@other$contaminantsRemoved <- TRUE
    } else if (txt[2] == "reverse"){
        obj@experimentData@other$reverseRemoved <- TRUE
    }
    return(obj)
}


##' Returns the indices of the lines of \code{exprs()} table to delete w.r.t. 
##' the conditions on the number of missing values.
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
##' @return An vector of indices that correspond to the lines to keep.
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' mvFilterGetIndices(UPSpep25, "wholeMatrix", 2)
mvFilterGetIndices <- function(obj,type, th)
{
#Check parameters
paramtype<-c("none", "wholeMatrix", "allCond", "atLeastOneCond") 
if (sum(is.na(match(type, paramtype)==TRUE))>0){
    warning("Param type is not correct.")
    return (NULL)
}

paramth<-c(seq(0, nrow(Biobase::pData(obj)), 1))
if (sum(is.na(match(th, paramth)==TRUE))>0){
    warning("Param th is not correct.")
    return (NULL)
}

keepThat <- NULL

if (type == "none"){
    keepThat <- seq(1:nrow(Biobase::exprs(obj)))
} else if (type == "wholeMatrix"){
    keepThat <- which(apply(!is.na(Biobase::exprs(obj)), 1, sum) >= th)
} else if (type == "atLeastOneCond" || type == "allCond"){
    
    conditions <- unique(Biobase::pData(obj)$Label)
    nbCond <- length(conditions)
    keepThat <- NULL
    s <- matrix(rep(0, nrow(Biobase::exprs(obj))*nbCond),nrow=nrow(Biobase::exprs(obj)), 
                ncol=nbCond)
    
    for (c in 1:nbCond){
        ind <- which(Biobase::pData(obj)$Label == conditions[c])
        s[,c] <- (apply(!is.na(Biobase::exprs(obj)[,ind]), 1, sum) >= th)
    }
    
    
    if (type == "allCond") {
        keepThat <- which(rowSums(s) == nbCond)
    }
    else if (type == "atLeastOneCond") {
        keepThat <- which(rowSums(s) >= 1)
    }
}
return(keepThat)
}
