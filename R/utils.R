#' Returns a list for the two conditions where each slot is a vector of indices
#' for the samples in the \code{pData()} table
#' 
#' @title Gets the conditions indices
#' @param obj An object of class \code{\link{MSnSet}}.
#' @param cond1 A vector of Labels (a slot in the \code{pData()} table) for
#' the condition 1.
#' @param cond2 A vector of Labels (a slot in the \code{pData()} table) for
#' the condition 2.
#' @return A list with two slots \code{iCond1} and \code{iCond2} containing
#' respectively the indices of samples in the \code{pData()} table of the
#' dataset. 
#' @author Florence Combes, Samuel Wieczorek
#' @examples data(UPSprotx2)
#' getIndicesConditions(UPSprotx2, "Cut3", "WT")
getIndicesConditions <- function(obj, cond1, cond2){
    indCondition1 <- indCondition2 <- NULL
    
    for(i in 1:length(cond1)){
        indCondition1 <- c(indCondition1,
                            which(pData(obj)[,"Label"] == cond1[i]))
    }
    for(i in 1:length(cond2)){
        indCondition2 <- c(indCondition2,
                            which(pData(obj)[,"Label"] == cond2[i]))
    }
    
    return(list(iCond1 = indCondition1, iCond2 = indCondition2))
}


##' Selects colors for the plots in DAPAR based on the different conditions in
##' the dataset. The palette is derived from
##' the brewer palette "Dark2" (see \code{\link{RColorBrewer}}).
##' 
##' @title Palette for plots in DAPAR
##' @param obj An object of class \code{\link{MSnSet}}.
##' @return A palette designed for the data manipulated in DAPAR
##' @author Florence Combes, Samuel Wieczorek
##' @examples data(UPSprotx2)
##' getPaletteForLabels(UPSprotx2)
getPaletteForLabels <- function(obj){
    nColors <- 8
    col <- c(1:nColors)
    palette(brewer.pal(nColors,"Dark2"))
    
    Label <- pData(obj)[,"Label"]
    
    ## Define one color per label/condition
    col.boxplot <- NULL
    for (i in 1:length(Label)){
        col.boxplot[which(Label==unique(Label)[i])] <- col[i]
    }
    
    return (col.boxplot)
}
