##' Returns the number of empty lines in a matrix.
##' 
##' @title Returns the number of empty lines in the data
##' @param qData A matrix corresponding to the quantitative data.
##' @return An integer
##' @author Samuel Wieczorek
##' @examples
##' library(DAPARdata)
##' data(UPSpep25)
##' qData <- exprs(UPSpep25)
##' getNumberOfEmptyLines(qData)
getNumberOfEmptyLines <- function(qData){
  
  n <- sum(apply(is.na(as.matrix(qData)), 1, all))
  return (n)
}




##' Returns a list for the two conditions where each slot is a vector of indices
##' for the samples
##' 
##' @title Gets the conditions indices
##' @param labels A vector of strings containing the column "Label" of the \code{pData()}.
##' @param cond1 A vector of Labels (a slot in the \code{pData()} table) for
##' the condition 1.
##' @param cond2 A vector of Labels (a slot in the \code{pData()} table) for
##' the condition 2.
##' @return A list with two slots \code{iCond1} and \code{iCond2} containing
##' respectively the indices of samples in the \code{pData()} table of the
##' dataset. 
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' library(DAPARdata)
##' data(UPSpep25)
##' labels <- pData(UPSpep25)[,"Label"]
##' getIndicesConditions(labels, "25fmol", "10fmol")
getIndicesConditions <- function(labels, cond1, cond2){
  indCondition1 <- indCondition2 <- NULL
  
  for(i in 1:length(cond1)){
    indCondition1 <- c(indCondition1,
                       which(labels == cond1[i]))
  }
  for(i in 1:length(cond2)){
    indCondition2 <- c(indCondition2,
                       which(labels == cond2[i]))
  }
  
  return(list(iCond1 = indCondition1, iCond2 = indCondition2))
}



##' Selects colors for the plots in DAPAR based on the different conditions in
##' the dataset. The palette is derived from
##' the brewer palette "Dark2" (see \code{\link{RColorBrewer}}).
##' 
##' @title Palette for plots in DAPAR
##' @param labels A vector of labels (strings).
##' @return A palette designed for the data manipulated in DAPAR
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' library(DAPARdata)
##' data(UPSpep25)
##' labels <- pData(UPSpep25)[,"Label"]
##' getPaletteForLabels(labels)
getPaletteForLabels <- function(labels){
  nColors <- 8
  col <- c(1:nColors)
  palette(brewer.pal(nColors,"Dark2"))
  
  ## Define one color per label/condition
  col.boxplot <- NULL
  for (i in 1:length(labels)){
    col.boxplot[which(labels == unique(labels)[i])] <- col[i]
  }
  
  return (col.boxplot)
}

##' Selects colors for the plots in DAPAR based on the replicates in
##' the dataset. The palette is derived from
##' the brewer palette "Dark2" (see \code{\link{RColorBrewer}}).
##' 
##' @title Palette for plot the replicates in DAPAR
##' @param nColors The desired number of colors
##' @return A palette designed for the data manipulated in DAPAR
##' @author Samuel Wieczorek
##' @examples
##' library(DAPARdata)
##' data(UPSpep25)
##' n <- nrow(pData(UPSpep25))
##' getPaletteForLabels(5)
getPaletteForReplicates <- function(nColors){
  col <- c(1:nColors)
  getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))
  
  return(col)
}



