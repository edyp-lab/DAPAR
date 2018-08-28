##' Returns the contents of the slot processing of an object of class \code{MSnSet}
##' 
##' @title Returns the contains of the slot processing  of an object of 
##' class \code{MSnSet}
##' @param  obj An object (peptides) of class \code{MSnSet}.
##' @return The slot processing of obj@processingData
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' getProcessingInfo(Exp1_R25_pept)
getProcessingInfo <- function(obj){
return(obj@processingData@processing)
}

##' Returns the number of empty lines in a matrix.
##' 
##' @title Returns the number of empty lines in the data
##' @param qData A matrix corresponding to the quantitative data.
##' @return An integer
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' qData <- Biobase::exprs(Exp1_R25_pept)
##' getNumberOfEmptyLines(qData)
getNumberOfEmptyLines <- function(qData){
n <- sum(apply(is.na(as.matrix(qData)), 1, all))
return (n)
}

##' Similar to the function \code{is.na} but focused on the equality with the paramter 'type'.
##' 
##' @title Similar to the function \code{is.na} but focused on the equality with the paramter 'type'.
##' @param data A data.frame
##' @param type The value to search in the dataframe
##' @return A boolean dataframe 
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept
##' data <- Biobase::fData(obj)[,obj@experimentData@other$OriginOfValues]
##' is.OfType(data, "MEC")
is.OfType <- function(data, type){
  return (type == data)
}





##' Similar to the function \code{is.na} but focused on the equality with the missing 
##' values in the dataset (type 'POV' and 'MEC')
##' 
##' @title Similar to the function \code{is.na} but focused on the equality with the missing 
##' values in the dataset (type 'POV' and 'MEC')
##' @param data A data.frame
##' @return A boolean dataframe 
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj <- Exp1_R25_pept
##' data <- Biobase::fData(obj)[,obj@experimentData@other$OriginOfValues]
##' is.MV(data)
is.MV <- function(data){
  #MV=is.OfType(data, "MV")
  POV=is.OfType(data, "POV")
  MEC=is.OfType(data, "MEC")
  
  df <- POV | MEC
  return (df)
}

##' Returns the possible number of values in lines in a matrix.
##' 
##' @title Returns the possible number of values in lines in the data
##' @param obj An object of class \code{MSnSet}
##' @param type xxxxxxx
##' @return An integer
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' getListNbValuesInLines(Exp1_R25_pept)
getListNbValuesInLines <- function(obj, type="wholeMatrix"){
  if (is.null(obj)){return()}
  
  if(is.null(obj@experimentData@other$OriginOfValues)){
    ll <- seq(0,ncol(obj))
  }
  data <- Biobase::fData(obj)[,obj@experimentData@other$OriginOfValues]
  switch(type,
         wholeMatrix= {
           ll <- unique(ncol(data) - apply(is.MV(data), 1, sum))
           },
         allCond = {
                    tmp <- NULL
                    for (cond in unique(Biobase::pData(obj)$Condition)){
                     tmp <- c(tmp, length(which(Biobase::pData(obj)$Condition== cond)))
                  }
                  ll <- seq(0,min(tmp))
                  },
         atLeastOneCond = {
                   tmp <- NULL
                  for (cond in unique(Biobase::pData(obj)$Condition)){
                       tmp <- c(tmp, length(which(Biobase::pData(obj)$Condition== cond)))
                   }
                   ll <- seq(0,max(tmp))
                    }
         )
  
  return (sort(ll))
}



##' Returns a list for the two conditions where each slot is a vector of 
##' indices for the samples.
##' 
##' @title Gets the conditions indices.
##' @param conds A vector of strings containing the column "Condition" of 
##' the \code{pData()}.
##' @param cond1 A vector of Conditions (a slot in the \code{pData()} table) for
##' the condition 1.
##' @param cond2 A vector of Conditions (a slot in the \code{pData()} table) for
##' the condition 2.
##' @return A list with two slots \code{iCond1} and \code{iCond2} containing
##' respectively the indices of samples in the \code{pData()} table of the
##' dataset. 
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
##' getIndicesConditions(conds, "25fmol", "10fmol")
getIndicesConditions <- function(conds, cond1, cond2){
indCondition1 <- indCondition2 <- NULL

for(i in 1:length(cond1)){
    indCondition1 <- c(indCondition1,
                        which(conds == cond1[i]))
}
for(i in 1:length(cond2)){
    indCondition2 <- c(indCondition2,
                        which(conds == cond2[i]))
}

return(list(iCond1 = indCondition1, iCond2 = indCondition2))
}



##' Selects colors for the plots in \code{DAPAR} based on the different conditions in
##' the dataset. The palette is derived from
##' the brewer palette "Dark2" (see \code{\link{RColorBrewer}}).
##' 
##' @title Palette for plots in \code{DAPAR}
##' @param conds A vector of conditions (strings).
##' @return A palette designed for the data manipulated in \code{DAPAR}
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
##' getPaletteForConditions(conds)
getPaletteForConditions <- function(conds){
    nColors <- 8
    col <- c(1:nColors)
    palette(brewer.pal(nColors,"Dark2"))
    
    ## Define one color per condition
    col.boxplot <- NULL
    for (i in 1:length(conds)){
        col.boxplot[which(conds == unique(conds)[i])] <- col[i]
    }
    
    return (col.boxplot)
}


##' Selects colors for the highcharter plots in \code{DAPAR} based on the different conditions in
##' the dataset. The palette is derived from
##' the brewer palette "Dark2" (see \code{\link{RColorBrewer}}).
##' 
##' @title Palette for highcharter plots used in \code{DAPAR}
##' @param conds A vector of conditions (strings).
##' @return A palette designed for the data manipulated in \code{DAPAR}
##' @author Florence Combes, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
##' getPaletteForConditions_HC(conds)
getPaletteForConditions_HC <- function(conds){
  nColors <- 8
  #col <- c(1:nColors)
  #pal <- c('#7cb5ec', '#434348', '#90ed7d', '#f7a35c', '#8085e9', '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1')
  #pal <- c("#002F80", "#002F80","#002F80","#002F80","#F9AF38","#F9AF38","#F9AF38","#F9AF38")
  pal <- brewer.pal(nColors, "Dark2")
  ## Define one color per condition
  col.boxplot <- NULL
  for (i in 1:length(conds)){
    col.boxplot[which(conds == unique(conds)[i])] <- pal[i]
  }
  
  return (col.boxplot)
}

##' Selects colors for the highcharter plots in \code{DAPAR} based on the replicates in
##' the dataset. The palette is derived from the brewer palette "Dark2" (see \code{\link{RColorBrewer}}).
##' 
##' @title Palette for highcharter plot the replicates in DAPAR
##' @param nColors The desired number of colors
##' @return A palette designed for the data manipulated in \code{DAPAR}
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' n <- nrow(Biobase::pData(Exp1_R25_pept))
##' getPaletteForReplicates_HC(n)
getPaletteForReplicates_HC <- function(nColors){
  #col <- c(1:nColors)
  #pal <- c("#002F80", "#002F80","#002F80","#002F80","#F9AF38","#F9AF38","#F9AF38","#F9AF38")
  pal <- brewer.pal(nColors, "Dark2")
  return(pal)
}



##' Selects colors for the plots in \code{DAPAR} based on the replicates in
##' the dataset. The palette is derived from
##' the brewer palette "Dark2" (see \code{\link{RColorBrewer}}).
##' 
##' @title Palette for plot the replicates in \code{DAPAR}
##' @param nColors The desired number of colors
##' @return A palette designed for the data manipulated in \code{DAPAR}
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' n <- nrow(Biobase::pData(Exp1_R25_pept))
##' getPaletteForReplicates(n)
getPaletteForReplicates <- function(nColors){
#col <- c(1:nColors)
#getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))
#return(col)
pal <- brewer.pal(nColors, "Dark2")
return(pal)
}



##' Customise the contextual menu of highcharts plots.
##' 
##' @title Customised contextual menu of highcharts plots
##' @param hc A highcharter object
##' @param filename The filename under which the plot has to be saved
##' @return A contextual menu for highcharts plots
##' @author Samuel Wieczorek
##' @examples
##' library("highcharter")
##' hc <- highchart() 
##' hc_chart(hc,type = "line") 
##' hc_add_series(hc,data = c(29, 71, 40))
##' my_hc_ExportMenu(hc,filename='foo')
my_hc_ExportMenu <- function(hc, filename){
  hc_exporting(hc, enabled=TRUE,
               filename = filename,
               buttons= list(
                 contextButton= list(
                   menuItems= list('downloadPNG', 'downloadSVG','downloadPDF')
                 )
               )
  )
}



##' Customise the resetZoomButton of highcharts plots.
##' 
##' @title Customised resetZoomButton of highcharts plots
##' @param hc A highcharter object
##' @param chartType The type of the plot
##' @param zoomType The type of the zoom (one of "x", "y", "xy", "None")
##' @return A highchart plot
##' @author Samuel Wieczorek
##' @examples
##' library("highcharter")
##' hc <- highchart() 
##' hc_chart(hc,type = "line") 
##' hc_add_series(hc,data = c(29, 71, 40))
##' my_hc_ExportMenu(hc,filename='foo')
my_hc_chart <- function(hc,  chartType,zoomType="None"){
  hc %>% 
    hc_chart(type = chartType, 
           zoomType=zoomType,
           resetZoomButton= list(
             position = list(
               align= 'left',
               verticalAlign = 'top')
           ))
}



##' This function retrieves the indices of non-zero elements in sparse matrices
##' of class dgCMatrix from package Matrix. Thi
##' 
##' @title Retrieve the indices of non-zero elements in sparse matrices
##' @param x A sparse matrix of class dgCMatrix
##' @return A two-column matrix
##' @author Samuel Wieczorek
##' @cite This function is largely inspired from the package RINGO
##' @examples
##' library(Matrix)
##' mat <- Matrix(c(0,0,0,0,0,1,0,0,1,1,0,0,0,0,1),nrow=5, byrow=TRUE, sparse=TRUE)
##' res <- nonzero(mat)
nonzero <- function(x){
    ## function to get a two-column matrix containing the indices of the
    ### non-zero elements in a "dgCMatrix" class matrix
    
    stopifnot(inherits(x, "dgCMatrix"))
    if (all(x@p == 0))
        return(matrix(0, nrow=0, ncol=2,
                      dimnames=list(character(0), c("row","col"))))
    res <- cbind(x@i+1, rep(seq(dim(x)[2]), diff(x@p)))
    colnames(res) <- c("row", "col")
    res <- res[x@x != 0, , drop = FALSE]
    return(res)
}



##' This function overloads the brackets to select lines of dataframes in a MSnset file. It takes 
##' into account the slots experimentData / other / OriginOfValues
##' 
##' @title Selects lines of dataframes in a MSnset file
##' @param obj A MSnset object
##' @param lineIndices The indices of lines to be extracted
##' @return A MSnset object 
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' res <- tabOperator(Exp1_R25_pept, c(1:10))
# tabOperator <- function(obj, lineIndices){
#     
#     tmp <- obj[lineIndices]
#     if (!is.null(tmp@experimentData@other$OriginOfValues)){
#         tmp@experimentData@other$OriginOfValues <- tmp@experimentData@other$OriginOfValues[lineIndices,]
#     }
#     
#    
#     return(tmp)
# }
