

#' @title Extends a base-palette of the package RColorBrewer to n colors.
#' 
#' @description The colors in the returned palette are always in the same order
#' 
#' @param n The number of desired colors in the palette
#' 
#' @param base The name of the palette of the package RColorBrewer from which the extended palette is built.
#' Default value is 'Set1'.
#' 
#' @return A vector composed of n color code.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' ExtendPalette(12)
#' nPalette <- 15
#' par(mfrow=c(nPalette,1))
#' par(mar=c(0.5, 4.5, 0.5, 0.5))
#' for (i in 8:nPalette){
#'   palette <- ExtendPalette(n=i, base = 'Set1')
#'   barplot(1:length(palette), col=palette)
#'   print(palette)
#' }
#' 
#' @export
#' 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' 
#' 
ExtendPalette <- function(n = 0, base = "Set1"){
  palette <- NULL
  nMaxColors <- RColorBrewer::brewer.pal.info[base, 'maxcolors']
  limit <- nMaxColors*(nMaxColors-1)/2
  if (n > limit){
    stop('Number of colors exceed limit of ', limit, ' colors per palette.')
  }
  
  if(n > nMaxColors){
    palette <- RColorBrewer::brewer.pal(nMaxColors, base)
    allComb <- combn(palette, 2)
    
    for (i in 1:(n-nMaxColors)){
      palette <- c(palette, grDevices::colorRampPalette(allComb[,i])(3)[2])
    }
    
  } else {
    palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(nMaxColors, base))(n)
  }
  palette
}






#' @title Builds a complete color palette for the conditions given in argument
#' 
#' @description xxxx
#' 
#' @param conds The extended vector of samples conditions
#' 
#' @param base_palette The basic color (HEX code) used to build the complete palette. This vector have the same length as unique(conds).
#' Default base palette is 'Dark2' from the package RColorBrewer.
#' 
#' @return A vector composed of HEX color code for the conditions
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' conditions <- Biobase::pData(Exp1_R25_pept)$Condition
#' BuildPalette(conditions, c('AAAAAA', 'BBBBBB'))
#' 
#' @export
#' 
#' @importFrom RColorBrewer brewer.pal
#' 
#' 
BuildPalette <- function(conds, base_palette=NULL){
  
  if (!is.null(base_palette) && length(unique(conds)) != length(base_palette)){
    stop('The length of `conds` must be equal to the length of `base_palette`.')
  }
  
  palette <- NULL
  if (is.null(base_palette)){
    base_palette <- grDevices::colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(conds)))
  }
  
  for (i in 1:length(conds)){
    palette[i] <- base_palette[which(conds[i] == unique(conds))]
  }
  return(palette)
  
}      



#' Returns the contents of the slot processing of an object of class \code{MSnSet}
#' 
#' @title Returns the contains of the slot processing  of an object of 
#' class \code{MSnSet}
#' 
#' @param  obj An object (peptides) of class \code{MSnSet}.
#' 
#' @return The slot processing of obj@processingData
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' getProcessingInfo(Exp1_R25_pept)
#' 
#' @export
#'
getProcessingInfo <- function(obj){
return(obj@processingData@processing)
}

#' Returns the number of empty lines in a matrix.
#' 
#' @title Returns the number of empty lines in the data
#' 
#' @param qData A matrix corresponding to the quantitative data.
#' 
#' @return An integer
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' getNumberOfEmptyLines(qData)
#' 
#' @export
#'
#' @importFrom Biobase pData exprs fData
#'
getNumberOfEmptyLines <- function(qData){
n <- sum(apply(is.na(as.matrix(qData)), 1, all))
return (n)
}

#' Similar to the function \code{is.na} but focused on the equality with the paramter 'type'.
#' 
#' @title Similar to the function \code{is.na} but focused on the equality with the paramter 'type'.
#' 
#' @param data A data.frame
#' 
#' @param type The value to search in the dataframe
#' 
#' @return A boolean dataframe 
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept
#' data <- Biobase::fData(obj)[,obj@experimentData@other$OriginOfValues]
#' is.OfType(data, "MEC")
#' 
#' @export
#'
is.OfType <- function(data, type){
  return (type == data)
}





#' Similar to the function \code{is.na} but focused on the equality with the missing 
#' values in the dataset (type 'POV' and 'MEC')
#' 
#' @title Similar to the function \code{is.na} but focused on the equality with the missing 
#' values in the dataset (type 'POV' and 'MEC')
#' 
#' @param data A data.frame
#' 
#' @return A boolean dataframe 
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept
#' data <- Biobase::fData(obj)[,obj@experimentData@other$OriginOfValues]
#' is.MV(data)
#' 
#' @export
#'
is.MV <- function(data){
  #MV=is.OfType(data, "MV")
  POV=is.OfType(data, "POV")
  MEC=is.OfType(data, "MEC")
  isNA = is.na(data)
  df <- POV | MEC | isNA
  return (df)
}

#' Returns the possible number of values in lines in a matrix.
#' 
#' @title Returns the possible number of values in lines in the data
#' 
#' @param obj An object of class \code{MSnSet}
#' 
#' @param type xxxxxxx
#' 
#' @return An integer
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' getListNbValuesInLines(Exp1_R25_pept)
#' 
#' @export
#'
#' @importFrom Biobase pData exprs fData
#' 
getListNbValuesInLines <- function(obj, type="WholeMatrix"){
  if (is.null(obj)){return(NULL)}
  
  if(is.null(obj@experimentData@other$OriginOfValues)){
    ll <- seq(0,ncol(obj))
  }
  data <- Biobase::fData(obj)[,obj@experimentData@other$OriginOfValues]
  switch(type,
         WholeMatrix= {
           ll <- unique(ncol(data) - apply(is.MV(data), 1, sum))
           },
         AllCond = {
                    tmp <- NULL
                    for (cond in unique(Biobase::pData(obj)$Condition)){
                     tmp <- c(tmp, length(which(Biobase::pData(obj)$Condition == cond)))
                  }
                  ll <- seq(0,min(tmp))
                  },
         AtLeastOneCond = {
                   tmp <- NULL
                  for (cond in unique(Biobase::pData(obj)$Condition)){
                       tmp <- c(tmp, length(which(Biobase::pData(obj)$Condition == cond)))
                   }
                   ll <- seq(0,max(tmp))
                    }
         )
  
  return (sort(ll))
}



#' Returns a list for the two conditions where each slot is a vector of 
#' indices for the samples.
#' 
#' @title Gets the conditions indices.
#' 
#' @param conds A vector of strings containing the column "Condition" of 
#' the \code{pData()}.
#' 
#' @param cond1 A vector of Conditions (a slot in the \code{pData()} table) for
#' the condition 1.
#' 
#' @param cond2 A vector of Conditions (a slot in the \code{pData()} table) for
#' the condition 2.
#' 
#' @return A list with two slots \code{iCond1} and \code{iCond2} containing
#' respectively the indices of samples in the \code{pData()} table of the
#' dataset. 
#' 
#' @author Florence Combes, Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' conds <- Biobase::pData(Exp1_R25_pept)[,"Condition"]
#' getIndicesConditions(conds, "25fmol", "10fmol")
#' 
#' @export
#' 
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



#' Customise the contextual menu of highcharts plots.
#' 
#' @title Customised contextual menu of highcharts plots
#' 
#' @param hc A highcharter object
#' 
#' @param filename The filename under which the plot has to be saved
#' 
#' @return A contextual menu for highcharts plots
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library("highcharter")
#' hc <- highchart() 
#' hc_chart(hc,type = "line") 
#' hc_add_series(hc,data = c(29, 71, 40))
#' my_hc_ExportMenu(hc,filename='foo')
#' 
#' @export
#' 
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



#' Customise the resetZoomButton of highcharts plots.
#' 
#' @title Customised resetZoomButton of highcharts plots
#' 
#' @param hc A highcharter object
#' 
#' @param chartType The type of the plot
#' 
#' @param zoomType The type of the zoom (one of "x", "y", "xy", "None")
#' 
#' @return A highchart plot
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library("highcharter")
#' hc <- highchart() 
#' hc_chart(hc,type = "line") 
#' hc_add_series(hc,data = c(29, 71, 40))
#' my_hc_ExportMenu(hc,filename='foo')
#' 
#' @export
#' 
my_hc_chart <- function(hc,  chartType,zoomType="None"){
  hc %>% 
    hc_chart(type = chartType, 
           zoomType=zoomType,
           showAxes = TRUE,
           resetZoomButton= list(
             position = list(
               align= 'left',
               verticalAlign = 'top')
           ))
}



#' This function retrieves the indices of non-zero elements in sparse matrices
#' of class dgCMatrix from package Matrix. This function is largely inspired from 
#' the package \code{RINGO}
#' 
#' @title Retrieve the indices of non-zero elements in sparse matrices
#' 
#' @param x A sparse matrix of class dgCMatrix
#' 
#' @return A two-column matrix
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(Matrix)
#' mat <- Matrix(c(0,0,0,0,0,1,0,0,1,1,0,0,0,0,1),nrow=5, byrow=TRUE, sparse=TRUE)
#' res <- nonzero(mat)
#' 
#' @export
#' 
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




#' @title Customised contextual menu of highcharts plots
#' 
#' @param hc A highcharter object
#' 
#' @param filename The filename under which the plot has to be saved
#' 
#' @return A contextual menu for highcharts plots
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library("highcharter")
#' hc <- highchart() 
#' hc_chart(hc,type = "line") 
#' hc_add_series(hc,data = c(29, 71, 40))
#' dapar_hc_ExportMenu(hc,filename='foo')
#' 
#' @export
#' 
#' @importFrom highcharter hc_exporting
#' 
dapar_hc_ExportMenu <- function(hc, filename){
  hc_exporting(hc, enabled=TRUE,
               filename = filename,
               buttons= list(
                 contextButton= list(
                   menuItems= list('downloadPNG', 'downloadSVG','downloadPDF')
                 )
               )
  )
}






#' @title Customised resetZoomButton of highcharts plots
#' 
#' @param hc A highcharter object
#' 
#' @param chartType The type of the plot
#' 
#' @param zoomType The type of the zoom (one of "x", "y", "xy", "None")
#' 
#' @param width xxx
#' 
#' @param height xxx
#' 
#' @return A highchart plot
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library("highcharter")
#' hc <- highchart() 
#' hc <- dapar_hc_chart(hc, chartType='line', zoomType='x')
#' hc_add_series(hc, data = c(29, 71, 40))
#' 
#' @export
#' 
#' @importFrom highcharter hc_chart
#' 
dapar_hc_chart <- function(hc,  chartType, zoomType="None", width=0, height=0){
  hc %>% 
    hc_chart(type = chartType, 
             zoomType=zoomType,
             showAxes = TRUE,
             width = width,
             height = height,
             resetZoomButton= list(
               position = list(
                 align= 'left',
                 verticalAlign = 'top')
             ))

}
