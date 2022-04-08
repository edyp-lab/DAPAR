()

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
#' Xshared <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], "Protein_group_IDs", FALSE)
#' Xunique <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], "Protein_group_IDs", TRUE)
#' ll.X <- list(matWithSharedPeptides = Xshared, matWithUniquePeptides = Xunique)
#' Exp1_R25_pept <- SetMatAdj(Exp1_R25_pept, ll.X)
#' ll.X <- GetMatAdj(Exp1_R25_pept)
#' 
#' @export
#'
GetMatAdj <- function(obj){
  if (GetTypeofData(obj) != 'peptide')
    warning("This function is only available for a peptide dataset.")
  return(obj@experimentData@other$matAdj)
}


#' @title Record the adjacency matrices in a slot of the 
#' dataset  of class \code{MSnSet}
#' 
#' @param  obj An object (peptides) of class \code{MSnSet}.
#' 
#' @param X A list of two adjacency matrices
#' 
#' @return NA
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' Xshared <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], "Protein_group_IDs", FALSE)
#' Xunique <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], "Protein_group_IDs", TRUE)
#' ll.X <- list(matWithSharedPeptides = Xshared, matWithUniquePeptides = Xunique)
#' Exp1_R25_pept <- SetMatAdj(Exp1_R25_pept, ll.X)
#' 
#' @export
#'
SetMatAdj <- function(obj, X){
  if (GetTypeofData(obj) != 'peptide')
    stop("This function is only available for a peptide dataset. Abort.")
  obj@experimentData@other$matAdj <- X
  return(obj)
}


#' @title Returns the contains of the slot processing  of an object of 
#' class \code{MSnSet}
#' 
#' @param  obj An object (peptides) of class \code{MSnSet}.
#' 
#' @return A list of connected components
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' Xshared <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], "Protein_group_IDs", FALSE)
#' Xunique <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], "Protein_group_IDs", TRUE)
#' ll.X <- list(matWithSharedPeptides = Xshared, matWithUniquePeptides = Xunique)
#' Exp1_R25_pept <- SetMatAdj(Exp1_R25_pept, ll.X)
#' ll1 <- get.pep.prot.cc(GetMatAdj(Exp1_R25_pept)$matWithSharedPeptides)
#' ll2 <- DAPAR::get.pep.prot.cc(GetMatAdj(Exp1_R25_pept)$matWithUniquePeptides)
#' cc <- list(allPep = ll1, onlyUniquePep = ll2)
#' Exp1_R25_pept <- SetCC(Exp1_R25_pept, cc)
#' ll.cc <- GetCC(Exp1_R25_pept)
#' 
#' @export
#'
GetCC <- function(obj){
  if (GetTypeofData(obj) != 'peptide')
    warning("This function is only available for a peptide dataset.")
  return(obj@experimentData@other$CC)
}


#' @title Returns the connected components
#' 
#' @param  obj An object (peptides) of class \code{MSnSet}.
#' 
#' @param cc The connected components list
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' Xshared <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], "Protein_group_IDs", FALSE)
#' Xunique <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], "Protein_group_IDs", TRUE)
#' ll.X <- list(matWithSharedPeptides = Xshared, matWithUniquePeptides = Xunique)
#' Exp1_R25_pept <- SetMatAdj(Exp1_R25_pept, ll.X)
#' ll1 <- get.pep.prot.cc(GetMatAdj(Exp1_R25_pept)$matWithSharedPeptides)
#' ll2 <- DAPAR::get.pep.prot.cc(GetMatAdj(Exp1_R25_pept)$matWithUniquePeptides)
#' cc <- list(allPep = ll1, onlyUniquePep = ll2)
#' Exp1_R25_pept <- SetCC(Exp1_R25_pept, cc)
#' 
#' @export
#'
SetCC <- function(obj, cc){
  if (GetTypeofData(obj) != 'peptide')
    stop("This function is only available for a peptide dataset. Abort.")
  obj@experimentData@other$CC <- cc
  return(obj)
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
#'
getNumberOfEmptyLines <- function(qData){
  n <- sum(apply(is.na(as.matrix(qData)), 1, all))
  return (n)
}


#' @title xxxx
#' 
#' @description
#' xxxx
#' 
#' @param obj xxx
#' 
#' @export
#' 
GetTypeofData <- function(obj){
  if (!is.null(obj))
    obj@experimentData@other$typeOfData
  else
    NULL
}


#' @title xxxx
#' 
#' @description
#' xxxx
#' 
#' @param obj xxx
#' 
#' @export
#' 
GetKeyId <- function(obj)
{
  if (!is.null(obj))
    obj@experimentData@other$keyId
  else
    NULL
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
#' getListNbValuesInLines(Exp1_R25_pept, 'WholeMatrix')
#' 
#' @export
#'
#' 
getListNbValuesInLines <- function(obj, type){
  if(missing(obj))
    stop("'obj' is required")
  else if (is.null(obj))
    stop("'obj' is NULL.")
  
  if(missing(type))
    stop("'type' is required")
  else if (!(type %in% MetacellFilteringScope()))
    stop(paste0("'type' must be one of the following values: ", paste0(MetacellFilteringScope(), collapse= ' ')))
  
  data <- GetMetacell(obj)
  conds <- Biobase::pData(obj)$Condition
  
  ll <- switch(type,
               WholeLine = NULL,
               WholeMatrix= seq(0, ncol(data)),
               AllCond = seq(0, min(unlist(lapply(unique(conds), function(x) length(which(conds == x)))))),
               AtLeastOneCond = seq(0, min(unlist(lapply(unique(conds), function(x) length(which(conds == x))))))
  )
  
  return (ll)
}



#' Returns a list for the two conditions where each slot is a vector of 
#' indices for the samples.
#' 
#' @title Gets the conditions indices.
#' 
#' @param conds A vector of strings containing the column "Condition" of 
#' the \code{Biobase::pData()}.
#' 
#' @param cond1 A vector of Conditions (a slot in the \code{Biobase::pData()} table) for
#' the condition 1.
#' 
#' @param cond2 A vector of Conditions (a slot in the \code{Biobase::pData()} table) for
#' the condition 2.
#' 
#' @return A list with two slots \code{iCond1} and \code{iCond2} containing
#' respectively the indices of samples in the \code{Biobase::pData()} table of the
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




#' 
#' @title Retrieve the indices of non-zero elements in sparse matrices
#' 
#' @description 
#' This function retrieves the indices of non-zero elements in sparse matrices
#' of class dgCMatrix from package Matrix. This function is largely inspired 
#' from the package \code{RINGO}
#' 
#' @param x A sparse matrix of class dgCMatrix
#' 
#' @return A two-column matrix
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' library(Matrix)
#' mat <- Matrix(c(0,0,0,0,0,1,0,0,1,1,0,0,0,0,1),nrow=5, byrow=TRUE, 
#' sparse=TRUE)
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
#' @import highcharter
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
#' @import highcharter
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
