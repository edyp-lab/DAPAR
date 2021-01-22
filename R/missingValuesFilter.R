

#' Returns the percentage of missing values in the quantitative
#' data (\code{exprs()} table of the dataset).
#' 
#' @title Percentage of missing values
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @return A floating number
#' 
#' @author Florence Combes, Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' getPourcentageOfMV(Exp1_R25_pept)
#' 
#' @export
#' 
#' @importFrom Biobase exprs fData pData
#' 
getPourcentageOfMV <- function(obj){
  
  df <- data.frame(Biobase::exprs(obj))
  
  NA.count<-apply(df, 2,
                  function(x) length(which(is.na(data.frame(x))==TRUE)) )
  
  
  pourcentage <- 100 * round(sum(NA.count) /(nrow(df)* ncol(df)), digits=4)
  
  return(pourcentage)
}


#' Returns the number of lines, in a given column, where content matches 
#' the prefix.
#' 
#' @title Number of lines with prefix
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param name The name of a column.
#' 
#' @param prefix A string
#' 
#' @return An integer
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' getNumberOf(Exp1_R25_pept, "Potential_contaminant", "+")
#' 
#' @export
#' 
#' @importFrom Biobase exprs fData pData
#' 
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



#' This function removes lines in the dataset based on numerical conditions.
#' 
#' @title Removes lines in the dataset based on numerical conditions.
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param name The name of the column that correspond to the line to filter
#' 
#' @param value A number 
#' 
#' @param operator A string
#' 
#' @return An list of 2 items :
#' obj : an object of class \code{MSnSet} in which the lines have been deleted
#' deleted : an object of class \code{MSnSet} which contains the deleted lines 
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' NumericalFiltering(Exp1_R25_pept, 'A_Count', '6', '==')
#' 
#' @export
#' 
NumericalFiltering <- function(obj, name=NULL, value=NULL, operator=NULL){
  if ((is.null(name) || (name == ""))) {return(NULL)}
  
  deleted <- NULL
  ind <- NULL
  ind <- NumericalgetIndicesOfLinesToRemove(obj,name, value, operator)
  
  if (!is.null(ind) && (length(ind) > 0)){
    deleted <- obj[ind]
    
    obj <- deleteLinesFromIndices(obj, ind,
                                  paste("\"",
                                        length(ind),
                                        " lines were removed from dataset.\"",
                                        sep="")
    )
    
  }
  
  return(list(obj=obj, deleted=deleted))
}




#' This function returns the indice of the lines to delete, based on a 
#' prefix string
#' 
#' @title Get the indices of the lines to delete, based on a prefix string
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param name The name of the column that correspond to the data to filter
#' 
#' @param value xxxx
#' 
#' @param operator A xxxx
#' 
#' @return A vector of integers.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' NumericalgetIndicesOfLinesToRemove(Exp1_R25_pept, "A_Count", value="6", operator='==')
#' 
#' @export
#' 
#' @importFrom Biobase exprs fData pData
#' 
NumericalgetIndicesOfLinesToRemove <- function(obj, name=NULL, value=NULL, operator=NULL)
{
  if ((value == "") || is.null(value)|| (operator=="") || is.null(operator)) {
    # warning ("No change was made")
    return (NULL)}
  
  data <- Biobase::fData(obj)[,name]
  ind <- which(eval(parse(text=paste0("data", operator, value))))
  
  return(ind)
}



#' Plots a barplot of proportion of contaminants and reverse. Same as the function
#' \code{proportionConRev} but uses the package \code{highcharter}
#' 
#' @title Barplot of proportion of contaminants and reverse
#' 
#' @param nBoth The number of both contaminants and reverse identified in the dataset.
#' 
#' @param nCont The number of contaminants identified in the dataset.
#' 
#' @param nRev The number of reverse entities identified in the dataset.
#' 
#' @param lDataset The total length (number of rows) of the dataset
#' 
#' @return A barplot
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' proportionConRev_HC(10, 20, 100)
#' 
#' @export
#' 
proportionConRev_HC <- function(nBoth = 0, nCont=0, nRev=0, lDataset=0){
  if (is.null(nCont) && is.null(nBoth) && is.null(nRev) && is.null(lDataset)){return(NULL)}
  
  total <- nBoth + nCont + nRev + lDataset
  pctGood <- 100 * round(lDataset/total,  digits=4)
  pctBoth <- 100 * round(nBoth/total,  digits=4)
  pctContaminants <- 100 * round(nCont/total,  digits=4)
  pctReverse <- 100 * round(nRev/total,  digits=4)
  
  counts <- c(lDataset, nCont, nRev, nBoth)
  slices <- c(pctGood, pctContaminants, pctReverse ,pctBoth)
  lbls <- c("Quantitative data", "Contaminants", "Reverse", "Both contaminants & Reverse")
  #pct <- c(pctGood, pctContaminants, pctReverse  ,pctBoth)
  lbls <- paste(lbls, " (", counts, " lines)", sep="")
  
  mydata <- data.frame(test=c(pctGood, pctContaminants, pctReverse ,pctBoth))
  
  highchart() %>%
    my_hc_chart(chartType = "bar") %>%
    hc_yAxis(title = list(text = "Pourcentage")) %>%
    hc_xAxis(categories=lbls) %>%
    hc_legend(enabled = FALSE) %>%
    hc_plotOptions(column = list(
      dataLabels = list(enabled = TRUE),
      stacking = "normal",
      enableMouseTracking = FALSE)
    ) %>%
    hc_add_series(data  = mydata$test,
                  dataLabels = list(enabled = TRUE, format='{point.y}%'),
                  colorByPoint = TRUE) %>%
    my_hc_ExportMenu(filename = "contaminants")
  
  
}




#' This function removes lines in the dataset based on a prefix string.
#' 
#' @title Removes lines in the dataset based on a prefix string.
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param idLine2Delete The name of the column that correspond to the 
#' data to filter
#' 
#' @param prefix A character string that is the prefix to find in the data
#' @return An object of class \code{MSnSet}.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' removeLines(Exp1_R25_pept, "Potential_contaminant")
#' removeLines(Exp1_R25_pept, "Reverse")
#' 
#' @export
#' 
removeLines <- function(obj, idLine2Delete=NULL, prefix=NULL){
  if ((prefix == "") || is.null(prefix)) {
    #warning ("No change was made")
    return (obj)}
  t <- (prefix == substring(Biobase::fData(obj)[,idLine2Delete],1,nchar(prefix)))
  ind <- which( t== TRUE)
  obj <- obj[-ind ]
  
  return(obj)
}



#' This function removes lines in the dataset based on prefix strings (contaminants, reverse or both).
#' 
#' @title Removes lines in the dataset based on a prefix strings (contaminants, reverse or both).
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param idCont2Delete The name of the column that correspond to the 
#' contaminants to filter
#' 
#' @param prefix_Cont A character string that is the prefix for the contaminants to find in the data
#' 
#' @param idRev2Delete The name of the column that correspond to the 
#' reverse data to filter
#' 
#' @param prefix_Rev A character string that is the prefix for the reverse to find in the data
#' 
#' @return An list of 4 items :
#' obj : an object of class \code{MSnSet} in which the lines have been deleted
#' deleted.both : an object of class \code{MSnSet} which contains the deleted lines 
#' corresponding to both contaminants and reverse, 
#' deleted.contaminants : n object of class \code{MSnSet} which contains the deleted lines 
#' corresponding to contaminants, 
#' deleted.reverse : an object of class \code{MSnSet} which contains the deleted lines 
#' corresponding to reverse,
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' StringBasedFiltering(Exp1_R25_pept, 'Potential_contaminant', '+', 'Reverse', '+')
#' 
#' @export
#' 
StringBasedFiltering <- function(obj, 
                                 idCont2Delete=NULL, prefix_Cont=NULL, 
                                 idRev2Delete=NULL, prefix_Rev=NULL){
  
  deleted.both <- deleted.contaminants <- deleted.reverse <- NULL
  
  ##
  ##Search for both
  ##
  if ((!is.null(idCont2Delete) || (idCont2Delete != "")) &&
      (!is.null(idRev2Delete) || (idRev2Delete != ""))) {
    indContaminants <- indReverse <- indBoth <- NULL
    indContaminants <- getIndicesOfLinesToRemove(obj,idCont2Delete,  prefix_Cont)
    indReverse <- getIndicesOfLinesToRemove(obj, idRev2Delete, prefix_Rev)
    indBoth <- intersect(indContaminants, indReverse)
    
    if (!is.null(indBoth) && (length(indBoth) > 0)){
      deleted.both <- obj[indBoth]
      obj <- deleteLinesFromIndices(obj, indBoth,
                                    paste("\"",
                                          length(indBoth),
                                          " both contaminants and reverse were removed from dataset.\"",
                                          sep="")
      )
    }
  }
  
  ##
  ##Search for contaminants
  ##
  if ((!is.null(idCont2Delete) || (idCont2Delete != ""))) {
    indContaminants <- NULL
    indContaminants <- getIndicesOfLinesToRemove(obj,idCont2Delete,  prefix_Cont)
    
    if (!is.null(indContaminants) && (length(indContaminants) > 0)){
      deleted.contaminants <- obj[indContaminants]
      
      obj <- deleteLinesFromIndices(obj, indContaminants,
                                    paste("\"",
                                          length(indContaminants),
                                          " contaminants were removed from dataset.\"",
                                          sep="")
      )
      
    }
  }
  
  
  ##
  ## Search for reverse
  ##
  if ((!is.null(idRev2Delete) || (idRev2Delete != ""))) {
    indReverse <- getIndicesOfLinesToRemove(obj, idRev2Delete, prefix_Rev)
    
    if (!is.null(indReverse)){
      if (length(indReverse) > 0)  {
        deleted.reverse <- obj[indReverse]
        
        obj <- deleteLinesFromIndices(obj, indReverse,
                                      paste("\"",
                                            length(indReverse),
                                            " reverse were removed from dataset.\"",
                                            sep="")
        )
        
      }
    }
  }
  
  
  return(list(obj=obj,
              deleted.both=deleted.both,
              deleted.contaminants=deleted.contaminants,
              deleted.reverse=deleted.reverse))
}






#' This function removes lines in the dataset based on prefix strings.
#' 
#' @title Removes lines in the dataset based on a prefix strings.
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param cname The name of the column that correspond to the line to filter
#' 
#' @param tag A character string that is the prefix for the contaminants to find in the data
#' 
#' @return An list of 4 items :
#' obj : an object of class \code{MSnSet} in which the lines have been deleted
#' deleted : an object of class \code{MSnSet} which contains the deleted lines 
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' StringBasedFiltering2(Exp1_R25_pept, 'Potential_contaminant', '+')
#' 
#' @export
#' 
StringBasedFiltering2 <- function(obj, cname=NULL, tag=NULL){
  
  deleted <- NULL
  
  ##
  ##Search for contaminants
  ##
  if ((!is.null(cname) || (cname != ""))) {
    ind <- NULL
    ind <- getIndicesOfLinesToRemove(obj,cname,  tag)
    
    if (!is.null(ind) && (length(ind) > 0)){
      deleted <- obj[ind]
      
      obj <- deleteLinesFromIndices(obj, ind,
                                    paste("\"",
                                          length(ind),
                                          " contaminants were removed from dataset.\"",
                                          sep="")
      )
      
    }
  }
  
  return(list(obj=obj, deleted=deleted))
}




#' This function returns the indice of the lines to delete, based on a 
#' prefix string
#' 
#' @title Get the indices of the lines to delete, based on a prefix string
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param idLine2Delete The name of the column that correspond to the data 
#' to filter
#' 
#' @param prefix A character string that is the prefix to find in the data
#' 
#' @return A vector of integers.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' getIndicesOfLinesToRemove(Exp1_R25_pept, "Potential_contaminant", prefix="+")
#' 
#' @export
#' 
#' @importFrom Biobase exprs fData pData
#' 
getIndicesOfLinesToRemove <- function(obj, idLine2Delete=NULL, prefix=NULL)
{
  if ((prefix == "") || is.null(prefix)) {
    # warning ("No change was made")
    return (NULL)}
  t <- (prefix == substring(Biobase::fData(obj)[,idLine2Delete],1,nchar(prefix)))
  ind <- which( t== TRUE)
  return(ind)
}


#' Filters the lines of \code{exprs()} table with conditions on the number
#' of missing values.
#' The user chooses the minimum amount of intensities that is acceptable and
#' the filter delete lines that do not respect this condition.
#' The condition may be on the whole line or condition by condition.
#' 
#' The different methods are :
#' "WholeMatrix": given a threshold \code{th}, only the lines that contain
#' at least \code{th} values are kept.
#' "AllCond": given a threshold \code{th}, only the lines which contain
#' at least \code{th} values for each of the conditions are kept.
#' "AtLeastOneCond": given a threshold \code{th}, only the lines that contain
#' at least \code{th} values, and for at least one condition, are kept.
#' 
#' @title Filter lines in the matrix of intensities w.r.t. some criteria
#' 
#' @param obj An object of class \code{MSnSet} containing
#' quantitative data.
#' 
#' @param type Method used to choose the lines to delete.
#' Values are : "None", "WholeMatrix", "AllCond", "AtLeastOneCond"
#' 
#' @param th An integer value of the threshold
#' 
#' @param processText A string to be included in the \code{MSnSet}
#' object for log. 
#' 
#' @return An instance of class \code{MSnSet} that have been filtered.
#' 
#' @author Florence Combes, Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' mvFilter(Exp1_R25_pept, "WholeMatrix", 2)
#' 
#' @export
#' 
mvFilter <- function(obj,
                     type,
                     th,
                     processText=NULL )
{
  #Check parameters
  paramtype <- c("None", "WholeMatrix", "AllCond", "AtLeastOneCond")
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
  
  keepThat <- mvFilterGetIndices(obj, condition = type, threshold=th)
  
  obj <- obj[keepThat]
  
  obj@processingData@processing <-
    c(obj@processingData@processing, processText)
  return(obj)
}



#' Filters the lines of \code{exprs()} table with conditions on the number
#' of missing values.
#' The user chooses the minimum amount of intensities that is acceptable and
#' the filter delete lines that do not respect this condition.
#' The condition may be on the whole line or condition by condition.
#' 
#' The different methods are :
#' "WholeMatrix": given a threshold \code{th}, only the lines that contain
#' at least \code{th} values are kept.
#' "AllCond": given a threshold \code{th}, only the lines which contain
#' at least \code{th} values for each of the conditions are kept.
#' "AtLeastOneCond": given a threshold \code{th}, only the lines that contain
#' at least \code{th} values, and for at least one condition, are kept.
#' 
#' @title Filter lines in the matrix of intensities w.r.t. some criteria
#' 
#' @param obj An object of class \code{MSnSet} containing
#' quantitative data.
#' 
#' @param keepThat A vector of integers which are the indices of lines to 
#' keep.
#' 
#' @param processText A string to be included in the \code{MSnSet}
#' object for log. 
#' 
#' @return An instance of class \code{MSnSet} that have been filtered.
#' 
#' @author Florence Combes, Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' mvFilterFromIndices(Exp1_R25_pept, c(1:10))
#' 
#' @export
#' 
mvFilterFromIndices <- function(obj,
                                keepThat = NULL, 
                                processText="" )
{
  
  if (is.null(keepThat)) {return(obj)}
  obj <- obj[keepThat]
  
  # if (!is.null(obj@experimentData@other$OriginOfValues)){
  #     obj@experimentData@other$OriginOfValues <- obj@experimentData@other$OriginOfValues[keepThat,]
  # }
  obj@processingData@processing <-
    c(obj@processingData@processing, processText)
  
  return(obj)
}


#' Delete the lines of \code{exprs()} table identified by their indice.
#' 
#' @title Delete the lines in the matrix of intensities and the metadata table
#' given their indice.
#' 
#' @param obj An object of class \code{MSnSet} containing
#' quantitative data.
#' 
#' @param deleteThat A vector of integers which are the indices of lines to 
#' delete.
#' 
#' @param processText A string to be included in the \code{MSnSet}
#' object for log. 
#' 
#' @return An instance of class \code{MSnSet} that have been filtered.
#' 
#' @author Florence Combes, Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' deleteLinesFromIndices(Exp1_R25_pept, c(1:10))
#' 
#' @export
#' 
deleteLinesFromIndices <- function(obj,deleteThat=NULL, processText="" )
{
  
  if (is.null(deleteThat)) {return(obj)}
  obj <- obj[-deleteThat]
  
  obj@processingData@processing <-  c(obj@processingData@processing, processText)
  if (grepl("contaminants", processText)){obj@experimentData@other$contaminantsRemoved <- TRUE}
  if (grepl("reverse", processText)){obj@experimentData@other$reverseRemoved <- TRUE }
  return(obj)
}



#' Returns the indices of the lines of \code{exprs()} table to delete w.r.t. 
#' the conditions on the number of missing values.
#' The user chooses the minimum amount of intensities that is acceptable and
#' the filter delete lines that do not respect this condition.
#' The condition may be on the whole line or condition by condition.
#' 
#' The different methods are :
#' "WholeMatrix": given a threshold \code{th}, only the lines that contain
#' at least \code{th} values are kept.
#' "AllCond": given a threshold \code{th}, only the lines which contain
#' at least \code{th} values for each of the conditions are kept.
#' "AtLeastOneCond": given a threshold \code{th}, only the lines that contain
#' at least \code{th} values, and for at least one condition, are kept.
#' 
#' @title Filter lines in the matrix of intensities w.r.t. some criteria
#' 
#' @param obj An object of class \code{MSnSet} containing
#' quantitative data.
#' 
#' @param percent TRUE or FALSE. Default is FALSE..
#' 
#' @param condition Method used to choose the lines to delete.
#' Values are : "None", "EmptyLines", "WholeMatrix", "AllCond", "AtLeastOneCond"
#' 
#' @param threshold An integer value of the threshold if percent is FALSE. Otherwise, a floating
#' number between 0 and 1.
#' 
#' @return An vector of indices that correspond to the lines to keep.
#' 
#' @author Enora Fremy, Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' mvFilterGetIndices(Exp1_R25_pept, condition = "WholeMatrix", threshold=2)
#' mvFilterGetIndices(Exp1_R25_pept, condition = "EmptyLines")
#' mvFilterGetIndices(Exp1_R25_pept, condition = "WholeMatrix", percent=TRUE, threshold=0.5)
#' 
#' @export
#' 
mvFilterGetIndices <- function(obj,
                               percent = FALSE,
                               condition = 'WholeMatrix', 
                               threshold = NULL){
  #Check parameters
  paramtype<-c("None", "EmptyLines", "WholeMatrix", "AllCond", "AtLeastOneCond")
  if (!(condition %in% paramtype)){
    warning("Param `type` is not correct.")
    return (NULL)
  }
  
  if (condition != 'EmptyLines')
    if (!(percent %in% c(T, F))){
      warning("Param `type` is not correct.")
      return (NULL)
    } else {
      if (!isTRUE(percent)){
        paramth <- c(seq(0, nrow(Biobase::pData(obj)), 1))
        if (!(threshold %in% paramth)){
          warning(paste0("Param `threshold` is not correct. It must an integer greater than or equal to 0 and less or equal than ",
                         nrow(Biobase::pData(obj))))
          return (NULL)
        }
      } else {
        if (threshold < 0 || threshold > 1){
          warning("Param `threshold` is not correct. It must be greater than 0 and less than 1.")
          return (NULL)
        }
      }
    }
  
  keepThat <- NULL
  if (is.null(obj@experimentData@other$OriginOfValues)){
    data <- Biobase::exprs(obj)
    warning('The dataset contains no slot OriginOfValues in which to search for indices. The search will
            be proceeded in the intensities tab based on NA values')
  } else {
    data <- dplyr::select(Biobase::fData(obj),
                          obj@experimentData@other$OriginOfValues)
  }
  
  if (condition == "None") {
    keepThat <- seq(1:nrow(data))
  } else if (condition == "EmptyLines") {
    keepThat <- which(apply(!DAPAR::is.MV(data), 1, sum) >= 1)
  } else if (condition == "WholeMatrix") {
    if (isTRUE(percent)) {
      keepThat <- which(rowSums(!DAPAR::is.MV(data))/ncol(data) >= threshold) 
    } else {
      keepThat <- which(apply(!DAPAR::is.MV(data), 1, sum) >= threshold)
    }
  } else if (condition == "AtLeastOneCond" || condition == "AllCond") {
    
    conditions <- unique(Biobase::pData(obj)$Condition)
    nbCond <- length(conditions)
    keepThat <- NULL
    s <- matrix(rep(0, nrow(data)*nbCond),
                nrow=nrow(data),
                ncol=nbCond)
    
    if (isTRUE(percent)) {
      for (c in 1:nbCond) {
        ind <- which(Biobase::pData(obj)$Condition == conditions[c])
        s[,c] <- (rowSums(!DAPAR::is.MV(data[,ind]))/length(ind)) >= threshold
      }
    } else {
      for (c in 1:nbCond) {
        ind <- which(Biobase::pData(obj)$Condition == conditions[c])
        if (length(ind) == 1){
          s[,c] <- (!DAPAR::is.MV(data[,ind]) >= threshold) 
        }
        else {
          s[,c] <- (apply(!DAPAR::is.MV(data[,ind]), 1, sum)) >= threshold
        }
      }
    }
    
    switch(condition,
           AllCond = keepThat <- which(rowSums(s) == nbCond),
           AtLeastOneCond = keepThat <- which(rowSums(s) >= 1)
    )
  }
  
  return(keepThat)
}







#' Returns the indices of the lines of \code{exprs()} table to delete w.r.t. 
#' the conditions on the number of missing values.
#' The user chooses the minimum amount of intensities that is acceptable and
#' the filter delete lines that do not respect this condition.
#' The condition may be on the whole line or condition by condition.
#' 
#' The different methods are :
#' "WholeMatrix": given a threshold \code{th}, only the lines that contain
#' at least \code{th} values are kept.
#' "AllCond": given a threshold \code{th}, only the lines which contain
#' at least \code{th} values for each of the conditions are kept.
#' "AtLeastOneCond": given a threshold \code{th}, only the lines that contain
#' at least \code{th} values, and for at least one condition, are kept.
#' 
#' @title Filter lines in the matrix of intensities w.r.t. some criteria
#' 
#' @param obj An object of class \code{MSnSet} containing
#' quantitative data.
#' 
#' @param type Method used to choose the lines to delete.
#' Values are : "None", "EmptyLines", "WholeMatrix", "AllCond", "AtLeastOneCond"
#' 
#' @param th An integer value of the threshold
#' 
#' @return An vector of indices that correspond to the lines to keep.
#' 
#' @author Florence Combes, Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' mvFilterGetIndices(Exp1_R25_pept, "wholeMatrix", 2)
#' 
#' @export
#' 
mvFilterGetIndices_old <- function(obj,
                                   type, 
                                   th=NULL)
{
  #Check parameters
  paramtype<-c("None", "EmptyLines", "WholeMatrix", "AllCond", "AtLeastOneCond")
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
  if (is.null(obj@experimentData@other$OriginOfValues)){
    data <- Biobase::exprs(obj)
    warning('There is no slot OriginOfValues in which to search for indices. The search will
            be proceeded in the intensities tab based on NA values')
  } else {
    data <- dplyr::select(Biobase::fData(obj),obj@experimentData@other$OriginOfValues)
  }
  
  if (type == "None"){
    keepThat <- seq(1:nrow(data))
  } else if (type == "EmptyLines"){
    keepThat <- which(apply(!is.MV(data), 1, sum) >= 1)
  } else if (type == "WholeMatrix"){
    keepThat <- which(apply(!is.MV(data), 1, sum) >= th)
  } else if (type == "AtLeastOneCond" || type == "AllCond"){
    
    conditions <- unique(Biobase::pData(obj)$Condition)
    nbCond <- length(conditions)
    keepThat <- NULL
    s <- matrix(rep(0, nrow(data)*nbCond),nrow=nrow(data),
                ncol=nbCond)
    
    for (c in 1:nbCond){
      ind <- which(Biobase::pData(obj)$Condition == conditions[c])
      if (length(ind) == 1){
        s[,c] <- (!is.MV(data[,ind]) >= th)}
      else {
        s[,c] <- (apply(!is.MV(data[,ind]), 1, sum) >= th)
      }
    }
    
    
    if (type == "AllCond") {
      keepThat <- which(rowSums(s) == nbCond)
    }
    else if (type == "AtLeastOneCond") {
      keepThat <- which(rowSums(s) >= 1)
    }
  }
  return(keepThat)
}


#' Filter peptides/features after their identification method
#'
#' @description Identification xxx (output MaxQuant): By MS/MS, By matching, Missing Value
#' #################################################################
#' Remove lines in the data according to the proportion of missing
#' values. This proportion is calculated differently depending on whether we
#' want a certain proportion of missing values (NA) to remain on:
#' * the entire matrix, regardless of the conditions: the rows containing a
#' proportion of NA equal or below the threshold will be kept.
#' * all the conditions: the lines for which all the conditions have a NA
#' proportion equal to or less than the fixed proportion will be kept.
#' * at least one condition: the lines for which at least one condition is
#' equal to or less than the fixed proportion of NA will be kept.
#' #############################################################
#' 
#' @param obj  An object of class \code{MSnSet} containing quantitative data
#' and phenotype data.
#' 
#' @param mode character string. Four possibilities corresponding to the
#' description above: "None", "WholeMatrix", "AllCond" and "AtLeastOneCond".
#'
#' @param th integer between 0 and either nb samples if "WholeMatrix" or min(nb samples per condition) if "AllCond" or "AtLeastOneCond".
#'
#' @return the object given as input minus the lines without enough identification
#'
#' @author Enora Fremy
#'
#' @examples
#' utils::data(Exp1_R25_prot, package='DAPARdata')
#' filtered <- filt_Marianne(obj = Exp1_R25_prot, th=2,  mode = "AtLeastOneCond")
#'
#' @export
#'
#' @importFrom Biobase pData fData
#'

filt_Marianne <- function(obj) {
  
  # for (i in 81:86){
  #   print(as.data.frame(table((Biobase::fData(obj))[i])))
  # }
  fdata <- Biobase::fData(obj)[81:86]
  
  # calculer identifcation en fonction PSMs et abondance ? Normalement non parce que filtrage sur msnset
  # compter nb MS/MS pour chaque features (et selon cas dans chaque condition)
  # somme de chaque, somme de chaque dans chaque condition plus le maximum de chaque valeur
  # creer vecteurs contenant ces valeurs ? ou nouvelles colonnes ?
  
  # vector containing sum of 'By MS/MS'
  # "None", "WholeMatrix", "AllCond" and "AtLeastOneCond".
  new_vec_sum <- vector()
  for (i in 1:nrow(fdata)) {
    print(i)
    if (mode=='AllCond'|mode='AtLeastOneCond'){
      condition <- Biobase::pData(obj)['Condition']
      condition <- condition[,1]
      # sep by condition
      # count nb of replicates by condition
      # condition_nouvelle_boucle
      indices_replicates <- list()
      for (i in 1:length(unique(condition))) {
        #find how increment name in list in loop > assign?
        indices_replicates <- list(paste0('cond',i) = which(condition==unique(condition)[i]))
      }
      ###############################################################################
      intermediR <- as.data.frame(table(t(fdata[i,indices_replicates])))
      
      new_vec_max <- c(new_vec_max,intermediR[intermediR$Var1=='By MS/MS',2])
      ###############################################################################
    }
    else if (mode=='WholeMatrix' {
      intermediR <- as.data.frame(table(t(fdata[i,])))
      new_vec_sum <- c(new_vec_sum,intermediR[intermediR$Var1=='By MS/MS',2])
    }
  }
}

#' #' Filter missing values by proportion
#' #'
#' #' @description Remove lines in the data according to the proportion of missing
#' #' values. This proportion is calculated differently depending on whether we
#' #' want a certain proportion of missing values (NA) to remain on:
#' #' * the entire matrix, regardless of the conditions: the rows containing a
#' #' proportion of NA equal or below the threshold will be kept.
#' #' * all the conditions: the lines for which all the conditions have a NA
#' #' proportion equal to or less than the fixed proportion will be kept.
#' #' * at least one condition: the lines for which at least one condition is
#' #' equal to or less than the fixed proportion of NA will be kept.
#' #'
#' #' @param obj  An object of class \code{MSnSet} containing quantitative data
#' #' and phenotype data.
#' #' 
#' #' @param intensities_proportion float between 0 and 1 corresponding to the proportion
#' #' of intensities to keep in the lines.
#' #' 
#' #' @param mode character string. Four possibilities corresponding to the
#' #' description above: "None", "WholeMatrix", "AllCond" and "AtLeastOneCond".
#' #' 
#' #' @return the object given as input but with the lines not respecting the
#' #' proportion of NA requested in less.
#' #' 
#' #' @author Hélène Borges
#' #' 
#' #' @examples
#' #' utils::data(Exp1_R25_prot, package='DAPARdata')
#' #' filtered <- filterByProportion(obj = Exp1_R25_prot, intensities_proportion = 0.8, mode = "atLeastOneCond")
#' #' 
#' #' @export
#' #' 
#' #' @importFrom stringr str_glue
#' #' @importFrom Biobase exprs pData fData
#' #' @import dplyr
#' #' @importFrom tidyr pivot_longer
#' #' @importFrom methods is
#' #' 
#' filterByProportion <- function(obj, intensities_proportion, mode = "None"){
#'   # check if mode is valid
#'   if(!(mode %in% c("None","WholeMatrix", "AllCond", "AtLeastOneCond"))){
#'     stop(stringr::str_glue("Wrong mode: {mode} is not a valid string.
#'                      Please choose between 'None', 'WholeMatrix', 'AllCond' or 'AtLeastOneCond'.",
#'                            call. =FALSE))
#'   }
#'   # check if intensities_proportion is valid
#'   if(!methods::is(intensities_proportion, "numeric" )){
#'     stop(stringr::str_glue("Wrong parameter: intensities_proportion needs to be numeric"))
#'   }else if(!dplyr::between(intensities_proportion,0,1)){
#'     stop(stringr::str_glue("Wrong parameter: intensities_proportion must be between 0 and 1"))
#'   }
#'   
#'   print(stringr::str_glue("chosen proportion of intensities to be present: {intensities_proportion}"))
#'   print(stringr::str_glue("chosen mode: {mode}"))
#'   intensities <- Biobase::exprs(obj)
#'   sTab <- Biobase::pData(obj)
#'   sTab$Condition <- as.factor(sTab$Condition)
#'   
#'   intensities_t <- as.data.frame(t(intensities))
#'   intensities_t <- dplyr::bind_cols(intensities_t,
#'                                     condition = sTab$Condition,
#'                                     sample = rownames(intensities_t))
#'   tbl_intensities <- dplyr::as_tibble(intensities_t, rownames = NA)
#'   longer_intensities <- tbl_intensities %>%
#'     tidyr::pivot_longer(-c(condition,sample), names_to = "feature", values_to = "intensity")
#'   # group_by does not keep the initial order when it is not a factor so to keep
#'   # the protein order, we cheat by transforming feature into a factor.
#'   longer_intensities$feature <- factor(longer_intensities$feature,
#'                                        levels = unique(longer_intensities$feature))
#'   if(mode == "None"){
#'     to_keep <- obj
#'   }else if(mode == "WholeMatrix"){
#'     nb_samples <- ncol(intensities)
#'     threshold <- ceiling(nb_samples*intensities_proportion)
#'     print(stringr::str_glue("missing value threshold {threshold}"))
#'     # for each feature (protein / peptide) we count the number of intensities present
#'     feat_grp <- longer_intensities %>%
#'       dplyr::group_by(feature) %>%
#'       dplyr::summarise(non_na = sum(!is.na(intensity)))
#'     to_keep <- obj[which(feat_grp$non_na >= threshold),]
#'     
#'   }else if(mode == "AllCond" || mode == "AtLeastOneCond"){
#'     workforces <- longer_intensities %>%
#'       dplyr::group_by(feature, condition) %>%
#'       dplyr::count(condition)
#'     # the number of samples per condition
#'     workforces <- workforces$n[seq_len(length(levels(sTab$Condition)))]
#'     
#'     # for each condition of each feature, we count the number of intensities present
#'     feat_grp <- longer_intensities %>%
#'       dplyr::group_by(feature, condition) %>%
#'       dplyr::summarise(non_na = sum(!is.na(intensity)))
#'     # the threshold for each condition
#'     thresholds <- ceiling(workforces*intensities_proportion)
#'     print(stringr::str_glue("for condition {unique(levels(longer_intensities$condition))} number of samples is {workforces}, so missing value threshold is {thresholds} "))
#'     # For each feature, each condition is compared with its respective
#'     # threshold, we put 0 if the protein has a number of intensities lower than
#'     # the threshold of the corresponding condition, and 1 otherwise
#'     check_th <- feat_grp %>%
#'       dplyr::group_by(feature) %>%
#'       dplyr::mutate(non_na = dplyr::case_when(
#'         non_na < thresholds ~ 0,
#'         TRUE ~ 1
#'       )) %>%
#'       dplyr::ungroup()
#'     # if it is allCond then we must find the features for which all the conditions
#'     # respect the threshold
#'     if(mode == "AllCond"){
#'       all_cond_ok <- check_th %>%
#'         dplyr::group_by(feature) %>%
#'         dplyr::filter(all(non_na ==1)) %>%
#'         dplyr::ungroup() %>%
#'         as.data.frame()
#'       all_cond_ok$feature <- as.character(all_cond_ok$feature)
#'       to_keep <- obj[which(rownames(obj) %in% all_cond_ok$feature),]
#'     }else if(mode == "AtLeastOneCond"){
#'       # if it is atLeastOneCond then we must find the features for which at
#'       # least one condition that respects the threshold
#'       any_cond_ok <- check_th %>%
#'         dplyr::group_by(feature) %>%
#'         dplyr::filter(any(non_na ==1)) %>%
#'         dplyr::ungroup() %>%
#'         as.data.frame()
#'       any_cond_ok$feature <- as.character(any_cond_ok$feature)
#'       to_keep <- obj[which(rownames(obj) %in% any_cond_ok$feature),]
#'     }
#'   }
#'   print(stringr::str_glue("There were initially {nrow(intensities)} features.
#'                  After filtering out the missing values, {nrow(exprs(to_keep))} remain."))
#'   return(to_keep)
#' }
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #' Filter missing values by proportion
#' #'
#' #' @description Remove lines in the data according to the proportion of missing
#' #' values. This proportion is calculated differently depending on whether we
#' #' want a certain proportion of missing values (NA) to remain on:
#' #' * the entire matrix, regardless of the conditions: the rows containing a
#' #' proportion of NA equal or below the threshold will be kept.
#' #' * all the conditions: the lines for which all the conditions have a NA
#' #' proportion equal to or less than the fixed proportion will be kept.
#' #' * at least one condition: the lines for which at least one condition is
#' #' equal to or less than the fixed proportion of NA will be kept.
#' #'
#' #' @param obj  An object of class \code{MSnSet} containing quantitative data
#' #' and phenotype data.
#' #' 
#' #' @param intensities_proportion float between 0 and 1 corresponding to the proportion
#' #' of intensities to keep in the lines.
#' #' 
#' #' @param mode character string. Four possibilities corresponding to the
#' #' description above: "None", "WholeMatrix", "AllCond" and "AtLeastOneCond".
#' #' 
#' #' @return the object given as input but with the lines not respecting the
#' #' proportion of NA requested in less.
#' #' 
#' #' @author Hélène Borges, Samuel Wieczorek
#' #' 
#' #' @examples
#' #' utils::data(Exp1_R25_prot, package='DAPARdata')
#' #' filtered <- filterByProportion(obj = Exp1_R25_prot, intensities_proportion = 0.8, mode = "AtLeastOneCond")
#' #' 
#' #' @export
#' #' 
#' #' @importFrom stringr str_glue
#' #' @importFrom Biobase exprs pData fData
#' #' @import dplyr
#' #' @importFrom tidyr pivot_longer
#' #' @importFrom methods is
#' #' 
#' GetIndices_Filter_By_Proportion <- function(obj, intensities_proportion, mode = "None"){
#'   # check if mode is valid
#'   if(!(mode %in% c("None","WholeMatrix", "AllCond", "AtLeastOneCond"))){
#'     stop(stringr::str_glue("Wrong mode: {mode} is not a valid string.
#'                      Please choose between 'None', 'WholeMatrix', 'AllCond' or 'AtLeastOneCond'.",
#'                            call. =FALSE))
#'   }
#'   # check if intensities_proportion is valid
#'   if(!methods::is(intensities_proportion, "numeric" )){
#'     stop(stringr::str_glue("Wrong parameter: intensities_proportion needs to be numeric"))
#'   }else if(!dplyr::between(intensities_proportion,0,1)){
#'     stop(stringr::str_glue("Wrong parameter: intensities_proportion must be between 0 and 1"))
#'   }
#'   
#'   
#'   keepThat <- NULL
#'   if (is.null(obj@experimentData@other$OriginOfValues)){
#'     data <- Biobase::exprs(obj)
#'   } else {
#'     data <- dplyr::select(Biobase::fData(obj),obj@experimentData@other$OriginOfValues)
#'   }
#'   
#'   
#'   print(stringr::str_glue("chosen proportion of intensities to be present: {intensities_proportion}"))
#'   print(stringr::str_glue("chosen mode: {mode}"))
#'   data <- Biobase::exprs(obj)
#'   sTab <- Biobase::pData(obj)
#'   sTab$Condition <- as.factor(sTab$Condition)
#'   
#'   data_t <- as.data.frame(t(data))
#'   data_t <- dplyr::bind_cols(data_t,
#'                              condition = sTab$Condition,
#'                              sample = rownames(data_t))
#'   tbl_data <- dplyr::as_tibble(data_t, rownames = NA)
#'   longer_data <- tbl_data %>%
#'     tidyr::pivot_longer(-c(condition,sample), names_to = "feature", values_to = "intensity")
#'   # group_by does not keep the initial order when it is not a factor so to keep
#'   # the protein order, we cheat by transforming feature into a factor.
#'   longer_data$feature <- factor(longer_data$feature,
#'                                 levels = unique(longer_data$feature))
#'   if(mode == "None"){
#'     indices_to_keep <- 1:nrow(obj)
#'   }else if(mode == "WholeMatrix"){
#'     nb_samples <- ncol(data)
#'     threshold <- ceiling(nb_samples*data_proportion)
#'     print(stringr::str_glue("missing value threshold {threshold}"))
#'     # for each feature (protein / peptide) we count the number of data present
#'     feat_grp <- longer_data %>%
#'       dplyr::group_by(feature) %>%
#'       dplyr::summarise(non_na = sum(!is.na(intensity)))
#'     indices_to_keep <- which(feat_grp$non_na >= threshold)
#'     
#'   }else if(mode == "AllCond" || mode == "AtLeastOneCond"){
#'     workforces <- longer_data %>%
#'       dplyr::group_by(feature, condition) %>%
#'       dplyr::count(condition)
#'     # the number of samples per condition
#'     workforces <- workforces$n[seq_len(length(levels(sTab$Condition)))]
#'     
#'     # for each condition of each feature, we count the number of data present
#'     feat_grp <- longer_data %>%
#'       dplyr::group_by(feature, condition) %>%
#'       dplyr::summarise(non_na = sum(!is.na(intensity)))
#'     # the threshold for each condition
#'     thresholds <- ceiling(workforces*data_proportion)
#'     print(stringr::str_glue("for condition {unique(levels(longer_data$condition))} number of samples is {workforces}, so missing value threshold is {thresholds} "))
#'     # For each feature, each condition is compared with its respective
#'     # threshold, we put 0 if the protein has a number of data lower than
#'     # the threshold of the corresponding condition, and 1 otherwise
#'     check_th <- feat_grp %>%
#'       dplyr::group_by(feature) %>%
#'       dplyr::mutate(non_na = dplyr::case_when(
#'         non_na < thresholds ~ 0,
#'         TRUE ~ 1
#'       )) %>%
#'       dplyr::ungroup()
#'     # if it is allCond then we must find the features for which all the conditions
#'     # respect the threshold
#'     if(mode == "AllCond"){
#'       all_cond_ok <- check_th %>%
#'         dplyr::group_by(feature) %>%
#'         dplyr::filter(all(non_na ==1)) %>%
#'         dplyr::ungroup() %>%
#'         as.data.frame()
#'       all_cond_ok$feature <- as.character(all_cond_ok$feature)
#'       indices_to_keep <- which(rownames(obj) %in% all_cond_ok$feature)
#'     }else if(mode == "AtLeastOneCond"){
#'       # if it is atLeastOneCond then we must find the features for which at
#'       # least one condition that respects the threshold
#'       any_cond_ok <- check_th %>%
#'         dplyr::group_by(feature) %>%
#'         dplyr::filter(any(non_na ==1)) %>%
#'         dplyr::ungroup() %>%
#'         as.data.frame()
#'       any_cond_ok$feature <- as.character(any_cond_ok$feature)
#'       indices_to_keep <- which(rownames(obj) %in% any_cond_ok$feature)
#'     }
#'   }
#'   print(stringr::str_glue("There were initially {nrow(data)} features.
#'                  After filtering out the missing values, {nrow(exprs(to_keep))} remain."))
#'   return(indices_to_keep)
#' }
