


#' @title xxx
#' @description xxx
#' @param obj xxx
#' @export
Check_Dataset_Validity <- function(obj){
  if (is.null(obj))
    return(NULL)
  
  is.valid <- TRUE
  
  test <- Check_NbValues_In_Columns(exprs(obj))
  is.valid <- length(which(test == 0)) == 0
  is.valid <- length(which(test == 1)) == 0
  return(is.valid)
}



#' @title xxx
#' @description xxx
#' @param qdata xxx
#' @export
Check_NbValues_In_Columns <- function(qdata){
  
  unlist(lapply(1:ncol(qdata), function(x)
    length(which(!is.na(qdata[,x])))))
}