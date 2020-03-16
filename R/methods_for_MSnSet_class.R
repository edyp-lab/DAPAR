

setGeneric("parentProtId", function(obj,...) standardGeneric("parentProtId"))

#' @export
setMethod("parentProtId", "MSnSet", function(obj, withDimnames=TRUE) {
  out <- obj@experimentData@other$parentProtId
  out
})


setGeneric("parentProtId<-", function(obj, value) standardGeneric("parentProtId<-"))

#' @export
setReplaceMethod("parentProtId", "MSnSet", function(obj, value) {
  obj@experimentData@other$parentProtId <- value
  obj
})




setGeneric("DAPARVersion", function(obj,...) standardGeneric("DAPARVersion"))

#' Function to get the version of DAPAR used to create the current object
#' @title Get the version of DAPAR used.
#' @description qkljsdh qsodhqs qsidsq.
#' @param obj xxxx
#' @return xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("DAPARVersion", "MSnSet", function(obj) {
  out <- obj@experimentData@other$DAPAR_Version
  out
})

setGeneric("DAPARVersion<-", function(obj, value) standardGeneric("DAPARVersion<-"))

#' Function to set the version of DAPAR used to create the current object
#' @title Set the version of DAPAR used.
#' @description qsdh qiosdhqs qisdq.
#' @param obj xxxx
#' @param value xxx
#' @return xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export

setReplaceMethod("DAPARVersion", "MSnSet", function(obj, value) {
  obj@experimentData@other$DAPAR_Version <- value
  obj
})



setGeneric("ProstarVersion", function(obj,...) standardGeneric("ProstarVersion"))

#' Function to get the version of Prostar used to create the current object
#' @title Get the version of Prostar used.
#' @description sfgsdff sjkfhksqdsq .
#' @param obj xxxx
#' @param value xxx
#' @return xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("ProstarVersion", "MSnSet", function(obj, withDimnames=TRUE) {
  out <- obj@experimentData@other$Prostar_Version
  out
})


setGeneric("ProstarVersion<-", function(obj, value) standardGeneric("ProstarVersion<-"))

#' Function to set the version of Prostar used to create the current object
#' @title Set the version of Prostar used.
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export

setReplaceMethod("ProstarVersion", "MSnSet", function(obj, value) {
  obj@experimentData@other$Prostar_Version <- value
  obj
})





setGeneric("Params", function(obj,...) standardGeneric("Params"))

#' @export
setMethod("Params", "MSnSet", function(obj, withDimnames=TRUE) {
  out <- obj@experimentData@other$Params
  out
})


setGeneric("Params<-", function(obj, value) standardGeneric("Params<-"))

#' @export
setReplaceMethod("Params", "MSnSet", function(obj, value) {
  obj@experimentData@other$Params <- value
  obj
})







setGeneric("OriginOfValues", function(obj,...) standardGeneric("OriginOfValues"))

#' @export
setMethod("OriginOfValues", "MSnSet", function(obj, withDimnames=TRUE) {
  out <- obj@experimentData@other$OriginOfValues
  out
})
setGeneric("OriginOfValues<-", function(obj, value) standardGeneric("OriginOfValues<-"))

#' @export
setReplaceMethod("OriginOfValues", "MSnSet", function(obj, value) {
  obj@experimentData@other$OriginOfValues <- value
  obj
})





setGeneric("typeOfData", function(obj,...) standardGeneric("typeOfData"))

#' @export
setMethod("typeOfData", "MSnSet", function(obj, withDimnames=TRUE) {
  out <- obj@experimentData@other$typeOfData
  out
})

setGeneric("typeOfData<-", function(obj, value) standardGeneric("typeOfData<-"))

#' @export
setReplaceMethod("typeOfData", "MSnSet", function(obj, value) {
  obj@experimentData@other$typeOfData <- value
  obj
})


setGeneric("RawPValues", function(obj,...) standardGeneric("RawPValues"))

#' @export
setMethod("RawPValues", "MSnSet", function(obj, withDimnames=TRUE) {
  out <- obj@experimentData@other$RawPValues
  out
})

setGeneric("RawPValues<-", function(obj, value) standardGeneric("RawPValues<-"))

#' @export
setReplaceMethod("RawPValues", "MSnSet", function(obj, value) {
  obj@experimentData@other$RawPValues <- value
  obj
})



setGeneric("properties", function(obj,...) standardGeneric("properties"))

#' @export
setMethod("properties", "MSnSet", function(obj, withDimnames=TRUE) {
  out <- obj@experimentData@other
  out
})

