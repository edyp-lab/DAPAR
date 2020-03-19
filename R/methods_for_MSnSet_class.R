#' @title qfsqfds
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setGeneric("keyId", function(obj) standardGeneric("keyId"))


#' @title Get the keyId value for the object (peptide)
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return Th value of parentProtId which is a column of the dataset which contains xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("keyId", "MSnSet", function(obj) {
  out <- obj@experimentData@other$keyId
  out
})



#' @title qfsqfds
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setGeneric("keyId<-", function(obj, value) standardGeneric("keyId<-"))

#' @title Set the keyId value for the object (peptide)
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return xxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setReplaceMethod("keyId", "MSnSet", function(obj, value) {
  obj@experimentData@other$keyId <- value
  obj
})






#' @title qfsqfds
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setGeneric("parentProtId", function(obj) standardGeneric("parentProtId"))


#' @title Get the parentProtId value for the object (peptide)
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return Th value of parentProtId which is a column of the dataset which contains xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("parentProtId", "MSnSet", function(obj) {
  out <- obj@experimentData@other$parentProtId
  out
})

#' @title qfsqfds
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setGeneric("parentProtId<-", function(obj, value) standardGeneric("parentProtId<-"))

#' @title Set the parentProtId value for the object (peptide)
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return xxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setReplaceMethod("parentProtId", "MSnSet", function(obj, value) {
  obj@experimentData@other$parentProtId <- value
  obj
})



#' @title qfsqfds
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setGeneric("DAPARVersion", function(obj) standardGeneric("DAPARVersion"))

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

#' @title qfsqfds
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
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


#' @title qfsqfds
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setGeneric("ProstarVersion", function(obj) standardGeneric("ProstarVersion"))

#' Function to get the version of Prostar used to create the current object
#' @title Get the version of Prostar used.
#' @description sfgsdff sjkfhksqdsq .
#' @param obj xxxx
#' @return xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("ProstarVersion", "MSnSet", function(obj) {
  out <- obj@experimentData@other$Prostar_Version
  out
})


#' @title qfsqfds
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
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




#' @export
setGeneric("Params", function(obj,...) standardGeneric("Params"))

#' @title Set the Params value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return The value of Params which is a column of the dataset which contains xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("Params", "MSnSet", function(obj) {
  out <- obj@experimentData@other$Params
  out
})

#' @export
setGeneric("Params<-", function(obj, value) standardGeneric("Params<-"))

#' @title Set the Params value for the object (peptide)
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return Th value of Params which is a column of the dataset which contains xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setReplaceMethod("Params", "MSnSet", function(obj, value) {
  obj@experimentData@other$Params <- value
  obj
})






#' @export
setGeneric("OriginOfValues", function(obj,...) standardGeneric("OriginOfValues"))

#' @title Get the OriginOfValues value for the object (peptide)
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return The value of OriginOfValues which is a column of the dataset which contains xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("OriginOfValues", "MSnSet", function(obj) {
  out <- obj@experimentData@other$OriginOfValues
  out
})

#' @export
setGeneric("OriginOfValues<-", function(obj, value) standardGeneric("OriginOfValues<-"))



#' @title Set the OriginOfValues value for the object (peptide)
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return The value of OriginOfValues which is a column of the dataset which contains xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setReplaceMethod("OriginOfValues", "MSnSet", function(obj, value) {
  obj@experimentData@other$OriginOfValues <- value
  obj
})




#' @export
setGeneric("typeOfData", function(obj,...) standardGeneric("typeOfData"))

#' @title Get the typeOfData value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return The value of typeOfData which is a column of the dataset which contains xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("typeOfData", "MSnSet", function(obj) {
  out <- obj@experimentData@other$typeOfData
  out
})

#' @export
setGeneric("typeOfData<-", function(obj, value) standardGeneric("typeOfData<-"))

#' @title Set the typeOfData value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return The value of typeOfData which is a column of the dataset which contains xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setReplaceMethod("typeOfData", "MSnSet", function(obj, value) {
  obj@experimentData@other$typeOfData <- value
  obj
})


#' @export
setGeneric("RawPValues", function(obj,...) standardGeneric("RawPValues"))

#' @title Get the RawPValues value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return The value of RawPValues which is a column of the dataset which contains xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("RawPValues", "MSnSet", function(obj) {
  out <- obj@experimentData@other$RawPValues
  out
})

#' @export
setGeneric("RawPValues<-", function(obj, value) standardGeneric("RawPValues<-"))

#' @title Set the RawPValues value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return The value of RawPValues which is a column of the dataset which contains xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setReplaceMethod("RawPValues", "MSnSet", function(obj, value) {
  obj@experimentData@other$RawPValues <- value
  obj
})


#' @export
setGeneric("properties", function(obj,...) standardGeneric("properties"))

#' @title Get the properties value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return The value of properties which is a column of the dataset which contains xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("properties", "MSnSet", function(obj) {
  out <- obj@experimentData@other
  out
})

