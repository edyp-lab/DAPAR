
#' @exportClass PipelinePeptide
.PipelinePeptide <- setClass("PipelinePeptide",
                             contains = "PipelineTemplate",
                             representation = representation(
                               proteinID = "character",
                               matAdj = "list",
                               CC = "list"
                             ),
                             prototype = prototype()
)



#' @export PipelinePeptide
PipelinePeptide <- function(proteinID = character(),
                            ...)
{
  if (is.null(proteinID) || length(proteinID)==0){
    warning("The proteinID is not set.")
    return(NULL)
  }
  
  template <- PipelineTemplate(...)
  matAdj <- BuildListAdjacencyMatrices(template[['original']], proteinID)
  cc <- ComputeConnexComposants(matAdj)
  obj <- new("PipelinePeptide", template, proteinID = proteinID, matAdj = matAdj, CC = cc)
  
  obj
}



###########################
# Getter functions



#' @export
setMethod("show", "PipelinePeptide", function(object) {
  cat(
    "The proteinID is: ", proteinID(object), " \n",
    sep=""
  )
  cat('-----------------------------------------\n')
  callNextMethod()
  
})





#' @export
setGeneric("proteinID", function(x) standardGeneric("proteinID"))

#' @export
setMethod("proteinID", "PipelinePeptide", function(x) {
  out <- x@proteinID
  out
})



#' @export
setGeneric("CC", function(x) standardGeneric("CC"))

#' @export
setMethod("CC", "PipelinePeptide", function(x) {
  out <- x@CC
  out
})



#' @export
setGeneric("matAdj", function(x) standardGeneric("matAdj"))

#' @export
setMethod("matAdj", "PipelinePeptide", function(x) {
  out <- x@matAdj
  out
})





#.validPipelinePeptide <- function(object) {
 # .checkProteinID(object)
#}

#S4Vectors::setValidity2("PipelinePeptide", .validPipelinePeptide)
