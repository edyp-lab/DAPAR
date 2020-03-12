#' @import BiocGenerics MultiAssayExperiment S4Vectors methods


#' @title PipelineTemplate class
#' @description
#' The \code{PipelineTemplate} class inherits from the \code{PipelineTemplate} and serves as a template for
#' instanciate pipeline classes such as proteinPipeline, proteinPipeline, etc. It can be used to manage results of
#' diverse assays on a collection of specimen. Currently,  the class can handle
#' assays that are organized instances of
#' \code{\linkS4class{MSnSet}},
#' \code{matrix}.
#'
#'
#' @slot analysis A character vector that is the name of the MS study
#' @slot pipelineType A character vector that indicates the type of data that are managed in this instance
#' @slot processes xxx
#'
#' @param ... Additional arguments for supporting functions. See details.
#'
#' @return A \code{PipelineTemplate} object
#'
#' @examples
#' example("PipelineTemp")
#' library(DAPARdata)
#' data('Exp1_R25_pept')
#' obj <- Exp1_R25_pept
#' samples <- Biobase::pData(obj)
#' mae <- PipelinePeptide(analysis= 'test',pipelineType = 'peptide', dataType = 'peptide',
#' processes='original',experiments=list(original=obj), 
#' colData=samples)


#' @exportClass PipelinePeptide
.PipelinePeptide <- setClass("PipelinePeptide",
                             contains = "PipelineTemplate",
                             representation = representation(
                               matAdj = "list",
                               CC = "list"
                             ),
                             prototype = prototype()
)



#' @export PipelinePeptide
PipelinePeptide <- function(...)
{
  
  
  template <- PipelineTemplate(...)
  
  parentProtId <- parentProtId(template[['original']])
  if (!is.null(parentProtId)) {
    matAdj <- BuildListAdjacencyMatrices(template[['original']], parentProtId)
    cc <- ComputeConnexComposants(matAdj)
  } else { matAdj <- cc <- NULL}
  
  obj <- new("PipelinePeptide", 
             template, 
             matAdj = matAdj, 
             CC = cc)
  
  obj
}



###########################
# Getter functions



#' @export
setMethod("show", "PipelinePeptide", function(object) {
  cat(
    "The proteinID is: ", parentProtId(object), " \n",
    sep=""
  )
  cat('-----------------------------------------\n')
  callNextMethod()
  
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
