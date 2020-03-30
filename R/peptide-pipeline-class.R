#' @import BiocGenerics MultiAssayExperiment S4Vectors methods


#' @title PipelinePeptide class
#' @description xxxxx
#'
#'
#' @slot matAdj A character vector that is the name of the MS study
#' @slot CC A character vector that indicates the type of data that are managed in this instance
#'
#' @return A \code{PipelinePeptide} object
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
#' @export
#' @exportClass PipelinePeptide
.PipelinePeptide <- setClass("PipelinePeptide",
                             contains = "PipelineTemplate",
                             representation = representation(
                               matAdj = "list",
                               CC = "list"
                             ),
                             prototype = prototype()
)



#' @title sdfdsfs
#' @description sfklsjhf qsjdhsqd.
#' @param ... xxxx
#' @return xxxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
PipelinePeptide <- function(...)
{
  template <- DAPAR::PipelineTemplate(...)
  
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





#' @export
setGeneric("CC", function(obj) standardGeneric("CC"))

#' @title sdfdsfs
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return xxxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("CC", "PipelinePeptide", function(obj) {
  out <- xobjCC
  out
})



#' @export
setGeneric("matAdj", function(obj) standardGeneric("matAdj"))

#' @title Get 
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return xxxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("matAdj", "PipelinePeptide", function(obj) {
  out <- obj@matAdj
  out
})



#' @export
setGeneric("pipelineType", function(obj) standardGeneric("pipelineType"))


#' @export
setMethod("pipelineType", "PipelinePeptide", function(obj) {
  callNextMethod()
})




#.validPipelinePeptide <- function(object) {
 # .checkProteinID(object)
#}

#S4Vectors::setValidity2("PipelinePeptide", .validPipelinePeptide)
