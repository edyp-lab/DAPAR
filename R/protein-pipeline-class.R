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
#' @slot PairwiseComparisons A \code{list} to store the result of hypothesis tests.
#' @slot indexNA A xxxxxs
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
#' data('Exp1_R25_prot')
#' samples <- Biobase::pData(Exp1_R25_prot)
#' mae <- PipelineProtein(analysis= 'test',pipelineType = 'protein', dataType = 'protein',
#' processes='original',experiments=list(original=Exp1_R25_prot), 
#' colData=samples)


#' @exportClass PipelineProtein
.PipelineProtein <- setClass("PipelineProtein",
          contains = "PipelineTemplate", 
          representation = representation(),
          prototype = prototype()
         )



#' @export PipelineProtein
PipelineProtein <- function(...)
{
  template <- PipelineTemplate(...)
  obj <- new("PipelineProtein", template)
  
  obj
}
