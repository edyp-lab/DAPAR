#' @import BiocGenerics MultiAssayExperiment S4Vectors methods


#' @title PipelineProtein class
#' @description xx x 
#'
#'
#' @param ... Additional arguments for supporting functions. See details.
#'
#' @return A \code{PipelineProtein} object
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
