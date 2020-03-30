#' @import BiocGenerics MultiAssayExperiment S4Vectors methods


#' @title PipelineTemplate class
#' @description
#' The \code{PipelineTemplate} class inherits from the \code{MultiAssayExperiment} and serves as a template for
#' instanciate pipeline classes such as proteinPipeline, proteinPipeline, etc. It can be used to manage results of
#' diverse assays on a collection of specimen. Currently,  the class can handle
#' assays that are organized instances of
#' \code{\linkS4class{MSnSet}},
#' \code{matrix}.
#'
#'
#' @slot PairwiseComparisons A \code{list} to store the result of hypothesis tests.
#' @slot indexNA A xxxxxs
#' @slot dataType xxxx.
#' @slot analysis A character vector that is the name of the MS study
#' @slot pipelineType A character vector that indicates the type of data that are managed in this instance
#' @slot processes xxx
#' @slot version xxx
#'
#' @param ... Additional arguments for supporting functions. See details.
#'
#' @return A \code{PipelineTemplate} object
#'
#' @examples
#' example("PipelineTemp")
#' library(DAPARdata)
#' data('Exp1_R25_prot')
#' obj <- Exp1_R25_prot
#' samples <- Biobase::pData(Exp1_R25_prot)
#' mae <- PipelineTemplate(analysis= 'test',pipelineType = 'protein', dataType = 'protein', 
#' processes='original',experiments=list(original=Exp1_R25_prot), 
#' colData=samples)

#' @import MultiAssayExperiment
#' @importFrom utils installed.packages
#' @exportClass PipelineTemplate

.PipelineTemplate <- setClass("PipelineTemplate",
                          slots= list(
                            PairwiseComparisons = "list",
                            analysis = "character",
                            dataType = "character",
                            pipelineType = "character",
                            processes = "character",
                            version = "character"
                          ),
                          contains=c("MultiAssayExperiment"),

)


#' Function to create a new instance of the class \code{PipelineTemplate}
#' @title xxxxx.
#' @param PairwiseComparisons xxxx
#' @param analysis A string that is the name of the object (ie the name of the MS study)
#' @param dataType A string that is the type of data
#' @param pipelineType the type of the pipeline used in Prostar
#' @param processes a character vector that indicates the list of the processes which composed the pipeline
#' @param version The version of Prostar used to create this object
#' @param ... sdfkjdshfhdsfk
#' @return An instance of class \code{PipelineTemplate}.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
PipelineTemplate <- function(
  PairwiseComparisons = list(),
  analysis = character(),
  dataType = character(),
  pipelineType = character(),
  processes = character(),
  version = character(),
  ...)
  {
  mae <- MultiAssayExperiment(...)
  
  if (length(experiments(mae)) != length(c('original',processes))){
    warning('The number of experiment dataset must be equal to the number of processes in the pipeline.')
    return(NULL)
  }
  
  
  ## on configure le dataset avant de le stocker dans mae
  experiments(mae)[[1]] <-.ConfigureDataset(experiments(mae)[[1]])
  .version = if (is.na(utils::installed.packages()["Prostar"])) 'NA' else utils::installed.packages()["Prostar",'Version']
  obj <- new("PipelineTemplate",mae,
                     PairwiseComparisons = PairwiseComparisons,
                    analysis=analysis, 
                    dataType =dataType,
                    pipelineType=pipelineType,
                    version = .version,
                    processes=c('original',processes)
  )
  validObject(obj)
  obj
  }


#' Function to do some modifications of a MSnSet class dataset to be compatible with Prostar
#' @title xxxxx.
#' @param x An object of class \code{PipelineTemplate}
#' @return An instance of class \code{PipelineTemplate}.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
.ConfigureDataset <- function(x) {
  if(class(x) != "MSnSet"){
    warning("The class of object is not MSnSet")
    return(NULL)
  }
  out <- addOriginOfValue(x)
  out <- ReplaceDotsByUnderscore(out)
  #tmp <- setIndexNA(tmp,which(is.na(data)))
  out
}





## Create getter methods for 1D data structures

#' @export
setGeneric("analysis", function(obj) standardGeneric("analysis"))


#' @export
setGeneric("pipelineType", function(obj) standardGeneric("pipelineType"))

#' @export
setGeneric("dataType", function(obj) standardGeneric("dataType"))


#' @export
setGeneric("processes", function(obj) standardGeneric("processes"))


#' @export
setGeneric("version", function(obj) standardGeneric("version"))


#' @title Get the version slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return The value of the slot version which is a column of the dataset which contains xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("version", "PipelineTemplate", function(obj) {
  out <- obj@version
  out
})




#' @title Get the dataType slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return The value of the dataType version which is a column of the dataset which contains xxxx.
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("dataType", "PipelineTemplate", function(obj) {
  out <- obj@dataType
  out
})

#' @title Get the analysis slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return The value of the slot called analysis 
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("analysis", "PipelineTemplate", function(obj) {
  out <- obj@analysis
  out
})

#' @title Get the pipelineType slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return The value of the slot called pipelineType 
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("pipelineType", "PipelineTemplate", function(obj) {
  out <- obj@pipelineType
  out
})

#' @title Get the processes slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return The value of the slot called processes 
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("processes", "PipelineTemplate", function(obj) {
  out <- obj@processes
  out
})




#' @export
setGeneric("PairwiseComps", function(obj) standardGeneric("PairwiseComps"))


##
## Create getter methods for 2D data structures
##
#' @title Get the PairwiseComps slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return The value of the slot called PairwiseComps 
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("PairwiseComps", "PipelineTemplate", function(obj) {
  out <- obj@PairwiseComparisons
  out
})





#' @export
setGeneric("designMap", function(obj) standardGeneric("designMap"))


##
## Create getter methods for 2D data structures
##
#' @title Get the designMap slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @return The value of the slot called designMap 
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("designMap", "PipelineTemplate", function(obj) {
  out <- Biobase::pData(MultiAssayExperiment::experiments(obj)[['original']])
  out
})







# For MultiAssayExperiment slots
# The getter methods defined in MultiAssayExperiment can be directly used to retrieve 
# data from slots in the base class. These should generally not require any re-defining 
# for a derived class. However, if it is necessary, the methods should use callNextMethod 
# internally. This will call the method for the base MultiAssayExperiment class, the output
# of which can be modified as required.


# We use setValidity2 to define a validity function for PipelineTemplate. We use the previously 
# defined getter functions to retrieve the slot values rather than using @. This is generally 
# a good idea to keep the interface separate from the implementation
# (This protects against changes to the slot names, and simplifies development when the 
#   storage mode differs from the conceptual meaning of the data, e.g., for efficiency 
#   purposes.)
# We also set withDimnames=FALSE in our getter calls, as consistent naming is not necessary 
# for internal functions.



# Creating a show method
# The default show method will only display information about the MultiAssayExperiment slots. 
# We can augment it to display some relevant aspects of the custom slots. This is done by 
# calling the base show method before printing additional fields as necessary.


#' @title show function override
#' @description sfklsjhf qsjdhsqd.
#' @param object xxxx
#' @return The value of the slot called designMap 
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("show", "PipelineTemplate", function(object) {
  cat(
    "The name of the pipeline is ", pipelineType(object), " \n",
    "The pipeline is composed of the following processes: ", processes(object), " \n",
    "The name of the analysis is ", analysis(object), " \n",
    "Built under Prostar version ", version(object), " \n",
     sep=""
  )
  cat('-----------------------------------------\n')
  callNextMethod()
  
})

# Creating setter methods
# 3.6.1 For 1D data structures
# We define some setter methods for the custom slots containing the 1D structures. 
# Again, this usually requires the creation of new generics.


#' @export
setGeneric("pipelineType<-", function(obj, value) standardGeneric("pipelineType<-"))


#' @export
setGeneric("analysis<-", function(obj, value) standardGeneric("analysis<-"))


#' @export
setGeneric("processes<-", function(obj, value) standardGeneric("processes<-"))


# We define the class-specific methods for these generics. Note that use of validObject 
# to ensure that the assigned input is still valid.





#' @title Set the pipelineType slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return xxxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setReplaceMethod("pipelineType", "PipelineTemplate", function(obj, value) {
  x@pipelineType <- value
  validObject(x)
  x
})

#' @title Set the processes slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return xxxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setReplaceMethod("processes", "PipelineTemplate", function(obj, value) {
  obj@processes <- value
  validObject(obj)
  obj
})

#' @title Set the analysis slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return xxxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setReplaceMethod("analysis", "PipelineTemplate", function(obj, value) {
  obj@analysis <- value
  validObject(obj)
  obj
})



# For 2D data structures
# We repeat this process for the 2D structures.


#' @export
setGeneric("PairwiseComps<-", function(obj, ..., value) standardGeneric("PairwiseComps<-"))

#' @title Set the PairwiseComps slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param value xxx
#' @return xxxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setReplaceMethod("PairwiseComps", "PipelineTemplate", function(obj, value) {
  obj@PairwiseComparisons <- value
  validObject(obj)
  obj
})



#' @export
setGeneric("rmDatasetByIndice", function(obj, ind) standardGeneric("rmDatasetByIndice"))

#' @title Set the rmDatasetByIndice slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param ind xxx
#' @return xxxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("rmDatasetByIndice", "PipelineTemplate", function(obj, ind) {
  #mae <- callNextMethod()
  experiments(obj) <- experiments(obj)[-ind] 
  validObject(obj)
  obj
})


#' @export
setGeneric("rmDatasetByName", function(obj, name) standardGeneric("rmDatasetByName"))

#' @title Set the rmDatasetByName slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param name xxx
#' @return xxxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("rmDatasetByName", "PipelineTemplate", function(obj, name) {
  #mae <- callNextMethod()
  experiments(obj) <- within(experiments(obj), rm(name)) 
  validObject(obj)
  obj
})

#' @export
setGeneric("addDataset", function(obj, name, dataset) standardGeneric("addDataset"))

#' @title Set the addDataset slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param name xxx
#' @param dataset xxx
#' @return xxxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("addDataset", "PipelineTemplate", function(obj, name, dataset) {
  #mae <- callNextMethod()
  ds <- list()
  ds[[name]]<- dataset
  obj <- c(obj, ds)
  validObject(obj)
  obj
})



#' @export
setGeneric("updateDataset<-", function(obj, name, value) standardGeneric("updateDataset<-"))


#' @title Set the updateDataset slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param name xxx
#' @param value xxx
#' @return xxxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setReplaceMethod("updateDataset", "PipelineTemplate", function(obj, name, value) {
  #mae <- callNextMethod()
  .checkIfAnalysisExists(obj, name)
  
  experiments(obj)[[name]] <- value
  validObject(obj)
  obj
})


#' @export
setGeneric("dataset", function(obj, name) standardGeneric("dataset"))

#' @title Get a dataset slot value for the object
#' @description sfklsjhf qsjdhsqd.
#' @param obj xxxx
#' @param name xxx
#' @return xxxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
setMethod("dataset", "PipelineTemplate", function(obj, name) {
  #mae <- callNextMethod()
  .checkIfAnalysisExists(obj, name)
  out <- experiments(obj)[[name]]
  out
})





#' @title Check if xxx
#' @description sfklsjhf qsjdhsqd.
#' @param object xxxx
#' @param name xxx
#' @return xxxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
.checkIfAnalysisExists <- function(object, name)
{
  if (!(name %in% names(assays(object)))){
    warning("The dataset called name was not found")
    return(NULL)
  }
  
  if (class(name) != "character"){
    warning("The name parameter must be a string.")
    return(NULL)
  }
}


#' @title Check if xxx
#' @description sfklsjhf qsjdhsqd.
#' @param object xxxx
#' @return xxxx
#' @author Samuel Wieczorek
#' @examples 
#' utils::data(Exp1_R25_prot)
#' @export
.validPipelineTemplate <- function(object) {
  #if (length(experiments(object)) != 0L) {
  
  # }
}



S4Vectors::setValidity2("PipelineTemplate", .validPipelineTemplate)
















