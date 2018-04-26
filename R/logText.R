##' Build the text information to be saved after a process on an object of class \code{MSnSet}
##' 
##' @title  Build the text information to be saved
##' @param  name The name of the process in Prostar
##' @param l.params A list of parameters related to the process of the dataset
##' @param ... Parameter for the function \code{getTextForImputation}
##' @return A string
##' @author Samuel Wieczorek
##' @examples
##' buildLogText("Original", list(filename="foo.MSnset"))
buildLogText <- function(name, l.params,...){
    
    txt <- NULL
    switch(name, 
           Original = {txt <- getTextForNewDataset(l.params)},
           Filtering={txt <- getTextForFiltering(l.params)},
           Normalization = {txt <- getTextForNormalization(l.params)},
           Imputation = {txt <- getTextForImputation(l.params, ...)},
           Aggregation = {txt <- getTextForAggregation(l.params)},
           anaDiff ={txt <- getTextForAnaDiff(l.params)},
           GOAnalysis = {txt <- getTextForGOAnalysis(l.params)}
    )
    
    return (txt)
}


##' Build the text information for a new dataset
##' 
##' @title  Build the text information for a new dataset
##' @param l.params A list of parameters related to the process of the dataset
##' @return A string
##' @author Samuel Wieczorek
##' @examples
##' getTextForNewDataset(list(filename="foo.MSnset"))
getTextForNewDataset <- function(l.params){
    txt <- as.character(tags$li(paste("Open dataset: ",l.params$filename)))
    return (txt)
}


##' Build the text information for the filtering process
##' 
##' @title  Build the text information for the filtering process
##' @param l.params A list of parameters related to the process of the dataset
##' @return A string
##' @author Samuel Wieczorek
##' @examples
##' getTextForFiltering(list(mvFilterType="wholeMatrix",mvThNA=3))
getTextForFiltering <- function(l.params){ 
    # str(l.params) = list(mvFilterType ,
    #                 mvThNA,
    #                 stringFilter.df)
    
    txt <- NULL
    
    if (l.params$mvFilterType != "None"){
        txt <- paste(txt, as.character(tags$li(paste("MV filter with ", l.params$mvFilterType, ", and minimal nb of values per lines = ", l.params$mvThNA))))
    }
    if (!is.null(l.params$stringFilter.df) && nrow(l.params$stringFilter.df) > 1){
        ll <- l.params$stringFilter.df$Filter
        txt <- paste(txt, as.character(tags$li(paste("Text filtering based on", paste(ll[-1], collapse=", ")))))
    }
    
    return (txt)
}




##' Build the text information for the Normalization process
##' 
##' @title  Build the text information for the Normalization process
##' @param l.params A list of parameters related to the process of the dataset
##' @return A string
##' @author Samuel Wieczorek
##' @examples
##' getTextForNormalization(list(method="Sum by columns"))
getTextForNormalization <- function(l.params){ 
    
    # str(l.params) = list(method ,
    #                 type,
    #                 varReduction,
    #                 quantile,
    #                 otherQuantile)
    
    
    txt <- NULL
    if (l.params$method=="None") return (NULL)
    else if (l.params$method == "Global quantile alignment"){
        txt <- paste(as.character(tags$li(paste(txt,l.params$method))))
    }
    
    else if  (l.params$method == "Sum by columns"){
        txt <- paste(txt,as.character(tags$li(paste(l.params$method, " - ", l.params$type))))
    }
    
    else if  (l.params$method == "Mean Centering"){
        
        if ( isTRUE(l.params$varReduction )){
            txt <- paste(txt,as.character(tags$li(paste(l.params$method, " - ", l.params$type, " with variance reduction"))))
        } else {
            txt <- paste(txt,as.character(tags$li(paste(l.params$method, " - ", l.params$type))))
        }
    }
    
    else if  (l.params$method == "Quantile Centering"){
        txt <- paste(txt,as.character(tags$li(l.params$method, " - ", l.params$type)))
        quant <- ifelse (l.params$quantile == "Other",l.params$otherQuantile, l.params$quantile)
        txt <- paste(txt,as.character(tags$li("with quantile =", quant)))
    }
    
    return (txt)
}



##' Build the text information for the Imputation process
##' 
##' @title  Build the text information for the Imputation process
##' @param l.params A list of parameters related to the process of the dataset
##' @param level The type of data (peptide or protein)
##' @return A string
##' @author Samuel Wieczorek
##' @examples
##' getTextForImputation(list(POV_algorithm="slsa",MEC_algorithm="fixedValue", MEC_fixedValue = 0), level="protein")
getTextForImputation <- function(l.params, level){ 
    
    txt <- NULL
    
    
    switch(level, 
           
           protein = {
               # l.params <- list(
            #    POV_algorithm,
            #    POV_detQuant_quantile,
            #    POV_detQuant_factor,
            #    POV_KNN_n,
            #    MEC_algorithm,
            #    MEC_detQuant_quantile,
            #    MEC_detQuant_factor,
            #    MEC_fixedValue)
    
    
    
            if (!is.null(l.params$POV_algorithm)){
                if (l.params$POV_algorithm=="None") {
                    return (NULL)
            } else {
                txt <- paste(txt,as.character(tags$li(paste("POV imputed with ", l.params$POV_algorithm))))
            }
        
            switch(l.params$POV_algorithm,
               slsa = {},
               detQuantile = {txt <- paste(txt,as.character(tags$li(paste("quantile= ", l.params$POV_detQuant_quantile, 
                                           ", factor=",l.params$POV_detQuant_factor))))
               },
               KNN = {txt <- paste(txt,as.character(tags$li(paste("neighbors= ", l.params$POV_KNN_n))))}
                )
            }
    
            if (!is.null(l.params$MEC_algorithm)){
                if (l.params$MEC_algorithm=="None") {
                return (NULL)
            } else {
                txt <- paste(txt,as.character(tags$li(paste("MEC imputed with ", l.params$MEC_algorithm))))
            }
        
            switch(l.params$MEC_algorithm,
                    detQuantile = {txt <- paste(txt,as.character(tags$li(paste("quantile= ", l.params$MEC_detQuant_quantile, 
                                           ", factor=",l.params$MEC_detQuant_factor))))
                                    },
                    fixedValue = {txt <- paste(txt,as.character(tags$li(paste("fixed value= ", l.params$MEC_fixedValue))))}
                )
        
                }
           },
    peptide = {
        # l.params <- list(pepLevel_algorithm
        #                  pepLevel_basicAlgorithm,
        #                  pepLevel_detQuantile,
        #                  pepLevel_detQuant_factor,
        #                  pepLevel_imp4p_nbiter,
        #                  pepLevel_imp4p_withLapala,
        #                  pepLevel_imp4p_qmin,
        #                  pepLevel_imp4pLAPALA_distrib)
        
        if (!is.null(l.params$pepLevel_algorithm)){
            
         
            if (l.params$pepLevel_algorithm=="None") {
                return (NULL)
            } else if (l.params$pepLevel_algorithm=="imp4p"){
                txt <- paste(txt,as.character(tags$li(paste("Missing values imputed with ", l.params$pepLevel_algorithm,
                                                            ", nb iter = ", l.params$pepLevel_imp4p_nbiter ))))
                if (!is.null(l.params$pepLevel_imp4p_withLapala) && isTRUE(l.params$pepLevel_imp4p_withLapala)){
                    txt <- paste(txt,as.character(tags$li(paste("Imputation of MEC: "))))
                    txt <- paste(txt,as.character(tags$li(paste("Upper bound: ", l.params$pepLevel_imp4p_qmin))))
                    txt <- paste(txt,as.character(tags$li(paste("DIstribution: ", l.params$pepLevel_imp4pLAPALA_distrib))))
                }
            }
            else if(!is.null(l.params$pepLevel_basicAlgorithm) && (l.params$pepLevel_algorithm=="Basic methods")){
                txt <- paste(txt,as.character(tags$li(paste("Dataset imputed with ", l.params$pepLevel_basicAlgorithm))))
                switch (l.params$pepLevel_basicAlgorithm,
                    detQuantile = {
                      txt <- paste(txt,as.character(tags$li(paste("Quantile = ", l.params$pepLevel_detQuantile,
                                                                  ", Factor =", l.params$pepLevel_detQuant_factor))))
                                   },
                    KNN = {},
                    MLE = {}
                )
            }
        }      
    }
    )
    
    return (txt)
    
}


##' Builds the text information for the Aggregation process
##' 
##' @title  Build the text information for the Aggregation process
##' @param l.params A list of parameters related to the process of the dataset
##' @return A string
##' @author Samuel Wieczorek
##' @examples
##' getTextForAggregation(list(POV_algorithm="slsa",MEC_algorithm="fixedValue", MEC_fixedValue = 0))
getTextForAggregation <- function(l.params){ 
    
    # l.params <- list(withSharedPeptides,
    #                  agregMethod,
    #                  proteinId,
    #                  topN
    #                 )
    
    txt <- NULL
    if (is.null(l.params$proteinId) || (l.params$proteinId =="None")) { return (NULL)}
    
    txt <- paste(txt,as.character(tags$li(paste("proteinId:", l.params$proteinId, " "))))
    if (!is.null(l.params$agregMethod) && (l.params$agregMethod =="none")) {
        txt <- paste(txt,as.character(tags$li(paste("method:", l.params$agregMethod," "))))
        if (l.params$agregMethod =="sum on top n"){
            txt <- paste(txt,as.character(tags$li(paste("n=", l.params$topN," "))))
        }
    }
    
    if (!is.null(l.params$withSharedPeptides) ){
        if(isTRUE(l.params$withSharedPeptides)){
        txt <- paste(txt,as.character(tags$li(paste("with shared peptides"))))
    } else {
        txt <- paste(txt,as.character(tags$li(paste("without shared peptides"))))
    }
    }
    return (txt)
    
}


##' Build the text information for the differential Analysis process
##' 
##' @title  Build the text information for the Aggregation process
##' @param l.params A list of parameters related to the process of the dataset
##' @return A string
##' @author Samuel Wieczorek
##' @examples
##' getTextForAnaDiff(list(design="OnevsOne",method="Limma"))
getTextForAnaDiff <- function(l.params){ 
    
    # l.params <- list(design,
    #                  method,
    #                  ttest_options,
    #                  th_logFC,
    #                  AllPairwiseCompNames,
    #                  comp,
    #                  th_pval,
    #                  calibMethod,
    #                  fdr,
    #                  swapVolcano,
    #                  filterType,
    #                  filter_th_NA,
    #                  numValCalibMethod)
    
    
    
    txt <- NULL                 
    
    if (is.null(l.params$design) || (l.params$design =="None")) { return (NULL)}
    
    txt <- paste(txt, as.character(tags$li(paste("Design:", l.params$design)))) 
    txt <- paste(txt, as.character(tags$li(paste("Method:", l.params$method))))
    
    if (!is.null(l.params$ttest_options)){
      txt <- paste(txt, as.character(tags$li(paste("t-test:", l.params$ttest_options))))
    } 
    
    txt <- paste(txt, as.character(tags$li(paste("Threshold logFC = ", l.params$th_logFC))))
    
    if (!is.null(l.params$comp)){
      txt <- paste(txt, as.character(tags$li(paste("Comparison: ", l.params$comp))))
        
    }
    
    if (!is.null(l.params$filterType) && (l.params$filterType != "None")){
      txt <- paste(txt, as.character(tags$li(paste("Filter type: ", l.params$filterType, 
                                                   ", min nb values / lines: ", l.params$filter_th_NA))))
    }
    
    if (!is.null(l.params$swapVolcano) && (isTRUE(l.params$swapVolcano))){
      txt <- paste(txt, as.character(tags$li("Swap volcano")))
    }
    
    if (!is.null(l.params$calibMethod) ){
        if (!is.null(l.params$numValCalibMethod)){
          txt <- paste(txt, as.character(tags$li(paste("Calibration with ", l.params$calibMethod, ", num value = ", l.params$numValCalibMethod))))
        } else {
          txt <- paste(txt, as.character(tags$li(paste("Calibration with ", l.params$calibMethod))))
        }
    }
    
    if (!is.null(l.params$th_pval)){
      txt <- paste(txt, as.character(tags$li(paste("Threshold pvalue = ", l.params$th_pval))))
    }
    
    return (txt)
}


##' Build the text information for the Aggregation process
##' 
##' @title  Build the text information for the Aggregation process
##' @param l.params A list of parameters related to the process of the dataset
##' @return A string
##' @author Samuel Wieczorek
##' @examples
##' getTextForGOAnalysis(list())
getTextForGOAnalysis <-  function(l.params){
  { 
      if (is.null(l.params$whichGO2Save)){return(NULL)}
      switch(l.params$whichGO2Save,
           Both =
              {
                txt <- paste(txt, as.character(tags$li(paste(textGOParams,", GO grouping for level(s):",
                           input$GO_level))))
                txt <- paste(txt, as.character(tags$li(paste("GO enrichment with",
                           ", adj p-value cutoff = ", input$pvalueCutoff,
                           ", universe =", input$universe))))
              },
           Enrichment ={
             txt <- paste(txt, as.character(tags$li(paste(textGOParams, " GO enrichment with",
                           ", adj p-value cutoff = ", input$pvalueCutoff,
                           ", universe =", input$universe, sep= " "))))
              },
           Classification= {
             txt <- paste(txt, as.character(tags$li(paste(textGOParams,", GO grouping for level(s):",
                           input$GO_level,sep=" "))))
            }
  )
  }
}