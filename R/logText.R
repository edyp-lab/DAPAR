
#' @title  Build the text information for a new dataset
#'
#' @param l.params A list of parameters related to the process of the dataset
#'
#' @return A string
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' getTextForNewDataset(list(filename = "foo.msnset"))
#'
#' @export
#'
getTextForNewDataset <- function(l.params) {
    if (is.null(l.params) || length(l.params) == 0) {
        return(NULL)
    }

    txt <- tags$ul(as.character(tags$li(paste("Open dataset: ", 
        l.params$filename))))
    return(txt)
}



#' @title  Build the text information for the filtering process
#'
#' @param l.params A list of parameters related to the process of the dataset
#'
#' @return A string
#'
#' @author Samuel Wieczorek
#'
#'
#' @examples 
#' getTextForFiltering(list(filename = "foo.msnset"))
#' 
#' @export
#' @importFrom utils str
#'
getTextForFiltering <- function(l.params) {
    if (is.null(l.params) || length(l.params) == 0) {
        return(NULL)
    }


    txt <- "<ul>"

    if (!is.null(l.params$metacellFilter) && 
            nrow(l.params$metacellFilter) > 1) {
        ll <- l.params$metacellFilter$query
        txt <- paste(txt, "<li>Metacell filtering: ", 
            paste(ll[-1], collapse = ", "), "</li>")
    }

    if (!is.null(l.params$stringFilter.df) && 
            nrow(l.params$stringFilter.df) > 1) {
        ll <- l.params$stringFilter.df$Filter
        txt <- paste(txt, "<li>Text filtering: ", 
            paste(ll[-1], collapse = ", "), "</li>")
    }

    if (!is.null(l.params$numericFilter.df) && 
            nrow(l.params$numericFilter.df) > 1) {
        ll <- l.params$numericFilter.df$Filter
        txt <- paste(txt, "<li>Numerical filtering: ", 
            paste(ll[-1], collapse = ", "), "</li>")
    }
    txt <- paste(txt, "</ul>")
    return(txt)
}





#' @title  Build the text information for the Normalization process
#' 
#' @description 
#' The items of the parameter list for the normalisation is:
#' * method,
#' * type,
#' * varReduction,
#' * quantile,
#  * spanLOESS

#'
#' @param l.params A list of parameters related to the process of the dataset
#'
#' @return A string
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' getTextForNormalization(list(method = "SumByColumns"))
#'
#' @export
#'
getTextForNormalization <- function(l.params) {


    if (is.null(l.params) || length(l.params) == 0) {
        return(NULL)
    }


    txt <- "<ul>"

    txt <- paste(txt, "<li>Norm. method: ", l.params$method, "</li>")
    if (l.params$method != "GlobalQuantileAlignment") {
        txt <- paste(txt, "<li>Application: ", l.params$type, "</li>")
    }

    switch(l.params$method,
        GlobalQuantileAlignment = { },
        SumByColumns = {},
        MeanCentering = {
            txt <- paste(txt, "<li>Variance reduction: ", 
                l.params$varReduction, "</li>")
        },
        QuantileCentering = {
            txt <- paste(txt, "<li>Quantile: ", l.params$quantile, "</li>")
        },
        LOESS = {
            txt <- paste(txt, "<li>Span: ", l.params$spanLOESS, "</li>")
        },
        vsn = {}
    )
    txt <- paste(txt, "</ul>")
    return(txt)
}




#' @title  Build the text information for the peptide Imputation process
#' 
#' @description 
#' * pepLevel_algorithm,
#' * pepLevel_basicAlgorithm,
#' * pepLevel_detQuantile,
#' * pepLevel_detQuant_factor,
#' * pepLevel_imp4p_nbiter,
#' * pepLevel_imp4p_withLapala,
#' * pepLevel_imp4p_qmin,
#' * pepLevel_imp4pLAPALA_distrib
#'
#' @param l.params A list of parameters related to the process of the dataset
#'
#' @return A string
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' params <- list()
#' getTextForpeptideImputation(params)
#'
#' @export
#'
getTextForpeptideImputation <- function(l.params) {


    if (is.null(l.params) || length(l.params) == 0) {
        return(NULL)
    }


    txt <- "<ul>"


    if (l.params$pepLevel_algorithm == "imp4p") {
        txt <- paste(txt, "<li>Algorithm: ", 
            l.params$pepLevel_algorithm, "</li>")
        txt <- paste(txt, "<li>Number of iterations: ", 
            l.params$pepLevel_imp4p_nbiter, "</li>")
        txt <- paste(txt, "<li>MEC imputation: ", 
            l.params$pepLevel_imp4p_withLapala, "</li>")
        if (l.params$pepLevel_imp4p_withLapala) {
            txt <- paste(txt, "<li>Upper lapala bound: ", 
                l.params$pepLevel_imp4p_qmin, "</li>")
            txt <- paste(txt, "<li>Distribution: ", 
                l.params$pepLevel_imp4pLAPALA_distrib, "</li>")
        }
    } else {
        txt <- paste(txt, "<li>Algorithm: ", 
            l.params$pepLevel_basicAlgorithm, "</li>")
        if (l.params$pepLevel_basicAlgorithm == "detQuantile") {
            txt <- paste(txt, "<li>Quantile: ", 
                l.params$pepLevel_detQuantile, "</li>")
            txt <- paste(txt, "<li>Factor: ", 
                l.params$pepLevel_detQuant_factor, "</li>")
        }
        if (l.params$pepLevel_basicAlgorithm == "KNN") {
            txt <- paste(txt, "<li>Nb neighnors: ", 
                l.params$pepLevel_KNN_n, "</li>")
        }
    }

    txt <- paste(txt, "</ul>")
    return(txt)
}


#' @title  Build the text information for the protein Imputation process
#' 
#' @description 
#' * POV_algorithm,
#' * POV_detQuant_quantile,
#' * POV_detQuant_factor,
#' * POV_KNN_n,
#' * MEC_algorithm,
#' * MEC_detQuant_quantile,
#' * MEC_detQuant_factor,
#' * MEC_fixedValue
#'
#' @param l.params A list of parameters related to the process of the dataset
#'
#' @return A string
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' params <- list()
#' getTextForproteinImputation(params)
#'
#' @export
#'
getTextForproteinImputation <- function(l.params) {
    if (is.null(l.params) || length(l.params) == 0) {
        return(NULL)
    }

    txt <- "<ul>"
    txt <- paste(txt, "<li>POV imputation: ", 
        l.params$POV_algorithm, "</li>")
    if (l.params$POV_algorithm == "detQuantile") {
        txt <- paste(txt, "<li>Quantile: ", 
            l.params$POV_detQuant_quantile, "</li>")
        txt <- paste(txt, "<li>Factor: ", 
            l.params$POV_detQuant_factor, "</li>")
    }
    if (l.params$POV_algorithm == "KNN") {
        txt <- paste(txt, "<li>N = ", l.params$POV_KNN_n, "</li>")
    }


    txt <- paste(txt, "<li>MEC imputation: ", l.params$MEC_algorithm, "</li>")
    if (l.params$MEC_algorithm == "detQuantile") {
        txt <- paste(txt, "<li>Quantile: ", 
            l.params$MEC_detQuant_quantile, "</li>")
        txt <- paste(txt, "<li>Factor: ", 
            l.params$MEC_detQuant_factor, "</li>")
    } else if (l.params$MEC_algorithm == "fixedValue") {
        txt <- paste(txt, "<li>Fixed value: ", 
            l.params$MEC_fixedValue, "</li>")
    }

    txt <- paste(txt, "</ul>")
    return(txt)
}



#' @title  Build the text information for the Aggregation process
#' 
#' @description 
#' * includeSharedPeptides,
#' * operator,
#' * considerPeptides,
#' * proteinId,
#' * topN
#'
#' @param l.params A list of parameters related to the process of the dataset
#'
#' @return A string
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' params <- list()
#' getTextForAggregation(params)
#'
#' @export
#'
getTextForAggregation <- function(l.params) {
    if (is.null(l.params) || length(l.params) == 0) {
        return(NULL)
    }

    txt <- "<ul>"

    txt <- paste(txt, "<li>Protein IDs: ", l.params$proteinId, "</li>")
    txt <- paste(txt, "<li>Include shared peptides: ", 
        l.params$includeSharedPeptides, "</li>")
    txt <- paste(txt, "<li>Which peptides to consider: ", 
        l.params$considerPeptides, "</li>")
    txt <- paste(txt, "<li>Operator: ", l.params$operator, "</li>")

    if (l.params$considerPeptides == "onlyN") {
        txt <- paste(txt, "<li>N (most abundant peptides) =", 
            l.params$topN, "</li>")
    }
    txt <- paste(txt, "</ul>")

    return(txt)
}




#' @title  Build the text information for the hypothesis test process
#' 
#' @description 
#' * design,
#' * method,
#' * ttest_options,
#' * th_logFC,
#' * AllPairwiseCompNames = list(
#' *    logFC,
#' *    P_Value)
#'
#' @param l.params A list of parameters related to the process of the dataset
#'
#' @return A string
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' params <- list(design = "OnevsOne", method = "limma")
#' getTextForHypothesisTest(params)
#'
#' @export
#'
getTextForHypothesisTest <- function(l.params) {


    if (is.null(l.params) || length(l.params) == 0) {
        return(NULL)
    }

    txt <- "<ul>"
    txt <- paste(txt, "<li>Constrast: ", l.params$design, "</li>")
    if (l.params$method == "ttests") {
        txt <- paste(txt, "<li>Test: ", l.params$ttest_options, "</li>")
    } else {
        txt <- paste(txt, "<li>Test: ", l.params$method, "</li>")
    }

    txt <- paste(txt, "<li>logFC threshold: ", l.params$th_logFC, "</li>")
    txt <- paste(txt, "</ul>")

    return(txt)
}



#' @title  Build the text information for the Aggregation process
#' 
#' @description 
#' * Condition1
#' * Condition2
#' * Comparison
#' * filterType
#' * filter_th_NA
#' * calibMethod
#' * numValCalibMethod
#' * th_pval
#' * FDR
#' * NbSelected
#'
#' @param l.params A list of parameters related to the process of the dataset
#'
#' @return A string
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' getTextForAnaDiff(list(design = "OnevsOne", method = "Limma"))
#'
#' @export
#'
getTextForAnaDiff <- function(l.params) {


    if (is.null(l.params) || length(l.params) == 0) {
        return(NULL)
    }
    txt <- "<ul>"
    txt <- paste(txt, "<li>The comparison is ", gsub("_", " ", 
        l.params$Comparison, fixed = TRUE), "</li>")
    txt <- paste(txt, "<li>The conditions are ", gsub("_", " ", 
        l.params$Condition1, fixed = TRUE), " and ", gsub("_", " ", 
            l.params$Condition2, fixed = TRUE), "</li>")

    if (!is.null(l.params$filterType) && (l.params$filterType != "None")) {
        txt <- paste(
            txt, "<li>The filter used is ", l.params$filterType,
            "with min nb values / lines: ", l.params$filter_th_NA, "</li>"
        )
    }


    if (!is.null(l.params$calibMethod)) {
        if (!is.null(l.params$numValCalibMethod)) {
            txt <- paste(txt, "<li>The calibration method is ", 
                l.params$calibMethod, ", with num value = ", 
                l.params$numValCalibMethod, "</li>")
        } else {
            txt <- paste(txt, "<li>The calibration method is ", 
                l.params$calibMethod, "</li>")
        }
    }

    if (!is.null(l.params$th_pval)) {
        txt <- paste(txt, "<li>The pvalue threshold is ", 
            l.params$th_pval, "</li>")
    }


    if (!is.null(l.params$FDR)) {
        txt <- paste(txt, "<li>FDR = ", l.params$FDR, "</li>")
    }
    txt <- paste(txt, "</ul>")

    return(txt)
}



#' @title  Build the text information for the Aggregation process
#'
#' @param l.params A list of parameters related to the process of the dataset
#'
#' @return A string
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' getTextForGOAnalysis(list())
#'
#' @export
#'
getTextForGOAnalysis <- function(l.params) {
    if (is.null(l.params) || length(l.params) == 0) {
        return(NULL)
    }

    if (is.null(l.params$whichGO2Save)) {
        return(NULL)
    }
    switch(l.params$whichGO2Save,
        Both = {
            txt <- paste(txt, as.character(tags$li(paste(
                textGOParams, ", GO grouping for level(s):",
                input$GO_level
            ))))
            txt <- paste(txt, as.character(tags$li(paste(
                "GO enrichment with",
                ", adj p-value cutoff = ", input$pvalueCutoff,
                ", universe =", input$universe
            ))))
        },
        Enrichment = {
            txt <- paste(txt, as.character(
                tags$li(paste(textGOParams, " GO enrichment with",
                ", adj p-value cutoff = ", input$pvalueCutoff,
                ", universe =", input$universe,
                sep = " "
            ))))
        },
        Classification = {
            txt <- paste(txt, as.character(
                tags$li(paste(textGOParams, ", GO grouping for level(s):",
                input$GO_level,
                sep = " "
            ))))
        }
    )
}
