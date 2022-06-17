



#' This function is a wrappper to the function adjust.p from the
#' cp4p package. It returns the FDR corresponding to the p-values of the
#' differential analysis.
#' The FDR is computed with the function \code{p.adjust}\{stats\}..
#'
#' @title Computes the FDR corresponding to the p-values of the
#' differential analysis using
#'
#' @param logFC The result (logFC values) of the differential analysis processed
#' by \code{\link{limmaCompleteTest}}
#'
#' @param pval The result (p-values) of the differential analysis processed
#' by \code{\link{limmaCompleteTest}}
#'
#' @param threshold_PVal The threshold on p-pvalue to
#' distinguish between differential and non-differential data
#'
#' @param threshold_LogFC The threshold on log(Fold Change) to
#' distinguish between differential and non-differential data
#'
#' @param pi0Method The parameter pi0.method of the method adjust.p
#' in the package \code{cp4p}
#'
#' @return The computed FDR value (floating number)
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_prot, package = "DAPARdata")
#' obj <- Exp1_R25_prot[1:1000]
#' level <- obj@experimentData@other$typeOfData
#' metacell.mask <- match.metacell(GetMetacell(obj), "missing", level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj <- MetaCellFiltering(obj, indices, cmd = "delete")
#' qData <- Biobase::exprs(obj$new)
#' sTab <- Biobase::pData(obj$new)
#' limma <- limmaCompleteTest(qData, sTab)
#' diffAnaComputeFDR(limma$logFC[, 1], limma$P_Value[, 1])
#'
#' @export
#'
diffAnaComputeFDR <- function(logFC,
                              pval, threshold_PVal = 0,
                              threshold_LogFC = 0,
                              pi0Method = 1) {
    # require(cp4p)
    if (is.null(logFC) || is.null(pval)) {
        return()
    }

    upItems <- which(abs(logFC) >= threshold_LogFC)

    selectedItems <- pval[upItems]

    padj <- cp4p::adjust.p(selectedItems, pi0Method)

    items <- which(-log10(padj$adjp[, 1]) >= threshold_PVal)

    BH.fdr <- max(padj$adjp[items, 2])

    return(BH.fdr)
}



#' This method returns a list of the statistical tests performed with DAPAR and
#' recorded in an object of class \code{MSnSet}.
#'
#' @title Returns list that contains a list of the statistical tests performed
#' with DAPAR and recorded in an object of class \code{MSnSet}.
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @return A list of two slots: logFC and P_Value
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_prot, package = "DAPARdata")
#' obj <- Exp1_R25_prot[1:1000]
#' level <- GetTypeofData(obj)
#' metacell.mask <- match.metacell(GetMetacell(obj), "missing", level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj <- MetaCellFiltering(obj, indices, cmd = "delete")
#' qData <- Biobase::exprs(obj$new)
#' sTab <- Biobase::pData(obj$new)
#' allComp <- limmaCompleteTest(qData, sTab)
#' data <- list(logFC = allComp$logFC[1], P_Value = allComp$P_Value[1])
#' obj$new <- diffAnaSave(obj$new, allComp, data)
#' ll <- Get_AllComparisons(obj$new)
#'
#' @export
#'
Get_AllComparisons <- function(obj) {
    logFC_KEY <- "_logFC"
    pvalue_KEY <- "_pval"

    ####### SAVE ALL THEPAIRWISE COMPARISON RESULTS
    res_AllPairwiseComparisons <- NULL

    # If there are already pVal values, then do no compute them
    if (length(grep(logFC_KEY, names(Biobase::fData(obj)))) > 0) {
        res_AllPairwiseComparisons <- list(
            logFC = dplyr::select(Biobase::fData(obj), 
                grep(logFC_KEY, names(Biobase::fData(obj)))),
            P_Value = dplyr::select(Biobase::fData(obj), 
                grep(pvalue_KEY, names(Biobase::fData(obj))))
        )
    }

    return(res_AllPairwiseComparisons)
}




#' This method returns a class \code{MSnSet} object with the results
#' of differential analysis.
#'
#' @title Returns a \code{MSnSet} object with the results of
#' the differential analysis performed with \code{\link{limma}} package.
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @param allComp A list of two items which is the result of the function
#' wrapper.limmaCompleteTest or xxxx
#'
#' @param data The result of the differential analysis processed
#' by \code{\link{limmaCompleteTest}}
#'
#' @param th_pval xxx
#'
#' @param th_logFC xxx
#'
#' @return A MSnSet
#'
#' @author Alexia Dorffer, Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_prot, package = "DAPARdata")
#' obj <- Exp1_R25_prot[1:1000]
#' level <- obj@experimentData@other$typeOfData
#' metacell.mask <- match.metacell(GetMetacell(obj), "missing", level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj <- MetaCellFiltering(obj, indices, cmd = "delete")
#' qData <- Biobase::exprs(obj$new)
#' sTab <- Biobase::pData(obj$new)
#' allComp <- limmaCompleteTest(qData, sTab)
#' data <- list(logFC = allComp$logFC[1], P_Value = allComp$P_Value[1])
#' diffAnaSave(obj$new, allComp, data)
#'
#' @export
#'
#'
diffAnaSave <- function(obj,
                        allComp,
                        data = NULL,
                        th_pval = 0,
                        th_logFC = 0) {
    if (is.null(allComp)) {
        warning("The analysis has not been completed. Maybe there
            are some missing values in the dataset. If so, please impute before
            running differential analysis")
        return(NULL)
    }


    ####### SAVE ALL THE PAIRWISE COMPARISON RESULTS

    .fc <- as.data.frame(allComp$logFC)
    .pval <- as.data.frame(allComp$P_Value)
    cnames <- c(colnames(allComp$logFC), colnames(allComp$P_Value))
    ind <- which(colnames(Biobase::fData(obj)) %in% cnames)
    if (length(ind) > 0) {
        Biobase::fData(obj) <- Biobase::fData(obj)[, -ind]
    }

    for (i in 1:ncol(.fc)) {
        Biobase::fData(obj) <- cbind(Biobase::fData(obj), .fc[, i], .pval[, i])
        coln <- colnames(Biobase::fData(obj))
        colnames(Biobase::fData(obj))[(length(coln) - 1):length(coln)] <- 
            c(colnames(allComp$logFC)[i], colnames(allComp$P_Value)[i])
    }

    text <- paste("Null hypothesis test")
    obj@processingData@processing <- c(obj@processingData@processing, text)
    # Save parameters

    obj@experimentData@other$RawPValues <- TRUE

    #### SAVE A COMPARISON ANALYSIS IF EXISTS
    if (!(is.null(data$logFC) && is.null(data$P_Value))) {
        Biobase::fData(obj)$P_Value <- data$P_Value
        Biobase::fData(obj)$logFC <- data$logFC
        Biobase::fData(obj)$Significant <- 0

        ## setSignificant info
        x <- data$logFC
        y <- -log10(data$P_Value)

        ipval <- which(y >= th_pval)
        ilogfc <- which(abs(x) >= th_logFC)
        Biobase::fData(obj)[intersect(ipval, ilogfc), ]$Significant <- 1
    }


    return(obj)
}


#' Returns a MSnSet object with only proteins significant after
#' differential analysis.
#'
#' @title Returns a MSnSet object with only proteins significant after
#' differential analysis.
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @return A MSnSet
#'
#' @author Alexia Dorffer
#'
#' @examples
#' utils::data(Exp1_R25_prot, package = "DAPARdata")
#' obj <- Exp1_R25_prot[1:1000]
#' level <- obj@experimentData@other$typeOfData
#' metacell.mask <- match.metacell(GetMetacell(obj), "missing", level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj <- MetaCellFiltering(obj, indices, cmd = "delete")
#' qData <- Biobase::exprs(obj$new)
#' sTab <- Biobase::pData(obj$new)
#' allComp <- limmaCompleteTest(qData, sTab)
#' data <- list(logFC = allComp$logFC[1], P_Value = allComp$P_Value[1])
#' obj$new <- diffAnaSave(obj$new, allComp, data)
#' signif <- diffAnaGetSignificant(obj$new)
#'
#' @export
#'
#'
diffAnaGetSignificant <- function(obj) {
    if (is.null(obj)) {
        warning("The dataset contains no data")
        return(NULL)
    }
    if (!("Significant" %in% colnames(Biobase::fData(obj)))) {
        warning("Please Set Significant data before")
        return(NULL)
    }
    temp <- obj
    signif <- which(Biobase::fData(temp)$Significant == TRUE)
    return(temp[signif, ])
}



#' This function is a wrapper to the calibration.plot method of the
#' \code{cp4p} package for use with \code{MSnSet} objects.
#'
#' @title Performs a calibration plot on an \code{MSnSet} object,
#' calling the \code{cp4p} package functions.
#'
#' @param vPVal A dataframe that contains quantitative data.
#'
#' @param pi0Method A vector of the conditions (one condition per sample).
#'
#' @return A plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_prot, package = "DAPARdata")
#' obj <- Exp1_R25_prot[1:1000]
#' level <- obj@experimentData@other$typeOfData
#' metacell.mask <- match.metacell(GetMetacell(obj), "missing", level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj <- MetaCellFiltering(obj, indices, cmd = "delete")
#' qData <- Biobase::exprs(obj$new)
#' sTab <- Biobase::pData(obj$new)
#' limma <- limmaCompleteTest(qData, sTab)
#' wrapperCalibrationPlot(limma$P_Value[, 1])
#'
#' @export
#'
wrapperCalibrationPlot <- function(vPVal, pi0Method = "pounds") {
    # require(cp4p)
    if (is.null(vPVal)) {
        return(NULL)
    }

    p <- cp4p::calibration.plot(vPVal, pi0.method = pi0Method)

    return(p)
}



#' This function plots a histogram ov p-values
#'
#' @title Plots a histogram ov p-values
#'
#' @param pval_ll xxx
#'
#' @param bins xxx
#'
#' @param pi0 xxx
#'
#' @return A plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_prot, package = "DAPARdata")
#' obj <- Exp1_R25_prot[1:1000]
#' level <- obj@experimentData@other$typeOfData
#' metacell.mask <- match.metacell(GetMetacell(obj), "missing", level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj <- MetaCellFiltering(obj, indices, cmd = "delete")
#' qData <- Biobase::exprs(obj$new)
#' sTab <- Biobase::pData(obj$new)
#' allComp <- limmaCompleteTest(qData, sTab)
#' histPValue_HC(allComp$P_Value[1])
#'
#' @export
#'
histPValue_HC <- function(pval_ll, bins = 80, pi0 = 1) {
    h <- hist(sort(unlist(pval_ll)), freq = F, breaks = bins)

    serieInf <- sapply(h$density, function(x) min(pi0, x))
    serieSup <- sapply(h$density, function(x) max(0, x - pi0))

    hc <- highchart() %>%
        hc_chart(type = "column") %>%
        hc_add_series(data = serieSup, name = "p-value density") %>%
        hc_add_series(data = serieInf, name = "p-value density") %>%
        hc_title(text = "P-value histogram") %>%
        hc_legend(enabled = FALSE) %>%
        hc_colors(c("green", "red")) %>%
        hc_xAxis(title = list(text = "P-value"), categories = h$breaks) %>%
        hc_yAxis(
            title = list(text = "Density"),
            plotLines = list(
                list(
                    color = "blue", 
                    width = 2, 
                    value = pi0, 
                    zIndex = 5)
                )
        ) %>%
        hc_tooltip(
            headerFormat = "",
            pointFormat = "<b> {series.name} </b>: {point.y} ",
            valueDecimals = 2
        ) %>%
        my_hc_ExportMenu(filename = "histPVal") %>%
        hc_plotOptions(
            column = list(
                groupPadding = 0,
                pointPadding = 0,
                borderWidth = 0
            ),
            series = list(
                stacking = "normal",
                animation = list(duration = 100),
                connectNulls = TRUE,
                marker = list(enabled = FALSE)
            )
        ) %>%
        hc_add_annotation(
            labelOptions = list(
                backgroundColor = "transparent",
                verticalAlign = "top",
                y = -30,
                borderWidth = 0,
                x = 20,
                style = list(
                    fontSize = "1.5em",
                    color = "blue"
                )
            ),
            labels = list(
                list(
                    point = list(
                        xAxis = 0,
                        yAxis = 0,
                        x = 80,
                        y = pi0
                    ),
                    text = paste0("pi0=", pi0)
                )
            )
        )
    return(hc)
}
