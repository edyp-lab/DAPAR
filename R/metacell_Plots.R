

#' @title Bar plot of missing values per lines using highcharter
#' 
#' @description 
#' This method plots a bar plot which represents the distribution of the
#' number of missing values (NA) per lines (ie proteins).
#'
#' @param obj xxx.
#' @param pattern xxx
#' @param detailed 'value' or 'percent'
#' @param indLegend The indice of the column name's in \code{Biobase::pData()} 
#' tab
#' @param showValues A logical that indicates whether numeric values should be
#' drawn above the bars.
#' @return A bar plot
#' @author Florence Combes, Samuel Wieczorek
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept[seq_len(10), ]
#' metacellPerLinesHisto_HC(obj, pattern = "Missing")
#'
#' @export
#'
metacellPerLinesHisto_HC <- function(obj,
    pattern,
    detailed = FALSE,
    indLegend = "auto",
    showValues = FALSE) {
    if (missing(obj)) {
        stop("'obj' is missing.")
    } else if (is.null(obj)) {
        stop("'obj' is NULL. Abort...")
    }
    if (missing(pattern)) {
        stop("'pattern' is missing.")
    } else if (pattern %in% c("", "None")) {
        warning("'pattern' is empty.")
        return(NULL)
    }

    qData <- Biobase::exprs(obj)
    samplesData <- Biobase::pData(obj)

    if (identical(indLegend, "auto")) {
        indLegend <- seq.int(from = 2, to = length(colnames(samplesData)))
    }


    for (j in seq_len(length(colnames(qData)))) {
        noms <- NULL
        for (i in seq_len(length(indLegend))) {
            noms <- paste(noms, samplesData[j, indLegend[i]], sep = " ")
        }
        colnames(qData)[j] <- noms
    }


    mask <- match.metacell(GetMetacell(obj),
        pattern = pattern,
        level = obj@experimentData@other$typeOfData
    )
    NbNAPerRow <- rowSums(mask)

    nb.col <- dim(qData)[2]
    nb.na <- NbNAPerRow
    temp <- table(NbNAPerRow)
    nb.na2barplot <- rep(0, ncol(qData))

    for (i in seq_len(length(temp))) {
        nb.na2barplot[as.integer(names(temp)[i])] <- temp[i]
    }


    df <- data.frame(
        y = nb.na2barplot,
        y_percent = round(100 * nb.na2barplot / dim(qData)[1], digits = 2)
    )

    myColors <- rep("lightgrey", nrow(df))

    h1 <- highchart() %>%
        hc_title(text = paste0("Nb of lines with x '", pattern, "' tags")) %>%
        hc_add_series(data = df, type = "column", colorByPoint = TRUE) %>%
        hc_colors(myColors) %>%
        hc_plotOptions(
            column = list(stacking = "normal"),
            animation = list(duration = 100)
        ) %>%
        hc_legend(enabled = FALSE) %>%
        hc_xAxis(categories = row.names(df), 
            title = list(
                text = paste0("Nb of '", pattern, "' tags in a line")
                )
            ) %>%
        my_hc_ExportMenu(filename = "missingValuesPlot1") %>%
        hc_tooltip(
            enabled = TRUE,
            headerFormat = "",
            pointFormat = paste0("{point.y} lines<br>
                ({point.y_percent}% of all lines)")
        )

    return(h1)
}






#' @title Bar plot of missing values per lines and per condition
#' 
#' @description 
#' This method plots a bar plot which represents the distribution of the
#' number of missing values (NA) per lines (ie proteins) and per conditions.
#'
#' @param obj xxx
#'
#' @param pattern xxx
#'
#' @param indLegend The indice of the column name's in \code{Biobase::pData()} 
#' tab
#'
#' @param showValues A logical that indicates wether numeric values should be
#' drawn above the bars.
#'
#' @param pal xxx
#'
#' @return A bar plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept
#' pal <- ExtendPalette(length(unique(Biobase::pData(obj)$Condition)), "Dark2")
#' metacellPerLinesHistoPerCondition_HC(obj, "Missing", pal = pal)
#' metacellPerLinesHistoPerCondition_HC(obj, "Quantified")
#'
#' @export
#'
metacellPerLinesHistoPerCondition_HC <- function(obj,
    pattern,
    indLegend = "auto",
    showValues = FALSE,
    pal = NULL) {
    if (missing(obj)) {
        stop("'obj' is missing.")
    } else if (is.null(obj)) {
        stop("'obj' is NULL. Abort...")
    }
    if (missing(pattern)) {
        stop("'pattern' is missing.")
    } else if (pattern %in% c("", "None")) {
        warning("'pattern' is empty.")
        return(NULL)
    }

    qData <- Biobase::exprs(obj)
    samplesData <- Biobase::pData(obj)
    conds <- samplesData$Condition
    u_conds <- unique(conds)
    nbConditions <- length(u_conds)

    myColors <- NULL
    if (is.null(pal)) {
        warning("Color palette set to default.")
        myColors <- GetColorsForConditions(conds, 
            ExtendPalette(length(unique(conds))))
    } else {
        if (length(pal) != length(u_conds)) {
            warning("The color palette has not the same dimension as the 
                number of samples")
            myColors <- GetColorsForConditions(conds, 
                ExtendPalette(length(unique(conds))))
        } else {
            myColors <- pal
        }
    }

    if (identical(indLegend, "auto")) {
        indLegend <- seq.int(from = 2, to = length(colnames(samplesData)))
    }


    ncolMatrix <- max(unlist(lapply(
        u_conds,
        function(x) {
            length(which(conds == x))
        }
    )))


    mask <- match.metacell(GetMetacell(obj),
        pattern = pattern,
        level = obj@experimentData@other$typeOfData
    )
    ll.df <- list()
    for (i in u_conds)
    {
        df <- as.data.frame(matrix(rep(0, 2 * (1 + nbConditions)),
            nrow = 1 + nbConditions,
            dimnames = list(
                seq(seq.int(from=0, to=(nbConditions))),
                c("y", "y_percent")
            )
        ))
        rownames(df) <- seq.int(from = 0, to = (nrow(df) - 1))
        ll.df[[i]] <- df
        nSample <- length(which(conds == i))
        t <- NULL
        if (nSample == 1) {
            t <- table(as.integer(mask[, which(conds == i)]))
        } else {
            t <- table(rowSums(mask[, which(conds == i)]))
        }

        df[as.integer(names(t)) + 1, "y"] <- t
        df[as.integer(names(t)) + 1, "y_percent"] <- round(100 * t / nrow(obj),
            digits = 2)
        ll.df[[i]] <- df
    }

    h1 <- highchart() %>%
        hc_title(text = paste0("Nb of lines containing x '", 
            pattern, "' tags (condition-wise)")) %>%
        my_hc_chart(chartType = "column") %>%
        hc_plotOptions(
            column = list(stacking = ""),
            dataLabels = list(enabled = FALSE),
            animation = list(duration = 100)
        ) %>%
        hc_colors(unique(myColors)) %>%
        hc_legend(enabled = FALSE) %>%
        hc_xAxis(categories = seq.int(from=0, to=ncolMatrix), 
            title = list(text = paste0("Nb of '", pattern, 
                "' tags in each line (condition-wise)"))) %>%
        my_hc_ExportMenu(filename = "missingValuesPlot_2") %>%
        hc_tooltip(
            headerFormat = "",
            pointFormat = "{point.y} lines<br>({point.y_percent}% of all lines)"
        )

    for (i in seq_len(nbConditions)) {
        h1 <- h1 %>% hc_add_series(data = ll.df[[u_conds[i]]])
    }

    return(h1)
}





#' @title Histogram of missing values
#' @description 
#' #' This method plots a histogram of missing values. Same as the function
#' \code{mvHisto} but uses the package \code{highcharter}
#' 
#' @param obj xxx
#' @param pattern xxx
#' @param indLegend The indices of the column name's in \code{Biobase::pData()}
#'  tab
#' @param showValues A logical that indicates wether numeric values should be
#' drawn above the bars.
#' @param pal xxx
#' @return A histogram
#' @author Florence Combes, Samuel Wieczorek
#'
#' @import highcharter
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept
#' pattern <- "Missing POV"
#' pal <- ExtendPalette(2, "Dark2")
#' metacellHisto_HC(obj, pattern, showValues = TRUE, pal = pal)
#'
#' @export
#'
metacellHisto_HC <- function(obj,
    pattern,
    indLegend = "auto",
    showValues = FALSE,
    pal = NULL) {
    if (missing(obj)) {
        stop("'obj' is missing.")
    } else if (is.null(obj)) {
        stop("'obj' is NULL. Abort...")
    }
    if (missing(pattern)) {
        stop("'pattern' is missing.")
    } else if (pattern %in% c("", "None")) {
        warning("'pattern' is empty.")
        return(NULL)
    }

    qData <- Biobase::exprs(obj)
    samplesData <- Biobase::pData(obj)
    conds <- samplesData[, "Condition"]

    myColors <- NULL
    if (is.null(pal)) {
        warning("Color palette set to default.")
        myColors <- GetColorsForConditions(conds, 
            ExtendPalette(length(unique(conds))))
    } else {
        if (length(pal) != length(unique(conds))) {
            warning("The color palette has not the same dimension as the 
                number of samples")
            myColors <- GetColorsForConditions(conds, 
                ExtendPalette(length(unique(conds))))
        } else {
            myColors <- GetColorsForConditions(conds, pal)
        }
    }

    if (identical(indLegend, "auto")) {
        indLegend <- seq.int(from=2, to = length(colnames(samplesData)))
    }



    mask <- match.metacell(GetMetacell(obj),
        pattern = pattern,
        level = obj@experimentData@other$typeOfData
    )

    NbNAPerCol <- colSums(mask)

    df <- data.frame(
        y = NbNAPerCol,
        y_percent = round(100 * NbNAPerCol / nrow(mask), digits = 2)
    )



    h1 <- highchart() %>%
        my_hc_chart(chartType = "column") %>%
        hc_title(text = paste0("Nb of '", pattern, "' tags by replicate")) %>%
        hc_add_series(df, type = "column", colorByPoint = TRUE) %>%
        hc_colors(myColors) %>%
        hc_plotOptions(
            column = list(stacking = "normal"),
            animation = list(duration = 100)
        ) %>%
        hc_legend(enabled = FALSE) %>%
        hc_xAxis(categories = conds, title = list(text = "Replicates")) %>%
        my_hc_ExportMenu(filename = "missingValuesPlot_3") %>%
        hc_tooltip(
            headerFormat = "",
            pointFormat = "{point.y} lines<br>({point.y_percent}% of all lines)"
        )

    return(h1)
}





#' @title Heatmap of missing values from a \code{MSnSet} object
#' @description 
#' #' Plots a heatmap of the quantitative data. Each column represent one of
#' the conditions in the object of class \code{MSnSet} and
#' the color is proportional to the mean of intensity for each line of
#' the dataset.
#' The lines have been sorted in order to vizualize easily the different
#' number of missing values. A white square is plotted for missing values.
#' 
#' @param obj An object of class \code{MSnSet}.
#'
#' @param pattern xxx
#'
#' @return A heatmap
#' @author Alexia Dorffer
#' @examples
#' data(Exp1_R25_prot, package="DAPARdata")
#' obj <- Exp1_R25_prot[seq_len(1000)]
#' level <- 'protein'
#' metacell.mask <- match.metacell(GetMetacell(obj), "Missing", level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
#' obj <- MetaCellFiltering(obj, indices, cmd = "delete")
#' wrapper.mvImage(obj$new)
#'
#' @export
#'
#'
wrapper.mvImage <- function(obj, pattern = "Missing MEC") {
    if (missing(obj)) {
        stop("'obj' is required.")
    } else if (is.null(obj)) {
        warning("'obj' is NULL. Return NULL.")
        return(NULL)
    }
    qData <- Biobase::exprs(obj)
    conds <- Biobase::pData(obj)[, "Condition"]
    metac <- Biobase::fData(obj)[, obj@experimentData@other$names_metacell]
    level <- obj@experimentData@other$typeOfData
    indices <- which(apply(match.metacell(metac, pattern, level), 1, sum) > 0)

    if (length(indices) == 0) {
        warning("The dataset contains no Missing value on Entire Condition. 
            So this plot is not available.")
        return(NULL)
    } else if (length(indices) == 1) {
        warning("The dataset contains only one Missing value on Entire 
        Condition. Currently, Prostar does not handle such dataset to build 
        the plot. As it has no side-effects on the results, you can continue 
        your imputation.")
        return(NULL)
    }

    mvImage(qData[indices, ], conds)
}





#' @title Heatmap of missing values
#' @description 
#' #' Plots a heatmap of the quantitative data. Each column represent one of
#' the conditions in the object of class \code{MSnSet} and
#' the color is proportional to the mean of intensity for each line of
#' the dataset.
#' The lines have been sorted in order to vizualize easily the different
#' number of missing values. A white square is plotted for missing values.
#' 
#' @param qData A dataframe that contains quantitative data.
#' @param conds A vector of the conditions (one condition per sample).
#' @return A heatmap
#' @author Samuel Wieczorek, Thomas Burger
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' conds <- Biobase::pData(Exp1_R25_pept)[, "Condition"]
#' mvImage(qData, conds)
#'
#' @export
#'
#'
mvImage <- function(qData, conds) {

    pkgs.require(c('grDevices', 'stats'))
    
    ### build indices of conditions
    indCond <- list()
    ConditionNames <- unique(conds)
    for (i in ConditionNames) {
        indCond <- append(indCond, list(which(i == conds)))
    }
    indCond <- stats::setNames(indCond, as.list(c("cond1", "cond2")))

    nNA1 <- apply(as.matrix(qData[, indCond$cond1]), 1, 
        function(x) sum(is.na(x)))
    nNA2 <- apply(as.matrix(qData[, indCond$cond2]), 1, 
        function(x) sum(is.na(x)))
    o <- order(((nNA1 + 1)^2) / (nNA2 + 1))
    exprso <- qData[o, ]

    for (i in seq_len(nrow(exprso))) {
        k <- order(exprso[i, indCond$cond1])
        exprso[i, rev(indCond$cond1)] <- exprso[i, k]
        .temp <- mean(exprso[i, rev(indCond$cond1)], na.rm = TRUE)
        exprso[i, which(!is.na(exprso[i, indCond$cond1]))] <- .temp

        k <- order(exprso[i, indCond$cond2])
        exprso[i, indCond$cond2] <- exprso[i, k + length(indCond$cond1)]
        .temp <- mean(exprso[i, indCond$cond2], na.rm = TRUE)
        exprso[i, length(indCond$cond1) +
            which(!is.na(exprso[i, indCond$cond2]))] <- .temp
    }


    heatmapForMissingValues(exprso,
        col = grDevices::colorRampPalette(c("yellow", "red"))(100),
        key = TRUE,
        srtCol = 0,
        labCol = conds,
        ylab = "Peptides / proteins",
        main = "MEC heatmap"
    )

    # heatmap_HC(exprso,col = colfunc(100),labCol=conds)
}



#' @title Distribution of Observed values with respect to intensity values
#' 
#' @description 
#' This method shows density plots which represents the repartition of
#' Partial Observed Values for each replicate in the dataset.
#' The colors correspond to the different conditions (slot Condition in in the
#' dataset of class \code{MSnSet}).
#' The x-axis represent the mean of intensity for one condition and one
#' entity in the dataset (i. e. a protein)
#' whereas the y-axis count the number of observed values for this entity
#' and the considered condition.
#'
#' @param obj xxx
#'
#' @param pal The different colors for conditions
#'
#' @param pattern xxx
#'
#' @param typeofMV xxx
#'
#' @param title The title of the plot
#'
#' @import highcharter
#'
#' @return Density plots
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept[seq_len(100)]
#' conds <- Biobase::pData(obj)$Condition
#' pal <- ExtendPalette(length(unique(conds)), "Dark2")
#' hc_mvTypePlot2(obj, pattern = "Missing MEC", title = "POV distribution", 
#' pal = pal)
#'
#' @import highcharter
#'
#' @export
#'
hc_mvTypePlot2 <- function(obj,
    pal = NULL,
    pattern,
    typeofMV = NULL,
    title = NULL) {
    pkgs.require('stats')
    
    conds <- Biobase::pData(obj)[, "Condition"]
    qData <- Biobase::exprs(obj)
    myColors <- NULL
    if (is.null(pal)) {
        warning("Color palette set to default.")
        pal <- ExtendPalette(length(unique(conds)))
    } else {
        if (length(pal) != length(unique(conds))) {
            warning("The color palette has not the same dimension as the 
                number of samples")
            pal <- ExtendPalette(length(unique(conds)))
        }
    }

    conditions <- conds
    mTemp <- nbNA <- nbValues <- matrix(
        rep(0, nrow(qData) * length(unique(conditions))
            ),
        nrow = nrow(qData),
        dimnames = list(NULL, unique(conditions))
    )
    dataCond <- data.frame()
    ymax <- 0
    series <- list()
    myColors <- NULL
    j <- 1

    level <- obj@experimentData@other$typeOfData

    for (iCond in unique(conditions)) {
        if (length(which(conditions == iCond)) == 1) {
            mTemp[, iCond] <- qData[, which(conditions == iCond)]
            nbNA[, iCond] <- as.integer(
                match.metacell(GetMetacell(obj)[, which(conditions == iCond)],
                pattern = pattern,
                level = level
            ))
            
            .op1 <- length(which(conditions == iCond))
            .op2 <- nbNA[, iCond]
            nbValues[, iCond] <- .op1 - .op2
        } else {
            .qcond <- which(conditions == iCond)
            mTemp[, iCond] <- apply(qData[, .qcond], 1, mean, na.rm = TRUE)
            nbNA[, iCond] <- rowSums(
                match.metacell(GetMetacell(obj)[, .qcond],
                    pattern = pattern,
                    level = level
            ))
            nbValues[, iCond] <- length(.qcond) - nbNA[, iCond]
        }


        for (i in seq_len(length(which(conditions == iCond)))) {
            data <- mTemp[which(nbValues[, iCond] == i), iCond]
            tmp <- NULL
            if (length(data) >= 2) {
                tmp <- stats::density(mTemp[which(nbValues[, iCond] == i), 
                    iCond])
                tmp$y <- tmp$y + i
                if (max(tmp$y) > ymax) {
                    ymax <- max(tmp$y)
                }
            }
            series[[j]] <- tmp
            myColors <- c(myColors, pal[which(unique(conditions) == iCond)])
            j <- j + 1
        }
    }


    hc <- highchart(type = "chart") %>%
        hc_title(text = title) %>%
        my_hc_chart(chartType = "spline", zoomType = "xy") %>%
        hc_legend(
            align = "left", verticalAlign = "top",
            layout = "vertical"
        ) %>%
        hc_xAxis(title = list(text = "Mean of intensities")) %>%
        hc_yAxis(
            title = list(text = "Number of quantity values per condition"),
            tickInterval = 0.5
        ) %>%
        hc_tooltip(
            headerFormat = "",
            pointFormat = "<b> {series.name} </b>: {point.y} ",
            valueDecimals = 2
        ) %>%
        my_hc_ExportMenu(filename = paste0(pattern, "_distribution")) %>%
        hc_plotOptions(
            series = list(
                showInLegend = TRUE,
                animation = list(
                    duration = 100
                ),
                connectNulls = TRUE,
                marker = list(
                    enabled = FALSE
                )
            )
        )

    for (i in seq_len(length(series))) {
        hc <- hc_add_series(hc,
            data = list_parse(data.frame(cbind(
                x = series[[i]]$x,
                y = series[[i]]$y
            ))),
            showInLegend = FALSE,
            color = myColors[i],
            name = conds[i]
        )
    }

    # add three empty series for the legend entries. Change color and marker 
    # symbol
    for (c in seq_len(length(unique(conds)))) {
        hc <- hc_add_series(hc,
            data = data.frame(),
            name = unique(conds)[c],
            color = pal[c],
            marker = list(symbol = "circle"),
            type = "line"
        )
    }

    hc
    return(hc)
}
