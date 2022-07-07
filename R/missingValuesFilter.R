#' @title Percentage of missing values
#'
#' @description
#' Returns the percentage of missing values in the quantitative
#' data (\code{Biobase::exprs()} table of the dataset).
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @return A floating number
#'
#' @author Florence Combes, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' getPourcentageOfMV(Exp1_R25_pept[seq_len(100), ])
#'
#' @export
#'
#'
getPourcentageOfMV <- function(obj) {
    df <- data.frame(Biobase::exprs(obj))

    NA.count <- apply(
        df, 2,
        function(x) length(which(is.na(data.frame(x)) == TRUE))
    )


    pourcentage <- 100 * round(sum(NA.count) / (nrow(df) * ncol(df)), 
        digits = 4)

    return(pourcentage)
}


#' @title Number of lines with prefix
#' 
#' @description 
#' Returns the number of lines, in a given column, where content matches
#' the prefix.
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @param name The name of a column.
#'
#' @param prefix A string
#'
#' @return An integer
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' getNumberOf(Exp1_R25_pept[seq_len(100)], "Potential_contaminant", "+")
#'
#' @export
#'
#'
getNumberOf <- function(obj, name = NULL, prefix = NULL) {
    if (is.null(name) || is.null(prefix) || (name == "") || (prefix == "")) {
        return(0)
    }
    if (!(is.null(name) || !is.null(name == "")) &&
        (is.null(prefix) || (prefix == ""))) {
        return(0)
    }

    if (nchar(prefix) > 0) {
        count <- length(
            which(substr(Biobase::fData(obj)[, name], 0, 1) == prefix))
    } else {
        count <- 0
    }

    return(count)
}




#' @title Removes lines in the dataset based on numerical conditions.
#' 
#' @description 
#' This function removes lines in the dataset based on numerical conditions.
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @param name The name of the column that correspond to the line to filter
#'
#' @param value A number
#'
#' @param operator A string
#'
#' @return An list of 2 items :
#' * obj : an object of class \code{MSnSet} in which the lines have been 
#' deleted,
#' * deleted : an object of class \code{MSnSet} which contains the deleted lines
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' NumericalFiltering(Exp1_R25_pept[seq_len(100)], "A_Count", "6", "==")
#'
#' @export
#'
NumericalFiltering <- function(
        obj, 
    name = NULL, 
    value = NULL, 
    operator = NULL) {
    if ((is.null(name) || (name == ""))) {
        return(NULL)
    }

    deleted <- NULL
    ind <- NULL
    ind <- NumericalgetIndicesOfLinesToRemove(obj, name, value, operator)

    if (!is.null(ind) && (length(ind) > 0)) {
        deleted <- obj[ind]

        obj <- deleteLinesFromIndices(
            obj, ind,
            paste("\"",
                length(ind),
                " lines were removed from dataset.\"",
                sep = ""
            )
        )
    }

    return(
        list(
            obj = obj, 
            deleted = deleted
            )
        )
}





#'
#' @title Get the indices of the lines to delete, based on a prefix string
#' 
#' @description 
#' This function returns the indices of the lines to delete, based on a
#' prefix string
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @param name The name of the column that correspond to the data to filter
#'
#' @param value xxxx
#'
#' @param operator A xxxx
#'
#' @return A vector of integers.
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' NumericalgetIndicesOfLinesToRemove(Exp1_R25_pept[seq_len(100)], "A_Count",
#' value = "6", operator = "==")
#'
#' @export
#'
#'
NumericalgetIndicesOfLinesToRemove <- function(
    obj, 
    name = NULL, 
    value = NULL, 
    operator = NULL
    ) {
    if ((value == "") || is.null(value) || 
            (operator == "") || is.null(operator)) {
        # warning ("No change was made")
        return(NULL)
    }

    data <- Biobase::fData(obj)[, name]
    ind <- which(eval(parse(text = paste0("data", operator, value))))

    return(ind)
}




#' @title Barplot of proportion of contaminants and reverse
#' 
#' @description 
#' Plots a barplot of proportion of contaminants and reverse. Same as the
#' function \code{proportionConRev} but uses the package \code{highcharter}
#'
#'
#' @param nBoth The number of both contaminants and reverse identified in
#' the dataset.
#'
#' @param nCont The number of contaminants identified in the dataset.
#'
#' @param nRev The number of reverse entities identified in the dataset.
#'
#' @param lDataset The total length (number of rows) of the dataset
#'
#' @return A barplot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' proportionConRev_HC(10, 20, 100)
#'
#' @export
#'
proportionConRev_HC <- function(nBoth = 0, nCont = 0, nRev = 0, lDataset = 0) {
    if (is.null(nCont) && is.null(nBoth) && 
            is.null(nRev) && is.null(lDataset)) {
        return(NULL)
    }

    total <- nBoth + nCont + nRev + lDataset
    pctGood <- 100 * round(lDataset / total, digits = 4)
    pctBoth <- 100 * round(nBoth / total, digits = 4)
    pctContaminants <- 100 * round(nCont / total, digits = 4)
    pctReverse <- 100 * round(nRev / total, digits = 4)

    counts <- c(lDataset, nCont, nRev, nBoth)
    slices <- c(pctGood, pctContaminants, pctReverse, pctBoth)
    lbls <- c("Quantitative data", "Contaminants", 
        "Reverse", "Both contaminants & Reverse")
    # pct <- c(pctGood, pctContaminants, pctReverse  ,pctBoth)
    lbls <- paste(lbls, " (", counts, " lines)", sep = "")

    mydata <- data.frame(
        test = c(pctGood, pctContaminants, pctReverse, pctBoth)
        )

    highchart() %>%
        my_hc_chart(chartType = "bar") %>%
        hc_yAxis(title = list(text = "Pourcentage")) %>%
        hc_xAxis(categories = lbls) %>%
        hc_legend(enabled = FALSE) %>%
        hc_plotOptions(column = list(
            dataLabels = list(enabled = TRUE),
            stacking = "normal",
            enableMouseTracking = FALSE
        )) %>%
        hc_add_series(
            data = mydata$test,
            dataLabels = list(enabled = TRUE, format = "{point.y}%"),
            colorByPoint = TRUE
        ) %>%
        my_hc_ExportMenu(filename = "contaminants")
}




#' @title Removes lines in the dataset based on a prefix string.
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @param idLine2Delete The name of the column that correspond to the
#' data to filter
#'
#' @param prefix A character string that is the prefix to find in the data
#' @return An object of class \code{MSnSet}.
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' removeLines(Exp1_R25_pept[seq_len(100)], "Potential_contaminant")
#' removeLines(Exp1_R25_pept[seq_len(100)], "Reverse")
#'
#' @export
#'
removeLines <- function(obj, idLine2Delete = NULL, prefix = NULL) {
    if ((prefix == "") || is.null(prefix)) {
        # warning ("No change was made")
        return(obj)
    }
    t <- (prefix == substring(Biobase::fData(obj)[, idLine2Delete], 1, 
        nchar(prefix)))
    ind <- which(t == TRUE)
    obj <- obj[-ind]

    return(obj)
}




#' @title Removes lines in the dataset based on a prefix strings (contaminants,
#' reverse or both).
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @param idCont2Delete The name of the column that correspond to the
#' contaminants to filter
#'
#' @param prefix_Cont A character string that is the prefix for the
#' contaminants to find in the data
#'
#' @param idRev2Delete The name of the column that correspond to the
#' reverse data to filter
#'
#' @param prefix_Rev A character string that is the prefix for the reverse to
#' find in the data
#'
#' @return An list of 4 items :
#' * obj : an object of class \code{MSnSet} in which the lines have been deleted
#' * deleted.both : an object of class \code{MSnSet} which contains the deleted
#' lines corresponding to both contaminants and reverse,
#' * deleted.contaminants : n object of class \code{MSnSet} which contains the
#' deleted lines corresponding to contaminants,
#' * deleted.reverse : an object of class \code{MSnSet} which contains the
#' deleted lines corresponding to reverse,
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' StringBasedFiltering(
#' Exp1_R25_pept[seq_len(100)], "Potential_contaminant", "+", "Reverse", "+")
#'
#' @export
#'
StringBasedFiltering <- function(obj,
    idCont2Delete = NULL, 
    prefix_Cont = NULL,
    idRev2Delete = NULL, 
    prefix_Rev = NULL
    ) {
    deleted.both <- deleted.contaminants <- deleted.reverse <- NULL

    ##
    ## Search for both
    ##
    if ((!is.null(idCont2Delete) || (idCont2Delete != "")) &&
        (!is.null(idRev2Delete) || (idRev2Delete != ""))) {
        indContaminants <- indReverse <- indBoth <- NULL
        indContaminants <- getIndicesOfLinesToRemove(obj, 
            idCont2Delete, 
            prefix_Cont)
        indReverse <- getIndicesOfLinesToRemove(obj, idRev2Delete, prefix_Rev)
        indBoth <- intersect(indContaminants, indReverse)

        if (!is.null(indBoth) && (length(indBoth) > 0)) {
            deleted.both <- obj[indBoth]
            obj <- deleteLinesFromIndices(
                obj, indBoth,
                paste("\"",
                    length(indBoth),
                    " both contaminants and reverse were removed from 
                    dataset.\"",
                    sep = ""
                )
            )
        }
    }

    ##
    ## Search for contaminants
    ##
    if ((!is.null(idCont2Delete) || (idCont2Delete != ""))) {
        indContaminants <- NULL
        indContaminants <- getIndicesOfLinesToRemove(obj, 
            idCont2Delete, 
            prefix_Cont)

        if (!is.null(indContaminants) && (length(indContaminants) > 0)) {
            deleted.contaminants <- obj[indContaminants]

            obj <- deleteLinesFromIndices(
                obj, indContaminants,
                paste("\"",
                    length(indContaminants),
                    " contaminants were removed from dataset.\"",
                    sep = ""
                )
            )
        }
    }


    ##
    ## Search for reverse
    ##
    if ((!is.null(idRev2Delete) || (idRev2Delete != ""))) {
        indReverse <- getIndicesOfLinesToRemove(obj, idRev2Delete, prefix_Rev)

        if (!is.null(indReverse)) {
            if (length(indReverse) > 0) {
                deleted.reverse <- obj[indReverse]

                obj <- deleteLinesFromIndices(
                    obj, indReverse,
                    paste("\"",
                        length(indReverse),
                        " reverse were removed from dataset.\"",
                        sep = ""
                    )
                )
            }
        }
    }


    return(list(
        obj = obj,
        deleted.both = deleted.both,
        deleted.contaminants = deleted.contaminants,
        deleted.reverse = deleted.reverse
    ))
}



#' @title Removes lines in the dataset based on a prefix strings.
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @param cname The name of the column that correspond to the line to filter
#'
#' @param tag A character string that is the prefix for the contaminants to
#' find in the data
#'
#' @return An list of 4 items :
#' * obj : an object of class \code{MSnSet} in which the lines have been deleted
#' * deleted : an object of class \code{MSnSet} which contains the deleted lines
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj.filter <- StringBasedFiltering2(Exp1_R25_pept[seq_len(100)], 
#' "Potential_contaminant", "+")
#'
#' @export
#'
StringBasedFiltering2 <- function(obj, cname = NULL, tag = NULL) {
    deleted <- NULL

    ##
    ## Search for contaminants
    ##
    if ((!is.null(cname) || (cname != ""))) {
        ind <- NULL
        ind <- getIndicesOfLinesToRemove(obj, cname, tag)

        if (!is.null(ind) && (length(ind) > 0)) {
            deleted <- obj[ind]

            obj <- deleteLinesFromIndices(
                obj, ind,
                paste("\"",
                    length(ind),
                    " contaminants were removed from dataset.\"",
                    sep = ""
                )
            )
        }
    }

    return(list(obj = obj, deleted = deleted))
}





#' @title Get the indices of the lines to delete, based on a prefix string
#'
#' @param obj An object of class \code{MSnSet}.
#'
#' @param idLine2Delete The name of the column that correspond to the data
#' to filter
#'
#' @param prefix A character string that is the prefix to find in the data
#'
#' @return A vector of integers.
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' ind <- getIndicesOfLinesToRemove(Exp1_R25_pept[seq_len(100)], 
#' "Potential_contaminant",
#'     prefix = "+"
#' )
#'
#' @export
#'
#'
getIndicesOfLinesToRemove <- function(
        obj, 
    idLine2Delete = NULL, 
    prefix = NULL) {
    if ((prefix == "") || is.null(prefix)) {
        # warning ("No change was made")
        return(NULL)
    }
    t <- (prefix == substring(
        Biobase::fData(obj)[, idLine2Delete], 1, nchar(prefix)))
    ind <- which(t == TRUE)
    return(ind)
}





#' @title Filter lines in the matrix of intensities w.r.t. some criteria
#' 
#' @description 
#' #' Filters the lines of \code{Biobase::exprs()} table with conditions on the 
#' number of missing values.
#' The user chooses the minimum amount of intensities that is acceptable and
#' the filter delete lines that do not respect this condition.
#' The condition may be on the whole line or condition by condition.
#'
#' The different methods are :
#' "WholeMatrix": given a threshold \code{th}, only the lines that contain
#' at least \code{th} values are kept.
#' "AllCond": given a threshold \code{th}, only the lines which contain
#' at least \code{th} values for each of the conditions are kept.
#' "AtLeastOneCond": given a threshold \code{th}, only the lines that contain
#' at least \code{th} values, and for at least one condition, are kept.
#'
#'
#' @param obj An object of class \code{MSnSet} containing
#' quantitative data.
#'
#' @param indices A vector of integers which are the indices of lines to
#' keep.
#'
#' @param cmd xxxx. Available values are: 'delete', 'keep'.
#'
#' @param processText A string to be included in the \code{MSnSet}
#' object for log.
#'
#' @return An instance of class \code{MSnSet} that have been filtered.
#'
#' @author Florence Combes, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept[seq_len(100)]
#' level <- 'peptide'
#' metacell.mask <- match.metacell(GetMetacell(obj), "missing", level)
#' indices <- GetIndices_WholeLine(metacell.mask)
#' obj.filter <- MetaCellFiltering(obj, indices, "delete")
#'
#' @export
#'
MetaCellFiltering <- function(obj,
    indices,
    cmd,
    processText = "") {
    if (missing(obj)) {
        stop("'obj' is required;")
    }
    if (missing(indices)) {
        stop("'indices' is required;")
    }
    if (missing(cmd)) {
        stop("'cmd' is required;")
    } else if (!(cmd %in% c("delete", "keep"))) {
        stop("'cmd' must be one of the following values: `delete` or `keep`.")
    }



    if (is.null(indices)) {
        warning("'indices' is NULL. No filtering will be process.")
        deleted <- obj[-c(seq_len(nrow(obj)))]
        new <- obj
    } else if (cmd == "delete") {
        deleted <- obj[indices]
        new <- obj[-indices]
    } else if (cmd == "keep") {
        deleted <- obj[-indices]
        new <- obj[indices]
    }

    new@processingData@processing <-
        c(new@processingData@processing, processText)

    return(list(
        new = new,
        deleted = deleted
    ))
}



#' @title Delete the lines in the matrix of intensities and the metadata table
#' given their indice.
#'
#' @param obj An object of class \code{MSnSet} containing
#' quantitative data.
#'
#' @param deleteThat A vector of integers which are the indices of lines to
#' delete.
#'
#' @param processText A string to be included in the \code{MSnSet}
#' object for log.
#'
#' @return An instance of class \code{MSnSet} that have been filtered.
#'
#' @author Florence Combes, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- deleteLinesFromIndices(Exp1_R25_pept[seq_len(100)], c(seq_len(10)))
#'
#' @export
#'
deleteLinesFromIndices <- function(obj, deleteThat = NULL, processText = "") {
    if (is.null(deleteThat)) {
        return(obj)
    }
    obj <- obj[-deleteThat]

    obj@processingData@processing <- c(obj@processingData@processing, 
        processText)
    if (grepl("contaminants", processText)) {
        obj@experimentData@other$contaminantsRemoved <- TRUE
    }
    if (grepl("reverse", processText)) {
        obj@experimentData@other$reverseRemoved <- TRUE
    }
    return(obj)
}



#' @title Delete the lines in the matrix of intensities and the metadata table
#' given their indice.
#'
#' @param obj An object of class \code{MSnSet} containing
#' quantitative data.
#'
#' @param level A vector of integers which are the indices of lines to
#' delete.
#'
#' @param pattern A string to be included in the \code{MSnSet}
#' object for log.
#'
#' @param type xxx
#'
#' @param percent xxx
#'
#' @param op xxx
#'
#' @param th xxx
#'
#' @return An instance of class \code{MSnSet} that have been filtered.
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept[seq_len(10), ]
#' level <- GetTypeofData(obj)
#' pattern <- "missing"
#' type <- "AllCond"
#' percent <- FALSE
#' op <- "=="
#' th <- 2
#' indices <- GetIndices_MetacellFiltering(obj, level, pattern, type, 
#' percent, op, th)
#'
#' @export
#'
GetIndices_MetacellFiltering <- function(
        obj, 
    level, 
    pattern, 
    type, 
    percent, 
    op, 
    th) {
    if (missing(obj)) {
        stop("'obj' is required.")
    }
    if (missing(level)) {
        stop("'level' is required.")
    }
    if (missing(pattern)) {
        stop("'pattern' is required.")
    }
    if (missing(type)) {
        tsop("'type' is required.")
    }
    if (missing(percent)) {
        stop("'percent' is required.")
    }
    if (missing(op)) {
        stop("'op' is required.")
    }
    if (missing(th)) {
        stop("'th' is required.")
    }


    indices <- NULL

    if (!(pattern %in% metacell.def(level)$node && 
            type != "None" && !is.null(type))) {
        warning("Either 'pattern' nor 'type' are equal to 'None'")
        return(NULL)
    }

    mask <- match.metacell(
        metadata = GetMetacell(obj),
        pattern = pattern,
        level = level
    )

    indices <- switch(type,
        WholeLine = GetIndices_WholeLine(metacell.mask = mask),
        WholeMatrix = GetIndices_WholeMatrix(
            metacell.mask = mask,
            op = op,
            percent = percent,
            th = th
        ),
        AllCond = GetIndices_BasedOnConditions(
            metacell.mask = mask,
            type = type,
            conds = Biobase::pData(obj)$Condition,
            percent = percent,
            op = op,
            th = th
        ),
        AtLeastOneCond = GetIndices_BasedOnConditions(
            metacell.mask = mask,
            type = type,
            conds = Biobase::pData(obj)$Condition,
            percent = percent,
            op = op,
            th = th
        )
    )

    return(indices)
}



#' @title
#' Lists the metacell scopes for filtering
#'
#' @export
#' 
#' @return xxx
#' 
#' @examples 
#' MetacellFilteringScope()
#'
MetacellFilteringScope <- function() {
    c("None", "WholeLine", "WholeMatrix", "AllCond", "AtLeastOneCond")
}



#' @title xxx
#'
#' @export
#' 
#' @return A `character()`
#' 
#' @examples 
#' SymFilteringOperators()
#'
SymFilteringOperators <- function() {
    c("<=", "<", ">=", ">", "==", "!=")
}


#' @title
#' Search lines which respects request on one or more conditions.
#'
#' @description
#' This function looks for the lines that respect the request in either all 
#' conditions or at least one condition.
#'
#' @param metacell.mask xxx
#'
#' @param op  String for operator to use. List of operators is available with 
#' SymFilteringOperators().
#'
#' @param percent A boolean to indicate whether the threshold represent an 
#' absolute value (percent = FALSE) or
#' a percentage (percent=TRUE).
#'
#' @param th A floating number which is in the interval [0, 1]
#'
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' level <- 'peptide'
#' pattern <- "missing"
#' metacell.mask <- match.metacell(metadata = GetMetacell(obj), 
#' pattern = pattern, level = level)
#' percent <- FALSE
#' th <- 3
#' op <- ">="
#' ind <- GetIndices_WholeMatrix(metacell.mask, op, percent, th)
#'
#' @export
#' 
#' @return xxx
#'
GetIndices_WholeMatrix <- function(metacell.mask,
    op = "==",
    percent = FALSE,
    th = 0) {

    # Check parameters
    if (missing(metacell.mask)) {
        stop("'metacell.mask' is required.")
    }
    if (isTRUE(percent)) {
        if (th < 0 || th > 1) {
            warning("With percent=TRUE, the threshold 'th' must be in the 
                interval [0, 1].")
            return(NULL)
        }
    } else {
        th.upbound <- ncol(metacell.mask)
        if (th > th.upbound) {
            warn.txt <- paste0(
                "Param `th` is not correct. It must be an integer greater 
                than or equal to 0 and less or equal than ",
                th.upbound
            )
            warning(warn.txt)
            return(NULL)
        }
    }

    if (!(op %in% SymFilteringOperators())) {
        warning(paste0(
            "'op' must be one of the following values: ",
            paste0(SymFilteringOperators(), collapse = " ")
        ))
        return(NULL)
    }

    indices <- NULL
    if (isTRUE(percent)) {
        inter <- rowSums(metacell.mask) / ncol(metacell.mask)
        indices <- which(eval(parse(text = paste0("inter", op, th))))
    } else {
        inter <- apply(metacell.mask, 1, sum)
        indices <- which(eval(parse(text = paste0("inter", op, th))))
    }


    if (length(indices) == 0) indices <- NULL
    return(indices)
}


#' @title
#' Search lines which respects query on all their elements.
#'
#' @description
#' This function looks for the lines where each element respect the query.
#'
#' @param metacell.mask xxx
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept[seq.int(from=20, to=30)]
#' level <- 'peptide'
#' pattern <- "missing POV"
#' metacell.mask <- match.metacell(metadata = GetMetacell(obj), 
#' pattern = pattern, level = level)
#' ind <- GetIndices_WholeLine(metacell.mask)
#'
#' @export
#' 
#' @return xxx
#'
GetIndices_WholeLine <- function(metacell.mask) {
    if (missing(metacell.mask)) {
        stop("'metacell.mask' is missing.")
    }

    indices <- which(rowSums(metacell.mask) == ncol(metacell.mask))
    if (length(indices) == 0) indices <- NULL
    return(indices)
}


#' @title
#' Search lines which respects request on one or more conditions.
#'
#' @description
#' This function looks for the lines that respect the request in either 
#' all conditions
#' or at least one condition.
#'
#' @param metacell.mask xxx
#'
#' @param type Available values are:
#' * 'AllCond' (the query is valid in all the conditions),
#' * 'AtLeatOneCond' (the query is valid in at leat one condition.
#'
#' @param conds xxx
#'
#' @param percent xxx
#'
#' @param op  String for operator to use. List of operators is available 
#' with SymFilteringOperators().
#'
#' @param th The theshold to apply
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' level <- GetTypeofData(obj)
#' pattern <- 'missing'
#' metacell.mask <- match.metacell(metadata=GetMetacell(obj), 
#' pattern=pattern, level=level)
#' type <- 'AllCond'
#' conds <- Biobase::pData(obj)$Condition
#' op <- '>='
#' th <- 0.5
#' percent <- TRUE
#' ind <- GetIndices_BasedOnConditions(metacell.mask, type, conds, 
#' percent, op, th)
#'
#' @return xxx
#'
#' @export
#'
GetIndices_BasedOnConditions <- function(metacell.mask,
    type,
    conds,
    percent,
    op,
    th) {

    # Check parameters
    if (missing(metacell.mask)) {
        stop("'metacell.mask' is missing.")
    }
    if (missing(conds)) {
        stop("'conds' is missing.")
    }
    if (missing(type)) {
        stop("'type' is missing.")
    } else if (!(type %in% c("AllCond", "AtLeastOneCond"))) {
        stop("'type' must be one of the following: AllCond, AtLeastOneCond.")
    }
    if (missing(percent)) {
        stop("'percent' is missing.")
    }
    if (missing(op)) {
        stop("'op' is missing.")
    }
    if (missing(th)) {
        stop("'th' is missing.")
    } else if (!(op %in% SymFilteringOperators())) {
        stop(paste0(
            "'op' must be one of the following values: ",
            paste0(SymFilteringOperators(), collapse = " ")
        ))
    }

    u_conds <- unique(conds)
    nbCond <- length(u_conds)

    if (isTRUE(percent)) {
        if (th < 0 || th > 1) {
            warning("With percent=TRUE, the threshold 'th' must be in the 
                interval [0, 1].")
            return(NULL)
        }
    } else {
        th.upbound <- min(unlist(lapply(u_conds, 
            function(x) length(which(conds == x)))))
        if (th > th.upbound) {
            warn.txt <- paste0(
                "Param `th` is not correct. It must be an integer greater 
                than or equal to 0 and less or equal than ",
                th.upbound
            )
            warning(warn.txt)
            return(NULL)
        }
    }

    indices <- NULL
    s <- matrix(rep(0, nrow(metacell.mask) * nbCond),
        nrow = nrow(metacell.mask),
        ncol = nbCond
    )

    for (c in seq_len(nbCond)) {
        ind.cond <- which(conds == u_conds[c])
        inter <- rowSums(metacell.mask[, ind.cond])
        if (isTRUE(percent)) {
            inter <- inter / length(ind.cond)
        }

        s[, c] <- eval(parse(text = paste0("inter", op, th)))
    }

    indices <- switch(type,
        AllCond = which(rowSums(s) == nbCond),
        AtLeastOneCond = which(rowSums(s) >= 1)
    )

    return(indices)
}
