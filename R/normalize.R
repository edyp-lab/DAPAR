#' @title List normalization methods with tracking option
#'
#' @param withTracking xxx
#'
#' @name normalizeMethods.dapar
#'
#' @export
#' 
#' @return xxx
#' 
#' @examples 
#' normalizeMethods.dapar()
#'
normalizeMethods.dapar <- function(withTracking = FALSE) {
    if (withTracking) {
        c("SumByColumns", "QuantileCentering", "MeanCentering")
    } else {
        c(
            "GlobalQuantileAlignment",
            "SumByColumns",
            "QuantileCentering",
            "MeanCentering",
            "LOESS",
            "vsn"
        )
    }
}




#' @title Normalisation
#' 
#' @description 
#' Provides several methods to normalize quantitative data from
#' a \code{MSnSet} object.
#' They are organized in six main families : GlobalQuantileAlignement,
#' sumByColumns, QuantileCentering, MeanCentering, LOESS, vsn
#' For the first family, there is no type.
#' For the five other families, two type categories are available :
#' "Overall" which means that the value for each protein
#' (ie line in the expression data tab) is computed over all the samples ;
#' "within conditions" which means that the value for each protein
#' (ie line in the \code{Biobase::exprs()} data tab) is computed condition
#' by condition.
#'
#'
#' @param obj An object of class \code{MSnSet}.
#' @param method One of the following : "GlobalQuantileAlignment" (for
#' normalizations of important magnitude), "SumByColumns", "QuantileCentering",
#' "Mean Centering", "LOESS" and "vsn".
#'
#' @param withTracking xxx
#'
#' @param ... xxx
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' conds <- Biobase::pData(Exp1_R25_pept)$Condition
#' obj <- wrapper.normalizeD(
#'     obj = Exp1_R25_pept, method = "QuantileCentering",
#'     conds = conds, type = "within conditions"
#' )
#'
#' @export
#' 
#' @return xxx
#'
wrapper.normalizeD <- function(obj, method, withTracking = FALSE, ...) {
    if (!(method %in% normalizeMethods.dapar(withTracking))) {
        stop("'method' is not correct")
    }

    conds <- Biobase::pData(obj)[, "Condition"]
    qData <- Biobase::exprs(obj)

    switch(method,
        GlobalQuantileAlignment = {
            Biobase::exprs(obj) <- GlobalQuantileAlignment(qData)
            },
        SumByColumns = {
            Biobase::exprs(obj) <- SumByColumns(qData, ...)
            },
        QuantileCentering = {
            Biobase::exprs(obj) <- QuantileCentering(qData, ...)
            },
        MeanCentering = {
            Biobase::exprs(obj) <- MeanCentering(qData, ...)
            },
        vsn = {
            Biobase::exprs(obj) <- vsn(qData, ...)
            },
        # data must be log-expressed.
        LOESS = {
            Biobase::exprs(obj) <- LOESS(qData, ...)
        }
    )

    return(obj)
}





#' @title Normalisation GlobalQuantileAlignement
#'
#' @param qData xxxx
#'
#' @return A normalized numeric matrix
#'
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier,
#' Enora Fremy
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' normalized <- GlobalQuantileAlignment(qData)
#'
#' @export
#'
#'
GlobalQuantileAlignment <- function(qData) {
    
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
        stop("Please install preprocessCore: 
            BiocManager::install('preprocessCore')")
    }
    e <- preprocessCore::normalize.quantiles(as.matrix(qData))
    return(e)
}


#' @title Normalisation SumByColumns
#'
#' @param qData xxxx
#'
#' @param conds xxx
#'
#' @param type  Available values are "overall" (shift all the
#' sample distributions at once) or "within conditions" (shift the sample
#' distributions within each condition at a time).
#'
#' @param subset.norm A vector of index indicating rows to be used for
#' normalization
#'
#' @return A normalized numeric matrix
#'
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier,
#' Enora Fremy
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' conds <- Biobase::pData(Exp1_R25_pept)$Condition
#' normalized <- SumByColumns(qData, conds,
#'     type = "within conditions",
#'     subset.norm = seq_len(10)
#' )
#'
#' @export
#'
SumByColumns <- function(qData,
    conds = NULL,
    type = NULL,
    subset.norm = NULL) {
    
    if (!requireNamespace("stats", quietly = TRUE)) {
        stop("Please install stats: BiocManager::install('stats')")
    }
    
    
    if (missing(conds)) {
        stop("'conds' is required")
    }
    if (missing(type)) {
        stop("'type' is required")
    }


    if (!(type %in% c("overall", "within conditions"))) {
        stop("'type' must have one of the following values: 
            'overall', 'within conditions'")
    }


    qData <- as.matrix(qData)

    e <- 2^qData

    if (is.null(subset.norm) || length(subset.norm) < 1) {
        subset.norm <- seq_len(nrow(qData))
    }

    if (type == "overall") {
        if (length(subset.norm) == 1) {
            sum_cols <- e[subset.norm, ]
        } else {
            sum_cols <- colSums(e[subset.norm, ], na.rm = TRUE)
        }

        for (i in seq_len(nrow(e))) {
            e[i, ] <- (e[i, ] / sum_cols) * (stats::median(sum_cols))
        }
    } else if (type == "within conditions") {
        for (l in unique(conds)) {
            indices <- which(conds == l)

            if (length(subset.norm) == 1) {
                sum_cols <- e[subset.norm, indices]
            } else {
                sum_cols <- colSums(e[subset.norm, indices], na.rm = TRUE)
            }

            for (i in seq_len(nrow(e))) {
                e[i, indices] <- (e[i, indices] / sum_cols) * 
                    stats::median(sum_cols)
            }
        }
    }
    e <- log2(e)
    return(e)
}


#' @title Normalisation QuantileCentering
#'
#' @param qData xxx
#'
#' @param conds xxx
#'
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each
#' condition at a time).
#'
#' @param subset.norm A vector of index indicating rows to be used for
#' normalization
#'
#' @param quantile A float that corresponds to the quantile used to
#' align the data.
#'
#' @return A normalized numeric matrix
#'
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier,
#' Enora Fremy
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept
#' conds <- Biobase::pData(Exp1_R25_pept)$Condition
#' normalized <- QuantileCentering(Biobase::exprs(obj), conds,
#' type = "within conditions", subset.norm = seq_len(10)
#' )
#'
#' @export
#'
QuantileCentering <- function(qData,
    conds = NULL,
    type = "overall",
    subset.norm = NULL,
    quantile = 0.15) {
    
    if (!requireNamespace("stats", quietly = TRUE)) {
        stop("Please install stats: BiocManager::install('stats')")
    }
    if (missing(conds)) {
        stop("'conds' is required")
    }
    if (missing(type)) {
        stop("'type' is required")
    }


    if (!(type %in% c("overall", "within conditions"))) {
        stop("'type' must have one of the following values: 'overall', 
            'within conditions'")
    }


    qData <- as.matrix(qData)

    if (is.null(subset.norm) || length(subset.norm) < 1) {
        subset.norm <- seq_len(nrow(qData))
    }

    q <- function(x) {
        stats::quantile(x, probs = quantile, na.rm = TRUE)
    }


    if (length(subset.norm) == 1) {
        quantileOverSamples <- qData[subset.norm, ]
    } else {
        quantileOverSamples <- apply(qData[subset.norm, ], 2, q)
    }


    if (type == "overall") {
        cOverall <- q(quantileOverSamples)
        qData <- sweep(qData, 2, quantileOverSamples)
        qData <- qData + cOverall
    } else if (type == "within conditions") {
        qData <- sweep(qData, 2, quantileOverSamples)
        cCond <- NULL
        for (l in unique(conds)) {
            indices <- which(conds == l)
            cCond[l] <- q(quantileOverSamples[indices])
            qData[, indices] <- qData[, indices] + cCond[l]
        }
    }
    return(qData)
}


#' @title Normalisation MeanCentering
#'
#' @param qData xxx
#'
#' @param conds xxx
#'
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each
#' condition at a time).
#'
#' @param subset.norm A vector of index indicating rows to be used for
#' normalization
#'
#' @param scaling A boolean that indicates if the variance of the data have to
#' be forced to unit (variance reduction) or not.
#'
#' @return A normalized numeric matrix
#'
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier,
#' Enora Fremy
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' conds <- Biobase::pData(Exp1_R25_pept)$Condition
#' normalized <- MeanCentering(qData, conds, type = "overall")
#'
#' @export
#'
MeanCentering <- function(qData,
    conds,
    type = "overall",
    subset.norm = NULL,
    scaling = FALSE) {
    if (missing(conds)) {
        stop("'conds' is required")
    }

    qData <- as.matrix(qData)

    if (is.null(subset.norm) || length(subset.norm) < 1) {
        subset.norm <- seq_len(nrow(qData))
    }

    if (length(subset.norm) == 1) {
        meanOverSamples <- qData[subset.norm, ]
    } else {
        meanOverSamples <- apply(qData[subset.norm, ], 2, mean, na.rm = TRUE)
    }

    if (type == "overall") {
        cOverall <- mean(meanOverSamples)
        qData <- sweep(qData, 2, meanOverSamples)
        if (scaling) {
            qData <- scale(qData, center = FALSE, scale = TRUE)
            attr(qData, "scaled:scale") <- NULL
        }
        qData <- qData + cOverall
    } else if (type == "within conditions") {
        .temp <- sweep(qData, 2, meanOverSamples)
        if (scaling) {
            qData <- scale(qData, center = FALSE, scale = TRUE)
            attr(qData, "scaled:scale") <- NULL
        }
        cCond <- NULL
        for (l in unique(conds)) {
            indices <- which(conds == l)
            cCond[l] <- mean(meanOverSamples[indices])
            qData[, indices] <- .temp[, indices] + cCond[l]
        }
    }
    return(qData)
}


#' @title Normalisation vsn
#'
#' @param qData A numeric matrix.
#'
#' @param conds xxx
#'
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each condition
#' at a time).
#'
#' @return A normalized numeric matrix
#'
#' @author Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' conds <- Biobase::pData(Exp1_R25_pept)$Condition
#' normalized <- vsn(qData, conds, type = "overall")
#'
#' @export
#'
#'
vsn <- function(qData, conds, type = NULL) {
    if (!requireNamespace("vsn", quietly = TRUE)) {
        stop("Please install vsn: BiocManager::install('vsn')")
    }

    if (missing(conds)) {
        stop("'conds' is required")
    }

    if (type == "overall") {
        vsn.fit <- vsn::vsnMatrix(2^(qData))
        qData <- vsn::predict(vsn.fit, 2^(qData))
    } else if (type == "within conditions") {
        for (l in unique(conds)) {
            indices <- which(conds == l)
            vsn.fit <- vsn::vsnMatrix(2^(qData[, indices]))
            qData[, indices] <- vsn::predict(vsn.fit, 2^(qData[, indices]))
        }
    }
    return(qData)
}


#' @title Normalisation LOESS
#'
#' @param qData A numeric matrix.
#'
#' @param conds xxx
#'
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each
#' condition at a time).
#'
#' @param span xxx
#'
#' @return A normalized numeric matrix
#'
#' @author Thomas Burger, Helene Borges, Anais Courtier, Enora Fremy
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' qData <- Biobase::exprs(Exp1_R25_pept)
#' conds <- Biobase::pData(Exp1_R25_pept)$Condition
#' normalized <- LOESS(qData, conds, type = "overall")
#'
#' @export
#'
LOESS <- function(qData, conds, type = "overall", span = 0.7) {
    
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Please install limma: BiocManager::install('limma')")
    }
    
    if (missing(conds)) {
        stop("'conds' is required")
    }

    if (type == "overall") {
        qData <- limma::normalizeCyclicLoess(
            x = qData, 
            method = "fast", 
            span = span
            )
    } else if (type == "within conditions") {
        for (l in unique(conds)) {
            indices <- which(conds == l)
            qData[, indices] <- limma::normalizeCyclicLoess(
                x = qData[, indices],
                method = "fast",
                span = span
            )
        }
    }
    return(qData)
}
