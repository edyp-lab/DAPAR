
#' @title Metadata vocabulary for entities
#'
#' @description
#' This function gives the vocabulary used for the metadata of each entity in
#' each condition.
#' Peptide-level vocabulary
#'
#' |-- 1.0 Quantitative Value
#' |    |
#' |    |-- 1.1 Identified (color 4, white)
#' |    |
#' |    |-- 1.2 Recovered (color 3, lightgrey)
#' |
#' |-- 2.0 Missing value (no color)
#' |    |
#' |    |-- 2.1 Missing POV (color 1)
#' |    |
#' |    |-- 2.2 Missing MEC (color 2)
#' |
#' |-- 3.0 Imputed value
#' |    |
#' |    |-- 3.1 Imputed POV (color 1)
#' |    |
#' |    |-- 3.2 Imputed MEC (color 2)
#'
#'
#'
#' Protein-level vocabulary:
#'
#' |-- 1.0 Quantitative Value
#' |    |
#' |    |-- 1.1 Identified (color 4, white)
#' |    |
#' |    |-- 1.2 Recovered (color 3, lightgrey)
#' |
#' |-- 2.0 Missing value
#' |    |
#' |    |-- 2.1 Missing POV (color 1)
#' |    |
#' |    |-- 2.2 Missing MEC (color 2)
#' |
#' |-- 3.0 Imputed value
#' |    |
#' |    |-- 3.1 Imputed POV (color 1)
#' |    |
#' |    |-- 3.2 Imputed MEC (color 2)
#' |
#' |-- 4.0 Combined value (color 3bis, light-lightgrey)
#'
#'
#' @param level A string designing the type of entity/pipeline.
#' Available values are: `peptide`, `protein`
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @export
#'
metacell.def <- function(level) {
    if (missing(level)) {
        stop("'level' is required.")
    }

    def <- switch(level,
        peptide = {
            node <- c(
                "all",
                "quanti",
                "identified",
                "recovered",
                "missing",
                "missing POV",
                "missing MEC",
                "imputed",
                "imputed POV",
                "imputed MEC"
            )
            parent <- c(
                "",
                "all",
                "quanti",
                "quanti",
                "all",
                "missing",
                "missing",
                "all",
                "imputed",
                "imputed"
            )
            data.frame(
                node = node,
                parent = parent
            )
        },
        protein = {
            node <- c(
                "all",
                "quanti",
                "identified",
                "recovered",
                "missing",
                "missing POV",
                "missing MEC",
                "imputed",
                "imputed POV",
                "imputed MEC",
                "combined"
            )
            parent <- c(
                "",
                "all",
                "quanti",
                "quanti",
                "all",
                "missing",
                "missing",
                "all",
                "imputed",
                "imputed",
                "all"
            )

            data.frame(
                node = node,
                parent = parent
            )
        }
    )


    colors <- list(
        "all" = "white",
        "missing" = "white",
        "missing POV" = "lightblue",
        "missing MEC" = "orange",
        "quanti" = "white",
        "recovered" = "lightgrey",
        "identified" = "white",
        "combined" = "red",
        "imputed" = "white",
        "imputed POV" = "#0040FF",
        "imputed MEC" = "#DF7401"
    )

    def <- cbind(def, color = rep("white", nrow(def)))

    for (n in 1:nrow(def)) {
        def[n, "color"] <- colors[[def[n, "node"]]]
    }

    return(def)
}



#' @title Sets the MEC tag in the metacell
#'
#' @description
#' This function is based on the metacell dataframe to look for either missing
#' values (used to update an initial dataset) or imputed values (used when
#' post processing protein metacell after aggregation)
#'
#' @param conds xxx
#'
#' @param df An object of class \code{MSnSet}
#'
#' @param level Type of entity/pipeline
#'
#' @return An instance of class \code{MSnSet}.
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_pept, package = "DAPARdata")
#' obj <- Exp1_R25_pept[1:10]
#' cols.for.ident <- obj@experimentData@other$names_metacell
#' conds <- Biobase::pData(obj)$Condition
#' df <- Biobase::fData(obj)[, cols.for.ident]
#' df <- Set_POV_MEC_tags(conds, df, level = "peptide")
#'
#' @export
#'
#'
Set_POV_MEC_tags <- function(conds, df, level) {
    u_conds <- unique(conds)

    for (i in 1:length(u_conds)) {
        ind.samples <- which(conds == u_conds[i])

        ind.imputed <- match.metacell(df[, ind.samples], "imputed", level)
        ind.missing <- match.metacell(df[, ind.samples], "missing", level)
        ind.missing.pov <- ind.missing & 
            rowSums(ind.missing) < length(ind.samples) & 
            rowSums(ind.missing) > 0
        ind.missing.mec <- ind.missing & 
            rowSums(ind.missing) == length(ind.samples)

        ind.imputed.pov <- ind.imputed & 
            rowSums(ind.imputed) < length(ind.samples) & 
            rowSums(ind.imputed) > 0
        ind.imputed.mec <- ind.imputed & 
            rowSums(ind.imputed) == length(ind.samples)

        df[, ind.samples][ind.imputed.mec] <- "imputed MEC"
        df[, ind.samples][ind.missing.mec] <- "missing MEC"
        df[, ind.samples][ind.imputed.pov] <- "imputed POV"
        df[, ind.samples][ind.missing.pov] <- "missing POV"
    }
    return(df)
}


#' @title xxxx
#'
#' @description
#' xxxxxx
#'
#' @param from xxx
#'
#' @param level xxx
#'
#' @param qdata An object of class \code{MSnSet}
#'
#' @param conds xxx
#'
#' @param df A list of integer xxxxxxx
#'
#' @return xxxxx
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' file <- system.file("extdata", "Exp1_R25_pept.txt", package = "DAPARdata")
#' data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#' metadataFile <- system.file("extdata", "samples_Exp1_R25.txt",
#'     package = "DAPARdata"
#' )
#' metadata <- read.table(metadataFile,
#'     header = TRUE, sep = "\t", as.is = TRUE,
#'     stringsAsFactors = FALSE
#' )
#' conds <- metadata$Condition
#' qdata <- data[, 56:61]
#' df <- data[, 43:48]
#' df <- BuildMetaCell(
#'     from = "maxquant", level = "peptide", qdata = qdata,
#'     conds = conds, df = df
#' )
#' df <- BuildMetaCell(
#'     from = "proline", level = "peptide", qdata = qdata,
#'     conds = conds, df = df
#' )
#'
#' @export
#'
#'
BuildMetaCell <- function(from,
                          level,
                          qdata = NULL,
                          conds = NULL,
                          df = NULL) {
    if (missing(from)) {
        stop("'from' is required.")
    }
    if (missing(level)) {
        stop("'level' is required.")
    }
    if (is.null(qdata)) {
        stop("'qdata' is required.")
    }
    if (is.null(conds)) {
        stop("'conds' is required.")
    }


    if (is.null(df)) {
        df <- Metacell_generic(qdata, conds, level)
    } else {
        switch(from,
            maxquant = df <- Metacell_maxquant(qdata, conds, df, level),
            proline = df <- Metacell_proline(qdata, conds, df, level)
        )
    }

    return(df)
}





#' @title Sets the metacell dataframe
#'
#' @description
#' In the quantitative columns, a missing value is identified by no value rather
#' than a value equal to 0.
#' Conversion rules
#' Quanti			Tag
#' NA or 0		NA
#'
#'
#' @param qdata An object of class \code{MSnSet}
#'
#' @param conds xxx
#'
#' @param level xxx
#'
#' @return xxxxx
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' file <- system.file("extdata", "Exp1_R25_pept.txt", package = "DAPARdata")
#' data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#' metadataFile <- system.file("extdata", "samples_Exp1_R25.txt",
#'     package = "DAPARdata"
#' )
#' metadata <- read.table(metadataFile,
#'     header = TRUE, sep = "\t", as.is = TRUE,
#'     stringsAsFactors = FALSE
#' )
#' conds <- metadata$Condition
#' qdata <- data[1:100, 56:61]
#' df <- data[1:100, 43:48]
#' df <- Metacell_generic(qdata, conds, level = "peptide")
#'
#' @export
#'
#'
Metacell_generic <- function(qdata, conds, level) {
    if (missing(qdata)) {
        stop("'qdata' is required")
    }
    if (missing(conds)) {
        stop("'conds' is required.")
    }
    if (missing(level)) {
        stop("'level' is required.")
    }

    df <- data.frame(matrix(rep("quanti", nrow(qdata) * ncol(qdata)),
        nrow = nrow(qdata),
        ncol = ncol(qdata)
    ),
    stringsAsFactors = FALSE
    )

    # Rule 1
    qdata[qdata == 0] <- NA
    df[is.na(qdata)] <- "missing"
    df <- Set_POV_MEC_tags(conds, df, level)

    colnames(df) <- paste0("metacell_", colnames(qdata))
    colnames(df) <- gsub(".", "_", colnames(df), fixed = TRUE)

    return(df)
}





#' @title Sets the metacell dataframe
#'
#' @description
#' In the quantitative columns, a missing value is identified by no value rather
#' than a value equal to 0.
#' Conversion rules
#' Initial conversion rules for maxquant
#' |--------------|-----------------|-----|
#' | Quanti       |    PSM count    | Tag |
#' |--------------|-----------------|-----|
#' |  == 0 | N.A.	|   whatever 			| 2.0 |
#' |  > 0		  		|    > 0		     	| 1.1 |
#' |  > 0		  		|    == 0	      	| 1.2 |
#' |  > 0		  		|   unknown col   | 1.0 |
#' |--------------|-----------------|-----|
#'
#' @param qdata An object of class \code{MSnSet}
#'
#' @param conds xxx
#'
#' @param df A list of integer xxxxxxx
#'
#' @param level xxx
#'
#' @return xxxxx
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' \dontrun{
#' file <- system.file("extdata", "Exp1_R25_pept.txt", package = "DAPARdata")
#' data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#' metadataFile <- system.file("extdata", "samples_Exp1_R25.txt",
#'     package = "DAPARdata"
#' )
#' metadata <- read.table(metadataFile,
#'     header = TRUE, sep = "\t", as.is = TRUE,
#'     stringsAsFactors = FALSE
#' )
#' conds <- metadata$Condition
#' qdata <- data[1:100, 56:61]
#' df <- data[1:100, 43:48]
#' df <- Metacell_proline(qdata, conds, df, level = "peptide")
#' }
#'
#' @export
#'
#'
Metacell_proline <- function(qdata, conds, df, level = NULL) {
    if (missing(qdata)) {
        stop("'qdata' is required")
    }
    if (missing(conds)) {
        stop("'conds' is required.")
    }
    if (missing(level)) {
        stop("'level' is required.")
    }


    if (is.null(df)) {
        df <- data.frame(matrix(rep("quanti", nrow(qdata) * ncol(qdata)),
            nrow = nrow(qdata),
            ncol = ncol(qdata)
        ),
        stringsAsFactors = FALSE
        )
    }

    # Rule 1
    df[is.na(qdata)] <- "missing"
    df <- Set_POV_MEC_tags(conds, df, level)

    # Rule 2
    df[df > 0 & qdata > 0] <- "identified"

    # Rule 3
    df[df == 0 & qdata > 0] <- "recovered"

    colnames(df) <- paste0("metacell_", colnames(qdata))
    colnames(df) <- gsub(".", "_", colnames(df), fixed = TRUE)

    return(df)
}

#' @title Sets the metacell dataframe
#'
#' @description
#' Initial conversion rules for maxquant
#' |------------|-----------------------|--------|
#' | Quanti     |     Identification    |    Tag |
#' |------------|-----------------------|--------|
#' |  == 0			|       whatever 				|    2.0 |
#' |  > 0				|       'By MS/MS'			|    1.1 |
#' |  > 0				|      'By matching'		|    1.2 |
#' |  > 0				|       unknown col			|    1.0 |
#' |------------|-----------------------|--------|
#'
#' @param qdata An object of class \code{MSnSet}
#'
#' @param conds xxx
#'
#' @param df A list of integer xxxxxxx
#'
#' @param level xxx
#'
#' @return xxxxx
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' file <- system.file("extdata", "Exp1_R25_pept.txt", package = "DAPARdata")
#' data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#' metadataFile <- system.file("extdata", "samples_Exp1_R25.txt",
#'     package = "DAPARdata"
#' )
#' metadata <- read.table(metadataFile,
#'     header = TRUE, sep = "\t", as.is = TRUE,
#'     stringsAsFactors = FALSE
#' )
#' conds <- metadata$Condition
#' qdata <- data[1:10, 56:61]
#' df <- data[1:10, 43:48]
#' df2 <- Metacell_maxquant(qdata, conds, df, level = "peptide")
#'
#' @export
#'
#'
Metacell_maxquant <- function(qdata, conds, df, level = NULL) {
    if (missing(qdata)) {
        stop("'qdata' is required")
    }
    if (missing(conds)) {
        stop("'conds' is required.")
    }
    if (missing(level)) {
        stop("'level' is required.")
    }


    if (is.null(df)) {
        df <- data.frame(matrix(rep(
            "quanti",
            nrow(qdata) * ncol(qdata)
        ),
        nrow = nrow(qdata),
        ncol = ncol(qdata)
        ),
        stringsAsFactors = FALSE
        )
    }


    # Rule 1
    qdata[qdata == 0] <- NA

    # Rule 2
    df[df == "byms/ms"] <- "identified"

    # Rule 3
    df[df == "bymatching"] <- "recovered"

    # Add details for NA values
    df[is.na(qdata)] <- "missing"
    df <- Set_POV_MEC_tags(conds, df, level)

    colnames(df) <- paste0("metacell_", colnames(qdata))
    colnames(df) <- gsub(".", "_", colnames(df), fixed = TRUE)

    return(df)
}






#' Similar to the function \code{is.na} but focused on the equality with
#' the paramter 'type'.
#'
#' @title Similar to the function \code{is.na} but focused on the equality
#' with the paramter 'type'.
#'
#' @param metadata A data.frame
#'
#' @param pattern The value to search in the dataframe
#'
#' @param level xxx
#'
#' @return A boolean dataframe
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_pept, package = "DAPARdata")
#' obj <- Exp1_R25_pept[1:10, ]
#' metadata <- GetMetacell(obj)
#' m <- match.metacell(metadata, pattern = "missing", level = "peptide")
#'
#' @export
#'
match.metacell <- function(metadata, pattern, level) {
    if (missing(metadata)) {
        stop("'metadata' is required")
    }
    if (missing(pattern)) {
        stop("'pattern' is required.")
    }
    if (missing(level)) {
        stop("'level' is required.")
    }


    if (!(pattern %in% metacell.def(level)$node)) {
        stop(paste0(
            "'pattern' is not correct. Available values are: ",
            paste0(metacell.def(level)$node, collapse = " ")
        ))
    }

    ll.res <- lapply(
        search.metacell.tags(pattern = pattern, level),
        function(x) {
            metadata == x
        }
    )

    res <- NULL
    for (i in 1:length(ll.res)) {
        if (i == 1) {
            res <- ll.res[[1]]
        } else {
            res <- res | ll.res[[i]]
        }
    }

    return(res)
}


#' @title xxxx
#'
#' @description
#' xxxx
#'
#' @param obj xxxx
#'
#' @export
#'
GetMetacell <- function(obj) {
    value <- Biobase::fData(obj)[, obj@experimentData@other$names_metacell]
    if (is.null(value)) {
        warning(" The metacell dataframe does not exist. Returns NULL.")
        return(NULL)
    } else {
        return(value)
    }
}

#' @title
#' Update metacell after imputation
#'
#' @description
#' Update the metacell information of missing values that were imputed
#'
#' @param obj xxx
#'
#' @param method xxx
#'
#' @param na.type xxx
#'
#' @author Samuel Wieczorek
#'
#' @export
#'
UpdateMetacell <- function(obj, method = "", na.type) {
    if (missing(obj)) {
        stop("'obj' is required.")
    }
    if (missing(na.type)) {
        values <- unname(search.metacell.tags(
            "missing",
            obj@experimentData@other$typeOfData
        ))
        stop(
            "'na.type' is required. Available values are: ",
            paste0(values, collapse = " ")
        )
    }

    level <- obj@experimentData@other$typeOfData
    ind <- match.metacell(
        metadata = Biobase::fData(obj)[, obj@experimentData@other$names_metacell],
        pattern = na.type,
        level = level
    ) & !is.na(Biobase::exprs(obj))

    Biobase::fData(obj)[, obj@experimentData@other$names_metacell][ind] <- 
        gsub("missing",
        "imputed",
        Biobase::fData(obj)[, obj@experimentData@other$names_metacell][ind],
        fixed = TRUE
    )
    return(obj)
}


#' @title
#' Search pattern in metacell vocabulary
#'
#' @description
#' Gives all the tags of the metadata vocabulary containing the pattern
#' (parent and all its children).
#'
#' @param pattern The string to search.
#'
#' @param level The available levels are : names()
#'
#' @param depth xxx
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' search.metacell.tags("missing POV", "peptide")
#' search.metacell.tags("quanti", "peptide")
#'
#' @export
#'

search.metacell.tags <- function(pattern, level, depth = "1") {
    if (missing(pattern)) {
        stop("'pattern' is required.")
    } else if (!(pattern %in% metacell.def(level)$node)) {
        stop(paste0("'pattern' must be one of the following: ", 
            paste0(metacell.def(level)$node, collapse = " ")))
    }

    if (missing(level)) {
        stop("'level' is required.")
    }
    if (!(depth %in% c("0", "1", "*"))) {
        stop("'depth' must be one of the following: 0, 1 or *")
    }

    .ind <- which(metacell.def(level)$parent == pattern)
    tags <- NULL
    tags <- switch(depth,
        "0" = pattern,
        "1" = c(pattern, 
            metacell.def(level)$node[.ind]),
        "*" = {
            if (length(metacell.def(level)$node[.ind]) == 0) {
                search.metacell.tags(pattern, level, "0")
            } else {
                c(pattern, unlist(lapply(
                    metacell.def(level)$node[.ind],
                    function(x) {
                        search.metacell.tags(x, level, depth)
                    }
                )))
            }
        }
    )

    return(tags)
}
