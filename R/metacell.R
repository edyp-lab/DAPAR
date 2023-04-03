
#' @title Metadata vocabulary for entities
#'
#' @description
#' This function gives the vocabulary used for the metadata of each entity in
#' each condition.
#' 
#' Peptide-level vocabulary
#'
#' |-- 'Any'
#' |    |
#' |    |-- 1.0 'Quantified'
#' |    |    |
#' |    |    |-- 1.1 "Quant. by direct id" (color 4, white)
#' |    |    |
#' |    |    |-- 1.2 "Quant. by recovery" (color 3, lightgrey)
#' |    |
#' |    |-- 2.0 "Missing" (no color)
#' |    |    |
#' |    |    |-- 2.1 "Missing POV" (color 1)
#' |    |    |
#' |    |    |-- 2.2 'Missing MEC' (color 2)
#' |    |
#' |    |-- 3.0 'Imputed'
#' |    |    |
#' |    |    |-- 3.1 'Imputed POV' (color 1)
#' |    |    |
#' |    |    |-- 3.2 'Imputed MEC' (color 2)
#'
#'
#'
#' Protein-level vocabulary:
#' |-- 'Any'
#' |    |
#' |    |-- 1.0 'Quantified'
#' |    |    |
#' |    |    |-- 1.1 "Quant. by direct id" (color 4, white)
#' |    |    |
#' |    |    |-- 1.2 "Quant. by recovery" (color 3, lightgrey)
#' |    |
#' |    |-- 2.0 "Missing"
#' |    |    |
#' |    |    |-- 2.1 "Missing POV" (color 1)
#' |    |    |
#' |    |    |-- 2.2 'Missing MEC' (color 2)
#' |    |
#' |    |-- 3.0 'Imputed'
#' |    |    |
#' |    |    |-- 3.1 'Imputed POV' (color 1)
#' |    |    |
#' |    |    |-- 3.2 'Imputed MEC' (color 2)
#' |    |
#' |    |-- 4.0 'Combined tags' (color 3bis, lightgrey)
#'
#'
#' @param level A string designing the type of entity/pipeline.
#' Available values are: `peptide`, `protein`
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @export
#' 
#' @return xxx
#' 
#' @examples 
#' metacell.def('protein')
#' metacell.def('peptide')
#' 
#'
metacell.def <- function(level) {
    if (missing(level)) {
        stop("'level' is required.")
    }
    


    def <- switch(level,
        peptide = {
            node <- c(
                "Any",
                "Quantified",
                "Quant. by direct id",
                "Quant. by recovery",
                "Missing",
                "Missing POV",
                "Missing MEC",
                "Imputed",
                "Imputed POV",
                "Imputed MEC"
            )
            parent <- c(
                "",
                "Any",
                "Quantified",
                "Quantified",
                "Any",
                "Missing",
                "Missing",
                "Any",
                "Imputed",
                "Imputed"
            )
            data.frame(
                node = node,
                parent = parent
            )
        },
        protein = {
            node <- c(
                "Any",
                "Quantified",
                "Quant. by direct id",
                "Quant. by recovery",
                "Missing",
                "Missing POV",
                "Missing MEC",
                "Imputed",
                "Imputed POV",
                "Imputed MEC",
                "Combined tags"
            )
            parent <- c(
                "",
                "Any",
                "Quantified",
                "Quantified",
                "Any",
                "Missing",
                "Missing",
                "Any",
                "Imputed",
                "Imputed",
                "Any"
            )

            data.frame(
                node = node,
                parent = parent
            )
        }
    )


    colors <- list(
        "Any" = "white",
        "Missing" = "#CF8205",
        "Missing POV" = "#E5A947",
        "Missing MEC" = "#F1CA8A",
        "Quantified" = "#0A31D0",
        "Quant. by recovery" = "#B9C4F2",
        "Quant. by direct id" = "#6178D9",
        "Combined tags" = "#1E8E05",
        "Imputed" = "#A40C0C",
        "Imputed POV" = "#E34343",
        "Imputed MEC" = "#F59898"
    )

    def <- cbind(def, color = rep("white", nrow(def)))

    for (n in seq_len(nrow(def))) {
        def[n, "color"] <- colors[[def[n, "node"]]]
    }

    return(def)
}


#' @title Parent name of a node
#' @description xxx
#' @param level xxx
#' @param node xxx
#' 
#' #' @examples 
#' Parent('protein', 'Missing')
#' Parent('protein', 'Missing POV')
#' Parent('protein', c('Missing POV', 'Missing MEC'))
#' Parent('protein', c('Missing', 'Missing POV', 'Missing MEC'))
#' 
#' 
#' @export
Parent <- function(level, node=NULL){
    parents <- NULL
    tags <- metacell.def(level)
    
    if (!is.null(node) && length(node) > 0){
      for (p in node){
        ind <- match(p, tags$node)
        if (length(ind) > 0)
          parents <- unique(c(parents, tags$parent[ind]))
      }
    }
    
    
    return(parents)
}

#' @title Names of all chidren of a node
#' @description xxx
#' @param level xxx
#' @param parent xxx
#' 
#' @examples 
#' Children('protein', 'Missing')
#' Children('protein', 'Missing POV')
#' Children('protein', c('Missing POV', 'Missing MEC'))
#' Children('protein', c('Missing', 'Missing POV', 'Missing MEC'))
#' @export
Children <- function(level, parent = NULL){
  childrens <- NULL
  tags <- metacell.def(level)
  if (!is.null(parent) && length(parent) > 0){
    for (p in parent){
      ind <- grepl(p, tags$parent)
      if (length(ind) > 0)
        childrens <- unique(c(childrens, tags$node[ind]))
    }
  }
  return(childrens)
}

#' @title xxxx
#' @description xxx
#' @param obj xxx
#' @export
GetUniqueTags <- function(obj){
    df <- Biobase::fData(obj)[, obj@experimentData@other$names_metacell]
    tmp <- sapply(colnames(df), function(x) unique(df[,x]))
    ll <- unique(as.vector(tmp))
    return(ll)
}

#' @title List of metacell tags
#'
#' @description
#' This function gives the list of metacell tags available in DAPAR.
#' 
#' - onlyPresent: In this case, the function gives the tags found in a dataset.
#' In addition, and w.r.t to the hierarchy of tags, if all leaves of a node are
#' present, then the tag corresponding to this node is added.
#'
#' @param level xxx
#'
#' @param obj An object of class \code{MSnSet}
#'
#' @param onlyPresent A boolean that indicates if one wants a list with only the tags
#' present in the dataset.
#' 
#' @param all A boolean that indicates if one wants the whole list
#'
#' @return A vector of tags..
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept
#' GetMetacellTags(level="peptide")
#' GetMetacellTags(level="peptide", obj, onlyPresent=TRUE)
#'
#' @export
#'
#'
GetMetacellTags <- function(level = NULL,
                            obj = NULL,
                            onlyPresent = FALSE, 
                            all = FALSE) {
    
    if (!onlyPresent && !all){
        if (!is.null(level) && is.null(obj))
            all <- TRUE
        if (is.null(level) && !is.null(obj))
            onlyPresent <- TRUE
        if ((!is.null(level) && !is.null(obj)) || (is.null(level) && is.null(obj)))
                stop("At least, one of level or obj must be defined")
    }
       
    if (onlyPresent && all)
        stop("Only of 'onlyPresent' or 'all' must be TRUE")
    
    if (is.null(level) && all)
        stop("level must be defined")
    if (is.null(level) && !all)
            all <- TRUE
    
    if (is.null(obj) && onlyPresent)
        stop("`obj` must be defined")
    
    #browser()
    ll <- NULL
    if(onlyPresent) {
        ll <- unique(unlist(GetUniqueTags(obj)))
        # Check if parent must be added
        test <- match (Children(level, 'Any'), ll)
        if (length(test) == length(Children(level, 'Any')) && !all(is.na(test)))
            ll <- c(ll, 'Any')
        
        test <- match (Children(level, 'Quantified'), ll)
        if (length(test) == length(Children(level, 'Quantified')) && !all(is.na(test)))
                ll <- c(ll, 'Quantified')
        
        test <- match (Children(level, 'Missing'), ll)
        if (length(test) == length(Children(level, 'Missing')) && !all(is.na(test)))
            ll <- c(ll, 'Missing')
        
        test <- match (Children(level, 'Imputed'), ll)
        if (length(test) == length(Children(level, 'Imputed')) && !all(is.na(test)))
            ll <- c(ll, 'Imputed')
        
        test <- match (Children(level, 'Combined tags'), ll)
        if (length(test) == length(Children(level, 'Combined tags')) && !all(is.na(test)))
            ll <- c(ll, 'Combined tags')
        
        
        
    } else if (all) {
      ll <- metacell.def(level)$node[-which(metacell.def(level)$node =='Any')]
      #  ll <- metacell.def(level)
    }
    
    return(ll)
    
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
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' cols.for.ident <- c("metacell_Intensity_C_R1", "metacell_Intensity_C_R2",
#' "metacell_Intensity_C_R3", "metacell_Intensity_D_R1",
#' "metacell_Intensity_D_R2", "metacell_Intensity_D_R3")
#' conds <- Biobase::pData(obj)$Condition
#' df <- Biobase::fData(obj)[, cols.for.ident]
#' df <- Set_POV_MEC_tags(conds, df, level = "peptide")
#'
#' @export
#'
#'
Set_POV_MEC_tags <- function(conds, df, level) {
    u_conds <- unique(conds)

    for (i in seq_len(length(u_conds))) {
        ind.samples <- which(conds == u_conds[i])

        ind.imputed <- match.metacell(df[, ind.samples], "Imputed", level)
        ind.missing <- match.metacell(df[, ind.samples], "Missing", level)
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

        df[, ind.samples][ind.imputed.mec] <- "Imputed MEC"
        df[, ind.samples][ind.missing.mec] <- "Missing MEC"
        df[, ind.samples][ind.imputed.pov] <- "Imputed POV"
        df[, ind.samples][ind.missing.pov] <- "Missing POV"
    }
    return(df)
}



#' @title The set of softwares available
#' 
#' @examples 
#' GetSoftAvailables()
#' @export

GetSoftAvailables <- function(){
    
    
    library(DAPAR)
    
    funcs <- ls('package:DAPAR')
    funcs <- funcs[grep('Metacell_', funcs)]
    funcs <- strsplit(funcs, 'Metacell_')
    funcs <- unlist(lapply(funcs, function(x) x[[2]]))
    funcs <- funcs[-which(funcs=='generic')]
    
    return(funcs)
}

#' @title Builds cells metadata
#'
#' @description
#' This function the cells metadata info base on the origin of identification
#' for entities.
#' There are actually two different type of origin which are managed by DAPAR:
#' - "Maxquant-like" info which is represented by strings/tags,
#' - Proline-like where the info which is used is an integer
#'
#' @param from A string which is the name of the software from which the data
#' are. Available values are 'maxquant', 'proline' and 'DIA-NN'
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
#' qdata <- data[, seq.int(from = 56, to = 61)]
#' df <- data[, seq.int(from = 43, to = 48)]
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
    if (!(from %in% GetSoftAvailables()))
        stop("'from' must be one of the following")
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
            proline = df <- Metacell_proline(qdata, conds, df, level),
          'DIA-NN' = df <- Metacell_proline(qdata, conds, df, level)
        )
    }

    return(df)
}





#' @title Sets the metacell dataframe for dataset without information about the
#' origin of identification
#'
#' @description
#' In the quantitative columns, a missing value is identified by no value rather
#' than a value equal to 0.
#' Conversion rules
#' QuantiTag
#' NA or 0 NA
#' The only information detected with this function are about missing values (
#' MEC and POV).
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
#' qdata <- data[seq_len(100), seq.int(from = 56, to = 61)]
#' df <- data[seq_len(100), seq.int(from = 43, to = 48)]
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

    df <- data.frame(
        matrix(rep("Quantified", nrow(qdata) * ncol(qdata)),
            nrow = nrow(qdata),
            ncol = ncol(qdata)
            ),
        stringsAsFactors = FALSE
        )

    # Rule 1
    qdata[qdata == 0] <- NA
    df[is.na(qdata)] <- "Missing"
    df <- Set_POV_MEC_tags(conds, df, level)

    colnames(df) <- paste0("metacell_", colnames(qdata))
    colnames(df) <- gsub(".", "_", colnames(df), fixed = TRUE)

    return(df)
}







#' @title Sets the metacell dataframe for datasets which are from Dia-NN software
#'
#' @description
#' Actually, this function uses the generic function to generate metacell info
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
#' qdata <- data[seq_len(100), seq.int(from = 56, to = 61)]
#' df <- data[seq_len(100), seq.int(from = 43, to = 48)]
#' df <- Metacell_DIA_NN(qdata, conds, df, level = "peptide")
#' 
#'
#' @export
#'
#'
Metacell_DIA_NN <- function(qdata, conds, df, level = NULL) {
    if (missing(qdata)) {
        stop("'qdata' is required")
    }
    if (missing(conds)) {
        stop("'conds' is required.")
    }
    if (missing(level)) {
        stop("'level' is required.")
    }
    
    
    return(df)
}






#' @title Sets the metacell dataframe for datasets which are from Proline software
#'
#' @description
#' In the quantitative columns, a missing value is identified by no value rather
#' than a value equal to 0.
#' 
#' In these datasets, the metacell info is computed from the 'PSM count' columns.
#' 
#' Conversion rules
#' Initial conversion rules for proline
#' |--------------|-----------------|-----|
#' | Quanti       |    PSM count    | Tag |
#' |--------------|-----------------|-----|
#' |  == 0 | N.A. |   whatever      | 2.0 |
#' |  > 0         |    > 0          | 1.1 |
#' |  > 0         |    == 0         | 1.2 |
#' |  > 0         |   unknown col   | 1.0 |
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
#' qdata <- data[seq_len(100), seq.int(from = 56, to = 61)]
#' df <- data[seq_len(100), seq.int(from = 43, to = 48)]
#' df <- Metacell_proline(qdata, conds, df, level = "peptide")
#' 
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
        df <- data.frame(matrix(rep("Quantified", nrow(qdata) * ncol(qdata)),
            nrow = nrow(qdata),
            ncol = ncol(qdata)
        ),
        stringsAsFactors = FALSE
        )
    }

    # Rule 1
    df[is.na(qdata)] <- "Missing"
    df <- Set_POV_MEC_tags(conds, df, level)

    # Rule 2
    df[df > 0 & qdata > 0] <- "Quant. by direct id"

    # Rule 3
    df[df == 0 & qdata > 0] <- "Quant. by recovery"

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
#' |  == 0      |       whatever        |    2.0 |
#' |  > 0       |       'By MS/MS'      |    1.1 |
#' |  > 0       |      'By matching'    |    1.2 |
#' |  > 0       |       unknown col     |    1.0 |
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
#' qdata <- data[seq_len(10), seq.int(from = 56, to = 61)]
#' df <- data[seq_len(10), seq.int(from = 43, to = 48)]
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
            "Quantified",
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
    df[df == "byms/ms"] <- "Quant. by direct id"

    # Rule 3
    df[df == "bymatching"] <- "Quant. by recovery"

    # Add details for NA values
    df[is.na(qdata)] <- "Missing"
    df <- Set_POV_MEC_tags(conds, df, level)

    colnames(df) <- paste0("metacell_", colnames(qdata))
    colnames(df) <- gsub(".", "_", colnames(df), fixed = TRUE)

    return(df)
}


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
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept[seq_len(10), ]
#' metadata <- GetMetacell(obj)
#' m <- match.metacell(metadata, pattern = "Missing", level = "peptide")
#' m <- match.metacell(metadata, pattern = NULL, level = "peptide")
#' m <- match.metacell(metadata, pattern = c('Missing', 'Missing POV'), level = "peptide")
#' @export
#'
match.metacell <- function(metadata, pattern = NULL, level) {
    if (missing(metadata))
        stop("'metadata' is required")

    if (missing(pattern))
        stop("'pattern' is required.")
  else if (is.null(pattern))
    return(NULL)

    if (missing(level))
        stop("'level' is required.")


    #is.subset <- pattern == intersect(pattern,  metacell.def(level)$node)
    if (sum(pattern == intersect(pattern,  metacell.def(level)$node)) !=  length(pattern)) {
        stop(paste0(
            "'pattern' is not correct. Available values are: ",
            paste0(metacell.def(level)$node, collapse = " ")
        ))
    }

    ll.res <- lapply(pattern, function(x) {metadata == x})

    res <- NULL
    for (i in seq_len(length(ll.res))) {
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
#' @return xxx
#' 
#' @examples 
#' NULL
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
#' Update the cells metadata tags after imputation
#'
#' @description
#' Update the metacell information of missing values that were imputed
#'
#' @param obj xxx
#'
#' @author Samuel Wieczorek
#'
#' @export
#' 
#' @return xxx
#' 
#' @examples 
#' data(Exp1_R25_pept, package="DAPARdata")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' obj.imp.pov <- wrapper.impute.KNN(obj, K = 3)
#'
UpdateMetacellAfterImputation <- function(obj) {
    if (missing(obj))
        stop("'obj' is required.")

    
  ind <- match.metacell(metadata = GetMetacell(obj),
                        pattern = c('Missing', 'Missing POV', 'Missing MEC'),
                        level = GetTypeofData(obj)) & !is.na(Biobase::exprs(obj))
  
  names.meta <- colnames(GetMetacell(obj))
  Biobase::fData(obj)[, names.meta][ind] <- gsub("Missing", "Imputed",
                                              Biobase::fData(obj)[, names.meta][ind],
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
#' search.metacell.tags("Missing POV", "peptide")
#' search.metacell.tags("Quantified", "peptide", depth = "0")
#'
#' @export
#' 
#' @return xxx
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
