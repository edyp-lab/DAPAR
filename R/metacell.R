
#' @title Metadata vocabulary for entities
#' 
#' @description
#' This function gives the vocabulary used for the metadata of each entity in
#' each condition.
#' Peptide-level vocabulary
#' 
#' ├── 1.0 Quantitative Value
#' |    |
#' │    |── 1.1 Direct
#' |    |
#' │    |── 1.2 Indirect
#' │
#' ├── 2.0 Missing value
#' |    |
#' │    |── 2.1 Missing POV
#' |    |
#' │    |── 2.2 Missing MEC
#' │
#' ├── 3.0 Imputed value
#' |    |
#' │    |── 3.1 Imputed POV
#' |    |
#' │    |── 3.2 Imputed MEC
#'        
#' Protein-level vocabulary: same as peptide-level with one more category (Combined Value)
#' 
#' @param level A string designing the type of entity/pipeline. Available values are:
#' `peptide`, `protein`
#' 
#' @author Thomas Burger, Samuel Wieczorek
#' 
#' @export
#' 
metacell.def <- function(level = NULL){
  if(is.null(level))
    stop("'level' is required.")
  
  switch(level,
         peptide = setNames(nm = c('quanti',
                                   'quanti_identified',
                                   'quanti_recovered',
                                   'missing',
                                   'missing_POV',
                                   'missing_MEC',
                                   'imputed',
                                   'imputed_POV',
                                   'imputed_MEC')) ,
         
         protein = setNames(nm = c('quanti',
                                   'quanti_identified',
                                   'quanti_recovered',
                                   'missing',
                                   'missing_POV',
                                   'missing_MEC',
                                   'imputed',
                                   'imputed_POV',
                                   'imputed_MEC',
                                   'combined')
         )
         
  )
}


#' Sets the MEC tag in the metacell
#' 
#' @title Sets the MEC tag in the metacell
#' 
#' @param qData xxx
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
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' cols.for.ident <- xxxxx
#' df <- Biobase::fData(obj)[, cols.for.ident]
#' setMEC(df, Exp1_R25_pept, level = 'peptide')
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#'  
setMEC <- function(qData, conds, df, level){
  
  conditions <- unique(conds)
  nbCond <- length(conditions)
  
  for (cond in 1:nbCond){
    ind <- which(conds == conditions[cond])
    
    if (length(ind) == 1)
      lNA <- which(is.na(qData[,ind]))
    else
      lNA <- which(apply(is.na(qData[,ind]), 1, sum)==length(ind))
    
    if (length(lNA) > 0)
      df[lNA, ind] <- metacell.def(level)['missing_MEC']
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
#' @param qData An object of class \code{MSnSet}
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
#' file <- system.file("extdata", "Exp1_R25_pept.txt", package="DAPARdata")
#' data <- read.table(file, header=TRUE, sep="\t",stringsAsFactors = FALSE)
#' metadataFile <- system.file("extdata", "samples_Exp1_R25.txt", package="DAPARdata")
#' metadata <- read.table(metadataFile, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
#' conds <- metadata$Condition
#' qData <- data[,56:61]
#' df <- data[ , 43:48]
#' df <- BuildMetaCell(from = 'maxquant', level='peptide', qData = qData, conds = conds, df = df)
#' df <- BuildMetaCell(from = 'proline', level='peptide', qData = qData, conds = conds, df = df)
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' 
BuildMetaCell <- function(from = NULL, level = NULL, qData = NULL, conds = NULL, df = NULL){
  #if (is.null(from))
  #  stop("'from' is required.")
  if (is.null(qData))
    stop("'qData' is required.")
  if (is.null(conds))
    stop("'conds' is required.")
  if (is.null(level))
    stop("'level' is required.")
  # if (is.null(df) && from != 'proline')
  #   stop("'df' is required.")
  
  if (is.null(df))
    df <- Metacell_generic(qData, conds, level)
  else
    switch(from,
           maxquant = df <- Metacell_maxquant(qData, conds, df, level),
           proline = df <- Metacell_proline(qData, conds, df, level)
    )
  
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
#' @param qData An object of class \code{MSnSet}
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
#' file <- system.file("extdata", "Exp1_R25_pept.txt", package="DAPARdata")
#' data <- read.table(file, header=TRUE, sep="\t",stringsAsFactors = FALSE)
#' metadataFile <- system.file("extdata", "samples_Exp1_R25.txt", package="DAPARdata")
#' metadata <- read.table(metadataFile, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
#' conds <- metadata$Condition
#' qData <- data[,56:61]
#' df <- data[ , 43:48]
#' df <- Metacell_generic(qData, conds)
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' 
Metacell_generic <- function(qData, conds, level){
  
  df <- data.frame(matrix(rep(metacell.def(level)['quanti'], nrow(qData)*ncol(qData)),
                          nrow = nrow(qData),
                          ncol = ncol(qData)),
                   stringsAsFactors = FALSE) 
  
  # Rule 1
  df[is.na(qData)] <-  metacell.def(level)['missing_POV']
  df[qData == 0] <-  metacell.def(level)['missing_POV']
  df <- setMEC(qData, conds, df)
  
  colnames(df) <- paste0("metacell_", colnames(qData))
  colnames(df) <- gsub(".", "_", colnames(df), fixed=TRUE)
  
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
#' @param qData An object of class \code{MSnSet}
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
#' file <- system.file("extdata", "Exp1_R25_pept.txt", package="DAPARdata")
#' data <- read.table(file, header=TRUE, sep="\t",stringsAsFactors = FALSE)
#' metadataFile <- system.file("extdata", "samples_Exp1_R25.txt", package="DAPARdata")
#' metadata <- read.table(metadataFile, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
#' conds <- metadata$Condition
#' qData <- data[,56:61]
#' df <- data[ , 43:48]
#' Metacell_proline(qData, conds, df, level = 'peptide')
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' 
Metacell_proline <- function(qData, conds, df, level=NULL){
  
  if (is.null(df))
    df <- data.frame(matrix(rep(metacell.def(level)['quanti'], nrow(qData)*ncol(qData)), 
                            nrow=nrow(qData),
                            ncol=ncol(qData)),
                     stringsAsFactors = FALSE) 
  
  # Rule 1
  df[is.na(qData)] <-  metacell.def(level)['missing_POV']
  df <- setMEC(qData, conds, df, level)
  
  # Rule 2
  df[df > 0 && qData > 0] <- metacell.def(level)['quanti_identified']
  
  # Rule 3
  df[df == 0 && qData > 0] <- metacell.def(level)['quanti_recovered']
  
  colnames(df) <- paste0("metacell_", colnames(qData))
  colnames(df) <- gsub(".", "_", colnames(df), fixed=TRUE)
  
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
#' @param qData An object of class \code{MSnSet}
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
#' file <- system.file("extdata", "Exp1_R2_prot.txt", package="DAPARdata")
#' data <- read.table(file, header=TRUE, sep="\t",stringsAsFactors = FALSE)
#' metadataFile <- system.file("extdata", "samples_Exp1_R25.txt", package="DAPARdata")
#' metadata <- read.table(metadataFile, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
#' conds <- metadata$Condition
#' qData <- data[,49:54]
#' df <- data[ , 36:41]
#' df <- Metacell_maxquant(qData, conds, df, level='protein')
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' 
Metacell_maxquant <- function(qData, conds, df, level=NULL){
  
  if (is.null(df))
    df <- data.frame(matrix(rep(metacell.def(level)['quanti'], nrow(qData)*ncol(qData)), 
                            nrow=nrow(qData),
                            ncol=ncol(qData)),
                     stringsAsFactors = FALSE) 
  
  
  # Rule 1
  df[qData == 0] <-  NA
  
  # Rule 2
  df[df=='By MS/MS'] <- metacell.def(level)['quanti_identified']
  
  # Rule 3
  df[df=='By matching'] <- metacell.def(level)['quanti_recovered']
  
  
  # Add details for NA values
  df[is.na(qData)] <-  metacell.def(level)['missing_POV']
  df <- setMEC(qData, conds, df, level)
  
  
  
  colnames(df) <- paste0("metacell_", colnames(qData))
  colnames(df) <- gsub(".", "_", colnames(df), fixed=TRUE)
  
  return(df)
}






#' Similar to the function \code{is.na} but focused on the equality with the paramter 'type'.
#'
#' @title Similar to the function \code{is.na} but focused on the equality with the paramter 'type'.
#'
#' @param data A data.frame
#'
#' @param type The value to search in the dataframe
#' 
#' @param level xxx
#'
#' @return A boolean dataframe
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept[1:10,]
#' data <- Biobase::fData(obj)[, obj@experimentData@other$names_metacell]
#' match.metacell(data, type="MV_MEC", level = 'peptide')
#'
#' @export
#'
match.metacell <- function(data, type, level=NULL){
  if (!(type %in% metacell.def(level)))
    stop(paste0("'type' is not correct. It must be one of the following: ", paste0(metacell.def(level), collapse = ' ')))
  
  ll.res <- lapply(unname(search.metacell.tags(type, level)), function(x){data==x})
  
  res <- NULL
  for (i in 1:length(ll.res))
    if (i==1){
      res <- ll.res[[1]]
    } else {
      res <- res | ll.res[[i]]
    }
  
  return(res)
}


GetMetacell <- function(obj){
  Biobase::fData(obj)[, obj@experimentData@other$names_metacell]
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
UpdateMetacell <- function(obj=NULL, method='', na.type=NULL){
  if (is.null(obj))
    stop("'obj' is required.")
  if (is.null(na.type))
    stop("'na.type' is required. Available values are: NA, POV, MEC.")
  
  level <- obj@experimentData@other$typeOfData
  ind <- match.metacell(Biobase::fData(obj)[, obj@experimentData@other$names_metacell], na.type) & Biobase::exprs(obj) > 0 & !is.na(Biobase::exprs(obj))
  Biobase::fData(obj)[, obj@experimentData@other$names_metacell][ind] <- paste0(metacell.def(level)['imputed'], '_', method)
  return(obj)
}


#' @title xxxx
#' 
#' @description
#' Gives tags containing pattern as parent and all its children
#' 
#' @param pattern xxx
#' 
#' @param level xxx
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' search.metacell.tags('POV', 'peptide')
#' search.metacell.tags('MV_POV', 'peptide')
#' search.metacell.tags('quanti', 'peptide')
#' 
#' @export
#' 
search.metacell.tags <- function(pattern, level){
  unlist(metacell.def(level)[unlist(lapply(metacell.def(level), function(x){length(grep(pattern, x))==1}))])
}
