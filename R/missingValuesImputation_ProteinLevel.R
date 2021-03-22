


#' This method finds the LAPALA in a dataset.
#'
#' @title Finds the LAPALA into a \code{MSnSet} object
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @return A data.frame that contains the indexes of LAPALA
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept[1:100]
#' lapala <- findMECBlock(obj)
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#'  
findMECBlock <- function(obj){
    
    conditions <- unique(Biobase::pData(obj)$Condition)
    nbCond <- length(conditions)
    
    s <- data.frame()
    
    for (cond in 1:nbCond){
        ind <- which(Biobase::pData(obj)$Condition == conditions[cond])
        lNA <- which(apply(is.na(Biobase::exprs(obj)[,ind]), 1, sum)==length(ind))
        if (length(lNA) > 0)
        {
            tmp <- data.frame(cond,which(apply(is.na(Biobase::exprs(obj)[,ind]), 1, sum)==length(ind)))
            names(tmp) <- c("Condition", "Line")
            s <- rbind(s,tmp)
            
        }
    }
    return(s)
}


#' @title Put back LAPALA into  a \code{MSnSet} object
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param MECIndex A data.frame that contains index of MEC (see findMECBlock) .
#' 
#' @return The object \code{obj} where LAPALA have been reintroduced
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept[1:100]
#' lapala <- findMECBlock(obj)
#' obj <- wrapper.impute.detQuant(obj, na.type='missing')
#' obj <- reIntroduceMEC(obj, lapala)
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' 
reIntroduceMEC <- function(obj, MECIndex){
    
    for (i in 1:nrow(MECIndex))
    {
        conditions <- unique(Biobase::pData(obj)$Condition)
        replicates <- which(Biobase::pData(obj)$Condition == conditions[MECIndex[i,"Condition"]])
        Biobase::exprs(obj)[MECIndex[i,"Line"], as.vector(replicates)] <- NA
    }
    return(obj)
}



#' This method is a wrapper for
#' objects of class \code{MSnSet} and imputes missing values with a fixed value.
#' This function imputes the missing values condition by condition.
#'
#' @title KNN missing values imputation from a \code{MSnSet} object
#' 
#' @description
#' Can impute only POV missing values
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param K the number of neighbors.
#' 
#' @param na.type A string which indicates the type of missing values to impute. 
#' Available values are: `missing_POV`.
#' 
#' @return The object \code{obj} which has been imputed
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj.imp.pov <- wrapper.impute.KNN(obj = Exp1_R25_pept[1:10], K=3, 
#' na.type = 'missing_POV')
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' 
wrapper.impute.KNN <- function(obj=NULL, K, na.type){
    #require(impute)
    if (missing(obj))
        stop("'obj' is required.")
    else if (is.null(obj))
        stop("'obj' is NULL")
    if (missing(na.type))
        stop("'na.type' is required. Available values are: 'missing_POV'.")
    else if (na.type != 'missing_POV')
        stop("Available value for na.type is: 'missing_POV'")

    data <- Biobase::exprs(obj)
    
    conditions <- unique(Biobase::pData(obj)$Condition)
    nbCond <- length(conditions)

    for (cond in 1:nbCond){
        ind <- which(Biobase::pData(obj)$Condition == conditions[cond])
        resKNN <- impute::impute.knn(Biobase::exprs(obj)[,ind] ,k = K, rowmax = 0.99, colmax = 0.99, maxp = 1500, rng.seed = sample(1:1000,1))
        Biobase::exprs(obj)[,ind] <- resKNN[[1]]
    }

    Biobase::exprs(obj)[Biobase::exprs(obj) == 0] <-NA
    obj <- UpdateMetacell(obj, 'knn', na.type) 
    
    return(obj)
}



#' This method is a wrapper to
#' objects of class \code{MSnSet} and imputes missing values with a fixed value.
#'
#' @title Missing values imputation from a \code{MSnSet} object
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param fixVal A float.
#' 
#' @param na.type A string which indicates the type of missing values to impute. 
#' Available values are: `NA` (for both POV and MEC), `POV`, `MEC`.
#' 
#' @return The object \code{obj} which has been imputed
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept[1:10,]
#' obj.imp.pov <- wrapper.impute.fixedValue(obj, 0.001, na.type = 'missing_POV')
#' obj.imp.mec <- wrapper.impute.fixedValue(obj, 0.001, na.type = 'missing_MEC')
#' obj.imp.na <- wrapper.impute.fixedValue(obj, 0.001, na.type = 'missing')
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' 
wrapper.impute.fixedValue <- function(obj, fixVal=0, na.type){
    if (missing(obj))
        stop("'obj' is required.")
    if(fixVal == 0)
        warning('Be aware that fixVal = 0. No imputation will be realize.')
    level <- obj@experimentData@other$typeOfData
    if (missing(na.type))
        stop(paste0("'na.type' is required. Available values are: ", paste0(metacell.def(level), collapse=' ')))
    else if (!(na.type %in% metacell.def(level)))
        stop(paste0("Available values for na.type are: ", paste0(metacell.def(level), collapse=' ')))
    
    ind.na.type <- match.metacell(Biobase::fData(obj)[, obj@experimentData@other$names_metacell], 
                                  na.type,
                                  level = obj@experimentData@other$typeOfData)
    Biobase::exprs(obj)[is.na(Biobase::exprs(obj)) & ind.na.type] <- fixVal
    obj <- UpdateMetacell(obj, 'fixedValue', na.type) 
    return (obj)
}



#' This method is a wrapper to the function \code{impute.pa} of the package
#' \code{imp4p} adapted to an object of class \code{MSnSet}.
#'
#' @title Imputation of peptides having no values in a biological condition.
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param q.min Same as the function \code{impute.pa} in the package 
#' \code{imp4p}
#' 
#' @param na.type A string which indicates the type of missing values to impute. 
#' Available values are: `NA` (for both POV and MEC).
#' 
#' @return The \code{exprs(obj)} matrix with imputed values instead of 
#' missing values.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept[1:10]
#' level <- obj@experimentData@other$typeOfData
#' metacell.mask <- match.metacell(GetMetacell(obj), 'missing', level)
#' indices <- GetIndices_WholeMatrix(metacell.mask, op='>=', th=1)
#' obj.imp.pov <- wrapper.impute.pa(obj, na.type = 'missing_POV')
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' @importFrom imp4p impute.pa
#' 
wrapper.impute.pa <- function(obj = NULL, q.min = 0.025, na.type){
    if (is.null(obj))
        stop("'obj' is required.")
    level <- obj@experimentData@other$typeOfData
    if (missing(na.type))
        stop(paste0("'na.type' is required. Available values are: ", paste0(metacell.def(level), collapse=' ')))
    else if (!(na.type %in% metacell.def(level)))
        stop(paste0("Available values for na.type are: ", paste0(metacell.def(level), collapse=' ')))
    
    cond <- as.factor(Biobase::pData(obj)$Condition)
    res <- impute.pa(Biobase::exprs(obj), conditions=cond, q.min = q.min, q.norm=3,  eps=0)
    Biobase::exprs(obj) <- res[["tab.imp"]]
    
    obj <- UpdateMetacell(obj, 'impute_pa', na.type) 
    
    return (obj)
}




#' This method is a wrapper of the function `impute.detQuant` for objects
#' of class \code{MSnSet} 
#' 
#' @title Wrapper of the function `impute.detQuant` for objects
#' of class \code{MSnSet}
#' 
#' @param obj An instance of class \code{MSnSet}
#' 
#' @param qval An expression set containing quantitative values of various 
#' replicates
#' 
#' @param factor A scaling factor to multiply the imputation value with 
#' 
#' @param na.type A string which indicates the type of missing values to impute. 
#' Available values are: `NA` (for both POV and MEC), `POV`, `MEC`.
#' 
#' @return An imputed instance of class \code{MSnSet}
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept[1:10]
#' obj.imp.pov <- wrapper.impute.detQuant(obj, na.type='missing_POV')
#' obj.imp.mec <- wrapper.impute.detQuant(obj, na.type='missing_MEC')
#' obj.imp.na <- wrapper.impute.detQuant(obj, na.type='missing')
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' 
wrapper.impute.detQuant <- function(obj, qval=0.025, factor=1, na.type){
    if (missing(obj))
        stop("'obj' is required.")
    availablePatterns <- unname(search.metacell.tags(pattern=na.type, 
                                              level=obj@experimentData@other$typeOfData))
    if (missing(na.type))
        stop(paste0("'na.type' is required. Available values are:.", 
                    paste0(availablePatterns, collapse=' '))
        )
    else if (!(na.type %in% search.metacell.tags(pattern=na.type, 
                                                 level=obj@experimentData@other$typeOfData)))
        stop(paste0("Available values for na.type are: ", 
                    paste0(availablePatterns, collapse=' '))
        )
    
    
    qdata <- Biobase::exprs(obj)
    values <- getQuantile4Imp(qdata, qval, factor)
    
   # Biobase::exprs(obj) <- impute.detQuant(qdata, values$shiftedImpVal, na.type)
    
    for(i in 1:ncol(qdata)){
        col <- qdata[,i]
        ind.na.type <- match.metacell(Biobase::fData(obj)[, obj@experimentData@other$names_metacell[i]], 
                                      pattern = na.type,
                                      level = obj@experimentData@other$typeOfData)
        
        col[ind.na.type] <- values$shiftedImpVal[i]
        qdata[,i] <- col
    }
    
    Biobase::exprs(obj) <- qdata
    msg <- "Missing values imputation using deterministic quantile"
    obj@processingData@processing <- c(obj@processingData@processing,msg)
    
    obj@experimentData@other$imputation.method <- "detQuantile"
    #browser()
    obj <- UpdateMetacell(obj = obj, method='detQuant', na.type=na.type) 
    
    return(obj)
}



#' @title Quantile imputation value definition
#' 
#' @description 
#' This method returns the q-th quantile of each column of an expression set,
#'  up to a scaling factor
#' 
#' @param qdata An expression set containing quantitative values of various 
#' replicates
#' 
#' @param qval The quantile used to define the imputation value
#' 
#' @param factor A scaling factor to multiply the imputation value with
#' 
#' @return A list of two vectors, respectively containing the imputation values 
#' and the rescaled imputation values
#' 
#' @author Thomas Burger
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' qdata <- Biobase::exprs(Exp1_R25_pept)
#' quant <- getQuantile4Imp(qdata) 
#' 
#' @export
#' 
getQuantile4Imp <- function(qdata, qval=0.025, factor=1){
    r1 <- apply(qdata, 2, quantile, qval, na.rm=TRUE)
    r2 <- r1*factor
    return(list(ImpVal = r1, shiftedImpVal = r2))
}


#' 
#' 
#' #' This method replaces each missing value by a given value
#' #'
#' #' @title Deterministic imputation
#' #' 
#' #' @param qdata An expression set containing quantitative or missing values
#' #' 
#' #' @param metadata xxx
#' #' 
#' #' @param values A vector with as many elements as the number of colums 
#' #' of qdata
#' #' 
#' #' @param na.type xxx
#' #' 
#' #' @return An imputed dataset
#' #' 
#' #' @author Thomas Burger, Samuel Wieczorek
#' #' 
#' #' @examples
#' #' utils::data(Exp1_R25_pept, package='DAPARdata')
#' #' qdata <- Biobase::exprs(Exp1_R25_pept)
#' #' meta <- 
#' #' values <- getQuantile4Imp(qdata)$shiftedImpVal
#' #' obj.imp.mec <- impute.detQuant(qdata, values, na.type = 'missing_POV') 
#' #' obj.imp.mec <- impute.detQuant(qdata, values, na.type = 'missing_MEC')
#' #' obj.imp.na <- impute.detQuant(qdata, values, na.type = 'missing') 
#' #' 
#' #' @export
#' #' 
#' impute.detQuant <- function(qdata, metadata, values, na.type=NULL){
#'     availablePatterns <- unname(search.metacell.tags(pattern=na.type, 
#'                                                      level=obj@experimentData@other$typeOfData))
#'     if (is.null(na.type))
#'         stop(paste0("'na.type' is required. Available values are:.", 
#'                     paste0(metacell.def(), collapse=' '))
#'         )
#'     else if (!(na.type %in% search.metacell.tags(pattern=na.type, 
#'                                                  level=obj@experimentData@other$typeOfData)))
#'         stop(paste0("Available values for na.type are: ", 
#'                     paste0(availablePatterns, collapse=' '))
#'         )
#'     
#'     if(missing(metadata))
#'         stop("`metadata` is missing.")
#'     
#'     
#'     #browser()
#'     for(i in 1:ncol(qdata)){
#'         col <- qdata[,i]
#'         ind.na.type <- match.metacell(Biobase::fData(obj)[, obj@experimentData@other$names_metacell[i]], 
#'                                       pattern = na.type,
#'                                       level = obj@experimentData@other$typeOfData)
#'         
#'         col[which(is.na(col) & ind.na.type)] <- values[i]
#'         qdata[,i] <- col
#'     }
#'     return(qdata)
#' }
#' 

#' This method is a wrapper to the function \code{impute.slsa} of the package
#' \code{imp4p} adapted to an object of class \code{MSnSet}.
#'
#' @title Imputation of peptides having no values in a biological condition.
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param na.type A string which indicates the type of missing values to impute. 
#' Available values are: `NA` (for both POV and MEC), `POV`, `MEC`.
#' 
#' @return The \code{exprs(obj)} matrix with imputed values instead of missing 
#' values.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept[1:10]
#' obj.slsa.pov <- wrapper.impute.slsa(obj, na.type = 'missing_POV')
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' 
wrapper.impute.slsa <- function(obj = NULL, na.type = NULL){
    if (is.null(obj))
        stop("'obj' is required.")
    if (is.null(na.type))
        stop("'na.type' is required. Available values are: 'missing_POV'.")
    else if (!(na.type %in% c('missing_POV')))
        stop("Available values for na.type are: 'missing_POV'.")
    
    
    MECIndex <- findMECBlock(obj)
    
    # sort conditions to be compliant with impute.slsa
    conds <- factor(Biobase::pData(obj)$Condition, levels=unique(Biobase::pData(obj)$Condition))
    sample.names.old <- Biobase::pData(obj)$Sample.name
    sTab <- Biobase::pData(obj)
    new.order <- unlist(lapply(split(sTab, conds), function(x) {x['Sample.name']}))
    qdata <- Biobase::exprs(obj)[,new.order]
    
    res <- imp4p::impute.slsa(qdata, 
                              conditions=conds, 
                              nknn=15, 
                              selec="all", 
                              weight=1,
                              ind.comp=1)
    
    #restore old order
    res <- res[,sample.names.old]
    
    Biobase::exprs(obj) <- res
    obj <- reIntroduceMEC(obj, MECIndex)
    
    obj <- UpdateMetacell(obj, 'slsa', na.type)
    return (obj)
}
