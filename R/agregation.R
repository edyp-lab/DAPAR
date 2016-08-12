##' This function computes the number of proteins that are only defined by 
##' specific peptides, shared peptides or a mixture of two. 
##' 
##' @title computes the number of proteins that are only defined by 
##' specific peptides, shared peptides or a mixture of two.
##' @param matUnique The adjacency matrix with only specific peptides.
##' @param matShared The adjacency matrix with both specific and shared peptides.
##' @return A list
##' @author Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' protID <- "Protein.group.IDs"
##' MShared <- BuildAdjacencyMatrix(UPSpep25, protID, FALSE)
##' MUnique <- BuildAdjacencyMatrix(UPSpep25, protID, TRUE)
##' getProteinsStats(MUnique,MShared)
getProteinsStats <- function(matUnique, matShared){
    if (is.null(matUnique) || is.null(matShared) ||
        !is.matrix(matUnique) || !is.matrix(matShared)){return(NULL)}
    
    t <- setdiff(union(rownames(matUnique), rownames(matShared)), intersect(rownames(matUnique), rownames(matShared)))
    sharedPeptides <- matShared[t,]
    sharedPeptides <- sharedPeptides[,-which(colSums(sharedPeptides)==0)]
    protOnlyUnique <- setdiff(union(colnames(sharedPeptides), colnames(matShared)), intersect(colnames(sharedPeptides), colnames(matShared)))
    
    protOnlyShared <- setdiff(union(colnames(matUnique), colnames(matShared)), intersect(colnames(matUnique), colnames(matShared)))
    a <- union(protOnlyUnique, protOnlyShared)
    b <- colnames(matShared)
    protMix <- setdiff(union(union(protOnlyUnique, protOnlyShared),colnames(matShared)),intersect(union(protOnlyUnique, protOnlyShared),colnames(matShared)))

    return (list(protOnlyUniquePep =protOnlyUnique,
                  protOnlySharedPep =protOnlyShared,
                  protMixPep = protMix))
}




##' This function creates a column for the protein dataset after agregation 
##' by using the previous peptide dataset.
##' 
##' @title creates a column for the protein dataset after agregation by
##'  using the previous peptide dataset.
##' @param peptideData A data.frame of meta data of peptides. It is the fData 
##' of the MSnset object.
##' @param matAdj The adjacency matrix used to agregate the peptides data.
##' @param columnName The name of the column in fData(peptides_MSnset) that 
##' the user wants to keep in the new protein data.frame.
##' @return A vector
##' @author Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' protID <- "Protein.group.IDs"
##' M <- BuildAdjacencyMatrix(UPSpep25, protID, FALSE)
##' data <- Biobase::fData(UPSpep25)
##' name <- "organism"
##' BuildColumnToProteinDataset(data, M, name )
BuildColumnToProteinDataset <- function(peptideData, matAdj, columnName){
nbProt <- ncol(matAdj)
newCol <- rep("", nbProt)

for (p in 1:nbProt){
    listeIndicePeptides <- which(matAdj[,p] == 1)
    listeData <- unique(peptideData[listeIndicePeptides,columnName])
    newCol[p] <- paste(listeData, collapse = ",")
}
return(newCol)
}



##' This function computes the number of peptides used to aggregate proteins.
##' 
##' @title Compute the number of peptides used to aggregate proteins
##' @param M A "valued" adjacency matrix in which lines and columns correspond 
##' respectively to peptides and proteins.
##' @return A vector of boolean which is the adjacency matrix 
##' but with NA values if they exist in the intensity matrix.
##' @author Alexia Dorffer
##' @examples
##' data(UPSpep25)
##' protID <- "Protein.group.IDs"
##' M <- BuildAdjacencyMatrix(UPSpep25, protID, FALSE)
##' CountPep(M)
CountPep <- function (M) {
    z <- M
    z[z!=0] <- 1
    return(z)
}


##' This function computes the intensity of proteins based on the sum of the 
##' intensities of their peptides.
##' 
##' @title Compute the intensity of proteins with the sum of the intensities
##' of their peptides.
##' @param matAdj An adjacency matrix in which lines and columns correspond 
##' respectively to peptides and proteins.
##' @param expr A matrix of intensities of peptides
##' @return A matrix of intensities of proteins
##' @author Alexia Dorffer
##' @examples
##' data(UPSpep25)
##' protID <- "Protein.group.IDs"
##' M <- BuildAdjacencyMatrix(UPSpep25, protID, FALSE)
##' SumPeptides(M, Biobase::exprs(UPSpep25))
SumPeptides <- function(matAdj, expr){
    # require(foreach)
    #register for use of parallel foreach
    cores <- detectCores()
    registerDoParallel(cores)
    j <- NULL
    ############### parallel #################
    t <- foreach (j=1:ncol(expr), .combine='cbind') %dopar% {  
            M <- matAdj * expr[rownames(matAdj),j]
            z <- CountPep(M)
            matrix(c(colSums(M, na.rm=TRUE),colSums(z, na.rm=TRUE)),ncol=1) 
            }
    Mp <-t[1:ncol(matAdj),]
    rownames(Mp) <- colnames(matAdj)
    colnames(Mp) <- colnames(expr)
    pep <-t[-c(1:ncol(matAdj)),]
    colnames(pep) <- paste("nb.pep.used.", colnames(expr), sep="")
    rownames(pep) <- colnames(matAdj)
    
    ############# original #####################
#   Mp <- matrix(c(rep(NA)), ncol=nrow(condition), nrow=ncol(matAdj), 
#    dimnames=c(list(colnames(matAdj), colnames(expr))))
#     pep <- Mp
#     for (j in 1:ncol(Mp)){  
#         M <- matAdj * expr[,j]
#         z <- CountPep(M)
#         pep[,j] <- colSums(z, na.rm=TRUE)
#         Mp[,j] <- colSums(M, na.rm=TRUE) 
#     }

    res <- list("idprot" = colnames(matAdj), "matfin"=Mp, "nbpep"=pep)
    return(res)


# return(Mp)
}


##' This function computes the intensity of proteins as the mean of the 
##' intensities of their peptides.
##' 
##' @title Compute the intensity of proteins as the mean of the intensities
##' of their peptides.
##' @param matAdj An adjacency matrix in which lines and columns correspond 
##' respectively to peptides and proteins.
##' @param expr A matrix of intensities of peptides
##' @return A matrix of intensities of proteins
##' @author Alexia Dorffer
##' @examples
##' data(UPSpep25)
##' protID <- "Protein.group.IDs"
##' matAdj <- BuildAdjacencyMatrix(UPSpep25, protID, FALSE)
##' MeanPeptides(matAdj, Biobase::exprs(UPSpep25))
MeanPeptides <- function(matAdj,expr){
    #require(foreach)
    ##register for use of parallel foreach
cores <- detectCores()
registerDoParallel(cores)

    ############## Parallel Version ################
    z <- j <- NULL
    t <- foreach (j=1:ncol(expr), .combine='cbind') %dopar% {  
            M <- matAdj * expr[rownames(matAdj),j]
            z <- CountPep(M)
            c(apply(M, 2, function(x) sum(x, na.rm=TRUE)/sum(!(x==0), 
                            na.rm=TRUE)),colSums(z, na.rm=TRUE))
            }
    Mp <-t[1:ncol(matAdj),]
    colnames(Mp) <- colnames(expr)
    rownames(Mp) <- colnames(matAdj)
    pep <-t[-c(1:ncol(matAdj)),]
    colnames(pep) <- paste("nb.pep.used.", colnames(expr), sep="")
    rownames(pep) <- colnames(matAdj)

    ####### Sequential version #########
##      pep <- Mp
##     for (j in 1:ncol(Mp)){  
##         M <- matAdj * expr[,j]
##         z <- CountPep(M)
##         pep[,j] <- colSums(z, na.rm=TRUE)
##         Mp[,j] <- apply(M, 2, function(x) sum(x, na.rm=T)/sum(!(x==0), 
##              na.rm=TRUE))
##     }  
    
    res <- list("idprot" = colnames(matAdj), "matfin"=Mp, "nbpep"=pep)
    return(res)
}

##' This function computes the intensity of proteins as the sum of the 
##' intensities of their n best peptides.
##' 
##' @title Compute the intensity of proteins as the sum of the 
##' intensities of their n best peptides.
##' @param matAdj An adjacency matrix in which lines and columns correspond 
##' respectively to peptides and proteins.
##' @param expr A matrix of intensities of peptides
##' @param n The maximum number of peptides used to aggregate a protein.
##' @return A matrix of intensities of proteins
##' @author Alexia Dorffer
##' @examples
##' data(UPSpep25)
##' protID <- "Protein.group.IDs"
##' matAdj <- BuildAdjacencyMatrix(UPSpep25, protID, FALSE)
##' TopnPeptides(matAdj, Biobase::exprs(UPSpep25), 3)
TopnPeptides<-function(matAdj,expr,n){
    #require(foreach)
    ##register for use of parallel foreach
cores <- detectCores()
registerDoParallel(cores)

    #Get the n indices of peptides with the best median
    med <- apply(expr[rownames(matAdj),], 1, median)
    xmed <- matAdj * med
    
    k <- apply(xmed, 2, function (x) topMaxUsingPartialSortIndices(x, n)) 
    j <- NULL
    t <- foreach (j=1:ncol(expr), .combine='cbind') %dopar% {   
        p <- mapply(function(x, y){x[y]}, 
                    as.data.frame(matAdj*expr[rownames(matAdj),j]), 
                    as.data.frame(k))
        z <- CountPep(p)
        c(colSums(p, na.rm=TRUE), colSums(z, na.rm=TRUE))
    }
    
    Mp <-t[1:ncol(matAdj),]
    colnames(Mp) <- colnames(expr)
    rownames(Mp) <- colnames(matAdj)
    pep <-t[-c(1:ncol(matAdj)),]
    colnames(pep) <- paste("nb.pep.used.", colnames(expr), sep="")
    rownames(pep) <- colnames(matAdj)
    
    
    
    ##       Mp <- matrix(c(rep(NA)), ncol=nrow(condition), nrow=ncol(M2), 
    ##dimnames=c(list(prot, condname)))

##     pep <- Mp
##     for (j in 1:ncol(Mp)){ 
##         med <- apply(expr, 1, median)
##         xmed <- M2 * med
##         #recupere les N indices des peptides ayant les meilleures medianes
##         k <- apply(xmed, 2,function (x) topMaxUsingPartialSortIndices(x, n))
##         #recupere les intensites correspondant aux indices
##         p <- matrix(expr[k,j], nrow=nrow(k), ncol=ncol(k)) 
##         z <- CountPep(p)
##         pep[,j] <- colSums(z, na.rm=TRUE)
##         Mp[,j] <- colSums(p)
##     }
    res <- list("idprot" = colnames(matAdj), "matfin"=Mp, "nbpep"=pep)

    ##           }
    
    return(res)  
}

##' Method to create a binary matrix with proteins in columns and peptides 
##' in lines on a MSnSet object (peptides)
##' 
##' @title Function matrix of appartenance group
##' @param obj.pep An object (peptides) of class \code{\link{MSnbase}}.
##' @param protID The name of proteins ID column 
##' @param unique A boolean to indicate whether only the unique peptides must 
##' be considered (TRUE) or if the shared peptides have to 
##' be integrated (FALSE).
##' @return A binary matrix  
##' @author Florence Combes, Samuel Wieczorek, Alexia Dorffer
##' @examples
##' data(UPSpep25) 
##' BuildAdjacencyMatrix(UPSpep25, "Protein.group.IDs")
BuildAdjacencyMatrix <- function(obj.pep, protID, unique=TRUE){

data <- Biobase::exprs(obj.pep)
    PG <- Biobase::fData(obj.pep)[,protID]
    PG.l <- strsplit(as.character(PG), split=";", fixed=TRUE)

    X <- matrix(c(rep(0)), 
                nrow = nrow(data), 
                ncol = length(unique(unlist(PG.l))),
                dimnames=list(rownames(data),unique(unlist(PG.l))))
    
    
    for (i in 1:nrow(data)){
    X[i, as.character(PG.l[[i]])] <- 1
    }

    
    
    if (unique == TRUE){
        X <- X[which(rowSums(X)==1),]
        X <- X[,which(colSums(X)>0)]
    }

return(X)
}

##' Method to return the indices of the n higher values in the vector
##' 
##' @title  Function to return the indices of the n higher values in the vector
##' @param  x A vector of numeric values
##' @param  n The number of values to be returned
##' @return A vector of the indices of the n highest values  
##' @author Alexia Dorffer
##' @examples topMaxUsingPartialSortIndices(c(1:10), 3)
topMaxUsingPartialSortIndices <- function(x, n) {
    v <- order(x,decreasing=TRUE)[1:n]
    return(v)
}




##' Method to agregate with a method peptides to proteins on
##' a MSnSet object (peptides)
##' 
##' @title Function agregate peptides to proteins 
##' @param obj.pep An object (peptides) of class \code{\link{MSnbase}}.
##' @param protID The name of proteins ID column 
##' @param method The method used to aggregate the peptides into proteins.
##' Values are "sum", "mean" or "sum on top n" : do the sum / mean of intensity
##' on all peptides belonging to proteins. Default is "sum"
##' @param matAdj An adjacency matrix
##' @param n The number of peptides considered for the aggregation.
##' @return An object of class \code{\link{MSnbase}} with proteins  
##' @author Alexia Dorffer, Samuel Wieczorek
##' @examples 
##' data(UPSpep25)
##' protID <- "Protein.group.IDs"
##' mat <- BuildAdjacencyMatrix(UPSpep25, protID, TRUE)
##' pepAgregate(UPSpep25, protID, "sum", mat)
pepAgregate <- function (obj.pep, protID, method="sum", matAdj=NULL, n=NULL){
    #Check the validity of parameters
    parammethod <- c("sum overall", "mean", "sum on top n") 
    if (sum(is.na(match(method, parammethod) == TRUE)) > 0){return (NULL)}
    if (is.null(matAdj)){warning("Adjacency matrix is missing.")
                            return (NULL)}

    if (!is.na(match(method, "sum on top n")) && is.null(n)){
        warning("With the top n method, the parameter n must not be NULL.")
    return (NULL)}
    
    condname <- Biobase::pData(obj.pep)$Experiment
    condition <- Biobase::pData(obj.pep)
    expr <- 2^(Biobase::exprs(obj.pep))

    if(method == "sum overall"){ res <- SumPeptides(matAdj, expr)}
    else if  (method == "mean"){ res <- MeanPeptides(matAdj, expr)}
    else if (method == "sum on top n"){ res <- TopnPeptides(matAdj, expr, n) }
    
    Mp <- res$matfin
    Mp[Mp == 0] <- NA
    Mp[is.nan(Mp)] <- NA
    Mp[is.infinite(Mp)] <-NA
    
    
    
    pep <- res$nbpep
    protId <- res$idprot
    fd <- data.frame(protId, pep)
    
    obj <- MSnSet(exprs = log2(Mp), fData = fd, pData = Biobase::pData(obj.pep))
    obj@experimentData@other  <- list(obj@experimentData@other,
                                    typeOfData ="protein")
    
    return(obj)
}



##' Method to create a plot with proteins and peptides on
##' a MSnSet object (peptides)
##' 
##' @title Function to create a histogram that shows the repartition of
##' peptides w.r.t. the proteins
##' @param mat An adjacency matrix.
##' @return A histogram  
##' @author Alexia Dorffer, Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' mat <- BuildAdjacencyMatrix(UPSpep25, "Protein.group.IDs")
##' GraphPepProt(mat)
GraphPepProt <- function(mat){
    if (is.null(mat)){return (NULL)} 

    #mat <- as.matrix(mat)
    t <- t(mat)
    t <- apply(mat, 2, sum, na.rm=TRUE)
    tab <- table(t)
    position <- seq(1, length(tab),by=3)
    label <- names(tab)

    #par(mar=c(6,4,4,8) + 0.1)#, mgp=c(3,0.5,0)
    barplot(tab, 
            xlim=c(1, length(tab)),
            xlab="Nb of peptides", ylab="Nb of proteins",
            names.arg=label, xaxp=c(1, length(tab), 3), las=1
            , col = "orange")

}




##' Method to create a binary matrix with proteins in columns and peptides 
##' in lines on a MSnSet object (peptides)
##' 
##' @title Function matrix of appartenance group
##' @param obj.pep An object (peptides) of class \code{\link{MSnbase}}.
##' @param protID The name of proteins ID column 
##' @param unique A boolean to indicate whether only the unique peptides must 
##' be considered (TRUE) or if the shared peptides have to 
##' be integrated (FALSE).
##' @return A binary matrix  
##' @author Florence Combes, Samuel Wieczorek, Alexia Dorffer
##' @examples
##' data(UPSpep25) 
##' BuildSparseAdjacencyMatrix(UPSpep25, "Protein.group.IDs", TRUE)
BuildSparseAdjacencyMatrix <- function(obj.pep, protID, unique=TRUE){
    
    data <- Biobase::exprs(obj.pep)
    PG <- Biobase::fData(obj.pep)[,protID]
    PG.l <- strsplit(as.character(PG), split=";", fixed=TRUE)
    
    Un1 <- unlist(PG.l)
    X<- sparseMatrix(i = rep(seq_along(PG.l), lengths(PG.l)),
                     j=as.integer(factor(Un1, levels = unique(Un1))),
                     x=1, 
                     dimnames=list(rownames(data),unique(Un1)))
    
    if (unique == TRUE){
        X <- X[which(rowSums(as.matrix(X))==1),]
        X <- X[,which(colSums(as.matrix(X))>0)]
    }
    
    return(X)
}



##' This function computes the intensity of proteins based on the sum of the 
##' intensities of their peptides.
##' 
##' @title Compute the intensity of proteins with the sum of the intensities
##' of their peptides.
##' @param matAdj An adjacency matrix in which lines and columns correspond 
##' respectively to peptides and proteins.
##' @param expr A matrix of intensities of peptides
##' @return A matrix of intensities of proteins
##' @author Alexia Dorffer
##' @examples
##' data(UPSpep25)
##' protID <- "Protein.group.IDs"
##' M <- BuildSparseAdjacencyMatrix(UPSpep25, protID, FALSE)
##' SumPeptidesSpeedUp(M, Biobase::exprs(UPSpep25))
SumPeptidesSpeedUp <- function(matAdj, expr){
    expr <- expr[rownames(matAdj),]
    Mp <- t(matAdj) %*% expr
    .temp <- expr
    .temp[!is.na(.temp)] <- 1
    .temp[is.na(.temp)] <- 0
    pep <- t(matAdj) %*% .temp
    
    colnames(pep) <- paste("nb.pep.used.", colnames(expr), sep="")
    rownames(pep) <- colnames(matAdj)
    
    res <- list("idprot" = colnames(matAdj), "matfin"=Mp, "nbpep"=pep)
    return(res)
    
}

##' This function computes the intensity of proteins as the mean of the 
##' intensities of their peptides.
##' 
##' @title Compute the intensity of proteins as the mean of the intensities
##' of their peptides.
##' @param matAdj An adjacency matrix in which lines and columns correspond 
##' respectively to peptides and proteins.
##' @param expr A matrix of intensities of peptides
##' @return A matrix of intensities of proteins
##' @author Alexia Dorffer
##' @examples
##' data(UPSpep25)
##' protID <- "Protein.group.IDs"
##' matAdj <- BuildSparseAdjacencyMatrix(UPSpep25, protID, FALSE)
##' MeanPeptidesSpeedUp(matAdj, Biobase::exprs(UPSpep25))
MeanPeptidesSpeedUp <- function(matAdj, expr){
    expr <- expr[rownames(matAdj),]
    Mp <- t(matAdj) %*% expr
    .temp <- expr
    .temp[!is.na(.temp)] <- 1
    .temp[is.na(.temp)] <- 0
    pep <- t(matAdj) %*% .temp
    Mp <- Mp / pep
    
    colnames(pep) <- paste("nb.pep.used.", colnames(expr), sep="")
    rownames(pep) <- colnames(matAdj)
    
    res <- list("idprot" = colnames(matAdj), "matfin"=Mp, "nbpep"=pep)
    return(res)
    
}


##' Method to agregate with a method peptides to proteins on
##' a MSnSet object (peptides)
##' 
##' @title Function agregate peptides to proteins 
##' @param obj.pep An object (peptides) of class \code{\link{MSnbase}}.
##' @param protID The name of proteins ID column 
##' @param method The method used to aggregate the peptides into proteins.
##' Values are "sum", "mean" or "sum on top n" : do the sum / mean of intensity
##' on all peptides belonging to proteins. Default is "sum"
##' @param matAdj An adjacency matrix
##' @param n The number of peptides considered for the aggregation.
##' @return An object of class \code{\link{MSnbase}} with proteins  
##' @author Alexia Dorffer, Samuel Wieczorek
##' @examples 
##' data(UPSpep25)
##' protID <- "Protein.group.IDs"
##' mat <- BuildSparseAdjacencyMatrix(UPSpep25, protID, TRUE)
##' pepAgregateSpeedUp(UPSpep25, protID, "sum overall", mat)
pepAgregateSpeedUp <- function (obj.pep, protID, method="sum overall", matAdj=NULL, n=NULL){
    #Check the validity of parameters
    parammethod <- c("sum overall", "mean", "sum on top n") 
    if (sum(is.na(match(method, parammethod) == TRUE)) > 0){return (NULL)}
    if (is.null(matAdj)){warning("Adjacency matrix is missing.")
        return (NULL)}
    
    if (!is.na(match(method, "sum on top n")) && is.null(n)){
        warning("With the top n method, the parameter n must not be NULL.")
        return (NULL)}
    
    condname <- Biobase::pData(obj.pep)$Experiment
    condition <- Biobase::pData(obj.pep)
    expr <- 2^(Biobase::exprs(obj.pep))
    
    if(method == "sum overall"){ res <- SumPeptidesSpeedUp(matAdj, expr)}
    else if  (method == "mean"){ res <- MeanPeptidesSpeedUp(matAdj, expr)}
    else if (method == "sum on top n"){ res <- TopnPeptides(matAdj, expr, n) }
    
    Mp <- as.matrix(res$matfin)
    Mp[Mp == 0] <- NA
    Mp[is.nan(Mp)] <- NA
    Mp[is.infinite(Mp)] <-NA
    
    
    
    pep <- as.matrix(res$nbpep)
    protId <- res$idprot
    fd <- data.frame(protId, pep)
    
    obj <- MSnSet(exprs = log2(Mp), fData = fd, pData = Biobase::pData(obj.pep))
    obj@experimentData@other  <- list(obj@experimentData@other,
                                      typeOfData ="protein")
    
    return(obj)
}

