##' This function computes the number of proteins that are only defined by 
##' specific peptides, shared peptides or a mixture of two. 
##' 
##' @title computes the number of proteins that are only defined by 
##' specific peptides, shared peptides or a mixture of two.
##' @param matUnique The adjacency matrix with only specific peptides.
##' @param matShared The adjacency matrix with both specific and 
##' shared peptides.
##' @return A list
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' protID <- "Protein.group.IDs"
##' MShared <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], protID, FALSE)
##' MUnique <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], protID, TRUE)
##' getProteinsStats(MUnique,MShared)
getProteinsStats <- function(matUnique, matShared){
    if (is.null(matUnique) || is.null(matShared) ||
        !is.matrix(matUnique) || !is.matrix(matShared)){return(NULL)}
    
    t <- setdiff(union(rownames(matUnique), rownames(matShared)), 
                 intersect(rownames(matUnique), rownames(matShared)))
    sharedPeptides <- matShared[t,]
    sharedPeptides <- sharedPeptides[,-which(colSums(sharedPeptides)==0)]
    protOnlyUnique <- setdiff(union(colnames(sharedPeptides), 
                                    colnames(matShared)), 
                              intersect(colnames(sharedPeptides), 
                                        colnames(matShared)))
    
    protOnlyShared <- setdiff(union(colnames(matUnique), 
                                    colnames(matShared)), 
                              intersect(colnames(matUnique), 
                                        colnames(matShared)))
    a <- union(protOnlyUnique, protOnlyShared)
    b <- colnames(matShared)
    protMix <- setdiff(union(union(protOnlyUnique, protOnlyShared),
                             colnames(matShared)),
                       intersect(union(protOnlyUnique, protOnlyShared),
                                 colnames(matShared)))

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
##' @param proteinNames The names of the protein in the new dataset (i.e. rownames)
##' @return A vector
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' protID <- "Protein.group.IDs"
##' obj.pep <- Exp1_R25_pept[1:1000]
##' M <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
##' data <- Biobase::fData(obj.pep)
##' protData <- DAPAR::aggregateMean(obj.pep, M)
##' name <- "Protein.group.IDs"
##' proteinNames <- rownames(Biobase::fData(protData))
##' BuildColumnToProteinDataset(data, M, name,proteinNames )
BuildColumnToProteinDataset <- function(peptideData, matAdj, columnName, proteinNames){
nbProt <- ncol(matAdj)
newCol <- rep("", nbProt)

#print(head(rownames(peptideData)))
i <- 1
for (p in proteinNames){
    listeIndicePeptides <- names(which(matAdj[,p] == 1))
    listeData <- unique(as.character(peptideData[listeIndicePeptides,columnName], ";"))
    newCol[i] <- paste0(listeData, collapse = ", ")
    i <- i +1
}
return(newCol)
}


##' This function creates a column for the protein dataset after agregation 
##' by using the previous peptide dataset. It is a parallel version of the function
##' \code{BuildColumnToProteinDataset}
##' 
##' @title creates a column for the protein dataset after agregation by
##'  using the previous peptide dataset.
##' @param peptideData A data.frame of meta data of peptides. It is the fData 
##' of the MSnset object.
##' @param matAdj The adjacency matrix used to agregate the peptides data.
##' @param columnName The name of the column in fData(peptides_MSnset) that 
##' the user wants to keep in the new protein data.frame.
##' @param proteinNames The names of the protein in the new dataset (i.e. rownames)
##' @return A vector
##' @author Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' protID <- "Protein.group.IDs"
##' obj.pep <- Exp1_R25_pept[1:1000]
##' M <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
##' data <- Biobase::fData(obj.pep)
##' protData <- DAPAR::aggregateSum(obj.pep, M)
##' name <- "Protein.group.IDs"
##' proteinNames <- rownames(Biobase::fData(protData))
##' BuildColumnToProteinDataset_par(data, M, name,proteinNames )
BuildColumnToProteinDataset_par <- function(peptideData, matAdj, columnName, proteinNames){
    doParallel::registerDoParallel()
    
    nbProt <- ncol(matAdj)
    newCol <- rep("", nbProt)
    i <- 1
    newCol <- foreach (i=1:length(proteinNames), .combine=rbind) %dopar% {
        listeIndicePeptides <- names(which(matAdj[,proteinNames[i]] == 1))
        listeData <- unique(as.character(peptideData[listeIndicePeptides,columnName], ";"))
        paste0(listeData, collapse = ", ")
    }
    return(as.vector(newCol))
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
##' library(DAPARdata)
##' data(Exp1_R25_pept)
##' protID <- "Protein.group.IDs"
##' M <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], protID, FALSE)
##' CountPep(M)
CountPep <- function (M) {
    z <- M
    z[z!=0] <- 1
    return(z)
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
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' mat <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], "Protein.group.IDs")
##' GraphPepProt(mat)
GraphPepProt <- function(mat){
    if (is.null(mat)){return (NULL)} 

    #mat <- as.matrix(mat)
    t <- t(mat)
    t <- apply(mat, 2, sum, na.rm=TRUE)
    tab <- table(t)
    position <- seq(1, length(tab),by=3)
    conds <- names(tab)

    #par(mar=c(6,4,4,8) + 0.1)#, mgp=c(3,0.5,0)
    barplot(tab, 
            xlim=c(1, length(tab)),
            xlab="Nb of peptides", 
            ylab="Nb of proteins",
            names.arg=conds, 
            xaxp=c(1, length(tab), 3), 
            las=1
            , col = "orange")

}




##' Method to create a binary matrix with proteins in columns and peptides 
##' in lines on a \code{MSnSet} object (peptides)
##' 
##' @title Function matrix of appartenance group
##' @param obj.pep An object (peptides) of class \code{MSnSet}.
##' @param protID The name of proteins ID column 
##' @param unique A boolean to indicate whether only the unique peptides must 
##' be considered (TRUE) or if the shared peptides have to 
##' be integrated (FALSE).
##' @return A binary matrix  
##' @author Florence Combes, Samuel Wieczorek, Alexia Dorffer
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept) 
##' BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], "Protein.group.IDs", TRUE)
BuildAdjacencyMatrix <- function(obj.pep, protID, unique=TRUE){
    
    data <- Biobase::exprs(obj.pep)
    PG <- Biobase::fData(obj.pep)[,protID]
    PG.l <- strsplit(as.character(PG), split=";", fixed=TRUE)
    
    Un1 <- unlist(PG.l)
    X<- sparseMatrix(i = rep(seq_along(PG.l), lengths(PG.l)),
                     j=as.integer(factor(Un1, levels = unique(Un1))),
                     x=1, 
                     dimnames=list(rownames(data),as.character(unique(Un1))))
    
    if (unique == TRUE){
        X[which(rowSums(as.matrix(X))>1),] <- 0
         }
    
    return(X)
}



##' This function computes the intensity of proteins based on the sum of the 
##' intensities of their peptides.
##' 
##' @title Compute the intensity of proteins with the sum of the intensities
##' of their peptides.
##' @param obj.pep A matrix of intensities of peptides
##' @param X An adjacency matrix in which lines and columns correspond 
##' respectively to peptides and proteins.
##' @param ... xxx
##' @return A matrix of intensities of proteins
##' @author Alexia Dorffer
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' protID <- "Protein.group.IDs"
##' obj.pep <- Exp1_R25_pept[1:1000]
##' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
##' DAPAR::aggregateSum(obj.pep, X)
aggregateSum <- function(obj.pep, X,...){
  qData <- 2^(Biobase::exprs(obj.pep))
  #qData <- expr[rownames(X),]
  Mp <- inner.sum(qData, X)
  obj.prot <- finalizeAggregation(obj.pep, qData, X, Mp, ...)
  return(obj.prot)
}



##' Method to xxxxx
##' 
##' @title xxxx 
##' @param obj.pep xxxxx
##' @param X xxxx
##' @param init.method xxxxx
##' @param method xxxxx
##' @param ... xxxx
##' @return xxxxx
##' @author Samuel Wieczorek
aggregateIterParallel <- function(obj.pep, X, init.method='sum', method='mean', ...){
  doParallel::registerDoParallel()
  ### a reproduire iterativement pour chaque condition
  # Initialisation: presque aucune dépendance à l'initialisation prendre "sum overall" et  matAdj = X par simplicité
  #X <- as.matrix(X)
  qData.pep <- 2^(Biobase::exprs(obj.pep))
  #qData.pep <- expr[rownames(X),]
  
  finalX <- matrix(rep(0,ncol(X)*ncol(obj.pep)), nrow=ncol(X))
  
  finalX <- foreach (cond=1:length(unique(Biobase::pData(obj.pep)$Condition)), 
                     .combine=cbind,.packages = "MSnbase") %dopar% {
    condsIndices <- which(Biobase::pData(obj.pep)$Condition == unique(Biobase::pData(obj.pep)$Condition)[cond])
    qData <- qData.pep[,condsIndices]
    #print(paste0("Condition ", cond))
    DAPAR::inner.aggregate.iter(qData, X, init.method, method)
   }
  
  
  finalX <- finalX[,colnames(Biobase::exprs(obj.pep))]
  
  obj.prot <- finalizeAggregation(obj.pep, qData, X, finalX, ...)
  return(obj.prot)
  
  #return(yprot)
}


##' Method to xxxxx
##' 
##' @title xxxx 
##' @param qData xxxxx
##' @param X xxxx
##' @param init.method xxx
##' @param method xxx
##' @param ... xxxx
##' @return xxxxx
##' @author Samuel Wieczorek
inner.aggregate.iter <- function(qData, X,...,init.method, method){
  yprot <- NULL
  switch(init.method,
         sum= yprot <- inner.sum(qData, X),
         mean= yprot <- inner.mean(qData, X)
  )
  conv <- 1
  
  while(conv > 10**(-10)){
    mean.prot <- rowMeans(as.matrix(yprot), na.rm = TRUE)
    mean.prot[is.na(mean.prot)] <- 0
    
    X.tmp <- mean.prot*X
    X.new <- X.tmp/rowSums(as.matrix(X.tmp), na.rm = TRUE)
    X.new[is.na(X.new)] <- 0
    
    method <- 'mean'
    # l'appel à la fonction ci-dessous dépend des paramètres choisis par l'utilisateur
    switch(method,
           mean = yprot <- inner.mean(qData, X.new),
           mean.topn = yprot <- DAPAR::aggregateTopn(qData,X.new,n=n, method=method)
    )
    
    mean.prot.new <- rowMeans(as.matrix(yprot), na.rm = TRUE)
    mean.prot.new[is.na(mean.prot.new)] <- 0
    
    conv <- mean(abs(mean.prot.new - mean.prot))
    print(paste0("conv : ", conv))
  }
  return(as.matrix(yprot))
}



##' Method to xxxxx
##' 
##' @title xxxx 
##' @param obj.pep xxxxx
##' @param X xxxx
##' @param init.method xxxxx
##' @param method xxxxx
##' @param ... xxxx
##' @return xxxxx
##' @author Samuel Wieczorek
##' @example 
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' protID <- "Protein.group.IDs"
##' X <- BuildAdjacencyMatrix(Exp1_R25_pept[1:1000], protID, FALSE)
##' aggregateIter(Exp1_R25_pept[1:1000],X=X)
aggregateIter <- function(obj.pep, X, init.method='sum', method='mean',  ...){
  
  ### a reproduire iterativement pour chaque condition
    # Initialisation: presque aucune dépendance à l'initialisation prendre "sum overall" et  matAdj = X par simplicité
    #X <- as.matrix(X)
  qData.pep <- 2^(Biobase::exprs(obj.pep))
    #qData.pep <- expr[rownames(X),]
  
    finalX <- matrix(rep(0,ncol(X)*ncol(obj.pep)), nrow=ncol(X))
    for (cond in unique(pData(obj.pep)$Condition)){
      condsIndices <- which(pData(obj.pep)$Condition == cond)
      qData <- qData.pep[,condsIndices]
      print(paste0("Condition ", cond))
    finalX[,condsIndices]  <- inner.aggregate.iter(qData, X, init.method, method)
     }
    
    obj.prot <- finalizeAggregation(obj.pep, qData, X, finalX, ...)
    return(obj.prot)
    
  #return(yprot)
}



##' Method to xxxxx
##' 
##' @title xxxx 
##' @param qData xxxxx
##' @param X xxxx
##' @return xxxxx
##' @author Samuel Wieczorek
GetNbPeptidesUsed <- function(X, qData){
   qData[!is.na(qData)] <- 1
  qData[is.na(qData)] <- 0
  pep <- t(X) %*% qData
  
  return(pep)
}

##' This function computes the intensity of proteins as the mean of the 
##' intensities of their peptides.
##' 
##' @title Compute the intensity of proteins as the mean of the intensities
##' of their peptides.
##' @param obj.pep A matrix of intensities of peptides
##' @param X An adjacency matrix in which lines and columns correspond 
##' respectively to peptides and proteins.
##' @param ... xxxxx
##' @return A matrix of intensities of proteins
##' @author Alexia Dorffer
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj.pep <- Exp1_R25_pept[1:1000]
##' protID <- "Protein.group.IDs"
##' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
##' aggregateMean(obj.pep, X)
aggregateMean <- function(obj.pep, X, ...){
  qData <- 2^(Biobase::exprs(obj.pep))
  #qData <- expr[rownames(X),]
  
  
  Mp <- inner.mean(qData, X)
   
    obj.prot <- finalizeAggregation(obj.pep, qData, X, Mp, ...)
    
    return(obj.prot)
    
}

##' Method to xxxxx
##' 
##' @title xxxx 
##' @param qData xxxxx
##' @param X xxxx
##' @return xxxxx
##' @author Samuel Wieczorek
inner.sum <- function(qData, X){
  qData[is.na(qData)] <- 0
  Mp <- t(X) %*% qData
  return(Mp)
}


##' Method to xxxxx
##' 
##' @title xxxx 
##' @param qData xxxxx
##' @param X xxxx
##' @return xxxxx
##' @author Samuel Wieczorek
inner.mean <- function(qData, X){
  Mp <- inner.sum(qData, X)
  Mp <- Mp / GetNbPeptidesUsed(X, qData)
  
  return(Mp)
  
}




##' Method to xxxxx
##' 
##' @title xxxx 
##' @param qData xxxxx
##' @param X xxxx
##' @param n xxxxx
##' @param method xxxxx
##' @return xxxxx
##' @author Samuel Wieczorek
inner.aggregate.topn <-function(qData,X, n, method='mean'){
  #qData <- expr[rownames(X),]
  #qData[is.na(qData)] <- 0
  
  med <- apply(qData, 1, median)
  xmed <- as(X * med, "dgCMatrix")
  for (c in 1:ncol(X)){
    v <- order(xmed[,c],decreasing=TRUE)[1:n]
    l <- v[which((xmed[,c])[v] != 0)]
    
    if (length(l) > 0){
      diff <- setdiff( which(X[,c] == 1), l)
      if (length(diff)) {X[diff,c] <- 0}
    }
  }
  
  Mp <- NULL
  switch(method,
         mean= Mp <- inner.mean(qData, X),
         sum= Mp <- inner.sum(qData, X)
  )
  
  return(Mp)
}

##' This function computes the intensity of proteins as the sum of the 
##' intensities of their n best peptides.
##' 
##' @title Compute the intensity of proteins as the sum of the 
##' intensities of their n best peptides.
##' @param obj.pep A matrix of intensities of peptides
##' @param X An adjacency matrix in which lines and columns correspond 
##' respectively to peptides and proteins.
##' @param n The maximum number of peptides used to aggregate a protein.
##' @param method xxx
##' @param ... xxx
##' @return A matrix of intensities of proteins
##' @author Alexia Dorffer, Samuel Wieczorek
##' @examples
##' require(DAPARdata)
##' data(Exp1_R25_pept)
##' obj.pep <- Exp1_R25_pept[1:1000]
##' protID <- "Protein.group.IDs"
##' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
##' DAPAR::aggregateTopn(obj.pep, X, n=3)
aggregateTopn <- function(obj.pep,X, n=10, method='mean', ...){
    
  #Get the indices of the n peptides with the best median
  qData <- 2^(Biobase::exprs(obj.pep))
  #qData <- expr[rownames(X),]
  
  finalX <- inner.aggregate.topn(qData, X, n, method='mean')
  
  obj.prot <- finalizeAggregation(obj.pep, qData, X, finalX, ...)
  return(obj.prot)
}




##' Method to xxxxx
##' 
##' @title xxxx 
##' @param obj.pep xxxxx
##' @param qData xxxx
##' @param X xxxxx
##' @param finalX xxxxx
##' @param lib.loc xxxx
##' @return xxxxx
##' @author Samuel Wieczorek
finalizeAggregation <- function(obj.pep, qData, X, finalX, lib.loc=NULL){
 
  finalX <- as.matrix(finalX)
  finalX[finalX==0] <- NA
  finalX[is.nan(finalX)] <- NA
  finalX[is.infinite(finalX)] <-NA
  
  
  pep <- as.matrix(GetNbPeptidesUsed(X, qData))
  colnames(pep) <- paste("nb.pep.used.", colnames(qData), sep="")
  rownames(pep) <- colnames(X)
  
   fd <- data.frame(colnames(X), pep)
  
  obj.prot <- MSnSet(exprs = log2(finalX), 
                fData = fd, 
                pData = Biobase::pData(obj.pep))
  obj.prot@experimentData@other  <- list(obj.prot@experimentData@other, typeOfData ="protein")
  obj.prot <- addOriginOfValue(obj.prot)
  obj.prot@experimentData@other$Prostar_Version <- installed.packages(lib.loc = lib.loc$Prostar.loc)["Prostar","Version"]
  obj.prot@experimentData@other$DAPAR_Version <- installed.packages(lib.loc = lib.loc$DAPAR.loc)["DAPAR","Version"]
  
  #obj.prot <- xxxx
  
  return (obj.prot)
}

