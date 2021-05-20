#' This function computes the number of proteins that are only defined by 
#' specific peptides, shared peptides or a mixture of two. 
#' 
#' @title Computes the number of proteins that are only defined by 
#' specific peptides, shared peptides or a mixture of two.
#' 
#' @param matShared The adjacency matrix with both specific and 
#' shared peptides.
#' 
#' @return A list
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' protID <- "Protein_group_IDs"
#' obj <- Exp1_R25_pept[1:20]
#' MShared <- BuildAdjacencyMatrix(obj, protID, FALSE)
#' getProteinsStats(matShared=MShared)
#' 
#' @export
#' 
getProteinsStats <- function(matShared){
  if(missing(matShared))
    stop("'matShared' is needed.")
  if (is.null(matShared))
    stop("'matShared' is NULL")
  
  nbPeptide <- 0
  
  ind.shared.Pep <- which(rowSums(as.matrix(matShared))>1)
  M.shared.Pep <- matShared[ind.shared.Pep,]
  if (length(ind.shared.Pep)==1){
    j <- which(as.matrix(M.shared.Pep)==0)
    M.shared.Pep <- M.shared.Pep[-j]
    pep.names.shared <- names(M.shared.Pep)
  } else {
    j <- which(colSums(as.matrix(M.shared.Pep))==0)
    M.shared.Pep <- M.shared.Pep[,-j]
    pep.names.shared <- colnames(M.shared.Pep)
  }
  
  
  ind.unique.Pep <- which(rowSums(as.matrix(matShared))==1)
  M.unique.Pep <- matShared[ind.unique.Pep,]
  if (length(ind.unique.Pep)==1){
    j <- which(as.matrix(M.unique.Pep)==0)
    M.unique.Pep <- M.unique.Pep[-j]
    pep.names.unique <- names(M.unique.Pep)
  } else {
    j <- which(colSums(as.matrix(M.unique.Pep))==0)
    M.unique.Pep <- M.unique.Pep[,-j]
    pep.names.unique <- colnames(M.unique.Pep)
  }
   

  
  protOnlyShared <- setdiff(pep.names.shared, intersect(pep.names.shared, pep.names.unique))
  protOnlyUnique <- setdiff(pep.names.unique, intersect(pep.names.shared, pep.names.unique))
  protMix <- intersect(pep.names.shared, pep.names.unique)
  
  
  
  return (list(nbPeptides = length(ind.unique.Pep) + length(ind.shared.Pep),
               nbSpecificPeptides = length(ind.unique.Pep),
               nbSharedPeptides = length(ind.shared.Pep),
               nbProt = length(protOnlyShared)+length(protOnlyUnique)+length(protMix),
               protOnlyUniquePep = protOnlyUnique,
               protOnlySharedPep = protOnlyShared,
               protMixPep = protMix))
}




#' This function creates a column for the protein dataset after aggregation 
#' by using the previous peptide dataset.
#' 
#' @title creates a column for the protein dataset after agregation by
#'  using the previous peptide dataset.
#'  
#' @param peptideData A data.frame of meta data of peptides. It is the fData 
#' of the MSnset object.
#' 
#' @param matAdj The adjacency matrix used to agregate the peptides data.
#' 
#' @param columnName The name of the column in fData(peptides_MSnset) that 
#' the user wants to keep in the new protein data.frame.
#' 
#' @param proteinNames The names of the protein in the new dataset 
#' (i.e. rownames)
#' 
#' @return A vector
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' protID <- "Protein_group_IDs"
#' obj.pep <- Exp1_R25_pept[1:10]
#' M <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' data <- Biobase::fData(obj.pep)
#' protData <- DAPAR::aggregateMean(obj.pep, M)
#' name <- "Protein_group_IDs"
#' proteinNames <- rownames(Biobase::fData(protData$obj.prot))
#' BuildColumnToProteinDataset(data, M, name, proteinNames )
#' 
#' @export
#' 
BuildColumnToProteinDataset <- function(peptideData, 
                                        matAdj, 
                                        columnName, 
                                        proteinNames){
  nbProt <- ncol(matAdj)
  newCol <- rep("", nbProt)
  i <- 1
  for (p in proteinNames){
    listeIndicePeptides <- names(which(matAdj[,p] == 1))
    listeData <- unique(as.character(peptideData[listeIndicePeptides,columnName], ";"))
    newCol[i] <- paste0(listeData, collapse = ", ")
    i <- i +1
  }
  return(newCol)
}


#' This function creates a column for the protein dataset after agregation 
#' by using the previous peptide dataset. It is a parallel version of 
#' the function \code{BuildColumnToProteinDataset}
#' 
#' @title creates a column for the protein dataset after agregation by
#'  using the previous peptide dataset.
#'  
#' @param peptideData A data.frame of meta data of peptides. It is the fData 
#' of the MSnset object.
#' 
#' @param matAdj The adjacency matrix used to agregate the peptides data.
#' 
#' @param columnName The name of the column in fData(peptides_MSnset) that 
#' the user wants to keep in the new protein data.frame.
#' 
#' @param proteinNames The names of the protein in the new dataset 
#' (i.e. rownames)
#' 
#' @return A vector
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' protID <- "Protein_group_IDs"
#' obj.pep <- Exp1_R25_pept[1:10]
#' M <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' data <- Biobase::fData(obj.pep)
#' protData <- DAPAR::aggregateSum(obj.pep, M)
#' name <- "Protein_group_IDs"
#' proteinNames <- rownames(Biobase::fData(protData$obj.prot))
#' BuildColumnToProteinDataset_par(data, M, name,proteinNames )
#' }
#' 
#' @export
#' 
#' @import doParallel 
#' @import foreach
#' 
BuildColumnToProteinDataset_par <- function(peptideData, 
                                            matAdj, 
                                            columnName, 
                                            proteinNames){
  
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  
  # https://cran.r-project.org/web/packages/policies.html
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores 
    num_workers <- 2L
  } else {
    # use all cores in devtools::test()
    num_workers <- parallel::detectCores()
  }
  
  doParallel::registerDoParallel(cores=num_workers)
  
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



#' This function computes the number of peptides used to aggregate proteins.
#' 
#' @title Compute the number of peptides used to aggregate proteins
#' 
#' @param M A "valued" adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' 
#' @return A vector of boolean which is the adjacency matrix 
#' but with NA values if they exist in the intensity matrix.
#' 
#' @author Alexia Dorffer
#' 
#' @examples
#' library(DAPARdata)
#' utils::data(Exp1_R25_pept)
#' protID <- "Protein_group_IDs"
#' M <- BuildAdjacencyMatrix(Exp1_R25_pept[1:10], protID, FALSE)
#' CountPep(M)
#' 
#' @export
#' 
CountPep <- function (M) {
  z <- M
  z[z!=0] <- 1
  return(z)
}


#' Method to compute the number of quantified peptides used for aggregating 
#' each protein
#' 
#' @title Computes the number of peptides used for aggregating each protein 
#' 
#' @param X An adjacency matrix
#' 
#' @param pepData A data.frame of quantitative data
#' 
#' @return A data.frame
#' 
#' @author Samuel Wieczorek
#' 
#' @export
#' 
GetNbPeptidesUsed <- function(X, pepData){
  pepData[!is.na(pepData)] <- 1
  pepData[is.na(pepData)] <- 0
  pep <- t(X) %*% pepData
  
  return(pep)
}


#' Method to compute the detailed number of quantified peptides used for 
#' aggregating each protein
#' 
#' @title Computes the detailed number of peptides used for aggregating 
#' each protein 
#' 
#' @param X An adjacency matrix
#' 
#' @param qdata.pep A data.frame of quantitative data
#' 
#' @return A list of two items
#' 
#' @author Samuel Wieczorek
#' 
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj.pep <- Exp1_R25_pept[1:10]
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' ll.n <- GetDetailedNbPeptidesUsed(X, Biobase::exprs(obj.pep))
#' 
#' @export
#' 
GetDetailedNbPeptidesUsed <- function(X, qdata.pep){
  X <- as.matrix(X)
  qdata.pep[!is.na(qdata.pep)] <- 1
  qdata.pep[is.na(qdata.pep)] <- 0
  
  mat <- splitAdjacencyMat(X)
  return(list(nShared=t(mat$Xshared) %*% qdata.pep, 
              nSpec=t(mat$Xspec) %*% qdata.pep))
  
}


#' Method to compute the detailed number of quantified peptides for each protein
#' 
#' @title Computes the detailed number of peptides for each protein 
#' 
#' @param X An adjacency matrix
#' 
#' @return A data.frame
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj.pep <- Exp1_R25_pept[1:10]
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' n <- GetDetailedNbPeptides(X)
#' 
#' @export
#' 
GetDetailedNbPeptides <- function(X){
  
  mat <- splitAdjacencyMat(as.matrix(X))
  
  
  return(list(nTotal = rowSums(t(as.matrix(X))),
              nShared=rowSums(t(mat$Xshared)), 
              nSpec=rowSums(t(mat$Xspec)))
  )
  
}



#' Method to create a plot with proteins and peptides on
#' a MSnSet object (peptides)
#' 
#' @title Function to create a histogram that shows the repartition of
#' peptides w.r.t. the proteins
#' 
#' @param mat An adjacency matrix.
#' 
#' @return A histogram  
#' 
#' @author Alexia Dorffer, Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' mat <- BuildAdjacencyMatrix(Exp1_R25_pept[1:10], "Protein_group_IDs")
#' GraphPepProt(mat)
#' 
#' @export
#' 
GraphPepProt <- function(mat){
  if (is.null(mat)){return (NULL)} 
  
  mat <- as.matrix(mat)
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




#' Method to create a binary matrix with proteins in columns and peptides 
#' in lines on a \code{MSnSet} object (peptides)
#' 
#' @title Function matrix of appartenance group
#' 
#' @param obj.pep An object (peptides) of class \code{MSnSet}.
#' 
#' @param protID The name of proteins ID column 
#' 
#' @param unique A boolean to indicate whether only the unique peptides must 
#' be considered (TRUE) or if the shared peptides have to 
#' be integrated (FALSE).
#' 
#' @return A binary matrix  
#' 
#' @author Florence Combes, Samuel Wieczorek, Alexia Dorffer
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' BuildAdjacencyMatrix(Exp1_R25_pept[1:10], "Protein_group_IDs", TRUE)
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' 
BuildAdjacencyMatrix <- function(obj.pep, protID, unique=TRUE){
  
  data <- Biobase::exprs(obj.pep)
  PG <- Biobase::fData(obj.pep)[,protID]
  PG.l <- strsplit(as.character(PG), split=";", fixed=TRUE)
  
  t <- table(data.frame(A=rep(seq_along(PG.l), lengths(PG.l)),
                        B=unlist(PG.l)
                        )
             )
  
  if (unique == TRUE){
    ll <- which(rowSums(t)>1)
    if (length(ll) > 0) {
      t[ll,] <- 0
    }
    
  }
  
  X <- Matrix::Matrix(t, sparse=T,
                      dimnames = list(rownames(obj.pep), colnames(t))
  )
  
  return(X)
}



#' This function computes the intensity of proteins based on the sum of the 
#' intensities of their peptides.
#' 
#' @title Compute the intensity of proteins with the sum of the intensities
#' of their peptides.
#' 
#' @param obj.pep A matrix of intensities of peptides
#' 
#' @param X An adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' 
#' @return A matrix of intensities of proteins
#' 
#' @author Alexia Dorffer
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' protID <- "Protein_group_IDs"
#' obj.pep <- Exp1_R25_pept[1:10]
#' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' ll.agg <- DAPAR::aggregateSum(obj.pep, X)
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' 
aggregateSum <- function(obj.pep, X){
  obj.prot <- NULL
  
  # Agregation of metacell data
  metacell <- AggregateMetacell(X = X, obj.pep = obj.pep)
  if (!is.null(metacell$issues))
    return(list(obj.prot = NULL,
                issues = metacell$issues)
    )
  else {
    # Agregation of quanti data
    pepData <- 2^(Biobase::exprs(obj.pep))
    protData <- inner.sum(pepData, X)
    # Build protein dataset
    obj.prot <- finalizeAggregation(obj.pep, pepData, protData, metacell$metacell, X)
    return(list(obj.prot = obj.prot,
                issues = NULL))
  }
}



#' Method to xxxxx
#' 
#' @title xxxx 
#' 
#' @param obj.pep xxxxx
#' 
#' @param X xxxx
#' 
#' @param init.method xxxxx
#' 
#' @param method xxxxx
#' 
#' @param n xxxx
#' 
#' @return xxxxx
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' protID <- "Protein_group_IDs"
#' obj.pep <- Exp1_R25_pept[1:10]
#' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' aggregateIterParallel(obj.pep, X)
#' }
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' @importFrom doParallel registerDoParallel 
#' @importFrom foreach foreach
#' 
aggregateIterParallel <- function(obj.pep, X, init.method='Sum', method='Mean', n=NULL){
  
  registerDoParallel()
  obj.prot <- NULL
  
  # Step 1: Agregation of metacell data
  metacell <- AggregateMetacell(X = X, obj.pep = obj.pep)
  if (!is.null(metacell$issues))
    return(list(obj.prot = NULL,
                issues = metacell$issues)
    )
  else { # Step 2 : Agregation of quantitative data
    qData.pep <- 2^(Biobase::exprs(obj.pep))
    protData <- matrix(rep(0,ncol(X)*nrow(X)), nrow=ncol(X))
  
    protData <- foreach(cond=1:length(unique(Biobase::pData(obj.pep)$Condition)), 
                        .combine=cbind, 
                        .packages = "MSnbase") %dopar% {
    
    condsIndices <- which(Biobase::pData(obj.pep)$Condition == unique(Biobase::pData(obj.pep)$Condition)[cond])
    qData <- qData.pep[,condsIndices]
    DAPAR::inner.aggregate.iter(qData, X, init.method, method, n)
    }
  
    protData <- protData[,colnames(Biobase::exprs(obj.pep))]
  
    # Step 3 : Build the protein dataset
    obj.prot <- DAPAR::finalizeAggregation(obj.pep, qData.pep, protData, metacell$metacell, X)
    return(list(obj.prot = obj.prot,
                issues = NULL))
    }
  
}


#' Method to xxxxx
#' 
#' @title xxxx 
#' 
#' @param pepData xxxxx
#' 
#' @param X xxxx
#' 
#' @param init.method xxx
#' 
#' @param method xxx
#' 
#' @param n xxxx
#' 
#' @return xxxxx
#' 
#' @author Samuel Wieczorek
#' 
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' protID <- "Protein_group_IDs"
#' obj.pep <- Exp1_R25_pept[1:10]
#' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' qdata.agg <- DAPAR::inner.aggregate.iter(exprs(obj.pep), X)
#' 
#' @export
#' 
inner.aggregate.iter <- function(pepData, 
                                 X,
                                 init.method='Sum', 
                                 method='Mean', 
                                 n=NULL){
  
  
  if (!(init.method %in% c("Sum", "Mean"))) {
    warning("Wrong parameter init.method")
    return(NULL)
  }
  
  if (!(method %in% c("onlyN", "Mean"))){
    warning("Wrong parameter method")
    return(NULL)
  }
  
  
  if (method=='onlyN' && is.null(n)){
    warning("Parameter n is null")
    return(NULL)
  }
  
  yprot <- NULL
  switch(init.method,
         Sum= yprot <- inner.sum(pepData, X),
         Mean= yprot <- inner.mean(pepData, X)
  )
  conv <- 1
  
  while(conv > 10**(-10)){
    mean.prot <- rowMeans(as.matrix(yprot), na.rm = TRUE)
    mean.prot[is.na(mean.prot)] <- 0
    
    X.tmp <- mean.prot*X
    X.new <- X.tmp/rowSums(as.matrix(X.tmp), na.rm = TRUE)
    X.new[is.na(X.new)] <- 0
    
    # l'appel ? la fonction ci-dessous d?pend des param?tres choisis par 
    # l'utilisateur
    switch(method,
           Mean = yprot <- inner.mean(pepData, X.new),
           onlyN = yprot <- inner.aggregate.topn(pepData,X.new,'Mean', n)
    )
    
    mean.prot.new <- rowMeans(as.matrix(yprot), na.rm = TRUE)
    mean.prot.new[is.na(mean.prot.new)] <- 0
    
    conv <- mean(abs(mean.prot.new - mean.prot))
  }
  return(as.matrix(yprot))
}



#' Method to xxxxx
#' 
#' @title xxxx 
#' 
#' @param obj.pep xxxxx
#' 
#' @param X xxxx
#' 
#' @param init.method xxxxx
#' 
#' @param method xxxxx
#' 
#' @param n xxxx
#' 
#' @return A protein object of class \code{MSnset}
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(Exp1_R25_pept[1:10], protID, FALSE)
#' ll.agg <- aggregateIter(Exp1_R25_pept[1:10],X=X)
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' 
aggregateIter <- function(obj.pep, 
                          X, 
                          init.method='Sum', 
                          method='Mean', 
                          n=NULL){
  
  obj.prot <- NULL
  
  # Step 1 : Agregation of metacell data
  metacell <- AggregateMetacell(X = X, obj.pep = obj.pep)
  
  if (!is.null(metacell$issues))
    return(list(obj.prot = NULL,
                issues = metacell$issues))
  else {
    # Step 2: Agregation of quantitative data
    # For each condition, reproduce iteratively
    # Initialisation: At initialization step, take "sum overall" and  matAdj = X
    # for simplicity. 
    # Note : X <- as.matrix(X)
    qData.pep <- 2^(Biobase::exprs(obj.pep))
    
    protData <- matrix(rep(0,ncol(X)*ncol(obj.pep)), nrow=ncol(X))
    
    for (cond in unique(Biobase::pData(obj.pep)$Condition)){
      condsIndices <- which(Biobase::pData(obj.pep)$Condition == cond)
      qData <- qData.pep[,condsIndices]
      protData[,condsIndices]  <- inner.aggregate.iter(qData, 
                                                       X, 
                                                       init.method, 
                                                       method, 
                                                       n)
    }
    
    # Step 3: Build the protein dataset
    obj.prot <- finalizeAggregation(obj.pep, qData.pep, protData, metacell$metacell, X)
    return(list(obj.prot = obj.prot,
                issues = NULL))
  }
}




#' This function computes the intensity of proteins as the mean of the 
#' intensities of their peptides.
#' 
#' @title Compute the intensity of proteins as the mean of the intensities
#' of their peptides.
#' 
#' @param obj.pep A peptide object of class \code{MSnset}
#' 
#' @param X An adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' 
#' @return A matrix of intensities of proteins
#' 
#' @author Alexia Dorffer
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj.pep <- Exp1_R25_pept[1:10]
#' obj.pep.imp <- wrapper.impute.detQuant(obj.pep, na.type='missing')
#' protID <- obj.pep@experimentData@other$proteinId
#' X <- BuildAdjacencyMatrix(obj.pep.imp, protID, FALSE)
#' ll.agg <- aggregateMean(obj.pep.imp, X)
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' 
aggregateMean <- function(obj.pep, X){
  obj.prot <- NULL
  # Agregation of metacell data
  cat("Aggregate metacell data...\n")
  metacell <- AggregateMetacell(X = X, obj.pep = obj.pep)
  if (!is.null(metacell$issues))
    return(list(obj.prot = NULL,
                issues = metacell$issues))
  else {
    # Step 2: Agregation of quantitative data
    cat("Computing quantitative data for proteins ...\n")
    pepData <- 2^(Biobase::exprs(obj.pep))
    protData <- inner.mean(pepData, as.matrix(X))
    
    # Step 3: Build protein dataset
    cat("Building the protein dataset...\n")
    obj.prot <- finalizeAggregation(obj.pep, pepData, protData, metacell$metacell, X)
    
    return(list(obj.prot = obj.prot,
                issues = NULL))
  }
}


#' Method to split an adjacency matrix into specific and shared
#' 
#' @title splits an adjacency matrix into specific and shared 
#' 
#' @param X An adjacency matrix
#' 
#' @return A list of two adjacency matrices
#' 
#' @author Samuel Wieczorek
#' 
#' @export
#' 
splitAdjacencyMat <- function(X){
  hasShared <- length( which(rowSums(X) > 1)) > 0
  hasSpec <- length( which(rowSums(X) == 1)) > 0
  
  
  if (hasShared && !hasSpec){
    tmpShared <- X
    tmpSpec <- X
    tmpSpec[which(rowSums(tmpSpec) > 1),] <- 0
  }
  else if (!hasShared && hasSpec){
    tmpSpec <- X
    tmpShared <- X
    tmpShared[which(rowSums(tmpShared) == 1),] <- 0
  }
  else if (hasShared && hasSpec){
    tmpSpec <- X
    tmpShared <- X
    tmpShared[which(rowSums(tmpShared) == 1),] <- 0
    tmpSpec[which(rowSums(tmpSpec) > 1),] <- 0
  } else {
    tmpSpec <- X
    tmpShared <- X
  }

  return (list(Xshared = tmpShared, Xspec = tmpSpec))
}



#' @title xxxx 
#' 
#' @param pepData A data.frame of quantitative data
#' 
#' @param X An adjacency matrix
#' 
#' @return A matrix
#' 
#' @author Samuel Wieczorek
#' 

inner.sum <- function(pepData, X){
  X <- as.matrix(X)
  pepData[is.na(pepData)] <- 0
  Mp <- t(as.matrix(X)) %*% pepData
  return(Mp)
}


#' @title xxxx 
#' 
#' @param pepData A data.frame of quantitative data
#' 
#' @param X An adjacency matrix
#' 
#' @return xxxxx
#' 
#' @author Samuel Wieczorek
#' 
inner.mean <- function(pepData, X){
  X <- as.matrix(X)
  Mp <- inner.sum(pepData, X)
  Mp <- Mp / GetNbPeptidesUsed(X, pepData)
  
  return(Mp)
  
}




#' @title xxxx 
#' 
#' @param pepData A data.frame of quantitative data
#' 
#' @param X An adjacency matrix
#' 
#' @param method xxxxx
#' 
#' @param n xxxxx
#' 
#' @return xxxxx
#' 
#' @author Samuel Wieczorek
#' 
inner.aggregate.topn <-function(pepData, X, method='Mean', n=10){
  
  X <- as.matrix(X)
  med <- apply(pepData, 1, median)
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
         Mean= Mp <- inner.mean(pepData, X),
         Sum= Mp <- inner.sum(pepData, X)
  )
  
  return(Mp)
}

#' This function computes the intensity of proteins as the sum of the 
#' intensities of their n best peptides.
#' 
#' @title Compute the intensity of proteins as the sum of the 
#' intensities of their n best peptides.
#' 
#' @param obj.pep A matrix of intensities of peptides
#' 
#' @param X An adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' 
#' @param method xxx
#' 
#' @param n The maximum number of peptides used to aggregate a protein.
#' 
#' @return A matrix of intensities of proteins
#' 
#' @author Alexia Dorffer, Samuel Wieczorek
#' 
#' @examples
#' \dontrun{
#' utils::data(Exp1_R25_pept, package='DAPARdata') 
#' obj.pep <- Exp1_R25_pept[1:10]
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' ll.agg <- DAPAR::aggregateTopn(obj.pep, X, n=3)
#' }
#' 
#' @export
#' 
#' @importFrom Biobase pData exprs fData
#' 
aggregateTopn <- function(obj.pep,X,  method='Mean', n=10){
  obj.prot <- NULL
  
  # Agregation of metacell data
  metacell <- AggregateMetacell(X = X, obj.pep = obj.pep)
  
  if (!is.null(metacell$issues))
    return(list(obj.prot = NULL,
                issues = metacell$issues))
  else {
    # Step 2 : Agregation of quantitative data
    pepData <- 2^(Biobase::exprs(obj.pep))
    protData <- inner.aggregate.topn(pepData, X, method=method, n)
    
    # Step 3: Build the protein dataset
    obj.prot <- finalizeAggregation(obj.pep, pepData, protData, metacell$metacell, X)
    
    return(list(obj.prot = obj.prot,
                issues = NULL))
  }
}




#' Method to finalize the aggregation process
#' 
#' @title Finalizes the aggregation process 
#' 
#' @param obj.pep A peptide object of class \code{MSnset}
#' 
#' @param pepData xxxx
#' 
#' @param X An adjacency matrix in which lines and columns correspond 
#' respectively to peptides and proteins.
#' 
#' @param protData xxxxx
#' 
#' @param protMetacell xxx
#' 
#' @return A protein object of class \code{MSnset}
#' 
#' @author Samuel Wieczorek
#' 
#' @export
#' 
#' @importFrom utils installed.packages
#' 
#' 
finalizeAggregation <- function(obj.pep, pepData, protData, protMetacell, X){
  
  if(missing(obj.pep))
    stop("'obj.pep' is missing")
  if(missing(pepData))
    stop("'pepData' is missing")
  if(missing(protData))
    stop("'protData' is missing")
  if(missing(protMetacell))
    stop("'protMetacell' is missing")
  if(missing(X))
    stop("'X' is missing")
 
  #(obj.pep, pepData, protData, metacell, X)
  
  protData <- as.matrix(protData)
  X <- as.matrix(X)
  protData[protData==0] <- NA
  protData[is.nan(protData)] <- NA
  protData[is.infinite(protData)] <-NA
  
  temp <- GetDetailedNbPeptidesUsed(X, pepData)
  
  pepSharedUsed <- as.matrix(temp$nShared)
  colnames(pepSharedUsed) <- paste("pepShared.used.", colnames(pepData), sep="")
  rownames(pepSharedUsed) <- colnames(X)
  
  pepSpecUsed <- as.matrix(temp$nSpec)
  colnames(pepSpecUsed) <- paste("pepSpec.used.", colnames(pepData), sep="")
  rownames(pepSpecUsed) <- colnames(X)
  
  pepTotalUsed <- as.matrix(GetNbPeptidesUsed(X, pepData))
  colnames(pepTotalUsed) <- paste("pepTotal.used.", colnames(pepData), sep="")
  rownames(pepTotalUsed) <- colnames(X)
  
  n <- GetDetailedNbPeptides(X)
  
  
  fd <- data.frame(nPepTotal = n$nTotal,
                   nPepShared = n$nShared, 
                   nPepSpec = n$nSpec, 
                   pepSpecUsed, 
                   pepSharedUsed, 
                   pepTotalUsed, 
                   protMetacell)
  rownames(fd) <- colnames(X)
  obj.prot <- MSnSet(exprs = log2(protData), 
                     fData = fd, 
                     pData = Biobase::pData(obj.pep))
  
  obj.prot@experimentData@other <- obj.pep@experimentData@other
  obj.prot@experimentData@other$typeOfData <- "protein"
  
  
  obj.prot@experimentData@other$Prostar_Version <- NA
  obj.prot@experimentData@other$DAPAR_Version <- NA
  tryCatch({
    find.package("Prostar")
    find.package("DAPAR")
    
    obj.prot@experimentData@other$Prostar_Version <- Biobase::package.version('Prostar')
    obj.prot@experimentData@other$DAPAR_Version <- Biobase::package.version('DAPAR')
  },
  error = function(e) {
    obj.prot@experimentData@other$Prostar_Version <- NA
    obj.prot@experimentData@other$DAPAR_Version <- NA
  }
  )
  
  return (obj.prot)
}




#' @title
#' Combine peptide metadata to build protein metadata
#' 
#' @description 
#' Agregation rules for the cells metadata of peptides. 
#' Please refer to the metacell vocabulary in `metacell.def()`
#' 
#' # Basic agreagtion
#' Agregation of non imputed values (2.X) with quantitative values 
#' (1.0, 1.X, 3.0, 3.X)
#' |----------------------------
#' Not possible
#' |----------------------------
#' 
#' Agregation of different types of missing values (among 2.1, 2.2)
#' |----------------------------
#' * Agregation of 2.1 peptides between each other gives a missing value 
#'   non imputed (2.0)
#' * Agreagtion of 2.2 peptides between each other givesa missing value 
#'   non imputed (2.0)
#' * Agregation of a mix of 2.1 and 2.2 gives a missing value non imputed (2.0)
#' |----------------------------
#' 
#' 
#' Agregation of a mix of quantitative values (among 1.0, 1.1, 1.2, 3.0, 3.X)
#' |----------------------------
#' * if the type of all the peptides to agregate is 1.0, 1.1 or 1.2, 
#'   then the final metadata is set the this tag
#' * if the set of metacell to agregate is a mix of 1.0, 1.1 or 1.2, 
#'   then the final metadata is set to 1.0
#' * if the set of metacell to agregate is a mix of 3.X and 3.0, 
#'   then the final metadata is set to 3.0
#' * if the set of metacell to agregate is a mix of 3.X and 3.0 and other (X.X),
#'   then the final metadata is set to 4.0
#' |----------------------------
#' 
#' # Post processing
#' Update metacell with POV/MEC status for the categories 2.0 and 3.0
#' TODO
#' 
#' @param met xxx
#' 
#' @param level xxx
#' 
#' @examples
#' \dontrun{
#' ll <- metacell.def('peptide')$node
#' for (i in 1:length(ll))
#' test <- lapply(combn(ll, i, simplify = FALSE), 
#' function(x) tag <- metacombine(x, 'peptide'))
#' }
#' 
#' 
metacombine <- function(met, level) {
  #browser()
  tag <- NULL
  if (length(met)==0)
    return('missing')
  
  u_met <- unique(met)
  
  ComputeNbTags <- function(tag)
    sum(unlist(lapply( search.metacell.tags(tag, level), function(x) length(grep(x, u_met)))))
  
  
  nb.tags <- lapply(metacell.def(level)$node, function(x) as.numeric(x %in% u_met))
  n.imputed <- ComputeNbTags('imputed')
  n.missing <- ComputeNbTags('missing')
  n.quanti <- ComputeNbTags('quanti')
 
  
  if(n.missing > 0 && (n.imputed > 0 || n.quanti > 0)) tag <- 'STOP'
   # stop("You try to combine missing values (2.X) with quantitative values (1.X or 3.X).")
  
  # sw : Agregation of a mix of 2.X gives a missing value non imputed (2.0)
  if (n.missing > 0 && n.quanti == 0 && n.imputed == 0) tag <- 'missing'
  
  
  # # Agregation of a mix of 2.1 and 2.2 gives a missing value non imputed (2.0)
  # if (length(u_met)== length(grep('missing_', u_met))) tag <- 'missing'
  # 
  # # Agreagtion of 2.2 peptides between each other givesa missing value non imputed (2.0)
  # if (length(u_met)==1 && 'missing_MEC' == u_met) tag <- 'missing'
  # 
  # # Agreagtion of 2.2 peptides between each other gives a missing value non imputed (2.0)
  # if (length(u_met)==1 && 'missing_POV' == u_met) tag <- 'missing'
  #     
  # # Agregation of 2.1 peptides between each other gives a missing value non imputed (2.0)
  # if (length(u_met)==1 && 'missing_MEC' == u_met) tag <- 'missing'
  
  # if the type of all the peptides to agregate is 1.0, 1.1 or 1.2, then the final
  # metadata is set the this tag
  if (length(u_met)==1 && u_met == 'quanti') tag <- 'quanti'
  if (length(u_met)==1 && u_met == 'identified') tag <- 'identified'
  if (length(u_met)==1 && u_met == 'recovered') tag <- 'recovered'
  
  
  # if the set of metacell to agregate is a mix of 1.0, 1.1 or 1.2, then the final
  # metadata is set to 1.0
  if (n.quanti > 1 && n.imputed == 0 && n.missing==0) tag <- 'quanti'
   

  # If the set of metacell to agregate is a mix of 3.X and 3.0, then the final
  # metadata is set to 3.0
  if (n.quanti == 0 && n.imputed > 0 && n.missing == 0) tag <- 'imputed'

  # If the set of metacell to agregate is a mix of 3.X and 3.0 and other (X.X), 
  # then the final metadata is set to 4.0
  if (n.quanti > 0 && n.imputed > 0 && n.missing == 0)
    tag <- 'combined'
  
  #print(paste0(paste0(u_met, collapse=' '), ' ---> ', tag))
  return(tag)
}


#' @title 
#' Symbolic product of matrices
#' 
#' @description 
#' Execute a product two matrices: the first is an adjacency one  while the
#' second if a simple dataframe
#' 
#' @param X An adjacency matrix between peptides and proteins
#' 
#' @param obj.pep A dataframe of the cell metadata for peptides
#' 
#' @return xxxx
#' 
#' @author Samuel Wieczorek
#' 
#' @examples 
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj.pep <- Exp1_R25_pept[1:100]
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' agg.meta <- AggregateMetacell(X, obj.pep)
#' 
#' @export
#' 
AggregateMetacell <- function(X, obj.pep){
  
  issues <- NULL
  meta = GetMetacell(obj.pep)
  level = obj.pep@experimentData@other$typeOfData
  rowcol <- function(meta.col, X.col) (meta.col)[X.col > 0]
  
  df <- data.frame(stringsAsFactors = TRUE)
  for (j in 1:ncol(meta))
    for(i in 1:ncol(X)){
      df[i, j] <- metacombine( rowcol(meta[,j], X[,i]), level)
    }

  df[df=='NA'] <- NA
  colnames(df) <- obj.pep@experimentData@other$names_metacell
  rownames(df) <- colnames(X)
  
  # Delete protein with only NA
  
  
  
  # Post processing of metacell to discover 'imputed POV', 'imputed MEC'
  conds <- Biobase::pData(obj.pep)$Condition
  df <- Set_POV_MEC_tags(conds, df, level)
  
  
  # Search for issues
  prot.ind <- unique(rownames(which(df == 'STOP', arr.ind = TRUE)))
  if (!is.null(prot.ind))
    issues <- setNames(
      lapply(prot.ind, 
             function(x) 
               rownames(X)[which(X[, which(colnames(X)==x)]==1)]
             ),
      prot.ind
      )
  
  list(metacell = df,
       issues = issues
  )
  
  
} 
