##' Builds an object of class \code{\link{MSnSet}} from a 
##' single tabulated-like file for quantitative and meta-data and a dataframe 
##' for the samples description. It differs from
##' the original \code{MSnSet} builder which requires three separated files 
##' tabulated-like quantitative proteomic data into a \code{MSnSet} object,
##' including metadata.
##' 
##' @title Creates an object of class \code{\link{MSnSet}} from text file
##' @param file The name of a tab-separated file that contains the data.
##' @param metadata A dataframe describing the samples (in lines).
##' @param indExpData A vector of string where each element is the name
##' of a column in designTable that have to be integrated in
##' the \code{fData()} table
##' of the \code{\link{MSnSet}} object.
##' @param indFData The name of column in \code{file} that will be the name of
##' rows for the \code{exprs()} and \code{fData()} tables
##' @param indiceID The indice of the column containing the ID of entities 
##' (peptides or proteins)
##' @param logData A boolean value to indicate if the data have to be
##' log-transformed (Default is FALSE)
##' @param replaceZeros A boolean value to indicate if the 0 and NaN values of
##' intensity have to be replaced by NA (Default is FALSE)
##' @param pep_prot_data A string that indicates whether the dataset is about 
##' peptides or proteins.
##' @return An instance of class \code{\link{MSnSet}}.
##' @author Florence Combes, Samuel Wieczorek
##' @examples 
##' exprsFile <- system.file("extdata", "UPSpep25.txt", package="DAPAR")
##' metadataFile <- system.file("extdata", "samples.txt", package="DAPAR")
##' metadata = read.table(metadataFile, header=TRUE, sep="\t", as.is=TRUE)
##' indExpData <- c(56:61)
##' indFData <- c(1:55,62:71)
##' indiceID <- 64
##' createMSnset(exprsFile, metadata,indExpData,  indFData, indiceID,
##' pep_prot_data = "peptide")
createMSnset <- function(file,metadata=NULL,indExpData,indFData,indiceID=NULL, 
                        logData=FALSE, replaceZeros=FALSE,
                        pep_prot_data=NULL){

if (!is.data.frame(file)){ #the variable is a path to a text file
data <- read.csv(file, header=TRUE, sep="\t", as.is=TRUE)
} else {data <- file}

    ## replace all blanks by a dot
    ##   cols <- gsub(" ","\\.",  colnames(data)[indExpData])
    ##   dotIndice <- regexpr(pattern = '.',cols, fixed=TRUE) [1]
    ##   pattern <- substr(cols,1,dotIndice)
    ##   cols <- sub(pattern[1], replacement="", cols)

    ##building exprs Data of MSnSet file
    Intensity <- data.matrix(data[,indExpData])
    colnames(Intensity) <- colnames(data)[indExpData]
    ##the name of lines are the same as the data of the first column
    if (is.null(indiceID)) {
    rownames(Intensity) <- rep(paste(pep_prot_data, "_", 1:nrow(Intensity), 
                                    sep=""))
    }else{rownames(Intensity) <- data[,indiceID]}

    ##building fData of MSnSet file
    fd <- data.frame( data[,indFData])
    if (is.null(indiceID)) {
    rownames(fd) <- rep(paste(pep_prot_data, "_", 1:nrow(fd), sep=""))
    }else{rownames(fd) <- data[,indiceID]}
    
    #rownames(fd) <- data[,indiceID]
    colnames(fd) <- colnames(data)[indFData]

    ##building pData of MSnSet file
    pd <- metadata
    rownames(pd) <- pd$Experiment

    ##Integrity tests
    if(identical(rownames(Intensity), rownames(fd))==FALSE)
        stop("Problem consistency between
            row names expression data and featureData")

    if(identical(colnames(Intensity), rownames(pd))==FALSE) 
        stop("Problem consistency between column names 
            in expression data and row names in phenoData")

    obj <- MSnSet(exprs = Intensity, fData = fd, pData = pd)

    if (logData) {
    Biobase::exprs(obj) <- log2(Biobase::exprs(obj))
        obj@processingData@processing <- 
            c(obj@processingData@processing, "Log2 tranformed data")
    }
    
    if (replaceZeros) {
    Biobase::exprs(obj)[Biobase::exprs(obj) == 0] <- NA
    Biobase::exprs(obj)[is.nan(Biobase::exprs(obj))] <- NA
    Biobase::exprs(obj)[is.infinite(Biobase::exprs(obj))] <-NA
        obj@processingData@processing <- c(obj@processingData@processing, 
                                            "All zeros were replaced by NA")
    }
    
    if (!is.null(pep_prot_data)) {
        obj@experimentData@other$typeOfData <- pep_prot_data
    }
    obj@experimentData@other$contaminantsRemoved <- FALSE
    obj@experimentData@other$reverseRemoved <- FALSE
    obj@experimentData@other$normalizationFamily <- NULL
    obj@experimentData@other$normalizationMethod <- NULL
    obj@experimentData@other$mvFilter.method <- NULL
    obj@experimentData@other$mvFilter.threshold <-NULL
    obj@experimentData@other$imputation.method <-NULL
    
    return(obj)
}


##' This function exports a \code{\link{MSnSet}} data object to a Excel file.
##' Each of the 
##' three data.frames in the \code{\link{MSnSet}} object (ie experimental data,
##' phenoData
##' and metaData are respectively integrated into separate sheets in
##' the Excel file).
##' 
##' @title This function exports a \code{\link{MSnSet}} object to a Excel file.
##' @param obj An object of class \code{\link{MSnSet}}.
##' @param filename A character string for the name of the Excel file.
##' @param id An integer to select in the fdata frame which column has to be
##' used as an index.
##' @return A Excel file
##' @author Samuel Wieczorek
##' @examples
##' data(UPSpep25)
##' writeMSnsetToExcel(UPSpep25, "foo", 1)
writeMSnsetToExcel <- function(obj, filename, id)
{
    name <- paste(filename, ".xls", sep="")
    file <- loadWorkbook(name, create=TRUE)
    createSheet(file, "Quantitative Data")
    createSheet(file, "Feature Meta Data")
    createSheet(file, "Samples Meta Data")
     #   temp.data <- cbind(Biobase::fData(obj)[,id],Biobase::exprs(obj))
        
    writeWorksheet(file, data = Biobase::exprs(obj),
                        sheet="Quantitative Data", header=TRUE)
    writeWorksheet(file, data = Biobase::fData(obj),
                        sheet="Feature Meta Data", header=TRUE)
    writeWorksheet(file, data = Biobase::pData(obj),
                        sheet="Samples Meta Data", header=TRUE)
    saveWorkbook(file)
    return(name)
}


