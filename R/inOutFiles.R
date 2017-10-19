##' Builds an object of class \code{MSnSet} from a 
##' single tabulated-like file for quantitative and meta-data and a dataframe 
##' for the samples description. It differs from
##' the original \code{MSnSet} builder which requires three separated files 
##' tabulated-like quantitative proteomic data into a \code{MSnSet} object,
##' including metadata.
##' 
##' @title Creates an object of class \code{MSnSet} from text file
##' @param file The name of a tab-separated file that contains the data.
##' @param metadata A dataframe describing the samples (in lines).
##' @param indExpData A vector of string where each element is the name
##' of a column in designTable that have to be integrated in
##' the \code{fData()} table
##' of the \code{MSnSet} object.
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
##' @return An instance of class \code{MSnSet}.
##' @author Florence Combes, Samuel Wieczorek
##' @examples 
##' require(DAPARdata)
##' require(Matrix)
##' exprsFile <- system.file("extdata", "Exp1_R25_pept.txt", package="DAPARdata")
##' metadataFile <- system.file("extdata", "samples_Exp1_R25.txt", package="DAPARdata")
##' metadata = read.table(metadataFile, header=TRUE, sep="\t", as.is=TRUE)
##' indExpData <- c(56:61)
##' indFData <- c(1:55,62:71)
##' indiceID <- 64
##' createMSnset(exprsFile, metadata,indExpData,  indFData, indiceID, pep_prot_data = "peptide")
createMSnset <- function(file,metadata=NULL,indExpData,indFData,indiceID=NULL, 
                        logData=FALSE, replaceZeros=FALSE,
                        pep_prot_data=NULL){

if (!is.data.frame(file)){ #the variable is a path to a text file
data <- read.table(file, header=TRUE, sep="\t",colClasses="character")
} else {data <- file}

    ## replace all blanks by a dot
    ##   cols <- gsub(" ","\\.",  colnames(data)[indExpData])
    ##   dotIndice <- regexpr(pattern = '.',cols, fixed=TRUE) [1]
    ##   pattern <- substr(cols,1,dotIndice)
    ##   cols <- sub(pattern[1], replacement="", cols)
    #intensities <- as.matrix(data[,indExpData])
    #intensities <- gsub(",", ".", intensities)
    ##building exprs Data of MSnSet file
   Intensity <- matrix(as.numeric(gsub(",", ".",as.matrix(data[,indExpData] )))
                       , ncol=length(indExpData)
                       , byrow=FALSE)
    
    #colnames(Intensity) <- colnames(data)[indExpData]
    colnames(Intensity) <- gsub(".", "_", colnames(data)[indExpData], fixed=TRUE)
    
    ##the name of lines are the same as the data of the first column
    if (is.null(indiceID)) {
    rownames(Intensity) <- rep(paste(pep_prot_data, "_", 1:nrow(Intensity), 
                                    sep=""))
    }else{rownames(Intensity) <- data[,indiceID]}

    ##building fData of MSnSet file
    fd <- data.frame( data[,indFData])
    if (is.null(indiceID)) {
        rownames(fd) <- rep(paste(pep_prot_data, "_", 1:nrow(fd), sep=""))
    }else{
        rownames(fd) <- data[,indiceID]
        }
    
    #rownames(fd) <- data[,indiceID]
    #colnames(fd) <- colnames(data)[indFData]
    colnames(fd) <- gsub(".", "_", colnames(data)[indFData], fixed=TRUE)
    

    ##building pData of MSnSet file
    if (!is.na(sum(match(metadata$Bio.Rep," ")))) 
        {metadata$Bio.Rep <-  as.factor(1:length(metadata$Bio.Rep))}
    if (!is.na(sum(match(metadata$Tech.Rep," ")))) 
        {metadata$Tech.Rep <-  as.factor(1:length(metadata$Tech.Rep))}
    if (!is.na(sum(match(metadata$Analyt.Rep," ")))) 
        {metadata$Analyt.Rep <-  as.factor(1:length(metadata$Analyt.Rep))}
    pd <- as.data.frame(metadata)
    #rownames(pd) <- pd$Experiment
    rownames(pd) <- gsub(".", "_", pd$Experiment, fixed=TRUE)
    
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
    
    if (is.null(obj@experimentData@other$isMissingValues)){
        obj@experimentData@other$isMissingValues <- 
            Matrix::Matrix(as.numeric(is.na(obj)),
                   nrow = nrow(obj), 
                   sparse=TRUE)
    }
    
    return(obj)
}


##' This function exports a \code{MSnSet} data object to a Excel file.
##' Each of the 
##' three data.frames in the \code{MSnSet} object (ie experimental data,
##' phenoData
##' and metaData are respectively integrated into separate sheets in
##' the Excel file).
##' 
##' @title This function exports a \code{MSnSet} object to a Excel file.
##' @param obj An object of class \code{MSnSet}.
##' @param filename A character string for the name of the Excel file.
##' @return A Excel file (.xlsx)
##' @author Samuel Wieczorek
##' @examples
##' \donttest{
##' Sys.setenv("R_ZIPCMD"= Sys.which("zip"))
##' require(DAPARdata)
##' data(Exp1_R2_pept)
##' obj <- Exp1_R2_pept[1:1000]
##' writeMSnsetToExcel(obj, "foo")
##' }
writeMSnsetToExcel <- function(obj, filename)
{
    #require(openxlsx)
    name <- paste(filename, ".xlsx", sep="")
    wb <- openxlsx::createWorkbook(name)
    n <- 1
    openxlsx::addWorksheet(wb, "Quantitative Data")
    openxlsx::writeData(wb, sheet=n, cbind(ID = rownames(Biobase::exprs(obj)),
                                 Biobase::exprs(obj)), rowNames = FALSE)
    #bodyStyleNumber <- createStyle(numFmt = "NUMBER")
    #addStyle(wb, sheet=1, bodyStyleNumber, rows = 2:nrow(Biobase::exprs(obj)), 
    #cols=2:ncol(Biobase::exprs(obj)),gridExpand = TRUE)
    
    openxlsx::addWorksheet(wb, "Samples Meta Data")
    n <- n +1
    openxlsx::writeData(wb, sheet=n, Biobase::pData(obj), rowNames = FALSE)
    n <- n +1
    if (dim(Biobase::fData(obj))[2] != 0){
        openxlsx::addWorksheet(wb, "Feature Meta Data")
        #numericCols <- which(sapply(Biobase::fData(obj), is.numeric))
        #Biobase::fData(obj)[,numericCols] <- 
        #format(Biobase::fData(obj)[,numericCols])
        
        openxlsx::writeData(wb, sheet=n, cbind(ID = rownames(Biobase::fData(obj)),
                                     Biobase::fData(obj)), rowNames = FALSE)
        #bodyStyleNumber <- createStyle(numFmt = "NUMBER")
        #addStyle(wb, sheet=3, bodyStyleNumber, 
        #rows = 2:nrow(Biobase::exprs(obj)), cols=numericCols, 
        #gridExpand = TRUE, stack=TRUE)
        
    }
    
    if (!is.null(obj@experimentData@other$GGO_analysis))
        {
        l <- length(obj@experimentData@other$GGO_analysis$ggo_res)
        for (i in 1:l){
            n <- n +1
            level <- as.numeric(obj@experimentData@other$GGO_analysis$levels[i])
            openxlsx::addWorksheet(wb, paste("Group GO - level ", level, sep=""))
            openxlsx::writeData(wb, sheet=n, obj@experimentData@other$GGO_analysis$ggo_res[[i]]$ggo_res@result)
            }
        }
    
    if (!is.null(obj@experimentData@other$EGO_analysis))
        {
        n <- n +1
        openxlsx::addWorksheet(wb, "Enrichment GO")
        openxlsx::writeData(wb, sheet=n, obj@experimentData@other$EGO_analysis$ego_res@result)
        
        }

    openxlsx::saveWorkbook(wb, name, overwrite=TRUE)
    return(name)
    
    
}

##' This function reads a sheet of an Excel file and put the data into a data.frame.
##' 
##' @title This function reads a sheet of an Excel file and put the data into a data.frame.
##' @param file The name of the Excel file.
##' @param extension ddddd
##' @param sheet The name of the sheet
##' @return A data.frame
##' @author Samuel Wieczorek
readExcel <- function(file, extension, sheet){
    data <- NULL
    if (extension=="xls") {
     data <- readxl::read_xls(file, sheet)
    }
    else if (extension=="xlsx") {
        data <- readxl::read_xlsx(file, sheet)
    }
    return(as.data.frame(data))

}


##' This function lists all the sheets of an Excel file.
##' 
##' @title This function returns the list of the sheets names in a Excel file.
##' @param file The name of the Excel file.
##' @return A vector
##' @author Samuel Wieczorek
listSheets <- function(file){
    #require(openxlsx)
    return(openxlsx::getSheetNames(file))
    
}