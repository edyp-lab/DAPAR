
#' Saves the parameters of a tool in the pipeline of Prostar
#' 
#' @title Saves the parameters of a tool in the pipeline of Prostar
#' 
#' @param obj An object of class \code{MSnSet}
#' 
#' @param name.dataset The name of the dataset
#' 
#' @param name The name of the tool. Available values are: "Norm, Imputation, 
#' anaDiff, GOAnalysis,Aggregation"
#' 
#' @param l.params A list that contains the parameters
#' 
#' @return An instance of class \code{MSnSet}.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples 
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' l.params=list(method="Global quantile alignment", type="overall")
#' saveParameters(Exp1_R25_pept, "Filtered.peptide", "Imputation",l.params)
#' 
#' @export
#' 
saveParameters <- function(obj,name.dataset=NULL,name=NULL,l.params=NULL){
  if ( is.null(name) || is.null(name.dataset)) {
    warning("No operation has been applied to the dataset.")
    return(obj)
  }
  tmp <- list()
  if(is.null(l.params)){
    tmp[[name]] <- list()
  } else {
    tmp[[name]] <- l.params
  }
  
  obj@experimentData@other$Params[[name.dataset]] <- tmp
  #obj@processingData@processing <- c(obj@processingData@processing , buildLogText(name, l.params, level=obj@experimentData@other$typeOfData))
  
  return(obj)
}



#' Builds an object of class \code{MSnSet} from a 
#' single tabulated-like file for quantitative and meta-data and a dataframe 
#' for the samples description. It differs from
#' the original \code{MSnSet} builder which requires three separated files 
#' tabulated-like quantitative proteomic data into a \code{MSnSet} object,
#' including metadata.
#' 
#' @title Creates an object of class \code{MSnSet} from text file
#' 
#' @param file The name of a tab-separated file that contains the data.
#' 
#' @param metadata A dataframe describing the samples (in lines).
#' 
#' @param indExpData A vector of string where each element is the name
#' of a column in designTable that have to be integrated in
#' the \code{fData()} table of the \code{MSnSet} object.
#' 
#' @param colnameForID The name of the column containing the ID of entities 
#' (peptides or proteins)
#' 
#' @param indexForMetacell xxxxxxxxxxx
#' 
#' @param logData A boolean value to indicate if the data have to be
#' log-transformed (Default is FALSE)
#' 
#' @param replaceZeros A boolean value to indicate if the 0 and NaN values of
#' intensity have to be replaced by NA (Default is FALSE)
#' 
#' @param pep_prot_data A string that indicates whether the dataset is about
#' 
#' @param proteinId xxxx
#' 
#' @param software xxx
#' 
#' @return An instance of class \code{MSnSet}.
#' 
#' @author Florence Combes, Samuel Wieczorek
#' 
#' @examples 
#' require(Matrix)
#' exprsFile <- system.file("extdata", "Exp1_R25_pept.txt", package="DAPARdata")
#' metadataFile <- system.file("extdata", "samples_Exp1_R25.txt", 
#' package="DAPARdata")
#' metadata = read.table(metadataFile, header=TRUE, sep="\t", as.is=TRUE)
#' indExpData <- c(56:61)
#' colnameForID <- 'id'
#' obj <- createMSnset(exprsFile, metadata,indExpData,  colnameForID, 
#' indexForMetacell = c(43:48), pep_prot_data = "peptide", software = 'maxquant')
#' 
#' 
#' exprsFile <- system.file("extdata", "Exp1_R25_pept.txt", package="DAPARdata")
#' metadataFile <- system.file("extdata", "samples_Exp1_R25.txt", package="DAPARdata")
#' metadata = read.table(metadataFile, header=TRUE, sep="\t", as.is=TRUE)
#' indExpData <- c(56:61)
#' colnameForID <- 'AutoID'
#' obj <- createMSnset(exprsFile, metadata, indExpData,  colnameForID, 
#' indexForMetacell = c(43:48), pep_prot_data = "peptide", software = 'maxquant')
#' 
#' 
#' @export
#' 
#' @importFrom MSnbase MSnSet
#' @importFrom utils read.table
#' 
createMSnset <- function(file,
                         metadata = NULL,
                         indExpData,
                         colnameForID = NULL,
                         indexForMetacell = NULL,
                         logData = FALSE, 
                         replaceZeros = FALSE,
                         pep_prot_data = NULL,
                         proteinId = NULL,
                         software = NULL){
  
  if (!is.data.frame(file)){ #the variable is a path to a text file
    data <- read.table(file, header=TRUE, sep="\t",stringsAsFactors = FALSE)
  } else {
    data <- file
    }
  
  colnames(data) <- gsub(".", "_", colnames(data), fixed=TRUE)
  colnameForID <- gsub(".", "_", colnameForID, fixed=TRUE)
  proteinId <- gsub(".", "_", proteinId, fixed=TRUE)
  colnames(data) <- gsub(" ", "_", colnames(data), fixed=TRUE)
  colnameForID <-  gsub(" ", "_", colnameForID, fixed=TRUE)
  proteinId <-  gsub(" ", "_", proteinId, fixed=TRUE)
   
  ##building exprs Data of MSnSet file
  Intensity <- matrix(as.numeric(gsub(",", ".",as.matrix(data[,indExpData] )))
                      , ncol=length(indExpData)
                      , byrow=FALSE)
  
  colnames(Intensity) <- gsub(".", "_", colnames(data)[indExpData], fixed=TRUE)
  rownames(Intensity) <- rownames(data)
 
  
  # Get teh metacell info
  metacell <- NULL
  if (!is.null(indexForMetacell)){
    metacell <- data[, indexForMetacell]
    metacell <- apply(metacell,2,tolower)
    metacell <- as.data.frame(apply(metacell,2, function(x) gsub(" ", '', x)),
                              stringsAsFactors = FALSE)
  }

  
  ##building fData of MSnSet file
  if(is.null(colnameForID))
    colnameForID <- 'AutoID'
  
  if (colnameForID == 'AutoID') {
    fd <- data.frame( data, 
                      AutoID = rep(paste(pep_prot_data, "_", 1:nrow(data), sep="")) , 
                      stringsAsFactors = FALSE)
    rownames(fd) <- paste(pep_prot_data, "_", 1:nrow(fd), sep="")
    rownames(Intensity) <- paste(pep_prot_data, "_",  1:nrow(Intensity), sep="")
  }else{
    fd <- data
    rownames(fd) <- data[ ,colnameForID]
    rownames(Intensity) <- data[ ,colnameForID]
  }
  
  colnames(fd) <- gsub(".", "_", colnames(fd), fixed=TRUE)
  
  pd <- as.data.frame(metadata, stringsAsFactors = FALSE)
  rownames(pd) <- gsub(".", "_", pd$Sample.name, fixed=TRUE)
  pd$Sample.name <- gsub(".", "_", pd$Sample.name, fixed=TRUE)
  
  ##Integrity tests
  if(identical(rownames(Intensity), rownames(fd))==FALSE)
    stop("Problem consistency between
             row names expression data and featureData")
  
  if(identical(colnames(Intensity), rownames(pd))==FALSE) 
    stop("Problem consistency between column names 
             in expression data and row names in phenoData")
  
  obj <- MSnSet(exprs = Intensity, fData = fd, pData = pd)
  
  
  
  if (replaceZeros) {
    exprs(obj)[exprs(obj) == 0] <- NA
    exprs(obj)[is.nan(exprs(obj))] <- NA
    exprs(obj)[is.infinite(exprs(obj))] <-NA
    obj@processingData@processing <- c(obj@processingData@processing, "All zeros were replaced by NA")
  }
  if (logData) {
    exprs(obj) <- log2(exprs(obj))
    obj@processingData@processing <- 
      c(obj@processingData@processing, "Data has been Log2 tranformed")
  }
  
  if (!is.null(pep_prot_data)) {
    obj@experimentData@other$typeOfData <- pep_prot_data
  }
  
  obj@experimentData@other$Prostar_Version <- NA
  tryCatch({
    find.package("Prostar")
    obj@experimentData@other$Prostar_Version <- package.version('Prostar')
  },
  error = function(e) obj.prot@experimentData@other$Prostar_Version <- NA
  )
  
  obj@experimentData@other$DAPAR_Version <- NA
  tryCatch({
    find.package("DAPAR")
    obj@experimentData@other$Prostar_Version <- package.version('DAPAR')
  },
  error = function(e) obj@experimentData@other$DAPAR_Version <- NA
  )
  
  obj@experimentData@other$proteinId <- proteinId
  obj@experimentData@other$keyId <- colnameForID
  
    obj@experimentData@other$RawPValues <- FALSE
  
 
  metacell <- BuildMetaCell(from = software,
                            level = pep_prot_data,
                            qdata = exprs(obj), 
                            conds = pData(obj)$Condition, 
                            df = metacell)
  
  fData(obj) <- cbind(fData(obj), 
                               metacell, 
                               deparse.level = 0)
  obj@experimentData@other$names_metacell <- colnames(metacell)
  
  
  return(obj)
}




#' @title This function exports a data.frame to a Excel file.
#' 
#' @param df An data.frame
#' 
#' @param tags xxx
#' 
#' @param colors xxx
#' 
#' @param tabname xxx
#' 
#' @param filename A character string for the name of the Excel file.
#' 
#' @return A Excel file (.xlsx)
#' 
#' @author Samuel Wieczorek
#' 
#' @export
#' 
#' @import openxlsx
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' df <- exprs(Exp1_R25_pept[1:100])
#' tags <- GetMetacell(Exp1_R25_pept[1:100])
#' colors <- list('missing POV' = "lightblue",
#'                'missing MEC' = "orange",
#'                'recovered' = "lightgrey",
#'                'identified' = "white",
#'                'combined' = "red")
#' write.excel(df, tags, colors, filename = 'toto')
write.excel <- function(df,
                        tags = NULL,
                        colors = NULL,
                        tabname = 'foo',
                        filename = NULL){
  
  if (is.null(filename))
    filename <- paste('data-', Sys.Date(), '.xlxs', sep='')
  else if(tools::file_ext(filename) != ""){
    if (tools::file_ext(filename) != "xlsx")
      stop("Filename extension must be equal to 'xlsx'. Abort...")
    else
      fname <- filename
  } else
    fname <- paste(filename, ".xlsx", sep="")
  
  unique.tags <- NULL
  if (!is.null(tags) && !is.null(colors)){
    unique.tags <- unique(as.vector(as.matrix(tags)))
    if (!isTRUE(sum(unique.tags %in% names(colors)) == length(unique.tags)))
      warning("The length of colors vector must be equal to the number of different tags. 
              As is it not the case, colors are ignored")
  }
  
  wb <- openxlsx::createWorkbook(fname)
  openxlsx::addWorksheet(wb, tabname)
  openxlsx::writeData(wb, sheet = 1, df, rowNames = FALSE)
  
  
  # Add colors w.r.t. tags
  if (!is.null(tags) && !is.null(colors))
    if (isTRUE(sum(unique.tags %in% names(colors)) == length(unique.tags))){
      lapply(1:length(colors), function(x){
        list.tags <- which(names(colors)[x]==tags, arr.ind=TRUE)
        openxlsx::addStyle(wb,
                           sheet = 1,
                           cols = list.tags[,"col"],
                           rows = list.tags[,"row"] + 1, 
                           style = openxlsx::createStyle(fgFill = colors[x])
        )
      })
    }
  
  openxlsx::saveWorkbook(wb, fname, overwrite=TRUE)
}


#' This function exports a \code{MSnSet} data object to a Excel file.
#' Each of the three data.frames in the \code{MSnSet} object (ie experimental 
#' data, phenoData and metaData are respectively integrated into separate sheets 
#' in the Excel file).
#' The colored cells in the experimental data correspond to the original 
#' missing values which have been imputed.
#' 
#' @title This function exports a \code{MSnSet} object to a Excel file.
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param filename A character string for the name of the Excel file.
#' 
#' @return A Excel file (.xlsx)
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' \donttest{
#' Sys.setenv("R_ZIPCMD"= Sys.which("zip"))
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R25_pept[1:10]
#' writeMSnsetToExcel(obj, "foo")
#' }
#' 
#' @export
#' 
#' @import openxlsx
#' 
writeMSnsetToExcel <- function(obj, filename)
{
  name <- paste(filename, ".xlsx", sep="")
  wb <- openxlsx::createWorkbook(name)
  n <- 1
  openxlsx::addWorksheet(wb, "Quantitative Data")
  openxlsx::writeData(wb, sheet=n, cbind(ID = rownames(exprs(obj)),
                                         exprs(obj)), rowNames = FALSE)
  
  
  # Add colors to quantitative table
  mc <- metacell.def(GetTypeofData(obj))
  colors <- as.list(setNames(mc$color, mc$node))
  tags <- cbind(keyId = rep('identified', nrow(obj)),
                GetMetacell(obj)
  )
  
  
  unique.tags <- NULL
  if (!is.null(tags) && !is.null(colors)){
    unique.tags <- unique(as.vector(as.matrix(tags)))
    if (!isTRUE(sum(unique.tags %in% names(colors)) == length(unique.tags)))
      warning("The length of colors vector must be equal to the number of different tags. 
              As is it not the case, colors are ignored")
    if (isTRUE(sum(unique.tags %in% names(colors)) == length(unique.tags))){
      lapply(1:length(colors), function(x){
        list.tags <- which(names(colors)[x]==tags, arr.ind=TRUE)
        openxlsx::addStyle(wb,
                           sheet = 1,
                           cols = list.tags[ ,"col"],
                           rows = list.tags[ ,"row"] + 1, 
                           style = openxlsx::createStyle(fgFill = colors[x])
        )
      })
    }
  }

  
  n <- 2
  openxlsx::addWorksheet(wb, "Samples Meta Data")
  openxlsx::writeData(wb, sheet = n, pData(obj), rowNames = FALSE)
  
  
  # Add colors for sample data sheet
  u_conds <- unique(pData(obj)$Condition)
  colors <- setNames(DAPAR::ExtendPalette(length(u_conds)),
                     u_conds)
  colors[['blank']] <- 'white'
  
  tags <- pData(obj)
  tags[,] <- 'blank'
  tags$Sample.name <- pData(obj)$Condition
  tags$Condition <- pData(obj)$Condition
  
  unique.tags <- NULL
  if (!is.null(tags) && !is.null(colors)){
    unique.tags <- unique(as.vector(as.matrix(tags)))
    if (!isTRUE(sum(unique.tags %in% names(colors)) == length(unique.tags)))
      warning("The length of colors vector must be equal to the number of different tags. 
              As is it not the case, colors are ignored")
    if (isTRUE(sum(unique.tags %in% names(colors)) == length(unique.tags))){
      lapply(1:length(colors), function(x){
        list.tags <- which(names(colors)[x]==tags, arr.ind=TRUE)
        openxlsx::addStyle(wb,
                           sheet = n,
                           cols = list.tags[ ,"col"],
                           rows = list.tags[ ,"row"] + 1, 
                           style = openxlsx::createStyle(fgFill = colors[x])
        )
      })
    }
  }

  
  ## Add feature Data sheet
   
  n <- 3
  if (dim(fData(obj))[2] != 0){
    openxlsx::addWorksheet(wb, "Feature Meta Data")
    openxlsx::writeData(wb, 
                        sheet = n, 
                        cbind(ID = rownames(fData(obj)),
                                           fData(obj)), rowNames = FALSE)
  }
  
  colors <- as.list(setNames(mc$color, mc$node))
  tags <- cbind(keyId = rep('identified', nrow(obj)),
                fData(obj)
                )
  
  tags[,] <- 'identified'
  tags[, 1 + which(colnames(fData(obj)) %in% obj@experimentData@other$names_metacell)] <- GetMetacell(obj)
  
  unique.tags <- NULL
  if (!is.null(tags) && !is.null(colors)){
    unique.tags <- unique(as.vector(as.matrix(tags)))
    if (!isTRUE(sum(unique.tags %in% names(colors)) == length(unique.tags)))
      warning("The length of colors vector must be equal to the number of different tags. 
              As is it not the case, colors are ignored")
    if (isTRUE(sum(unique.tags %in% names(colors)) == length(unique.tags))){
      lapply(1:length(colors), function(x){
        list.tags <- which(names(colors)[x]==tags, arr.ind=TRUE)
        openxlsx::addStyle(wb,
                           sheet = n,
                           cols = list.tags[ ,"col"],
                           rows = list.tags[ ,"row"] + 1, 
                           style = openxlsx::createStyle(fgFill = colors[x])
        )
      })
    }
  }
  
  # Add GO tab
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




#' @title This function reads a sheet of an Excel file and put the data 
#' into a data.frame.
#' 
#' @param file The name of the Excel file.
#' 
#' @param extension The extension of the file
#' 
#' @param sheet The name of the sheet
#' 
#' @return A data.frame
#' 
#' @author Samuel Wieczorek
#' 
#' @export
#' 
#' @importFrom readxl read_excel
#' 
readExcel <- function(file, extension, sheet){
  # data <- NULL
  # if (extension=="xls") {
  #     data <- readxl::read_xls(file, sheet)
  # }
  # else if (extension=="xlsx") {
  #     data <- readxl::read_xlsx(file, sheet)
  # }
  # return(as.data.frame(data,asIs=T))
  
  #options(digits=10)
  data <- NULL
  data <- readxl::read_excel(file, sheet)
  
  return(as.data.frame(data,asIs=T, stringsAsFactors=F))
  
}


#' This function lists all the sheets of an Excel file.
#' 
#' @title This function returns the list of the sheets names in a Excel file.
#' 
#' @param file The name of the Excel file.
#' 
#' @return A vector
#' 
#' @author Samuel Wieczorek
#' 
#' @export
#' 
#' @importFrom openxlsx getSheetNames
#' 
listSheets <- function(file){
  #require(openxlsx)
  return(openxlsx::getSheetNames(file))
  
}


#' @title Exports a MSnset dataset into a zip archive containing three 
#' zipped CSV files.
#' 
#' @param obj An object of class \code{MSnSet}.
#' 
#' @param fname The name of the archive file.
#' 
#' @return A compressed file
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' \donttest{
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' obj <- Exp1_R2_pept[1:1000]
#' writeMSnsetToCSV(obj, "foo")
#' }
#' 
#' @export
#' 
#' @importFrom utils write.csv zip
#'
writeMSnsetToCSV <- function(obj, fname){
  
  #fname <- paste(tempdir(),fname,  sep="/")
  write.csv(exprs(obj), paste(tempdir(), "exprs.csv", sep='/'))
  write.csv(fData(obj), paste(tempdir(), "fData.csv", sep='/'))
  write.csv(pData(obj), paste(tempdir(), "pData.csv", sep='/'))
  files <- c(paste(tempdir(), "exprs.csv", sep='/'),
             paste(tempdir(), "fData.csv", sep='/'),
             paste(tempdir(), "pData.csv", sep='/'))
  zip(fname, files, zip = Sys.getenv("R_ZIPCMD", "zip"))
  
  return(fname)
}



#' @title Similar to the function \code{rbind} but applies on two subsets of 
#' the same \code{MSnSet} object.
#' 
#' @param df1 An object (or subset of) of class \code{MSnSet}. May be NULL
#' 
#' @param df2 A subset of the same object as df1
#' 
#' @return An instance of class \code{MSnSet}.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' utils::data(Exp1_R25_pept, package='DAPARdata')
#' df1 <- Exp1_R25_pept[1:100]
#' df2 <- Exp1_R25_pept[200:250]
#' rbindMSnset(df1, df2)
#' 
#' @export
#' 
#' @importFrom MSnbase MSnSet
#' 
rbindMSnset <- function(df1=NULL, df2){
  
  if (is.null(df1)){
    obj <- df2
    return(obj)
  }
  if (is.null(df1) && is.null(df2)){return(NULL)}
  
  tmp.exprs <- rbind(exprs(df1), exprs(df2))
  tmp.fData <- rbind(fData(df1), fData(df2))
  tmp.pData <- pData(df1)
  
  obj <-  MSnSet(exprs = tmp.exprs, fData = tmp.fData, pData = tmp.pData)
  obj@protocolData <- df1@protocolData
  obj@experimentData <- df1@experimentData
  
  return(obj)
  
}
