data(Exp1_R25_pept, package="DAPARdata")
obj <- Exp1_R25_pept[seq_len(100)]
level <- 'peptide'

#' 
#' Delete lines which are entirely filled with any missing values ('Missing MEC' and 'Missing POV')
metacell.mask <- match.metacell(GetMetacell(obj), "Missing", level)
indices <- GetIndices_WholeLine(metacell.mask)
obj.filter <- MetaCellFiltering(obj, indices, "delete")


obj <- obj[1:10]

pattern <- "Quantified"
type <- "AtLeastOneCond"
percent <- FALSE
op <- ">="
th <- 3
indices <- GetIndices_MetacellFiltering(obj, level, pattern, type, percent, op, th)
obj <- MetaCellFiltering(obj, indices, "keep")$new
#fData(obj)[, obj@experimentData@other$names_metacell]

pattern <- "Quant. by direct id"
type <- "AtLeastOneCond"
percent <- FALSE
op <- ">="
th <- 3
indices <- GetIndices_MetacellFiltering(obj, level, pattern, type, percent, op, th)
obj <- MetaCellFiltering(obj, indices, "keep")$new
#fData(obj)[, obj@experimentData@other$names_metacell]
names.1 <- rownames(obj)


obj <- Exp1_R25_pept[seq_len(100)]
pattern <- "Quant. by direct id"
type <- "AtLeastOneCond"
percent <- FALSE
op <- ">="
th <- 3
indices <- GetIndices_MetacellFiltering(obj, level, pattern, type, percent, op, th)
obj <- MetaCellFiltering(obj, indices, "keep")$new
#fData(obj)[, obj@experimentData@other$names_metacell]

pattern <- "Quantified"
type <- "AtLeastOneCond"
percent <- FALSE
op <- ">="
th <- 3
indices <- GetIndices_MetacellFiltering(obj, level, pattern, type, percent, op, th)
obj <- MetaCellFiltering(obj, indices, "keep")$new
#fData(obj)[, obj@experimentData@other$names_metacell]
names.2 <- rownames(obj)