data(Exp1_R25_pept, package="DAPARdata")
protID <- "Protein_group_IDs"
obj.pep <- Exp1_R25_pept[seq_len(10)]
matAdj <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
peptideData <- Biobase::fData(obj.pep)
protData <- aggregateSum(obj.pep, M)
columnName <- "Protein_group_IDs"
proteinNames <- rownames(Biobase::fData(protData$obj.prot))
BuildColumnToProteinDataset_par(peptideData, matAdj, columnName, proteinNames)

