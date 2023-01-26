
testHT <- function(obj){
  level <- 'protein'
  metacell.mask <- match.metacell(GetMetacell(obj), "Missing", level)
  indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
  obj <- MetaCellFiltering(obj, indices, cmd = "delete")$new
  
   test.ttest <- compute_t_tests(obj)
  
  
  design <- 'OnevsOne'
  
  test.limma <- limmaCompleteTest(exprs(obj), pData(obj), design)
  
  cat('compute_t_tests(obj)\n')
  print(head(test.ttest$logFC))
  cat('limmaCompleteTest()\n')
  print(head(test.limma$logFC))
}



obj <- Exp1_R10_pept

testHT(obj)
