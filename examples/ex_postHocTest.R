data(Exp1_R25_prot, package="DAPARdata")
obj <- Exp1_R25_prot[seq_len(1000)]
level <- 'protein'
metacell.mask <- match.metacell(GetMetacell(obj), "Missing", level)
indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
obj <- MetaCellFiltering(obj, indices, cmd = "delete")
anova_tests <- t(apply(Biobase::exprs(obj$new), 1, classic1wayAnova,
    conditions = as.factor(Biobase::pData(obj$new)$Condition)
))
names(anova_tests) <- rownames(Biobase::exprs(obj$new))
pht <- postHocTest(aov_fits = anova_tests)