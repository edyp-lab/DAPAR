ll <- metacell.def("peptide")$node
for (i in seq_len(length(ll))) {
  test <- lapply(
    combn(ll, i, simplify = FALSE),
    function(x) tag <- metacombine(x, "peptide")
  )
}

metacombine(c('Quant. by direct id', 'Missing POV'), 'peptide')
