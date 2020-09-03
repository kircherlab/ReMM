df <- load(snakemake@input[[1]])

new_m <- cbind(m.subset,folds.subset)
new_m <- cbind(new_m,labels.subset)

write.table(new_m,snakemake@output[[1]],row.names = FALSE, col.names = FALSE, ,sep = "\t")
