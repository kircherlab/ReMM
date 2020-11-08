data <- str(load(snakemake@input[[1]]))
out <- str(load(snakemake@output[[1]]))
#data <-'output/features/annotated/test/SNVs.test.data.txt'
data <- read.delim(data, header = F, sep = "\t", dec = ".")
#out <-'output/features/annotated/SNVs.test.data.bin'

outfile <- file( out, "wb" )
for (i in 1:ncol(data)) {
    dataColumn = data[,i]
    writeBin( dataColumn, outfile )
}
#print(data)
