args = commandArgs(trailingOnly=TRUE)
infile <- args[1]

t <- read.table(infile,header=TRUE,row.names=1)
f <- as.numeric(as.character(t$classicFisher))
f[is.na(f)]<-0
p <- p.adjust(f,method='fdr')
t$classicFisher<-p

write.table(file=infile,t,sep="\t")
