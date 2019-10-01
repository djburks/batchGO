# Make sure these are installed through biocLite prior to running this script.

library(topGO)
library(org.At.tair.db)
library(Rgraphviz)

# Make arguments possible.

args = commandArgs(trailingOnly=TRUE)

# Prompt for cluster file.
# Should be separated as lists of genes per cluster on each row (default dump of MCL)

clusterfile <- args[1]

# Read in the file using readLine

conn <- file(clusterfile,open="r")
clusters <- readLines(conn)

# Prompt for universe list file.
# Should contain of all genes contained in the cluster file.

#unifile <- readline(prompt="What is the universe file name? ")

# Scan the universe file, and make a vector.

uni <- scan("Universe.txt",what=character(),sep="\t")

# Work through the clusters one at a time, generating topGO results.
# Unlist and split the lines into a character vector.

for (i in 1:length(clusters)) {
	clusterlist <- clusters[i]
	clusterlist <- unlist(strsplit(clusterlist,split="\t"))
	
	if (length(clusterlist) > 5000) next
	if (length(clusterlist) < 25) next

# Make a named vector between the cluster list and the universe.

	geneList <- factor(as.integer(uni %in% clusterlist))
	names(geneList) <- uni

# Make the TopGOdata Object

	GOdata <- new("topGOdata",ontology="BP",allGenes=geneList,annot=annFUN.org,mapping="org.At.tair.db")
	GOdata2 <- new("topGOdata",ontology="MF",allGenes=geneList,annot=annFUN.org,mapping="org.At.tair.db")
	GOdata3 <- new("topGOdata",ontology="CC",allGenes=geneList,annot=annFUN.org,mapping="org.At.tair.db")

# Perform Fisher test on the TopGOdata Object

	resultFisher <- runTest(GOdata,algorithm="classic",statistic="fisher")
	resultFisher2 <- runTest(GOdata2,algorithm="classic",statistic="fisher")
	resultFisher3 <- runTest(GOdata3,algorithm="classic",statistic="fisher")
	#resultKS <- runTest(GOdata,algorithm="classic",statistic="ks")
	#resultKS.elim <- runTest(GOdata,algorithm="elim",statistic="ks")

# Make the GO Enrichment Table

	#allRes <- GenTable(GOdata,classicFisher=resultFisher,classicKS=resultKS,elimKS=resultKS.elim,orderBy="elimKS",ranksOf="classicFisher",topNodes=length(score(resultFisher)))
	allRes <- GenTable(GOdata,classicFisher=resultFisher,orderBy="classicFisher",ranksOf="classicFisher",topNodes=length(score(resultFisher)))
	allRes2 <- GenTable(GOdata2,classicFisher=resultFisher2,orderBy="classicFisher",ranksOf="classicFisher",topNodes=length(score(resultFisher2)))
	allRes3 <- GenTable(GOdata3,classicFisher=resultFisher3,orderBy="classicFisher",ranksOf="classicFisher",topNodes=length(score(resultFisher3)))

# Print it to file
	
	outname <- paste("cluster",i,"BP","txt",sep=".")
	outname2 <- paste("cluster",i,"MF","txt",sep=".")
	outname3 <- paste("cluster",i,"CC","txt",sep=".")
	write.table(allRes,file=outname,sep="\t",col.names=NA)
	write.table(allRes2,file=outname2,sep="\t",col.names=NA)
	write.table(allRes3,file=outname3,sep="\t",col.names=NA)

# Print graphs using the 5 most significant GO terms of the Fisher test.

	graphname <- paste("cluster",i,"BP",sep="_")
	graphname2 <- paste("cluster",i,"MF",sep="_")
	graphname3 <- paste("cluster",i,"CC",sep="_")
	#printGraph(GOdata,resultFisher,firstSigNodes=5,fn.prefix=graphname,useInfo='all',pdfSW=FALSE)
	#printGraph(GOdata2,resultFisher2,firstSigNodes=5,fn.prefix=graphname2,useInfo='all',pdfSW=FALSE)
	#printGraph(GOdata3,resultFisher3,firstSigNodes=5,fn.prefix=graphname3,useInfo='all',pdfSW=FALSE)

# Close Loop

}

# Close File

close(conn)
