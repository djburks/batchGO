# batchGO
Scripts for Batch GO Enrichment and p-value adjustment using MCL Cluster Files

Requires gene clusters in MCL's standard output format: a tab-delimited list of genes in a cluster per line.
```
Rscript --vanilla batchGO.R input_file.txt
```
Exports BP, CC, and MF enrichment using the classicFisher test from topGO for each cluster, as well as optional .PS enrichment graphs.


The included fdrmaker.R script will convert (in-place) output from batchGO to use FDR-adjusted p-values.
```
Rscript --vanilla fdrmaker.R batchGO_output.txt
```

Default settings of script require a user-supplied Universe.txt file containing every gene in your enrichment universe.
Mapping is handled via org.At.tair.db for Arabidopsis, though this can be changed to any organism from available as an org.x db.

If you use this, make sure to cite the appropriate packages.

### topGO
Alexa A, Rahnenfuhrer J (2019). topGO: Enrichment Analysis for Gene Ontology. R package version 2.36.0. 
### Rgraphviz
Hansen KD, Gentry J, Long L, Gentleman R, Falcon S, Hahne F, Sarkar D (2019). Rgraphviz: Provides plotting capabilities for R graph objects. R package version 2.28.0.
### At.org.tair.db
Carlson M (2019). org.At.tair.db: Genome wide annotation for Arabidopsis. R package version 3.8.2. 
