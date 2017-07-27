#!/usr/bin/Rscript
ags <- commandArgs(TRUE)
fq<-ags[1]; genome<-ags[2]; num_clusters<-ags[3]; wd<-ags[4]; results_dir<-ags[5]; gid<-ags[6]; verbose<-ags[7]

source(paste(wd, "bin/piClusterBuster_source.R", sep=""))

## 1. proTRAC
setwd(paste(wd, "proTRAC/", sep=""))

fa <- resuffix2(fq, ".fa")

script <- paste(wd, "bin/fq2proTRAC.R", sep="")
system(paste("Rscript", script, fq, fa, sep=" ")) # Be sure to use full directory

map <- paste(gid,".map",sep="")

# Align reads using proTRAC
system(paste("perl sRNAmapper.pl -i", fa, "-g", genome, "-o", map, "-a best", sep=" "))

# Run proTRAC
system(paste("perl proTRAC_2.1.pl -genome", genome, "-map", map, "-pimax 34 -distr 1-90", sep=" "))

# Get output file
pro_dir <- list.files(path=".", pattern="map_")[1]
setwd(pro_dir)
file.copy("results.table", results_dir)
if(!isTRUE(verbose)) { setwd(".."); system(paste("rm -rf", pro_dir, sep=" ")) }

## 2. Post-proTRAC run: proTRAC2bed
setwd(results_dir)

r <- readLines("results.table")

# Sort by normalized piRNA counts
# Normalization: reads/genomic hits/(total mapped reads *10^6)
clusts <- r[grepl("Cluster", r)]
clusts <- clusts[grepl("Coordinates", clusts)]
spl <- strsplit(clusts, "Hits (normalized): ", fixed=TRUE)

temp <- NULL
for(i in 1:length(spl)) { temp <- append(temp, spl[[i]][2]) }

spl2 <- strsplit(temp, "\t")

temp <- NULL
for(i in 1:length(spl2)) { temp <- append(temp, spl2[[i]][1]) }

# Order top piRNA clusters
pinames <- paste("Cluster ", seq(length(temp)), sep="")
df1 <- data.frame(pinames, as.double(temp))
df1 <- df1[order(-df1[,2]),]
df1 <- df1[1:num_clusters,]

bed <- data.frame()
for(i in 1:length(df1[,1])) {
	# Isolate the cluster of interest
	temp <- clusts[grepl(paste(df1[i,1], "\t", sep=""), clusts)][1]
	
	# Retrieve scaffold or chromosome
	loc <- strsplit(temp, "Location: ")
	loc2 <- strsplit(loc[[1]][2], "\t")
	loc2 <- strsplit(as.character(loc2[[1]][1]), " ",)
	loc <- loc2[[1]][1]

	# Retrieve coordinates on scaffold or chromosome
	coord <- strsplit(temp, "Coordinates: ")
	coord2 <- strsplit(coord[[1]][2], "\t")
	coord2 <- coord2[[1]][1]
	coord <- strsplit(coord2, "-")
	beg <- coord[[1]][1]
	fin <- coord[[1]][2]

	# Retrieve Normalized Hits Value
	hits <- df1[i,2]

	curr <- data.frame(loc, beg, fin, df1[i,1], hits)
	bed <- rbind(bed, curr)
}
bed <- bed[order(-bed[,5]),]

outfile <- paste(results_dir, gid, "-proTRAC-piRNAclusters.bed", sep="")
write.table(bed, outfile, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
if(verbose == "FALSE") { 
	unlink("results.table") 
} else {
	system(paste("mv results.table ", gid, "-proTRAC-results.table", sep=""))
}
