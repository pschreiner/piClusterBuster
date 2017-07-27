#!/usr/bin/Rscript
#
# Returns whether or not lib1 piRNA clusters are in
# the lib2 piRNA clusters
#
# Example Run:
# Rscript proTRAC-ClusterComparison.R ~/lib1-proTRAC_dir/results.table lib2-proTRAC_dir/results.table

ags <- commandArgs(TRUE)
in1<-ags[1]; in2<-ags[2]
options(warn=-1)

suppressPackageStartupMessages(library(GenomicRanges))

## Post-proTRAC run: proTRAC2bed
in1 <- readLines(in1)
in2 <- readLines(in2)

# Sort by normalized piRNA counts
# Normalization: reads/genomic hits/(total mapped reads *10^6)
clusts1 <- in1[grepl("Cluster", in1)]
clusts1 <- clusts1[grepl("Coordinates", clusts1)]
bed1 <- GRanges()
for(i in 1:length(clusts1)) {
	# Retrieve scaffold or chromosome
	loc <- strsplit(clusts1, "Location: ")
	loc2 <- strsplit(loc[[i]][2], "\t")
	loc2 <- strsplit(as.character(loc2[[1]][1]), " ",)
	loc <- loc2[[1]][1]

	# Retrieve coordinates on scaffold or chromosome
	coord <- strsplit(clusts1, "Coordinates: ")
	coord2 <- strsplit(coord[[i]][2], "\t")
	coord2 <- coord2[[1]][1]
	coord <- strsplit(coord2, "-")
	beg <- coord[[1]][1]
	fin <- coord[[1]][2]

	curr <- GRanges(loc, IRanges(as.integer(beg), as.integer(fin)))
	bed1 <- c(bed1, curr)
}

clusts2 <- in2[grepl("Cluster", in2)]
clusts2 <- clusts2[grepl("Coordinates", clusts2)]
bed2 <- GRanges()
for(i in 1:length(clusts2)) {
        # Retrieve scaffold or chromosome
        loc <- strsplit(clusts2, "Location: ")
        loc2 <- strsplit(loc[[i]][2], "\t")
        loc2 <- strsplit(as.character(loc2[[1]][1]), " ",)
        loc <- loc2[[1]][1]

        # Retrieve coordinates on scaffold or chromosome
        coord <- strsplit(clusts2, "Coordinates: ")
        coord2 <- strsplit(coord[[i]][2], "\t")
        coord2 <- coord2[[1]][1]
        coord <- strsplit(coord2, "-")
        beg <- coord[[1]][1]
        fin <- coord[[1]][2]

        curr <- GRanges(loc, IRanges(as.integer(beg), as.integer(fin)))
        bed2 <- c(bed2, curr)
}
## If you just want to write to BED format
#write.table(bed, paste(results_dir, outfile, sep=""), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

ol <- NULL
for(i in 1:length(start(bed1))) {
	if(bed1[i] %over% bed2) { ol <- append(ol, 1) }
	else { ol <- append(ol, 0) }
}

tot_ol <- sum(ol)
perc_ol <- round((tot_ol / length(start(bed1)) ) * 100, digits=1)

cat("\n\n")
cat("piRNA Cluster Overlap: ", tot_ol, "/", length(start(bed1)), "\n")
cat("Percent of Overlapping piRNA Clusters: ", perc_ol, "%\n\n")
