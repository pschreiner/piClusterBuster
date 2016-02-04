#!/usr/bin/Rscript
ags = commandArgs(TRUE)
i=ags[1]; ref_gen=ags[2]; f_res=ags[3]; gn_set=ags[4]; te_set=ags[5]; ncbi_nt=ags[6]; base_dir=ags[7]; results_dir=ags[8]; qsub=ags[9];
source("~/scripts/piClusterBuster/piClusterBuster_source-v0.3.R")

setwd(results_dir)
dir.create(paste("picluster", i, sep=""))
setwd(paste("picluster", i, sep=""));

base = paste("picluster", i, sep="")
gn_blast = paste(base, "_Genelandscape.BLAST", sep="")
gn_filt =  paste(base, "_Genelandscape-FILTERED.BLAST", sep="")
te_blast = paste(base, "_TElandscape.BLAST", sep="")
te_filt =  paste(base, "_TElandscape-FILTERED.BLAST", sep="")
oth_blast = paste(base, "_Otherlandscape.BLAST", sep="")
oth_filt =  paste(base, "_Otherlandscape-FILTERED.BLAST", sep="")
gnte_blast = paste(base, "-GnTE.BLAST", sep="")
piseq = paste(base, ".fa", sep="")
unann_seqs = paste(base, "-UnannotatedSeqs.fa", sep="")

pi_loci = read.table(f_res)
suppressPackageStartupMessages(library(Biostrings))
ref_gen = readDNAStringSet(ref_gen)
loci2seq(ref_gen, pi_loci[i,], piseq)
piclust_size = pi_loci[i,3] - pi_loci[i,2]

run_blast(piseq, gn_set, gn_blast, base_dir, i, "Gene", qsub) # BLAST gene set
run_blast(piseq, te_set, te_blast, base_dir, i, "TE", qsub) # BLAST TE set

gn_summ = feature_summary(gn_filt, piclust_size) # Gene Summary
te_summ = feature_summary(te_filt, piclust_size) # TE Summary

if(ncbi_nt != gn_set)
{
	## Unannotated Summary
	system("cat *-FILTERED.BLAST >", gnte_blast, sep= " ") 	# merge all annotations 
	gnte_blast = read.table(gnte_blast)
	gnte_blast = gnte_blast[order(gnte_blast[,11]),]		

	# Differentiate annotated vs unannotated loci
	piclust_loci = GRanges(as.character(pi_loci[i,1]), IRanges(as.integer(as.character(pi_loci[i,2])),as.integer(as.character(pi_loci[i,3]))))

	# Convert from BLAST to Genome coordinates (based on the location of the piRNA cluster)
	col1=NULL; for(j in 1:length(gnte_blast[,1])) { col1 = append(col1, pi_loci[i,1]) }
	col2=NULL; for(j in 1:length(gnte_blast[,1])) { temp = (as.integer(as.character(pi_loci[i,2])) + as.integer(as.character(gnte_blast[j,9])) - 1); col2 = append(col2, temp)  } 
	col3=NULL; for(j in 1:length(gnte_blast[,1])) { temp = (as.integer(as.character(pi_loci[i,3])) + as.integer(as.character(gnte_blast[j,10])) - 1); col3 = append(col3, temp)  }
	annotated_loci = GRanges(as.character(col1), IRanges(as.integer(col2),as.integer(col3)))

	unann_loci = setdiff(piclust_loci, annotated_loci) # Extract regions not covered
	#heads=NULL; for(j in 1:length(unann_loci[,1])) { temp = paste(seqnames(unann_loci)[j], ":", start(unann_loci)[j], "..", end(unann_loci)[j], sep=""); heads = append(heads, temp) }
	loci2seq(ref_gen, unann_loci, unann_seqs)

	## Other Summary (NCBI results)
	run_blast(unann_seqs, ncbi_nt, oth_blast, base_dir, i, "Other", qsub) 

	oth_summ = feature_summary(oth_filt) 

} else { oth_summ = list(0,0,0,0,0,0,0,0,0) }

while( !exists("gn_summ") || !exists("te_summ") || !exists("oth_summ") ) { Sys.sleep(60) } # Program will wait until the necessary analyses have been performed

## Summarize Results
#### Create Nice Plots
occ_vals = create_piClust_Summary(piclust_size, gn_summ, te_summ, oth_summ, base)

# Print to temporary files values that need to be passed to genome summary
res = append(occ_vals, base)
res = append(res, piclust_size)

setwd("..")
temp_out = paste(base, "-TEMPORARY.xls", sep="")
write.table(res, temp_out, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
