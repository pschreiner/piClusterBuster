#!/usr/bin/Rscript
ags <- commandArgs(TRUE)
i<-ags[1];
ref_gen<-ags[2];
f_res<-ags[3];
gn_set<-ags[4];
te_set<-ags[5];
ncbi_nt<-ags[6];
base_dir<-ags[7];
results_dir<-ags[8];
qsub<-ags[9]; 
srun<-ags[10]
threads<-ags[11]; 
picall<-ags[12]; 
verbose<-ags[13]; 
gid<-ags[14]; 
fa<-ags[15]; 
go<-ags[16];

source(paste(base_dir, "bin/piClusterBuster_source.R", sep=""))

setwd(as.character(results_dir))

base <- paste("picluster", i, sep="")
dir.create(base)
setwd(base);

# Establish Filenames
gn_rm <- paste(results_dir, base, "/", base, "_Genelandscape-RM.tsv", sep="")
gn_filt1 <-  paste(results_dir, base, "/", base, "_Genelandscape-FILTERED_RM.tsv", sep="")
gn_cens <- paste(results_dir, base, "/", base, "_Genelandscape-CENSOR.tsv", sep="")
gn_filt2 <-  paste(results_dir, base, "/", base, "_Genelandscape-FILTERED_CENSOR.tsv", sep="")
gn_all <-  paste(results_dir, base, "/", base, "_Genelandscape-ALL.tsv", sep="")

te_rm <- paste(results_dir, base, "/", base, "_TElandscape-RM.tsv", sep="")
te_filt1 <-  paste(results_dir, base, "/", base, "_TElandscape-FILTERED_RM.tsv", sep="")
te_cens <- paste(results_dir, base, "/", base, "_TElandscape-CENSOR.tsv", sep="")
te_filt2 <-  paste(results_dir, base, "/", base, "_TElandscape-FILTERED_CENSOR.tsv", sep="")
te_all <-  paste(results_dir, base, "/", base, "_TElandscape-ALL.tsv", sep="")

oth_blast <- paste(results_dir, base, "/", base, "_Otherlandscape-BLAST.tsv", sep="")
oth_filt <-  paste(results_dir, base, "/", base, "_Otherlandscape-FILTERED_BLAST.tsv", sep="")
oth_all <- paste(results_dir, base, "/", base, "_Otherlandscape-ALL.tsv", sep="")

gnte1 <- paste(results_dir, base, "/", base, "-GnTE1-RM.tsv", sep="")
thru_te2 <- paste(results_dir, base, "/", base, "-RM_TEcens.tsv", sep="")
gnte_all <- paste(results_dir, base, "/", base, "-GnTE.tsv", sep="")

all_anno <- paste(results_dir, base, "/", base, "-AllAnnotation.tsv", sep="")
all_anno_fin <- paste(results_dir, base, "/", base, "-FinalAnnotation.tsv", sep="")

piseq <- paste(results_dir, base, "/", base, ".fa", sep="")
rm_out <- paste(piseq, ".out",sep="")
unann_seqs <- paste(results_dir, base, "/", base, "-UnannotatedSeqs.fa", sep="")

suppressPackageStartupMessages(library(Biostrings))
ref_gen <- readDNAStringSet(ref_gen)

pi_loci <- read.table(f_res, sep="\t")
pi_loci <- pi_loci[i,]
piclust_size <- pi_loci[1,3] - pi_loci[1,2] + 1

if(as.character(fa) != "FALSE") {
	library(seqinr)
	fst <- read.fasta(as.character(fa), as.string=TRUE, forceDNAtolower=FALSE)
	fst <- fst[as.integer(i)]
	write.fasta(as.list(as.character(fst)), as.character(attr(fst, "name")), piseq)  ## THIS LINE must be failing
} else { loci2seq(ref_gen, pi_loci, piseq) }

# TE & Gene Summary #
run_rm(piseq, te_set, te_rm, te_filt1, base_dir, results_dir, i, "TE", qsub, srun, threads)
while( !file.exists(te_filt1) ) { Sys.sleep(60) }
if(!isTRUE(verbose)) { unlink("*.o"); unlink("*.e"); unlink("*.id"); unlink("*.found"); unlink("*.ori.out"); unlink("*.tbl"); unlink("*.cat"); unlink("*.masked"); }

if(length(readLines(te_rm))>1) { 
	unann_loci <- getUnannotated(te_filt1, "BLAST", pi_loci)
} else { ## RM returned no hits
	library(GenomicRanges)
	unann_loci <- GRanges(as.character(pi_loci[1,1]), IRanges(as.integer(as.character(pi_loci[1,2])), as.integer(as.character(pi_loci[1,3]))))
}

if(class(unann_loci) == "GRanges") { # Sequence is available to be annotated
	loci2seq(ref_gen, unann_loci, unann_seqs)

	sw <- searchWorthy(unann_seqs, 11) # only perform search if greater than 11 nucleotides of sequence (minimum seed) is available
	if(isTRUE(sw)) {
		run_rm(unann_seqs, gn_set, gn_rm, gn_filt1, base_dir, results_dir, i, "Gene", qsub, srun, threads)
		
		while( !file.exists(gn_filt1) ) { Sys.sleep(60) }
		if(!isTRUE(verbose)) { unlink("*.o"); unlink("*.e"); unlink("*.id"); unlink("*.found"); unlink("*.ori.out"); unlink("*.tbl"); unlink("*.cat"); unlink("*.masked"); }
	
		system(paste("cat *-FILTERED_RM.tsv >", gnte1, sep=" ")) # merge RM annotations 
	} else { file.copy(te_filt1, gnte1) }
} else { file.copy(te_filt1, gnte1) }

unann_loci <- getUnannotated(gnte1, "BLAST", pi_loci)
if(class(unann_loci) == "GRanges") { # Sequence is available to be annotated
	loci2seq(ref_gen, unann_loci, unann_seqs)

	sw <- searchWorthy(unann_seqs, 11) # only perform search if greater than 11 nucleotides of sequence (minimum seed) is available
	if(isTRUE(sw)) {
		run_censor(unann_seqs, te_set, te_cens, te_filt2, base_dir, results_dir, i, "TE", "blastn", qsub, srun)

		while( !file.exists(te_filt2) ) { Sys.sleep(60) }

		system(paste("cat", gnte1, te_filt2, ">", thru_te2, sep=" "))

		if(!isTRUE(verbose)) { unlink(paste(unann_seqs, ".*", sep="")) }
	} else { file.copy(gnte1, thru_te2) }
} else { file.copy(gnte1, thru_te2) }

unann_loci <- getUnannotated(thru_te2, "BLAST", pi_loci)
if(class(unann_loci) == "GRanges") { # Sequence is available to be annotated
	loci2seq(ref_gen, unann_loci, unann_seqs)

	sw <- searchWorthy(unann_seqs, 11) # only perform search if greater than 11 nucleotides of sequence (minimum seed) is available
	if(isTRUE(sw)) {
		run_censor(unann_seqs, gn_set, gn_cens, gn_filt2, base_dir, results_dir, i, "Gene", "blastn", qsub, srun)

		while( !file.exists(gn_filt2) ) { Sys.sleep(60) }
 
		system(paste("cat", thru_te2, gn_filt2, ">", gnte_all, sep=" ")) 	# merge all annotations

		if(!isTRUE(verbose)) { unlink(paste(unann_seqs, ".*", sep="")) }
	} else { file.copy(thru_te2, gnte_all) }
} else { file.copy(thru_te2, gnte_all) }

# Other Summary #
	# Annotate Undefined Loci using BLAST against the NCBI database
unann_loci <- getUnannotated(gnte_all, "BLAST", pi_loci)
if(class(unann_loci) == "GRanges" && as.character(ncbi_nt) != "FALSE") { # Sequence is available to be annotated
	loci2seq(ref_gen, unann_loci, unann_seqs)

	sw <- searchWorthy(unann_seqs, 11) # only perform search if greater than 11 nucleotides of sequence (minimum seed) is available
	if(isTRUE(sw)) {
		run_blast(unann_seqs, ncbi_nt, oth_blast, oth_filt, base_dir, results_dir, i, "Other", "blastn", qsub, srun, threads)
		
		while( !file.exists(oth_filt) ) { Sys.sleep(60) }
		
		if(!isTRUE(verbose)) { unlink(oth_blast) }
	
		system(paste("cat", gnte_all, oth_filt, ">", all_anno, sep=" ")) 
	} else { file.copy(gnte_all,all_anno) }
} else { file.copy(gnte_all,all_anno) }

# Return final loci2seq to get unannotated loci after NCBI BLAST search
unann_loci <- getUnannotated(all_anno, "BLAST", pi_loci)
if(class(unann_loci) == "GRanges") { loci2seq(ref_gen, unann_loci, unann_seqs) }

ord <- read.table(all_anno)
temp<-NULL
for(ii in 1:length(ord[,1])) {
	if(as.character(ord[ii,3]) == "TE") { temp <- append(temp,1) 
	} else if(as.character(ord[ii,3]) == "Gene") {
		temp <- append(temp,2) 
	} else { temp <- append(temp,3) } }
ord <- cbind(ord, temp)
ord <- ord[order(ord$temp),]
write.table(ord[,1:12], all_anno, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

# Filter hits across features
system(paste("Rscript ", base_dir, "bin/filterFeatures.R", " ", all_anno, " ", all_anno_fin, " ", base_dir, sep=""))

if(!isTRUE(verbose)) { unlink("*.n*"); unlink("*.e*"); unlink("*.o*"); unlink("*.log"); unlink(all_anno) }
if(file.size(all_anno_fin)!=0) { all_anno_df <- read.table(all_anno_fin)
} else { all_anno_df <- data.frame("","","") }

te <- all_anno_df[all_anno_df[,3] == "TE",]
write.table(te, te_all, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
gn <- all_anno_df[all_anno_df[,3] == "Gene",]
write.table(gn, gn_all, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
oth <- all_anno_df[all_anno_df[,3] == "Other",]
write.table(oth, oth_all, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

# Gather values necessary for picluster summaries
if(file.info(te_all)$size != 0) { te_summ <- feature_summary(te_all, piclust_size) } else { te_summ <- list(0,0,0,0,0,0,0,0,0) }
if(file.info(gn_all)$size != 0) { gn_summ <- feature_summary(gn_all, piclust_size)  } else { gn_summ <- list(0,0,0,0,0,0,0,0,0) }
if(file.info(oth_all)$size != 0) { oth_summ <- feature_summary(oth_all, piclust_size) } else { oth_summ <- list(0,0,0,0,0,0,0,0,0) }

# Summarize Results #
occ_vals <- create_piClust_Summary(piclust_size, gn_summ, te_summ, oth_summ, base, all_anno_df, pi_loci, te, gn, gid, go)

# Print to temporary files values that need to be passed to genome summary
res <- append(occ_vals, base)
res <- append(res, piclust_size)
if(picall == "TRUE") { res <- append(res, round(pi_loci[1,5], digits=0)) }  # Get Normalized Read Count associated with piRNA cluster call

setwd("..")
temp_out <- paste(base, "-TEMPORARY.xls", sep="")
write.table(res, temp_out, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
