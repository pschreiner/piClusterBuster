processArgs <-
function(args)
{
# Convert the bash arguments provided by the user
# into usable R variables

	vars<-NULL; temp<-NULL; temp.fa<-NULL; temp.fq<-NULL; temp.bed<-NULL; temp.x<-NULL; temp.n<-NULL; temp.ncbidb<-NULL; temp.gndb<-NULL; temp.tedb<-NULL; temp.p<-NULL; temp.qsub<-NULL; temp.srun<-NULL; temp.d<-NULL; temp.gid<-NULL; temp.verbose<-NULL; temp.go<-NULL; temp.allsrna<-NULL;

	for(z in 1:length(args)) 
	{ 
		if(grepl("^-", args[z], perl=TRUE)) { mod <- gsub("-","", args[z], fixed=TRUE); assign(paste("temp.", mod, sep=""), args[z+1]) }  ### FIXXX: what if --qsub last arg?

		vars$fa <- temp.fa; vars$fq <- temp.fq; vars$bed<-temp.bed; vars$x<-temp.x; vars$n<-temp.n; vars$ncbidb<-temp.ncbidb; vars$gndb<-temp.gndb; vars$tedb<-temp.tedb; vars$p<-temp.p; vars$qsub<-temp.qsub; vars$srun<-temp.srun; vars$d<-temp.d; vars$gid<-temp.gid; vars$verbose<-temp.verbose; vars$go<- temp.go; vars$allsrna<-temp.allsrna;
	}
	return(vars)
}

# Load necessary libraries with libraries needed as the argument
install_unloadedLibs <-
function(req_libs, lib_source="http://www.bioconductor.org/biocLite.R")
{
# Installs the libraries that are required for the program
# as provided by a character vector (req_libs)

	source(lib_source)
	if(length(setdiff(req_libs, rownames(installed.packages()))) > 0) { biocLite(setdiff(req_libs, rownames(installed.packages())),) }

	return("Libraries Installed")
} 

bed2bam <-
function(bed, ref_gen)
{
# Converts a BED file to a BAM file.
# Returns:
# String specifying the location of the new
# BAM file in the same directory as the initial
# BED file

	library(systemPipeR)
	moduleload("bedtools")

	bam <- resuffix(bed, ".bam")
	system(paste("bedToBam -i", bed, "-g", ref_gen, ">", bam, sep=" "))

	return(bam)
}

get_ncbiDB <-
function(repository = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz")
{
# Downloads the NCBI nt database

	library(seqinr)
	get.ncbi(repository)
	return(cat("\n\nDownloading NCBI nt database.  This may take several hours.\n\n"))
}

get_teDB <-
function(repository = "http://www.girinst.org/server/RepBase/protected/RepBase20.12.fasta.tar.gz")
{
# Downloads the RepBase database

	wget(repository)
	return(cat("\n\nDownloading RepBase TE database.\n\n"))
}

decompress_file <-
function(file1)
{
# Identifies and decompresses the input file
# which should be specified by a character string

	if(grepl(".zip",file1, fixed=TRUE)) { system(paste("unzip", file1, sep=" ")) }
	else if(grepl(".gzip",file1, fixed=TRUE) || grepl(".gz",file1, fixed=TRUE)) { system(paste("gunzip", file1, sep=" ")) }
	else if(grepl(".tar.gz",file1, fixed=TRUE) || grepl(".tgz",file1, fixed=TRUE)) { system(paste("tar -xzf", file1, sep=" ")) }
	else { return("Decompression not neccessary or unrecognized file suffix.") }
	
	return(paste(file1,"was successfully decompressed.", sep=" "))
}

filterpiRNAs <-
function(sRNA_fastq, min_size)
{
# Filter sRNAs to potential piRNAs
# by length > 23 nucleotides

	library(ShortRead); library(seqinr);

	if(grepl(".fq", sRNA_fastq)) { mod<-gsub(".fq", "", sRNA_fastq, ignore.case=TRUE) } else if(grepl(".fastq", sRNA_fastq, ignore.case=TRUE)) { mod<-gsub(".fastq", "", sRNA_fastq, fixed=TRUE) }
	mod <- strsplit(mod,"/"); mod <- mod[[1]][length(mod[[1]])]
	out <- paste(mod, "-greater23.fa", sep="");

	if(file.exists(out)) { print(paste("File:", out, "already exists. Skipping...", sep=" "))
	} else { 
		fun <- function(x) { x[width(x)>=24] }
		filterFastq(sRNA_fastq, out, filter=fun, compress=FALSE)
	}
	return(out)
}
## Example Use:
## filterpiRNAs("~/piClusterBuster", "my.fastq")

resuffix <-
function(file1, new_suff)
{
# Changes suffix from the original file
# into the same name, but with the new
# (provided) suffix

	if(grepl(".fastq", file1, fixed=TRUE)) { temp <- gsub(".fastq", "", file1, fixed=TRUE); fout <- paste(temp, new_suff, sep="") }
	else if(grepl(".fasta", file1, fixed=TRUE)) { temp <- gsub(".fasta", "", file1, fixed=TRUE); fout <- paste(temp, new_suff, sep="") }
	else if(grepl(".fq", file1, fixed=TRUE)) { temp <- gsub(".fq", "", file1, fixed=TRUE); fout <- paste(temp, new_suff, sep="") }
	else if(grepl(".fa", file1, fixed=TRUE)) { temp <- gsub(".fa", "", file1, fixed=TRUE); fout <- paste(temp, new_suff, sep="") }
	else if(grepl(".sam", file1, fixed=TRUE)) { temp <- gsub(".sam", "", file1, fixed=TRUE); fout <- paste(temp, new_suff, sep="") }
	else if(grepl(".bam", file1, fixed=TRUE)) { temp <- gsub(".bam", "", file1, fixed=TRUE); fout <- paste(temp, new_suff, sep="") }
	else if(grepl(".bed", file1, fixed=TRUE)) { temp <- gsub(".bed", "", file1, fixed=TRUE); fout <- paste(temp, new_suff, sep="") }

	else( paste("Error in replacing input file suffix: ",  file1, sep=""))
	return(fout)
}
## Example Use:
## resuffix("my.fastq", ".bam")  
## Will return "my.bam"

resuffix2 <-
function(file1, new_suff)
{
# Same as resuffix except that it will trim
# and replace any suffix after the period
	temp <- strsplit(file1, ".", fixed=TRUE)	
	base <- temp[[1]][1]

	fout <- paste(base, new_suff, sep="") 
	return(fout)
}

run_bowtie2 <-
function(reference_genome, fasta, out)
{
# Runs Bowtie2 given the reference genome,
# fasta file of reads, and an output file name.
# All of which should be provided as character strings

	library(systemPipeR)
	moduleload("bowtie2")

	index <- paste(reference_genome, ".1.bt2", sep="")
	if (!file.exists(index)) { system(paste("bowtie2-build", reference_genome, reference_genome, sep=" ")) }
	system(paste("bowtie2 -f -x ", reference_genome, " -U ", fasta, " -S ", out, sep=""))
	return(cat("\nBowtie2 Run Complete\n"))
}

sam2bam <- 
function(infile, out)
{
# Conversion of a SAM file to a BAM file

	library(systemPipeR)
	moduleload("samtools")

	system(paste("samtools view -Sb ", infile, " > ", out, sep=""))
	return("SAM to BAM Conversion Complete")
}

percentTo100 <-
function(vals_round0, vals_round1, num_vals)
{
	if(sum(vals_round0)==99) {
		spl <- strsplit(as.character(vals_round1), ".", fixed=TRUE)
		high <- 0;
		for(ii in 1:num_vals) {
			if(is.na(spl[[ii]][2])) { next; }
			else if(spl[[ii]][2] > high) { 
				high <- spl[[ii]][2]
				fx <- ii
			}
		}
		vals_round0[fx] <- vals_round0[fx] + 1
	} else if(sum(vals_round0)==101) {
		spl <- strsplit(as.character(vals_round1), ".", fixed=TRUE)
		low <- 9;
		for(ii in 1:num_vals) {
			if(is.na(spl[[ii]][2])) { next; }
			else if(spl[[ii]][2] < low) { 
				low <- spl[[ii]][2]
				fx <- ii
			}
		}
		vals_round0[fx] <- vals_round0[fx] - 1
	}
	return(vals_round0)
}

run_fseq <-
function(wd, gid, bam, bed)
{
# Run the F-seq program to make piRNA cluster calls

	moduleload("bedtools")
	system(paste("bamToBed -i ", bam, " > ", bed, sep="")) 

	# Download and configure F-seq
	setwd("F-seq")
	system("ant")
	setwd("dist~/fseq/bin")
	system("chmod 755 fseq")

	system(paste("./fseq -v -d ", wd, "data/ -l 1000 -o ", wd, gid, "_results/ -of bed -t 15 ", bed, sep="")) ###  " > out.xls", after bed and before sep

	out <- paste(gid, "-piRNAcluster_loci.xls", sep="")
	system(paste("cat *.bed >", out, sep=" ")) 
	
	files <- list.files()
	for(z in 1:length(files))
	{
		if(!grepl(".xls", files[z])) { system(paste("rm ", files[z], sep="")) } # Delete all files produced by F-seq that are no longer needed 
	}

	return(out)
}

run_proTRAC <-
function(fq, genome_fasta, out, num_clusters, wd, results_dir)
{
# Run proTRAC

	script <- paste(wd, "bin/run_proTRAC.R", sep="")
	name <- paste(gid, "-", "proTRAC", sep="")

	system(paste("Rscript", script, fq, genome_fasta, out, num_clusters, wd, results_dir, verbose, sep=" "))

	cat("\nproTRAC Run Completed", "\n", sep=" ")
	return(out)
}

loci2seq <- 
function(ref_gen, loci, out)
{
# Given a reference genome, loci coordinates, and
# an output file name, extracts the corresponding
# sequences.

	library(Biostrings); library(GenomicRanges); library(seqinr)
	
	seqs<-NULL; heads<-NULL;
	names(ref_gen) <- paste(names(ref_gen)," ",sep="")
	if(class(loci)[1] != "GRanges")
	{
		for(i in 1:length(loci[,1]))
		{
			seqs <- append(seqs, DNAStringSet(ref_gen[grepl(paste("^", as.character(loci[i,1])," ",sep=""), names(ref_gen), ignore.case=TRUE, perl=TRUE),][[1]], as.integer(loci[i,2]), as.integer(loci[i,3])));  
			heads <- append(heads, paste(as.character(loci[i,1]), ":", as.integer(loci[i,2]), "..", as.integer(loci[i,3]), sep=""));
        	}
	}
	else
	{
		for(i in 1:length(start(loci)))
		{
			seqs <- append(seqs, DNAStringSet(ref_gen[grepl(paste("^", as.character(seqnames(loci)[i])," ",sep=""), names(ref_gen), ignore.case=TRUE, perl=TRUE),][[1]], start(loci)[i], end(loci)[i])); 
			heads <- append(heads, paste(seqnames(loci)[i], ":", start(loci)[i], "..", end(loci)[i], sep=""));
	 	}
	}
	write.fasta(as.list(as.character(seqs)), as.list(as.character(heads)), out)
	return("Genomic Loci Converted to Sequence")
}

run_blast <-
function(db, query, out, filt_file, wd, results_dir, i, content, task, qsub, srun, threads)
{
# Run BLAST

	script <- paste(wd, "bin/run_blast.sh", sep="")
	index <- paste(db, ".nhr", sep="")
	base <- paste("picluster", i, sep="")
	name <- paste(base, "-", content, "_BLAST", sep="")
	id<-NULL; id <- paste(base, "-", content, ".id", sep="")

	if(isTRUE(qsub)) { 
		system(paste("echo \"", script, db, query, task, out, filt_file, content, results_dir, base, threads, wd, "\" | qsub -l mem=16gb -N ", name, sep=" "))
	} else if(isTRUE(srun)) {
		j1 <- paste("--job-name=", name, sep="")
		system(paste("sbatch --mem-per-cpu=16G", j1, "\"", script, db, query, task, out, filt_file, content, results_dir, base, threads, wd, "\"", sep=" "))
	} else {
		system(paste(script, db, query, task, out, filt_file, content, results_dir, base, threads, wd, sep=" "))
	}

	return(cat("\nBLAST Search Submitted -  DB:", db, "Query:", query, "\n", sep=" "))
}
run_blast2 <-
function(wd, out, filt_file, content)
{
# Run BLAST

	script <- paste(wd, "bin/run_blast2.sh", sep="")

	system(paste(script, wd, out, filt_file, content, sep=" "))

	return(cat("\nBLAST Filtering...", "\n", sep=" "))
}

run_censor <-
function(db, query, out, filt_file, wd, results_dir, i, content, task, qsub, srun)
{
# Run CENSOR

	script <- paste(wd, "bin/run_censor.sh", sep="")
	index <- paste(db, ".nhr", sep="")
	base <- paste("picluster", i, sep="")
	name <- paste(base, "-", content, "_CENSOR", sep="")

	if(isTRUE(qsub)) { 
		system(paste("echo \"", script, db, query, task, out, filt_file, content, results_dir, base, wd, "\" | qsub -N ", name, sep=" "))
	} else if(isTRUE(srun)) {
		j1 <- paste("--job-name=", name, sep="")
		system(paste("sbatch", j1, "\"", script, db, query, task, out, filt_file, content, results_dir, base, wd, "\"", sep=" "))
	} else { 
		system(paste(script, db, query, task, out, filt_file, content, results_dir, base, wd, sep=" "))
	}

	return(cat("\nCENSOR Search Submitted -  DB:", db, "Query:", query, "\n", sep=" "))
}

run_rm <-
function(db, query, out, filt_file, wd, results_dir, i, content, qsub, srun, threads)
{
# Run RepeatMasker

	script <- paste(wd, "bin/run_rm.sh", sep="")
	index <- paste(db, ".nhr", sep="")
	base <- paste("picluster", i, sep="")
	name <- paste(base, "-", content, "_RM", sep="")
	id<-NULL; id <- paste(base, "-", content, ".id", sep="")

	if(isTRUE(qsub)) { 
		system(paste("echo \"", script, db, query, results_dir, base, threads, out, filt_file, content, wd, "\" | qsub -l mem=8gb -N ", name, sep=" ")) 
	} else if(isTRUE(srun)) {
		j1 <- paste("--job-name=", name, sep="")
		system(paste("sbatch --mem-per-cpu=8G", j1, "\"", script, db, query, results_dir, base, threads, out, filt_file, content, wd, "\"", sep=" "))
	} else { 
		system(paste(script, db, query, results_dir, base, threads, out, filt_file, content, wd, sep=" "))
	}

	return(cat("\nRepeatMasker Search Submitted -  DB:", db, "Query:", query, "\n", sep=" "))
}

run_gm <-
function(i, ref_gen, wd, results_dir, qsub, srun)
{
# Run Genemark

	script <- paste(wd, "run_gm.sh", sep="")
	base <- paste("picluster", i, sep="")
	fa <- paste(base, ".fa", sep="")
	name <- paste(base, "-", "Genemark", sep="")
	gmk <- paste(wd, "gm_key", sep="")
	mod <- paste(results_dir, "GeneMark_hmm_heuristic.mod", sep="")

	if(isTRUE(qsub)) { 
		system(paste("echo \"", script, mod, fa, ref_gen, wd, results_dir, base, "\" | qsub -N ", name, sep=" ")) 
	} else if(isTRUE(srun)) {
		j1 <- paste("--job-name=", name, sep="")
		system(paste("sbatch", j1, "\"", script, mod, fa, ref_gen, wd, results_dir, base, "\"", sep=" "))
	} else {
		system(paste(script, mod, fa, ref_gen, wd, results_dir, base, sep=" "))
	}

	return(cat("\nGenemark Run Submitted\n"))
}

searchWorthy <-
function(fa_file, min_length)
{
	library(seqinr)

	fa <- read.fasta(fa_file, forceDNAtolower=FALSE, as.string=TRUE)
	for(ii in 1:length(fa))	{
		if(nchar(as.character(fa[ii]))>min_length) { return(TRUE) }
	}
	return(FALSE)
}

blast2GRanges <-
function(blast_df)
{
# Convert a BLAST output file (in tabular format)
# to a GRanges object

	library(GenomicRanges)	

	## Eliminate negative widths from BLAST hit
	 beg<-NULL; fin<-NULL; strand<-NULL;
        for(i in 1:length(blast_df[,1]))
        {
            if(as.integer(as.character(blast_df[i,9])) > as.integer(as.character(blast_df[i,10])))
            {
                    beg <- append(beg, as.integer(blast_df[i,10]))
                    fin <- append(fin, as.integer(blast_df[i,9]))
                    strand <- append(strand, "-")
            }
            else
            {
                    beg <- append(beg, as.integer(blast_df[i,9]))
                    fin <- append(fin, as.integer(blast_df[i,10]))
                    strand <- append(strand, "+")
            }
        }	

	gr <- GRanges(as.character(blast_df[,1]),IRanges(beg, fin), strand=strand, sim=as.integer(as.character(blast_df[,11])), method=as.character(blast_df[,12]), feature=blast_df[,3])
	return(gr)
}

censor2GRanges <-
function(cens_df)
{
# Convert a CENSOR output file (in tabular format)
# to a GRanges object

	library(GenomicRanges)

	col1 <- as.character(cens_df[,4])
	col2 <- as.integer(as.character(cens_df[,2]))
	col3 <- as.integer(as.character(cens_df[,3]))
	col4 <- gsub("d","+",as.character(cens_df[,7]))
	col4 <- gsub("c","-",col4)
	align_len <- as.integer(as.character(cens_df[,3])) - as.integer(as.character(cens_df[,2]))

	gr <- GRanges(col1, IRanges(col2,col3), strand=col4, sim=cens_df[,8], method="CENSOR")
	return(gr)
}

rm2GRanges <-
function(rm_df)
{
# Convert a RepeatMasker (.out) output file 
# (in tabular format) to a GRanges object

	library(GenomicRanges)	

	sq <- strsplit(as.character(rm_df[,5]), ":")
	sq <- sq[[1]][1]
	beg <- as.integer(rm_df[,6])
	fin <- as.integer(rm_df[,7])
	sim <- 100 - rm_df[,2]
	te <- rm_df[,10]
	#strand <- gsub("C", "-", rm_df[,9])  ## orientation of where in the picluster it is annotated

	pos <- rm_df[,12:14]
	strand=NULL;
	for(ii in 1:length(pos[,1])) {
		if( grepl("(", as.character(pos[ii,1]), fixed=TRUE)  ) { strand <- append(strand, "-")	}
		else { strand <- append(strand, "+") }
	}

	gr <- GRanges(as.character(sq),IRanges(beg, fin), strand=strand, sim=sim, te=te, method="RM")

	return(gr)
}

gff2GRanges <-
function(gff_df)
{
# Convert a Genemark output file 
# (in tabular format) to a GRanges object

	library(GenomicRanges)	
	gff_df <- read.table(gff_df)
	gff_df <- gff_df[gff_df[,3]=="gene"]

	sq <- strsplit(as.character(gff_df[,1]), ":")
	sq <- sq[[1]][1]
	beg <- gff_df[,4]
	fin <- gff_df[,5]
	strand <- gff_df[,7]
	score <- gff_df[,6]
	gn <- gff_df[,9]
	gn <- gsub("ID=", "", gn)
	gr <- GRanges(as.character(sq), IRanges(beg, fin), strand=strand, sim=score, gene=gn, method="Genemark")

	return(gr)
}

filterAnnotation <-
function(gr)
{
# Filter picluster annotation to remove redundant annotation.

	gr <- gr[order(-elementMetadata(gr)[[1]],-width(gr)),]
	
	new_gr<-GRanges(); been_ann<-NULL; keep<-TRUE;
	for(zz in 1:length(start(gr)))
	{		
		if(isTRUE(keep)) {
			new_gr <- append(new_gr, gr[zz,])
			been_ann <- append(been_ann,start(gr)[zz]:end(gr)[zz])
			been_ann <- sort(been_ann)
		} # Store current line since first annotation in GRanges object take priority
		keep<-TRUE;
		if(zz+1 > length(start(gr))) { break; }
	
		temp <- start(gr)[zz+1]:end(gr)[zz+1]
		logic <- !temp %in% been_ann
		if(logic[1] == "FALSE") {
			if(logic[length(logic)] == "FALSE") { keep = FALSE; next;  # Full annotation is redundant, so ignore
			} else { # Start of the next annotation was already annotated
				for(ii in 1:length(logic)) { if(isTRUE(logic[ii])) { new_val <- temp[ii]; break; } }
				start(gr)[zz+1] <-  new_val  
			}
		} else if(logic[length(logic)] == "FALSE") { # Beginning of annotation is useful
			for(ii in length(logic):1) { if(isTRUE(logic[ii])) { new_val <- temp[ii]; break; } }
			end(gr)[zz+1] <-  new_val
		}
	}

	return(new_gr)
}

filterBLASThits <-
function(xls, out, wd, base, qsub, srun, filt)
{
# Filter BLAST hits to remove redundant or insignificant annotation

	script <- paste(wd, "filterBLASThits.R", sep="")
	submit <- paste("| qsub -l mem=6gb -N ", base, "-BLAST_Filter", sep="")
	submit2 <- paste("--mem-per-cpu=6G --job-name=", base, "-BLAST_Filter", sep="")
	filt <- paste("\"", filt,"\"\"",sep="")

	if(isTRUE(qsub)) { 
		system(paste("echo \"Rscript", script, xls, out, filt, submit, sep=" "))
	} else if(isTRUE(srun)) {
		system(paste("sbatch", submit2, "--wrap=\"Rscript", script, xls, out, filt, "\"", sep=" "))
	} else { system(paste("Rscript", script, xls, out, filt, sep=" ")) }

	return(paste(out, "file created.",sep=" "));
}
## Example Use:

filterCENSORhits <-
function(xls, out, wd, base, qsub, srun, filt)
{
# Filter CENSOR hits to remove redundant or insignificant annotation

	script <- paste(wd, "filterCENSORhits.R", sep="")
	submit <- paste("| qsub -N ", base, "-CENSOR_Filter", sep="")
	submit2 <- paste("--job-name=", base, "-CENSOR_Filter", sep="")
	filt <- paste("\"", filt,"\"\"",sep="")

	if(isTRUE(qsub)) { 
		system(paste("echo \"Rscript", script, xls, out, filt, submit, sep=" "))
	} else if(isTRUE(srun)) {
		system(paste("sbatch", submit2, "--wrap=\"Rscript", script, xls, out, filt, "\"", sep=" "))
	} else { system(paste("Rscript", script, xls, out, filt, sep=" ")) }

	return(paste(out, "file created.",sep=" "));
}
## Example Use:

filterRMhits <-
function(xls, out, wd, base, qsub, srun, filt)
{
# Filter RepeatMasker hits to remove redundant or insignificant annotation

	script <- paste(wd, "filterRMhits.R", sep="")
	submit <- paste("| qsub -N ", base, "-RM_Filter", sep="")
	submit2 <- paste("--job-name=", base, "-RM_Filter", sep="")
	filt <- paste("\"", filt,"\"\"",sep="")

	if(isTRUE(qsub)) { 
		system(paste("echo \"Rscript", script, xls, out, filt, submit, sep=" "))
	} else if(isTRUE(srun)) {
		system(paste("sbatch", submit2, "--wrap=\"Rscript", script, xls, out, filt, "\"", sep=" "))
	} else { system(paste("Rscript", script, xls, out, filt, sep=" ")) }

	return(paste(out, "file created.",sep=" "));
}
## Example Use:

filterGMhits <-
function(xls, out, wd, base, qsub, srun, filt)
{
# Filter Genemark hits to remove redundant or insignificant annotation

	script <- paste(wd, "filterGMhits.R", sep="")
	submit <- paste("| qsub -N ", base, "-GM_Filter", sep="")
	submit2 <- paste("--job-name=", base, "-GM_Filter", sep="")
	filt <- paste("\"", filt,"\"\"",sep="")

	if(isTRUE(qsub)) { 
		system(paste("echo \"Rscript", script, xls, out, filt, submit, sep=" "))
	} else if(isTRUE(srun)) {
		system(paste("sbatch", submit2, "--wrap=\"Rscript", script, xls, out, filt, "\"", sep=" "))
	} else { system(paste("Rscript", script, xls, out, filt, sep=" ")) }

	return(paste(out, "file created.",sep=" "));
}
## Example Use:

feature_summary <-
function(blast_file, piclust_size)
{
# Calculates a summary of the loci of interest and returns
# the number of sense, antisense, and total feature calls


	library(GenomicRanges)

	if(!file.exists(blast_file)) { return(list(0,0,0,0,0,0,0,0,0)) }
	feat <- read.table(blast_file)
	feat <- feat[order(-feat[,11],-feat[,4]),]
	
	gr <-blast2GRanges(feat)
	#gr <- censor2GRanges(feat)
	gr <- reduce(gr)
	
	tot_feat_calls <- length(start(gr))
	tot_feat_occ <- sum(width(gr))
	tot_feat_perc <- ( tot_feat_occ / piclust_size ) * 100
	tot_feat_perc <- round(tot_feat_perc, digits = 2)

	# Check Annotations on positive strand
	pos_gr <- gr[strand(gr)=="+",]
	pos_feat_calls <- length(start(pos_gr))
	pos_feat_occ <- sum(width(pos_gr))
	pos_feat_perc <- ( pos_feat_occ / piclust_size ) * 100
	pos_feat_perc <- round(pos_feat_perc, digits = 2)

	# Check Annotations on negative strand
	neg_gr <- gr[strand(gr)=="-",]
	neg_feat_calls <- length(start(neg_gr))
	neg_feat_occ <- sum(width(neg_gr))
	neg_feat_perc <- ( neg_feat_occ / piclust_size ) * 100
	neg_feat_perc <- round(neg_feat_perc, digits = 2)

	return(list(pos_feat_calls, pos_feat_occ, pos_feat_perc, neg_feat_calls, neg_feat_occ, neg_feat_perc, tot_feat_calls, tot_feat_occ, tot_feat_perc))
}
## Example Use;
## feature_summary("piclust15_TElandscape.BLAST")

getUnannotated <-
function(all_blast, form, pi_loci, prnt=FALSE, ...)
{
# Given the loci that have been annotated and the total
# loci to be analyzed, determines if there are still
# loci that have not been annotated

	library(GenomicRanges)

	piclust_loci <- GRanges(as.character(pi_loci[1,1]), IRanges(as.integer(as.character(pi_loci[1,2])),as.integer(as.character(pi_loci[1,3]))))

	if(class(all_blast) != "data.frame") {
		if(file.info(all_blast)$size==0) { return(piclust_loci) 
		} else { all_blast <- read.table(all_blast) }	
	} 

	# Convert from BLAST to Genome coordinates (based on the location of the piRNA cluster)
	col1<-NULL; col2<-NULL; col3<-NULL; 
	for(j in 1:length(all_blast[,1])) { 
		col1 <- append(col1, as.character(pi_loci[1,1]))

		if(form == "RM") {  #### FIXXXX: Speed this up with vectorization
			temp <- (as.integer(as.character(pi_loci[1,2])) + as.integer(as.character(all_blast[j,6])) - 1); col2 <- append(col2, temp)
			temp <- (as.integer(as.character(pi_loci[1,2])) + as.integer(as.character(all_blast[j,7])) - 1); col3 <- append(col3, temp)		
		} else if(form == "CENSOR") { #### FIXXXX: Speed this up with vectorization
			temp <- (as.integer(as.character(pi_loci[1,2])) + as.integer(as.character(all_blast[j,2])) - 1); col2 <- append(col2, temp)
			temp <- (as.integer(as.character(pi_loci[1,2])) + as.integer(as.character(all_blast[j,3])) - 1); col3 <- append(col3, temp)
		} else if(form == "BLAST") {
			if( as.integer(as.character(all_blast[j,9])) < as.integer(as.character(all_blast[j,10])) )
			{
				temp <- (as.integer(as.character(pi_loci[1,2])) + as.integer(as.character(all_blast[j,9])) - 1); col2 <- append(col2, temp)
				temp <- (as.integer(as.character(pi_loci[1,2])) + as.integer(as.character(all_blast[j,10])) - 1); col3 <- append(col3, temp) 
			} else {			
				temp <- (as.integer(as.character(pi_loci[1,2])) + as.integer(as.character(all_blast[j,10])) - 1); col2 <- append(col2, temp)
				temp <- (as.integer(as.character(pi_loci[1,2])) + as.integer(as.character(all_blast[j,9])) - 1); col3 <- append(col3, temp)  		
			}
		} else { stop("Unrecognized File Format passed to getUnannotated()") }
	}
	annotated_loci <- GRanges(as.character(col1), IRanges(as.integer(col2),as.integer(col3)))
	annotated_loci <- reduce(annotated_loci)
	annotated_loci <- annotated_loci[width(annotated_loci)>7,]

	# Differentiate annotated vs unannotated loci
	unann_loci <- setdiff(piclust_loci, annotated_loci) # Extract regions not covered
	if(length(start(unann_loci)) == 0) { return(0) }
	unann_loci <- reduce(unann_loci)

	return(unann_loci)
}

summarizeTEs <-
function(ids, base)
{
# Using the a vector of strings (ids) in the format:
# TEclass:TEname:Species, determine the
# number of times a TE class was characterized and the
# number of times a specific TE was characterized in 
# the picluster.

##### FIXXXX: add nucleotide occupancy plot?

	library(plyr)

	ids <- ids[!(grepl(")n", ids, fixed=TRUE))]
	if(length(ids) == 0) { 
		plot.new()
                text(x = 0.5, y = 0.5, paste("No TE Superfamilies Identified"), cex=1.6, col="black")		

		return("TE Class Summary Complete")
	}
	ids <- strsplit(as.character(ids),":",fixed=TRUE)
	ids <- t(as.data.frame(ids))	

	tab <- table(ids[,1])
	tab <- tab[!(grepl(")n",names(tab),fixed=TRUE))]  # Rid of Simple Repeats for this summary
	tab <- tab[order(-tab)]	

	opa <- par()

	par(mar=c(10.1, 6.1, 4.1, 4.1), xpd=TRUE)
	if(is.na(tab[1])) { 
		text(x = 0.5, y = 0.5, paste("No TE Superfamilies Characterized"), cex=1.6, col="black")
		return("TE Class Summary Complete")
	} else if(length(tab)<30) { barplot(as.integer(as.character(tab)), ylab="Number of TE Calls", las=2, ylim=c(0,max(round_any(as.integer(tab), 10, f=ceiling))), names.arg=rownames(tab))
	} else { barplot(as.integer(as.character(tab[1:30])), ylab="Number of TE Calls", las=2, ylim=c(0,max(round_any(as.integer(tab), 10, f=ceiling))), names.arg=rownames(tab)[1:30]) }

	par(opa)
	
	res<-data.frame();
	for(ii in 1:length(names(tab))) {
		curr <- ids[grepl(names(tab)[ii],ids[,1]),]
		# Account for the case in which there is only 1 hit
		if(class(curr) == "character") {
			curr_tab <- table(curr[2])
		} else {
			curr_tab <- table(curr[,2])
			curr_tab <- curr_tab[order(-curr_tab)]
		}

		if(length(curr_tab)<5) { res <- rbind(res,cbind(names(tab)[ii],names(curr_tab),as.integer(curr_tab))) }
		else { res <- rbind(res,cbind(names(tab)[ii],names(curr_tab)[1:5],as.integer(curr_tab)[1:5])) }
	}

	uniq <- res[!duplicated(as.character(res[,1])),1]
	uniq <- uniq[1:5]

	# If common Class I element, then cool color, and if common Class II element, then hot color, else gray
	c1_cols <- c("NonLTR","Gypsy","Jockey","LINE","SINE", "Copia", "BEL", "DIRS", "Endogenous Retrovirus", "ERV1", "ERV2", "ERV3", "Lentivirus", "ERV4", "Non-LTR Retrotransposon", "SINE1/7SL", "SINE2/tRNA", "SINE3/5S", "SINE4", "CRE", "NeSL", "R4", "R2", "L1", "RTE", "I", "CR1", "Rex1", "RandI", "Penelope", "Tx1", "RTEX", "Crack", "Nimb", "Proto1", "Proto2", "RTETP", "Hero", "L2", "Tad1", "Loa", "Ingi", "Outcast", "R1", "Daphne", "L2A", "L2B", "Ambal", "Vingi", "Kiri") #col=rainbow(20,start=0.25,end=0.5)
	c2_cols <- c("DNA", "Harbinger", "P", "Mariner/Tc1", "hAT", "MuDR", "EnSpm/CACTA", "piggyBac", "Merlin", "Transib", "Novosib", "Helitron", "Polinton", "Kolobok", "ISL2EU", "Crypton", "CryptonA", "CryptonF", "CryptonI", "CryptonS", "CryptonV", "Sola", "Sola1", "Sola2", "Sola3", "Zator", "Ginger1", "Ginger2/TDD", "Academ", "Zisupton", "IS3EU", "Dada") # col=heat.colors(10, alpha=0.65)
	# Other: col=gray.colors(5)
	
	c1s <- uniq[as.character(uniq) %in% c1_cols]
	cool <- c("darkblue", "steelblue", "purple", "violet", "springgreen")
	cool <- cool[1:length(c1s)]
	#cool <- rainbow(length(c1s), start=0.25, end=0.5)
	c2s <- uniq[as.character(uniq) %in% c2_cols]
	hot <- c("darkred", "red", "orange", "yellow", "lightgoldenrod")
	hot <- hot[1:length(c2s)]
	#hot <- heat.colors(length(c2s), alpha=0.65)
	grys <- uniq[!(as.character(uniq) %in% c1_cols) & !(as.character(uniq) %in% c2_cols)]
	neu <- gray.colors(length(grys))
	
	# Store TE class name, number of hits, and color to be printed in sorted data frame	
	if(length(c1s)>0) { 
		c1_df <- cbind(as.character(c1s),cool); colnames(c1_df) <- c("TEclass", "cid") 
	} else { c1_df <- data.frame() }
	if(length(c2s)>0) {
		c2_df <- cbind(as.character(c2s),hot); colnames(c2_df) <- c("TEclass", "cid")
	} else { c2_df <- data.frame() }
	if(length(grys)>0) {
		un_df <- cbind(as.character(grys),neu); colnames(un_df) <- c("TEclass", "cid")
	} else { un_df <- data.frame() }
	col_df <- rbind(c1_df, c2_df, un_df)

	## FIXXXXX: Only Top 5 TE Classes (and in order) -- check counts
	for(ii in 1:length(col_df[,1])) {
		curr <- res[grepl(as.character(col_df[ii,1]),as.character(res[,1])),]
		curr <- curr[!is.na(curr[,1]),]
		if(length(curr[,1]) == 0) { break; }
		
		report <- as.integer(as.character(curr[,3]));
		names(report) <- as.character(curr[,2])
		report <- report[order(-report)]

		# Center the bar for each TE Class in plot
		if(length(report) >= 5) {
			report <- report[1:5]	
		}
		else if(length(report) == 3) {
			report <- append(0, report); report <- append(report, 0);
		}
		else if(length(report) == 2) {
			report <- append(0, report); report <- append(report, 0); report <- append(report, 0);
		}
		else if(length(report) == 1) { report <- append(0, report); report <- append(0, report); report <- append(report, 0); report <- append(report, 0); }

		if(ii==1) { 
			par(mfrow=c(1,length(col_df[,1])), mar=c(12.1, 2.1, 8.1, 4), mgp=c(5,0,0), xpd=TRUE)
			barplot(as.integer(as.character(report)), beside=TRUE, col=col_df[ii,2], las=2, xlim=c(0,1), ylim=c(0,max(as.integer(as.character(res[,3])))), width=0.3, space=0.25, border=NA, axes=FALSE)
			mtext(side = 1, text = names(report), at=c(0.22, 0.58, 0.97, 1.35, 1.72), line = 0.5, las=2, cex=0.68)
			title(xlab=names(tab)[ii], line=8, cex.lab=1.25, font.lab=2, adj=1)
			#title(ylab="Number of TE Calls", cex.lab=1.5, line=0.7)
	axis(2, at=seq(0, max(as.integer(as.character(res[,3]))), by=1), labels=paste(seq(0, max(as.integer(as.character(res[,3]))), by=1), "   ", sep=""), line=0.25, las=1)

		} else {
			if(ii==4) { title(main="TE Calls", cex.main=1.75) }
			par(mar=c(12.1, 0.5, 8.1, 0.5), mgp=c(5, 0, 0), xpd=TRUE)
			barplot(as.integer(as.character(report)), beside=TRUE, col=col_df[ii,2], las=2, xlim=c(0,1), ylim=c(0,max(as.integer(as.character(res[,3])))), width=0.15, space=0.2, border=NA, axes=FALSE)
			mtext(side = 1, text = names(report), at=c(0.10, 0.28, 0.47, 0.655, 0.84), line = 0.5, las=2, cex=0.68)
			title(xlab=names(tab)[ii], line=8, cex.lab=1.25, font.lab=2)
		}
	}
	#axis(4, at=seq(0, max(as.integer(as.character(res[,3]))), by=1), labels=paste("   ", seq(0, max(as.integer(as.character(res[,3]))), by=1), sep=""), line=-1, las=1)

	par(opa)

	return("TE Class Summary Complete") #FIXXXX: pass along results?
}

summarizeGenes <-
function(gid, gn_ids)
{
	library(gProfileR)

	opa <- par()
	
	gp <- gprofiler(gn_ids, organism=gid, significant=F, src_filter="GO:BP")
	gp <- gp[order(gp$p.value, -gp$overlap.size, -gp$relative.depth),]
	gp <- gp[!duplicated(gp$intersection),] # remove several terms returned by the same genes
	gp <- gp[gp$term.name != "biological_process",] # remove "biological_process" term
	gp <- gp[order(-gp$overlap.size),]

	if(length(gp[,1])==0) { 
		plot.new()
                text(x = 0.5, y = 0.5, paste("No GO terms associated with Gene Hits"), cex=1.6, col="black")
		return("Gene Summary Complete")
	}

	par(mar=c(25.1, 4.1, 4.1, 4.1), xpd=TRUE)
	if(length(gp[,1]) < 20) { 
		barplot(gp$overlap.size, names.arg=as.character(gp$term.name), ylab="GO Term Frequency", las=2, cex.axis=0.75, cex.names=0.9, xlim=c(0,10), ylim=c(0,1))
	} else {
		barplot(gp$overlap.size[1:20], names.arg=as.character(gp$term.name)[1:20], ylab="GO Term Frequency", las=2, cex.axis=0.75, cex.names=0.9, xlim=c(0,10), ylim=c(0,1)) 
	}
	
	par(opa)

	#fa <- read.fasta(gn_set)
	gp <- gprofiler(gn_ids, organism=gid, significant=T, src_filter="GO:BP")
	if (length(gp[, 1]) > 1) {
		gp <- gp[order(gp$p.value, -gp$overlap.size, -gp$relative.depth),]
		gp <- gp[!duplicated(gp$intersection),]
		gp <- gp[gp$term.name != "biological_process",] # remove "biological_process" term

		par(mar=c(4.1, 25.1, 4.1, 4.1), xpd=TRUE)
		barplot(abs(log(gp$p.value)), main="Overrepresented GO terms", horiz=TRUE, names.arg=as.character(gp$term.name), las=2, axes=F, xlim=c(0,10), xlab="P Value")  # FIXXXX: this plot
		axis(1, at=c(0, abs(log(0.5)), abs(log(0.05)), abs(log(0.005)), 10), labels=c("", 0.5, 0.05, 0.005, ""))
	} else { # print text to plot or otherwise indicate no enriched genes
		plot.new()
		text(x = 0.5, y = 0.5, paste("No overrespresented GO terms in picluster"), cex=1.6, col="black")  
	}

	return("Gene Summary Complete")
}

squarrows <-
function(x1, y, x2, orient, color, arrowhead=10, isgr=FALSE)
{
	for(ii in 1:length(x1))
	{
		if( (x2[ii]-x1[ii]) < arrowhead) {  #Just print triangle
			if(orient[ii] == "+") { polygon(c(x1[ii],x2[ii],x1[ii]), y=c(-y,0,y), col=color, border=NA) }
			else if(orient[ii] == "-") { polygon(c(x2[ii],x1[ii],x2[ii]), y=c(-y,0,y), col=color, border=NA) }
			else { rect(x1, -y, x2, y, col=color, border=NA)  } # No stranded information provided. Rectangles used.
		} else {
			#print triangle in proper orientation
			if(orient[ii] == "+") { 
				polygon(c(x2[ii]-arrowhead, x2[ii], x2[ii]-arrowhead), y=c(-y,0,y), col=color, border=NA)
				rect(x1[ii], -y, x2[ii]-arrowhead, y, col=color, border=NA)  #fill rest with rectangle
			} else if(orient[ii] == "-") {
				polygon(c(x1[ii]+arrowhead, x1[ii], x1[ii]+arrowhead), y=c(-y,0,y), col=color, border=NA)
				rect(x1[ii]+arrowhead, -y, x2[ii], y, col=color, border=NA)  #fill rest with rectangle
			} else { rect(x1, -y, x2, y, col=color)  } # No stranded information provided. Rectangles used.	
		}
	}
}

plotStrandedCoverage <-
function(positive, negative, beg, fin, name, subt, all_anno_df, gr, xlab="Position", ylab="Coverage")
{
	library(plyr)

	xlim <- c(as.integer(as.character(beg)), as.integer(as.character(fin)))
	len <- fin - beg
	if(max(negative)==0) {
		y_low <- -1
	} else { 
		y_low <- round_any(-max(negative), 10, f=floor) 
	}
	if(max(positive)==0) { 
		y_high <- 1 
	} else {
		y_high <- round_any(max(positive), 10, f=ceiling)
	}
	ylim <- c(y_low, y_high)

	plot.new()
	plot.window(xlim=xlim, ylim=ylim)

	# Separate Annotations
	gn <- all_anno_df[all_anno_df[,3]=="Gene",]
	te <- all_anno_df[all_anno_df[,3]=="TE",]
	oth <- all_anno_df[all_anno_df[,3]=="Other",]

	if(length(gn[,1]) > 0) { gn_gr <- blast2GRanges(gn) } else { gn_gr <- GRanges("NA", IRanges(0,0)) }
	if(length(te[,1]) > 0) { te_gr <- blast2GRanges(te) } else { te_gr <- GRanges("NA", IRanges(0,0)) }
	if(length(oth[,1]) > 0) { oth_gr <- blast2GRanges(oth) } else { oth_gr <- GRanges("NA", IRanges(0,0)) }

	title(main=name)
	mtext(subt)	
	
	x_beg <- round_any(xlim[1], 10, f=floor)
	x_fin <- round_any(xlim[2], 10, f=ceiling)
	fb <- seq(xlim[1], xlim[2], (x_fin-x_beg)/5)
	mids <- round_any(fb[2:5], 100, f=ceiling)
	corr <- c(fb[1], mids, fin)
	
	if(!is.null(gr)) {
		title(main=name, xlab="Position",ylab="Coverage")
		axis(side=1,at=corr)
		axis(side=2,at=seq(ylim[1],ylim[2],as.integer(max(abs(ylim))/5)))

		lines(c(start(negative),length(negative)),
			- c(runValue(negative),tail(runValue(negative),1)),
			type="s",col="azure4")
		lines(c(start(positive),length(positive)),
			+ c(runValue(positive),tail(runValue(positive),1)),
			type="s",col="black")
	}

	# Annotate X axis
	abline(h=0,col="black")

	rec_height <- max(ylim) / 25; # worked ok at 150
	ahead <- len/100;
	orient <- NULL;

	if(length(gn[,1]) != 0) { 
		#rect(gn[,9]+beg, -rec_height, gn[,10]+beg, rec_height, col="blue") 
		squarrows(start(gn_gr)+beg, rec_height, end(gn_gr)+beg, as.character(strand(gn_gr)), col="blue", arrowhead=ahead)
	} # Representing Gene Calls
	if(length(te[,1]) != 0) { 
		#rect(te[,9]+beg, -rec_height, te[,10]+beg, rec_height, col="red") 
		squarrows(start(te_gr)+beg, rec_height, end(te_gr)+beg, as.character(strand(te_gr)), col="red", arrowhead=ahead)
	} # Representing TE Calls
	if(length(oth[,1]) != 0) { 
		#rect(oth[,9]+beg, -rec_height, oth[,10]+beg, rec_height, col="yellow") 
		squarrows(start(oth_gr)+beg, rec_height, end(oth_gr)+beg, as.character(strand(oth_gr)), col="yellow", arrowhead=ahead)
	} # Representing "Other" Calls
}
gr2StrandCoverage <-
function(gr, beg, fin, ptitle, num_reads, all_anno_df)
{
# Uses a GRanges object to print a stranded coverage plot

	library(GenomicRanges)
	
	if(!is.null(gr)) {
		chrom <- seqnames(gr)[1]

		pos <- gr[strand(gr)=="+",]
		neg <- gr[strand(gr)=="-",]
		pos_cov <- coverage(pos)
		neg_cov <- coverage(neg)

		plotStrandedCoverage(eval(parse(text=paste("pos_cov$'", chrom, "'", sep=""))), eval(parse(text=paste("neg_cov$'", chrom, "'", sep=""))), as.integer(as.character(beg)), as.integer(as.character(fin)), ptitle, paste("Number of Reads:", num_reads, sep=" "), all_anno_df, gr)
	} else {
		plotStrandedCoverage(0, 0, as.integer(as.character(beg)), as.integer(as.character(fin)), ptitle, "", all_anno_df, gr)
	}

	return(cat("Stranded Coverage plot printed\n"))
}

getReadDensity <-
function(all_anno_df, pi_loci, base, bam=NULL)
{
# Extracts reads from a BAM file provided by the user
# to count the number of reads in the region and create
# a stranded coverage plot
#
# Args:
# bam: Character string specifying the name of the BAM file
# pi_loci: Data frame specifying the chromosome, start and 
# end position to be used in the subset
# base: A name for the picluster under investigation
# all_anno: file name for annotation file
#
# Returns:
# The total number of reads that mapped to the region

	library(systemPipeR); library(GenomicRanges)
	
	num_reads <- ""
	if(!is.null(bam)) {
		what <- c("rname", "strand", "pos", "qwidth")
		which <- GRanges(as.character(pi_loci[1,1]),IRanges(as.integer(as.character(pi_loci[1,2])), as.integer(as.character(pi_loci[1,3]))))
		param <- ScanBamParam(which=which, what=what)
		bam <- readGappedReads(bam, param=param)
		gr <- granges(bam)
		num_reads <- length(start(gr))

		gr2StrandCoverage(gr, as.integer(as.character(pi_loci[1,2])), as.integer(as.character(pi_loci[1,3])), paste(base, " - ", as.character(seqnames(gr)[1]), ":", as.integer(as.character(pi_loci[1,2])), "..", as.integer(as.character(pi_loci[1,3])), sep=""), num_reads, all_anno_df)

	## if mapping directly to picluster (rather than whole genome)
	# which <- GRanges(as.character(pi_loci[1, 1]), IRanges(1,as.integer(as.character(pi_loci[1,3])) - as.integer(as.character(pi_loci[1,2])))) ### FIXXX: move above
	#st = start(gr) + as.integer(as.character(pi_loci[1,2]))
	#en = end(gr) + as.integer(as.character(pi_loci[1,2]))
	#gr = GRanges(as.character(pi_loci[1,1]), IRanges(st,en), strand=strand(gr))

	} else { 
		gr <- NULL 
		gr2StrandCoverage(gr, as.integer(as.character(pi_loci[1,2])), as.integer(as.character(pi_loci[1,3])), base, num_reads, all_anno_df)
	}

	return(num_reads)
}

create_piClust_Summary <-
function(cluster_length, gn_summ, te_summ, oth_summ, base, all_anno_df, pi_loci, te, gn, gid, go=FALSE)
{
# Calculates the piRNA cluster summary that resides within each picluster's respective directory
#
# Allows for all annotation to be maintained, while still calculating an accurate percentage of feature occupancy
# Sense columns are only considering other sense annotations as redundant, Antisense columns are only considering other antisense annotations as redundant, and total annotations are considering both sense and antisense with the possibility of overlapping annotation (e.g. a TE in a gene)
#
# Args: 
# cluster_length: size of the picluster
# XX_summ: are with regard to cluster contents in the form of list(calls, occupancy, percentage)
# base: name for current picluster
# all_anno: file name specifying the file that contains all current annotation
# pi_loci: provides the coordinates of the picluster under investigation
#
# Returns:
# picluster level text file and plots and passes on feature nucleotide occupancy values
# 1. Number of Feature Calls
# 2. Nucleotide Occupancy of each Feature
# 3. Strandedness of Feature Calls
# 4. Stranded Coverage Plot of Reads to picluster
# 5. TE Class and Subclass Summary
#
	library(GenomicRanges); library(qcc)

	## Calculate Unannotated equivalents
	## Only in nt Occupancy & Percentage
	unann_loci <- getUnannotated(all_anno_df, "BLAST", pi_loci)
	if(class(unann_loci) == "GRanges") { unann_occ <- sum(width(unann_loci)) } else { unann_occ<-0; }
	unann_perc <- round((unann_occ / cluster_length) * 100,2) 
 
	# Positive Strand piRNA Cluster Summaries
	pos_calls <- c(gn_summ[[1]], te_summ[[1]], oth_summ[[1]], 0)
	pos_occ <- c(gn_summ[[2]], te_summ[[2]], oth_summ[[2]], unann_occ)
	pos_perc <- c(gn_summ[[3]], te_summ[[3]], oth_summ[[3]], unann_perc)

	# Negative Strand piRNA Cluster Summaries
	neg_calls <- c(gn_summ[[4]], te_summ[[4]], oth_summ[[4]], 0)
	neg_occ <- c(gn_summ[[5]], te_summ[[5]], oth_summ[[5]], 0)
	neg_perc <- c(gn_summ[[6]], te_summ[[6]], oth_summ[[6]], 0)

	# Total piRNA Cluster Summaries 
	tot_calls <- c(gn_summ[[7]], te_summ[[7]], oth_summ[[7]], 0)
	tot_occ <- c(gn_summ[[8]], te_summ[[8]], oth_summ[[8]], unann_occ)
	tot_perc <- c(gn_summ[[9]], te_summ[[9]], oth_summ[[9]], unann_perc)
	
	pie_perc <- sum(tot_perc)
	pie_perc <- tot_perc / pie_perc

	## Print Data to text table	
	# piRNA Cluster Contents - Table Summary
	lbls2 <- c("Gene", "TE", "Other", "Unannotated")

	pos_tab <- data.frame(as.numeric(pos_calls), as.numeric(pos_occ), as.numeric(pos_perc))
	colnames(pos_tab) <- c("Sense (+) Feature Calls", "Sense (+) Nucleotide Occupancy", "Sense (+) Percent piRNA Cluster Occupancy"); rownames(pos_tab) <- lbls2;
	neg_tab <- data.frame(as.numeric(neg_calls), as.numeric(neg_occ), as.numeric(neg_perc))
	colnames(neg_tab) <- c("Antisense (-) Feature Calls", "Antisense (-) Nucleotide Occupancy", "Antisense (-) Percent piRNA Cluster Occupancy"); rownames(neg_tab) <- lbls2;
	tot_tab <- data.frame(as.numeric(tot_calls), as.numeric(tot_occ), as.numeric(tot_perc))
	colnames(tot_tab) <- c("Total Number of Feature Calls", "Total Nucleotide Occupancy", "Total Percent piRNA Cluster Occupancy"); rownames(tot_tab) <- lbls2;
	
	pi_table <- cbind(pos_tab, neg_tab, tot_tab); rownames(pi_table) <- lbls2
	write.table(pi_table, paste(base,"-Summary.xls",sep=""), quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

	# Print Charts to file
	lbls <- c("Gene", "TE", "Other")	

	pdf(paste(base, "-Summary.pdf", sep=""))

	rd_den <- getReadDensity(all_anno_df, pi_loci, base)   ### FIXXX: coverage option
	pareto <- pi_table[1:3,7]
	names(pareto) <- lbls

	# Keep the colors in the plot consistent
	assoc <- data.frame(pareto, c("blue","red","yellow"))
	cols <- as.character(assoc[order(-rank(assoc[,1], ties.method="first")), 2])
	pareto.chart(pareto, ylab="Number of Feature Calls", main=base, col=cols, las=1)

	lbls3 <- paste(lbls2, " (", pie_perc, "%)", sep="")

	pareto <- pi_table[,8]
	names(pareto) <- lbls2
	# Keep the colors in the plot consistent
	assoc <- data.frame(pareto, c("blue","red","yellow","gray"))
	cols <- as.character(assoc[order(-rank(assoc[,1], ties.method="first")),2])
	pareto.chart(pareto, ylab="Nucleotide Occupancy", main=paste(base, "Contents", sep=" "), col=cols, las=1)

	# Plot sense/antisense info
	str_plot <- c(sum(pos_occ[1:3]),sum(neg_occ[1:3]))
	names(str_plot) <- c("Sense", "Antisense")
	# Keep the colors in the plot consistent
	assoc <- data.frame(str_plot, c("azure4","azure2"))
	cols <- as.character(assoc[order(-rank(assoc[,1], ties.method="first")),2])

	pareto.chart(str_plot, ylab="Nucleotide Occupancy", main=paste(base, "Annotation Orientation", sep=" "), col=cols, las=2)

	#te <- read.table(te)
	suppressWarnings(summarizeTEs(as.character(te[,1]), base))
	if(isTRUE(go)) { suppressWarnings(summarizeGenes(gid, as.character(gn[,1]))) }
	
	dev.off()
		
	return(append(pos_occ, neg_occ))
}

create_genome_Summary <- 
function(gid, fseq_xls, genome_size, fin_out, picall)
{
# Calculates the genome level summaries (in the results directory)
#
# Args:
# gid: identifier of the piClusterBuster run
# fseq_xls: specifies a file containing the loci of interest
# genome_size: size of the reference genom
# fin_out: name of the resulting XLS summary file
# picall: indicates whether of not piclusters were called during the run of
#	  piClusterBuster
#
# Returns:
# 1. Comparison of Read density in each picluster FIXXXXXXXXXXXX
# 2. Degree of Confidence in each picluster call
# 3. Average Occupancy of each Feature
# 4. Average Occupancy on each Strand
# 5. Line plot comparing Percent Genome Occupancy of piclusters
#    across species
#
	library(GenomicRanges); library(plyr)

	gn_agg_pos<-0; te_agg_pos<-0; oth_agg_pos<-0; unann_agg_pos<-0;
	gn_agg_neg<-0; te_agg_neg<-0; oth_agg_neg<-0; unann_agg_neg<-0;
	loc_perc_strand<-list(); loc_perc_contents<-list();
	piSumm_df <- read.table(fin_out, header=TRUE, sep="\t")	

	pdf("piClusters-GenomeSummary.pdf")
	opar <- par()
	par(mar=c(8,6,4,4), oma=c(0,0,0,0), xpd=TRUE)

	barplot(as.integer(piSumm_df[,10])/1000, main="piCluster Size Comparison", names.arg=as.character(rownames(piSumm_df)), cex.names=0.75, ylim=c(0, round_any(max(as.integer(piSumm_df[,10])/1000), 50, f=ceiling))) # Barplot to compare piRNA cluster lengths	

	#barplot(as.integer(piSumm_df[,10])/1000, main="piCluster Size Comparison", names.arg=as.character(piSumm_df[,9]), las=2, ylim=c(0, round_any(max(as.integer(piSumm_df[,10])/1000), 50, f=ceiling))) # Barplot to compare piRNA cluster lengths	
	title(ylab="Nucleotide Occupancy (kbps)", line=4)	

	# Get number of reads in piRNA cluster region
	tab <- read.table(fseq_xls)
	if(isTRUE(picall)) {
		#barplot(as.integer(as.character(rds)), main="piCluster Density", ylab="Number of Reads", names.arg=as.character(piSumm_df[,9]),las=2) # Barplot to the number of reads associated with each piRNA cluster
	
		barplot(as.integer(piSumm_df[,11]), ylab="Normalized Reads to Picluster", names.arg=as.character(piSumm_df[,9]), las=2) # Barplot to compare scores associated with piRNA cluster calls
	
	## Implement Average F-Seq Score of piRNA cluster call when I have multiple species? ----> put results of this into a file as I get it?
		plot(mean(as.integer(piSumm_df[,11])), type="b", main="Average Degree of Confidence in piCluster Calls", ylab="Average piRNA Cluster F-seq Score", xaxt='n', xlab="")
        	axis(1, at=1,labels=gid, las=2)
	}

	for(z in 1:length(piSumm_df[,1])) 
	{
		pos_loc<-0; neg_loc<-0; gn_loc<-0; te_loc<-0; oth_loc<-0; unann_loc<-0; loc_tot<-0;
 
		piSumm_df[z,1] <- as.integer(as.character(piSumm_df[z,1]));
		piSumm_df[z,2] <- as.integer(as.character(piSumm_df[z,2]));
		piSumm_df[z,3] <- as.integer(as.character(piSumm_df[z,3]));
		piSumm_df[z,4] <- as.integer(as.character(piSumm_df[z,4]));
		piSumm_df[z,5] <- as.integer(as.character(piSumm_df[z,5]));
		piSumm_df[z,6] <- as.integer(as.character(piSumm_df[z,6]));
		piSumm_df[z,7] <- as.integer(as.character(piSumm_df[z,7]));
		piSumm_df[z,8] <- as.integer(as.character(piSumm_df[z,8]));
	
		##	Calculate Feature Occupancy of Positive (+) strand	##
		# Calculate Total Occupancy amongst piRNA clusters
		gn_agg_pos <- gn_agg_pos + as.integer(piSumm_df[z,1]);
		te_agg_pos <- te_agg_pos + as.integer(piSumm_df[z,2]);
		oth_agg_pos <- oth_agg_pos + as.integer(piSumm_df[z,3]);
		unann_agg_pos <- unann_agg_pos + as.integer(piSumm_df[z,4]);
		
		##	Calculate Feature Occupancy of Negative (-) strand	##
		# Calculate Total Occupancy amongst piRNA clusters
		gn_agg_neg <- gn_agg_neg + as.integer(piSumm_df[z,5]);
		te_agg_neg <- te_agg_neg + as.integer(piSumm_df[z,6]);
		oth_agg_neg <- oth_agg_neg + as.integer(piSumm_df[z,7]);
		unann_agg_neg <- unann_agg_neg + as.integer(piSumm_df[z,8]); 
		
		# Calculate (or retrieve) percentages of each piRNA cluster occupancy (to be printed in stacked bar plot)
		pos_loc <- as.integer(piSumm_df[z,1]) + as.integer(piSumm_df[z,2]) + as.integer(piSumm_df[z,3])
		neg_loc <- as.integer(piSumm_df[z,5]) + as.integer(piSumm_df[z,6]) + as.integer(piSumm_df[z,7])
		str_tot <- pos_loc + neg_loc

		gn_loc <- as.integer(piSumm_df[z,1]) + as.integer(piSumm_df[z,5])
		te_loc <- as.integer(piSumm_df[z,2]) + as.integer(piSumm_df[z,6])
		oth_loc <- as.integer(piSumm_df[z,3]) + as.integer(piSumm_df[z,7])
		unann_loc <- as.integer(piSumm_df[z,4])
		
		loc_tot <- gn_loc + te_loc + oth_loc + unann_loc

		loc_perc_strand[[z]] <- c(round((pos_loc/str_tot)*100,0), round((neg_loc/str_tot)*100,0))
		
		
		corr <- c(round((gn_loc/loc_tot)*100,0), round((te_loc/loc_tot)*100,0), round((oth_loc/loc_tot)*100,0), round((unann_loc/loc_tot)*100,0))
		corr0 <- round(corr/sum(corr)*100,0)
		corr1 <- round(corr/sum(corr)*100,1)

		corr <- percentTo100(corr0, corr1, 4)

		loc_perc_contents[[z]] <- corr
	}

	par(mar=c(8,6,4,8), oma=c(0,0,0,0), xpd=TRUE)
	
	# Make stacked bar plot of features (going up) and piclust #'s (going right)	
	lbls <- c("Gene", "TE", "Other", "Unannotated")
	barplot(as.matrix(as.data.frame(loc_perc_contents)), ylab="Percentage", col=c("blue","red","yellow","gray"), names.arg=as.character(piSumm_df[,9]), las=2)
	legend("bottomright", legend=lbls, fill=c("blue","red","yellow","gray"), inset=c(-0.3,0), bty="n")

	# Barplot to show % positive and negative strand of feature calls
	barplot(as.matrix(as.data.frame(loc_perc_strand)), ylab="Percentage", names.arg=as.character(piSumm_df[,9]), las=2)
	legend("bottomright", legend=c("Sense", "Antisense"), fill=c("azure4","azure2"), inset=c(-0.3,0), bty="n")
	
	# Calculate Total Genomic Feature Occupancy
	gn_agg_tot <- gn_agg_pos + gn_agg_neg
	te_agg_tot <- te_agg_pos + te_agg_neg
	oth_agg_tot <- oth_agg_pos + oth_agg_neg	
	unann_agg_tot <- unann_agg_pos
	tot <- gn_agg_tot + te_agg_tot + oth_agg_tot + unann_agg_tot

	# Calculate Average Genomic Feature Occupancy (to compare with other species)
	gn_agg_perc <- (gn_agg_tot/tot) * 100
	te_agg_perc <- (te_agg_tot/tot) * 100
	oth_agg_perc <- (oth_agg_tot/tot) * 100
	unann_agg_perc <- (unann_agg_tot/tot) * 100
	ave_occ <- c(gn_agg_perc, te_agg_perc, oth_agg_perc, unann_agg_perc)	
	ave_occ <- round(ave_occ,0)

	# Stacked Barplot to show aggregate average Feature Occupancy of Top piRNA Clusters
	par(mar=c(8.1, 8.1, 4.1, 8.1), xpd=TRUE)
	
	lbls2 <- paste(lbls, " (", ave_occ, "%)", sep="")
	pie(ave_occ, labels=lbls2, main="Average Feature Occupancy Amongst Top piRNA Clusters", col=c("blue","red","yellow","gray"), radius=0.6)  ###FIXXXX: Don't print 0's names?
		
	# Print Text Table
	gr <- GRanges(as.character(tab[,1]), IRanges(as.integer(tab[,2]),as.integer(tab[,3])))
	
	cluster_calls <- length(tab[,1]) # Comparison of Number of piRNA Cluster Calls
	wid <- sum(width(gr))
	gen_perc <- round((wid/genome_size) * 100,2)
	gen_perc2 <- paste(gen_perc, "%", sep="")

	gen_df <- data.frame(cluster_calls, wid, gen_perc2) 
	colnames(gen_df) <- c("Number of piRNA Cluster Calls", "piRNA Cluster Nucleotide Occupancy", "Percent of Genome")
	rownames(gen_df) <- gid
	write.table(gen_df, paste(gid, "-AggregateSummary.xls", sep=""), quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t") 
	
	# Line Plot of Percent Genome Occupancy
	####### Can append species that have been analyzed, when applicable
	if(!isTRUE(picall)) {
		plot(gen_perc, type="b", main="piRNA Cluster Percent Genome Occupancy", xaxt='n', xlab="", ylab="")
		axis(1, at=1,labels=gid, las=2)	
 	}

	dev.off()

	unlink("*-TEMPORARY*")
	####### FIXXX: Can append species that have been analyzed to file, when applicable
	####### Comparison of Nucleotide Occupancy of piRNA Clusters between species
	####### Comparison of Percent Genome Occupancy of piRNA Clusters between species

	####### Aggregate comparison between species (ex: top 10 clusters)?
	return(cat("\n\n", "** piClusterBuster Run Complete **", "\n\n"))
}

print_help <-
function(wd, error=NULL)
{
# Prints a help menu describing the parameters available 
# and how to run piClusterBuster
#
	system(paste("bash ", wd, "bin/usage.sh \"", error, "\"", sep="")); 
}
