#!/usr/bin/Rscript
source("./bin/piClusterBuster_source.R")

options(echo=TRUE)
ags <- commandArgs(TRUE)
args <- processArgs(ags)

wd<-paste(getwd(), "/", sep="")

####################################
##     Step 0: User Arguments     ##
####################################
## Setting Argument Defaults
if(is.null(args$x))	{ print_help(wd, "Reference genome was not indicated"); q(save="no") 
} else { ref_gen<-normalizePath(args$x) }

if(is.null(args$fa) && is.null(args$fq) && is.null(args$bed))	{ print_help(wd, "No FASTA, FASTQ, nor piRNA Cluster file was indicated"); q(save="no") }
if(!is.null(args$fa)) {
	fa <- normalizePath(args$fa);
	fq <- FALSE;
	bed <- FALSE;
        picall <- FALSE;
}
if(!is.null(args$fq)) {
	fa <- FALSE;
        fq <- normalizePath(args$fq);
        bed <- FALSE;
        picall <- TRUE;
}
if(!is.null(args$bed)) { 
	fa <- FALSE;
	fq <- FALSE;
	bed <- normalizePath(args$bed); 
	picall <- FALSE;
	pi_loci <- read.table(bed);
}

num_clusters <- 5 
if(!is.null(args$n))	{ num_clusters <- args$n }

ncbi_nt <- "FALSE"
if(!is.null(args$ncbidb)) { ncbi_nt <- normalizePath(args$ncbidb) }

if(is.null(args$gndb)) { print_help(wd, "No Gene database was indicated."); q(save="no");
} else	{ gn_set<-normalizePath(args$gndb) }

if(is.null(args$tedb)) { print_help(wd, "No TE database was indicated.  You can choose to have piClusterBuster retrieve the database by using -tedb \"NA\""); q(save="no")
} else	{ te_set<-normalizePath(args$tedb) }

if(is.null(args$qsub)) { qsub<-FALSE
} else	{ qsub<-TRUE }

if(is.null(args$srun)) { srun<-FALSE
} else  { srun<-TRUE; qsub<-FALSE }

if(is.null(args$d)) { setwd(wd)
} else	{ setwd(args$d) }

if(is.null(args$gid)) { gid<-"genome"
} else	{ gid<-args$gid }

if(is.null(args$go)) { go<-FALSE
} else  { go<-TRUE }

if(is.null(args$allsrna)) { allsrna<-FALSE
} else	{ allsrna<-TRUE }

if(is.null(args$verbose)) { verbose<-FALSE
} else { verbose<-TRUE  }

## Check for and Load Appropriate Libraries
install_unloadedLibs(c("Biostrings", "systemPipeR", "plyr", "doMC", "R.utils", "qcc", "plyr", "gProfileR", "seqinr"))
suppressPackageStartupMessages(library(Biostrings)); suppressPackageStartupMessages(library(seqinr)); suppressPackageStartupMessages(library(systemPipeR)); suppressPackageStartupMessages(library(foreach)); suppressPackageStartupMessages(library(doMC)); suppressPackageStartupMessages(library(R.utils)); suppressPackageStartupMessages(library(qcc)); suppressPackageStartupMessages(library(plyr)); suppressPackageStartupMessages(library(gProfileR));
source("./bin/piClusterBuster_source.R")

## Option to run in parallel (and haven't specified number of cores)
if(is.null(args$p)) { threads <- 1; registerDoMC(threads)
} else	{ threads<-args$p; registerDoMC(threads) }

## Provide update 3
cat(sprintf("\n", "R Libraries Loaded", "\n"))

##########################################
##	Step 1: Filter piRNAs and	##
## 	map piRNAs to reference genome	##
## 	(if BAM file is not provided)	##
##########################################
results_dir <- paste(wd, gid,"_results",sep=""); 
if(file.exists(results_dir)) { system(paste("mv -f", results_dir, paste(results_dir, "-OLD", sep=""), sep=" ")); dir.create(paste(results_dir, sep="")); setwd(results_dir)
} else { dir.create(paste(results_dir, sep="")); setwd(results_dir) }  
results_dir <- paste(results_dir,"/",sep="");

setwd(results_dir)

if(fq != "FALSE" && !isTRUE(allsrna)) { fq <- filterpiRNAs(fq, 24) } # Filter out non-piRNAs (by length)  FIXXXXXX: confirm this works?

##########################################
##	Step 2: Run proTRAC program	##
##########################################
if(isTRUE(picall)) {
	f_res <- run_proTRAC(fq, ref_gen, num_clusters, wd, results_dir, gid, verbose)
} else if(fa != "FALSE") {
	spl <- strsplit(as.character(fa), "/")
        f_res <- paste(results_dir, spl[[1]][length(spl[[1]])], sep="")	
	system(paste("cp", fa, f_res, sep=" "))

	f_res <- resuffix(f_res, ".bed")

	fst <- read.fasta(as.character(fa), as.string=TRUE, forceDNAtolower=FALSE)
	num_clusters <- length(fst)
	
	temp <- strsplit(as.character(attr(fst, "name")), ":", fixed=TRUE)   ###FIXXX: stop with error if not found
	one<-NULL; temp2<-NULL;
	for(ii in 1:length(temp)) { one <- append(one, as.character(temp[[ii]][1])); temp2 <- append(temp2, as.character(temp[[ii]][2])) }

	temp <- strsplit(temp2, "..", fixed=TRUE)
	two<-NULL; three<-NULL;
	for(ii in 1:length(temp)) { two <- append(two, as.character(temp[[ii]][1])); three <- append(three, as.character(temp[[ii]][2])) }
	pi_loci <- data.frame(one, two, three)

	write.table(pi_loci, f_res, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
} else {
	spl <- strsplit(bed,"/")
	f_res <- paste(results_dir, spl[[1]][length(spl[[1]])], sep="")
	file.copy(bed, f_res)
}

##########################################################
##	Step 3: Run Summary for each piRNA Cluster	##
##########################################################
for (i in 1:num_clusters)
{
	base <- paste("picluster", i, sep="")
	f2 <- paste(wd, "bin/annotate_piCluster.R", sep="")

	if(isTRUE(qsub)) { 
		system(paste("echo \"Rscript", f2, i, ref_gen, f_res, gn_set, te_set, ncbi_nt, wd, results_dir, qsub, srun, threads, picall, verbose, gid, fa, go, "\" | qsub -l mem=6gb -N", base, sep=" "))
	} else if(isTRUE(srun)) {
		j1 <- paste("--job-name=", base, sep="") 
		system(paste("sbatch --mem-per-cpu=6G", j1, "--wrap=\"Rscript", f2, i, ref_gen, f_res, gn_set, te_set, ncbi_nt, wd, results_dir, qsub, srun, threads, picall, verbose, gid, fa, go, "\"", sep=" "))
	} else {
		system(paste("Rscript", f2, i, ref_gen, f_res, gn_set, te_set, ncbi_nt, wd, results_dir, qsub, srun, threads, picall, verbose, gid, fa, go, sep=" "))
	}		
	cat(sprintf("\npiRNA Cluster", i, " submitted for processing\n"))
}

num_files <- 0; elapsed = 0;
while(num_files != num_clusters) # Allow system to sleep until the individual piRNA clusters have been analyzed
{
	num_files <- length(list.files(pattern="*-TEMPORARY.xls")) # Test to see how many files are done
	Sys.sleep(60)

	## Roughly every hour, remind user to check that jobs ran (or are running) properly
	elapsed = elapsed + 60;
	if(elapsed > 54001) { print("piClusterBuster is taking much longer than expected.  Be sure that the submitted jobs or processes have been completed successfully."); elapsed <- elapsed - 3601 }
}

#################################################
##  Step 4: Data Processing and Visualization  ##
#################################################
for(i in 1:num_clusters) 
{
	if(file.exists(paste("picluster", i, "-TEMPORARY.xls", sep=""))) { 
		assign(paste("piclust_temp", i, sep=""), read.table(paste("picluster", i, "-TEMPORARY.xls", sep=""))) 
	}
}
piSumm_df<-NULL;
temps <- ls()[grepl("piclust_temp",ls())]
for(i in 1:length(temps)) 
{ 
	if(is.null(piSumm_df)) { piSumm_df <- eval(parse(text = temps[i])) }
	else { piSumm_df <- cbind(piSumm_df, eval(parse(text = temps[i]))); }
}
piSumm_df <- t(piSumm_df)
fin_out <- "piClusters-GenomeSummary.xls";
if(length(piSumm_df[1,]) == 10) { colnames(piSumm_df) <- c("Pos_Gene", "Pos_TE", "Pos_Other", "Pos_Unann", "Neg_Gene", "Neg_TE", "Neg_Other", "Neg_Unann", "Cluster_Name", "Length") 
} else {
	colnames(piSumm_df) <- c("Pos_Gene", "Pos_TE", "Pos_Other", "Pos_Unann", "Neg_Gene", "Neg_TE", "Neg_Other", "Neg_Unann", "Cluster_Name", "Length", "Normalized_Reads")
}
piSumm_df <- cbind(piSumm_df,nchar(piSumm_df[,9]))  # Allows for 1,2,3 sorting
piSumm_df <- piSumm_df[order(piSumm_df[,length(piSumm_df[1,])],piSumm_df[,grep("Cluster_Name",colnames(piSumm_df))]),]
piSumm_df <- piSumm_df[,1:(length(piSumm_df[1,])-1)] # return all but column for sorting
write.table(piSumm_df, fin_out, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

ref_gen <- readDNAStringSet(ref_gen)
suppressWarnings(create_genome_Summary(gid, f_res, sum(as.numeric(width(ref_gen))), fin_out, picall))
if(!isTRUE(verbose)) { unlink("*.e*"); unlink("*.o*") }
