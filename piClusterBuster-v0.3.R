#!/usr/bin/Rscript
source("piClusterBuster_source-v0.3.R")

options(echo=TRUE)
ags = commandArgs(TRUE)
args = process_args(ags)

wd=paste(getwd(), "/", sep="")

####### Bowtie2 and BLAST deviations from default 
####### Word size, E-value

####################################
##     Step 0: User Arguments     ##
####################################
## Setting Argument Defaults
if(is.null(args$i) && is.null(args$b))	{ print_help("No FASTQ file or BAM file were indicated"); q(save="no") 
} else if (exists("args$b")) { bam_file=normalizePath(args$b) 
} else { fastq_file=normalizePath(args$i) }

if(is.null(args$x))	{ print_help("Reference genome not indicated"); q(save="no") 
} else { ref_gen=normalizePath(args$x) }

if(is.null(args$n))	{ num_clusters=10 
} else { num_clusters=args$n  }
#if(num_clusters == "all") { last = length(pi_loci[,1]) } else { last = as.integer(num_clusters) }

if(is.null(args$ncbidb)) { print_help("No NCBI nt database indicated.  You can choose to have piClusterBuster retrieve the database by using \"-ncbidb=\"NA\"\""); q(save="no")
} else	{ ncbi_nt = normalizePath(args$ncbidb) }
#else if(args$ncbidb == "NA") { get_ncbiDB(); x=getwd(); ncbi_nt = paste(x, "/nt.gz", sep="") } ## FIXXXXX: kill program

if(is.null(args$gn)) { gn_set = ncbi_nt
} else	{ gn_set=normalizePath(args$gn) }

if(is.null(args$te)) { print_help("No TE database indicated.  You can choose to have piClusterBuster retrieve the database by using \"-tedb=\"NA\"\""); q(save="no")
} else	{ te_set=normalizePath(args$te) }
#else if(args$te == "NA")  { get_teDB(); x=getwd(); te_set = paste(x, "/RepBase20.12.fasta.tar.gz", sep="") }  # FIXXXX: kill program   

if(is.null(args$qsub)) { qsub="FALSE"
} else	{ qsub="TRUE" }

if(is.null(args$d)) { setwd(wd)
} else	{ setwd(args$d) }

if(is.null(args$gid)) { gid="genome"
} else	{ gid=args$gid }

###### Add verbose option?

## Check for and Load Appropriate Libraries
install_unloadedLibs(c("Biostrings", "systemPipeR", "foreach", "doMC", "R.utils"))
suppressPackageStartupMessages(library(Biostrings)); suppressPackageStartupMessages(library(systemPipeR)); suppressPackageStartupMessages(library(foreach)); suppressPackageStartupMessages(library(doMC)); suppressPackageStartupMessages(library(R.utils));
source("piClusterBuster_source-v0.2.R")

# Option to run in parallel (and haven't specified number of cores)
if(is.null(args$p)) { num_cores = 1; registerDoMC(num_cores)
} else	{ num_cores=args$p; registerDoMC(num_cores) }

#FIXXXX: Provide some kind of update - bash script that reacts given a number (ex: "Libraries loaded")
cat(sprintf("\n", "R Libraries Loaded", "\n"))

##########################################
##	Step 1: Filter piRNAs and	##
## 	map piRNAs to reference genome	##
## 	(if BAM file is not provided)	##
##########################################
if(exists("fastq_file")) 
{ 
	fa_file = filterpiRNAs(fastq_file) # Filter out non-piRNAs (by length)

	b2_out = resuffix(fastq_file, "-piRNAs.sam")
	run_bowtie2(ref_gen, fa_file, b2_out)
 
	bam = resuffix(b2_out, ".bam")
	sam2bam(b2_out, bam)
} else if(grepl(".sam", bam_file)) { sam2bam(bam_file) }

##########################################
##	Step 2: Run F-seq program	##
##########################################

##### FIXXXX: If file exists, move it elsewhere
results_dir = paste(gid,"_results",sep=""); 
#results_dir2 = paste(gid, "_results-", Sys.Date(), sep="");
if(file.exists("results_dir")) { system(paste("mv", results_dir, paste(results_dir, "-OLD", sep=""), sep=" ")); dir.create(paste(results_dir, sep="")); setwd(results_dir) 
} else { dir.create(paste(results_dir, sep="")); setwd(results_dir) }

bed = resuffix(bam, ".bed")
#f_res = run_fseq(wd, gid, bam, bed) ###### FIXXXXXX: didn't confirm that it worked


###################################################################
#### DELETE: just for testing
system("cp ~/scripts/piClusterBuster/test.xls .") 
f_res = "~/scripts/piClusterBuster/test_results/test.xls"
###################################################################






pi_loci = read.table(f_res)
pi_loci = pi_loci[order(-pi_loci[,5]),]
write.table(pi_loci, f_res, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

##########################################################
##	Step 3: Run Summary for each piRNA Cluster	##
##########################################################
i=1;
#foreach (i:num_clusters) %dopar%
for (i in 1:num_clusters)
{ annotate_piCluster(i, ref_gen, f_res, gn_set, te_set, ncbi_nt, wd, results_dir, qsub) }

num_files = 0;
while(num_files != num_clusters) # Allow system to sleep until the individual piRNA clusters have been analyzed
{
	num_files = length(list.files(path="..",pattern="*-TEMPORARY.xls")) # Test to see how many files are done
	Sys.sleep(60) 
}
for(i in 1:num_clusters) { assign(paste("piclust_temp",i,sep=""), read.table(paste("picluster",i,"-TEMPORARY.xls",sep=""))) }

piSumm_df=NULL;
temps = ls()[grepl("piclust_temp",ls())]
for(i in 1:length(temps)) { piSumm_df = cbind(piSumm_df, eval(parse(text = temps[i]))); }
piSumm_df = t(piSumm_df); fin_out = "piCluster-GenomeSummary.xls";
colnames(piSumm_df) = c("Pos_Gene", "Neg_Gene", "Pos_TE", "Neg_TE", "Pos_Other", "Neg_Other", "Pos_Unann", "Neg_Unann", "Cluster Name", "Length")
piSumm_df = piSumm_df[order(piSumm_df[,9]),]
write.table(piSumm_df, fin_out, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

ref_gen = readDNAStringSet(ref_gen)
create_genome_Summary(gid, f_res, sum(width(ref_gen)), piSumm_df)
