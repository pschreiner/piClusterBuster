process_args <-
function(args)
{
	vars=NULL; temp=NULL;
	temp.i=NULL; temp.b=NULL; temp.x=NULL; temp.n=NULL; temp.ncbidb=NULL; temp.gn=NULL; temp.te=NULL; temp.p=NULL; temp.qsub=NULL; temp.d=NULL; temp.gid=NULL;
	for(z in 1:length(args)) 
	{ 
		if(grepl("^-", args[z], perl=TRUE)) { mod=gsub("-","", args[z], fixed=TRUE); assign(paste("temp.", mod, sep=""), args[z+1]); z = z+1 }

		vars$i = temp.i; vars$b=temp.b; vars$x=temp.x; vars$n=temp.n; vars$ncbidb=temp.ncbidb; vars$gn=temp.gn; vars$te=temp.te; vars$p=temp.p; vars$qsub=temp.qsub; vars$d=temp.d; vars$gid=temp.gid;	
	}
	return(vars)
}

# Load necessary libraries with libraries needed as the argument
install_unloadedLibs <-
function(req_libs, lib_source="http://www.bioconductor.org/biocLite.R")
{
	source(lib_source)
	if(length(setdiff(req_libs, rownames(installed.packages()))) > 0) { install.packages(setdiff(req_libs, rownames(installed.packages()))) }

	return("Libraries Installed")
} 

get_ncbiDB <-
function(repository = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz")
{
	library(seqinr)
	get.ncbi(repository)
	return("\n\nDownloading NCBI nt database.  This may take several hours.\n\n")
}

get_teDB <-
function(repository = "http://www.girinst.org/server/RepBase/protected/RepBase20.12.fasta.tar.gz")
{
	wget(repository)
	return("\n\nDownloading RepBase TE database.\n\n")
}

decompress_file <-
function(file1)
{
	if(grepl(".zip",file1, fixed=TRUE)) { system(paste("unzip", file1, sep=" ")) }
	else if(grepl(".gzip",file1, fixed=TRUE) || grepl(".gz",file1, fixed=TRUE)) { system(paste("gunzip", file1, sep=" ")) }
	else if(grepl(".tar.gz",file1, fixed=TRUE) || grepl(".tgz",file1, fixed=TRUE)) { system(paste("tar -xzf", file1, sep=" ")) }
	else { return("Decompression not neccessary or unrecognized file suffix.") }
	
	return(paste(file1,"was successfully decompressed.", sep=" "))
}
#if(grepl(".gz",file1, fixed=TRUE) || grepl(".tar.gz",file1, fixed=TRUE) || grepl(".tgz",file1, fixed=TRUE)  || grepl(".gzip",file1, fixed=TRUE) || grepl(".zip",file1, fixed=TRUE))

# Filter sRNA to potential piRNAs
filterpiRNAs <-
function(sRNA_fastq)
{
	library(ShortRead); library(seqinr);

	fq=NULL;	
	if(grepl(".fq", sRNA_fastq)) { mod=gsub(".fq", "", sRNA_fastq, ignore.case=TRUE) } else if(grepl(".fastq", sRNA_fastq, ignore.case=TRUE)) { mod=gsub(".fastq", "", sRNA_fastq, fixed=TRUE) }
	mod = strsplit(mod,"/"); mod = mod[[1]][length(mod[[1]])]
	out = paste("./data/", mod, "-greater23.fa", sep="");

	#if(file.exists(out)) { print(paste("File:", out, "already exists. Skipping...", sep=" ")) } 
	#else
	#{
	fq = readFastq(sRNA_fastq)
	fq = fq[width(fq)>23];
	write.fasta(as.list(as.character(sread(fq))), as.list(as.character(id(fq))), out, open="w")
	#}
	return(out)
}
## Example Use:
## filterpiRNAs("my.fastq")

# Changes suffix from the original file
# into the same name, but with the new
# (provided) suffix
resuffix <-
function(file1, new_suff)
{
	if(grepl(".fastq", file1, fixed=TRUE)) { temp = gsub(".fastq", "", file1, fixed=TRUE); fout = paste(temp, new_suff, sep="") }
	else if(grepl(".fasta", file1, fixed=TRUE)) { temp = gsub(".fasta", "", file1, fixed=TRUE); fout = paste(temp, new_suff, sep="") }
	else if(grepl(".fq", file1, fixed=TRUE)) { temp = gsub(".fq", "", file1, fixed=TRUE); fout = paste(temp, new_suff, sep="") }
	else if(grepl(".fa", file1, fixed=TRUE)) { temp = gsub(".fa", "", file1, fixed=TRUE); fout = paste(temp, new_suff, sep="") }
	else if(grepl(".sam", file1, fixed=TRUE)) { temp = gsub(".sam", "", file1, fixed=TRUE); fout = paste(temp, new_suff, sep="") }
	else if(grepl(".bam", file1, fixed=TRUE)) { temp = gsub(".bam", "", file1, fixed=TRUE); fout = paste(temp, new_suff, sep="") }

	else( paste("Error in replacing input file suffix: ",  file1, sep=""))
	return(fout)
}
## Example Use:
## resuffix("my.fastq", ".bam")  
## Will return my.bam

run_bowtie2 <-
function(reference_genome, fastq, out)
{
	library(systemPipeR)
	moduleload("bowtie2")

	index = paste(reference_genome, ".1.bt2", sep="")
	if (!file.exists(index)) { system(paste("bowtie2-build", reference_genome, reference_genome, sep=" ")) }
	system(paste("bowtie2 -f -x ", reference_genome, " -U ", fastq, " -S ", out, sep=""))
	return("\nBowtie2 Run Complete\n")
}

sam2bam <- 
function(infile, out)
{
	library(systemPipeR)
	moduleload("samtools")

	system(paste("samtools view -Sb ", infile, " > ", out, sep=""))
	return("SAM to BAM Conversion Complete")
}

run_fseq <-
function(wd, gid, bam, bed)
{
	moduleload("bedtools")
	system(paste("bamToBed -i ", bam, " > ", bed, sep="")) 
	system(paste("fseq -d ", wd, "data/ -l 1000 -o ", wd, gid, "_results/ -of bed -t 15 ", bed, sep=""))

	out = paste(gid, "-piRNAcluster_loci.xls", sep="")
	system(paste("cat * >", out, sep=" ")) 
	
	files = list.files()
	for(i in 1:length(files))
	{
		if(!grepl(".xls", files)) { system(paste("rm ", files[i], sep="")) } # Delete all files produced by F-seq that are no longer needed 
	}

	return(out)
}

loci2seq <- 
function(ref_gen, loci, out)
{
	library(Biostrings); library(GenomicRanges); library(seqinr)
	
	seqs=NULL; heads=NULL;	
	if(class(loci)[1] != "GRanges")
	{
		for(i in 1:length(loci[,1]))
		{
			seqs = append(seqs, DNAStringSet(ref_gen[grep(as.character(loci[i,1]), names(ref_gen)),][[1]], as.integer(loci[i,2]), as.integer(loci[i,3])));
			heads = append(heads, paste(as.character(loci[i,1]), ":", as.integer(loci[i,2]), "..", as.integer(loci[i,3]), sep=""));
        	}
	}
	else
	{
		for(i in 1:length(loci[,1]))
		{
			seqs = append(seqs, DNAStringSet(ref_gen[grep(as.character(seqnames(loci)[i]), names(ref_gen)),][[1]], start(loci)[i], end(loci)[i]));
       			heads = append(heads, paste(seqnames(loci)[i], ":", start(loci)[i], "..", end(loci)[i], sep=""));
	 	}
	}
	write.fasta(as.list(as.character(seqs)), as.list(as.character(heads)), out)
	return("Genomic Loci Converted to Sequence")
}

run_blast <-
function(db, query, out, wd, i, content, qsub="TRUE")
{		
	script = paste(wd, "run_blast.sh", sep="")
	index = paste(db, ".nhr", sep="")
	name = paste("picluster", i, "-", content, "BLAST", sep="")

	if(qsub == "TRUE") { system(paste("echo \"", script, db, query, out, "\" | qsub -N", name, sep=" ")) } 
	else 
	{ 
		library(systemPipeR); moduleload("ncbi-blast");
		if(!file.exists(index)) { system(paste("makeblastdb -dbtype nucl -in", db, "-out", db, sep=" ")) }
		system(paste("blastn -db", db, "-outfmt 6 -evalue 1e-3 -query", query, "-out", out, sep=" "))
	}
	return(paste("BLAST Search Complete -  DB:", db, "Query:", query, sep=" "))
}

#retrieveUnannotatedRegions <-
#function(direc,ref_genome,supercont,ostart,oend,size_unann,out)
#{
#	library(Biostrings); library(seqinr)
#	seqs=NULL; max = ostart;
#	setwd(direc)
#
#	system("cat *.BLAST > ALL.BLAST")
#	data = read.table("ALL.BLAST",fill=TRUE)
#	colnames(data) = c("qid","sid","perc_id","len","mism","gap","qstart","qend","sstart","send","evalue","bit")
#	data = data[order(data$sstart),]
#
#	data$sstart = data$sstart + ostart;
#	data$send = data$send + ostart;
#
#	scaff_ref <- readDNAStringSet(ref_genome)
#
#	for(i in 1:length(data[,1]))
#	{
#		# Account for annotations that are contained within other annotations
#		if( max < data$send[i] ) { max = data$send[i]}
#
#		# Ends at last Annotation 
#		if( (i+1) > length(data[,1]) ) {break;}
#		# Get sequences of unannotated loci above "size_unann" threshold
#		else if( (max+size_unann) <= data$sstart[i+1] )
#		{
#			temp = DNAStringSet(scaff_ref[grep(paste(supercont," ",sep=""),names(scaff_ref)),][[1]],(data$send[i]+1),(data$sstart[i+1]-1));
#			if(is.null(seqs)) {seqs = temp } else{ seqs = c(seqs,temp) }
#			temp = paste(supercont,":",(data$send[i]+1),"..",(data$sstart[i+1]-1),sep="");
#			names(seqs)[length(seqs)] = as.character(temp)	
#		}
#		# If not big enough gap in annotation, continue to search
#		else { next; }
#	}
#	write.fasta(as.list(as.character(seqs)),as.list(names(seqs)),out)
#
#	x = paste("** Retrieved Unannotated sequences. The results have been printed to: ", out, sep="");
#	return(x)
#}

filterBLASThits <-
function(xls, out, wd)
{
	script=paste(wd, "filterBLASThits.R", sep="")
	system(paste("Rscript", script, xls, out, wd, sep=" "))
	
	return(d);
}
## Example Use:
## filterBLASThits("./piclust12/intergenic_piclust12_wb-greater2kb-NCBInt.BLAST", "./piclust12/intergenic_piclust12_wb-greater2kb-NCBInt-FILTERED.BLAST")

## Summary of the feature of interest
feature_summary <-
function(blast_file, piclust_size)
{
	if(!exists("blast_file")) { return(list(0,0,0,0,0,0,0,0,0)) }
	feat = read.table(blast_file)
	feat = feat[order(feat[,11]),]

	beg=NULL; fin=NULL; strand=NULL;
	for(i in 1:length(feat[,1]))
	{
	    if(as.integer(feat[i,9]) > as.integer(feat[i,10]))
	    {
		    beg = append(beg, as.integer(feat[i,10]))
		    fin = append(fin, as.integer(feat[i,9]))
		    strand = append(strand, "-")
	    }
	    else
	    {
		    beg = append(beg, as.integer(feat[i,9]))
		    fin = append(fin, as.integer(feat[i,10]))
		    strand = append(strand, "+")
	    }
	}

	gr = GRanges(as.character(feat[,2]), IRanges(as.integer(beg),as.integer(fin)), strand=strand)
	gr = reduce(gr)
	tot_feat_calls = length(start(gr))
	tot_feat_occ = sum(width(gr))
	tot_feat_perc = ( tot_feat_occ / piclust_size ) * 100
	tot_feat_perc = round(tot_feat_perc, digits = 2)

	# Check Annotations on positive strand
	pos_gr = gr[strand(gr)=="+",]
	pos_feat_calls = length(start(pos_gr))
	pos_feat_occ = sum(width(pos_gr))
	pos_feat_perc = ( pos_feat_occ / piclust_size ) * 100
	pos_feat_perc = round(pos_feat_perc, digits = 2)

	# Check Annotations on negative strand
	neg_gr = gr[strand(gr)=="-",]
	neg_feat_calls = length(start(neg_gr))
	neg_feat_occ = sum(width(neg_gr))
	neg_feat_perc = ( neg_feat_occ / piclust_size ) * 100
	neg_feat_perc = round(neg_feat_perc, digits = 2)

	return(list(pos_feat_calls, pos_feat_occ, pos_feat_perc, neg_feat_calls, neg_feat_occ, neg_feat_perc, tot_feat_calls, tot_feat_occ, tot_feat_perc))
}
## Example Use;
## feature_summary("piclust15_TElandscape.BLAST")

### 	Arguments are with regard to cluster contents: 
###	list(calls, occupancy, percentage)
create_piClust_Summary <- 
function(cluster_length, gn_summ, te_summ, oth_summ, base)
{
	## Calculate Unannotated equivalents
	## Only in nt Occupancy & Percentage
	unann_occ = cluster_length - (gn_summ[[8]] + te_summ[[8]] + oth_summ[[8]])  #FIXXXXX: unannotated value can go into the negative here!!
	unann_perc = 1 - (gn_summ[[9]] + te_summ[[9]] + oth_summ[[9]])

	# Positive Strand piRNA Cluster Summaries
	pos_calls = c(gn_summ[[1]], te_summ[[1]], oth_summ[[1]], " ")
	pos_occ = c(gn_summ[[2]], te_summ[[2]], oth_summ[[2]], unann_occ)
	pos_perc = c(gn_summ[[3]], te_summ[[3]], oth_summ[[3]], unann_perc)

	# Negative Strand piRNA Cluster Summaries
	neg_calls = c(gn_summ[[4]], te_summ[[4]], oth_summ[[4]], " ")
	neg_occ = c(gn_summ[[5]], te_summ[[5]], oth_summ[[5]], unann_occ)
	neg_perc = c(gn_summ[[6]], te_summ[[6]], oth_summ[[6]], unann_perc)

	# Total piRNA Cluster Summaries 
	tot_calls = c(gn_summ[[7]], te_summ[[7]], oth_summ[[7]], " ")
	tot_occ = c(gn_summ[[8]], te_summ[[8]], oth_summ[[8]], unann_occ)
	tot_perc = c(gn_summ[[9]], te_summ[[9]], oth_summ[[9]], unann_perc)

	## Print Data to text table	
	# piRNA Cluster Contents - Table Summary
	lbls2 <- c("Gene", "TE", "Other", "Unannotated")

	pos_tab = data.frame(as.numeric(pos_calls), as.numeric(pos_occ), as.numeric(pos_perc))
	colnames(pos_tab) = c("Sense (+) Feature Calls", "Sense (+) Nucleotide Occupancy", "Sense (+) Percent piRNA Cluster Occupancy"); rownames(pos_tab) = lbls2;
	neg_tab = data.frame(as.numeric(neg_calls), as.numeric(neg_occ), as.numeric(neg_perc))
	colnames(neg_tab) = c("Antisense (-) Feature Calls", "Antisense (-) Nucleotide Occupancy", "Antisense (-) Percent piRNA Cluster Occupancy"); rownames(neg_tab) = lbls2;
	tot_tab = data.frame(as.numeric(tot_calls), as.numeric(tot_occ), as.numeric(tot_perc))
	colnames(tot_tab) = c("Total Number of Feature Calls", "Total Nucleotide Occupancy", "Total Percent piRNA Cluster Occupancy"); rownames(tot_tab) = lbls2;
	
	pi_table = cbind(pos_tab, neg_tab, tot_tab)
	write.table(pi_table, paste(base,"-Summary.xls",sep=""), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

	# Print Charts to file
	lbls <- c("Gene", "TE", "Other")	

	pdf(paste(base, "-Summary.pdf", sep=""))
	barplot(pi_table[1:3,7], main=base, ylab="Number of Feature Calls", names.arg=lbls)
	barplot(pi_table[,8], main=base, ylab="Nucleotide Occupancy", names.arg=lbls2)
	
	lbls3 <- paste(lbls2, " (", tot_perc, "%)", sep="")

	### FIXXXXX: have to fix margins
	pie(pi_table[,9], labels=lbls3, main=paste(base, "- Content Summary",sep=" "))
	dev.off()
		
	return(append(pos_occ, neg_occ))
}

annotate_piCluster <-
function(i, ref_gen, f_res, gn_set, te_set, ncbi_nt, wd, results_dir, qsub="TRUE")
{
	setwd(wd)
	
	base = paste("picluster", i, sep="")
	f1 = paste("\"", wd, "run_annotate_piCluster.sh", sep="");	
	f2 = paste(wd, "annotate_piCluster.R", sep="")
	 	
	if(qsub=="TRUE") { system(paste("echo", f1, i, ref_gen, f_res, gn_set, te_set, ncbi_nt, wd, results_dir, qsub, "\" | qsub -N", base, sep=" "))
	} else { system(paste("Rscript", f2, i, ref_gen, f_res, gn_set, te_set, ncbi_nt, wd, results_dir, qsub, sep=" ")) } 

	return(cat("\npiRNA Cluster", i, "submitted for processing\n"))
}

create_genome_Summary <- 
function(gid, fseq_xls, genome_size, piSumm_df)
{
	gn_agg_pos=0; te_agg_pos=0; oth_agg_pos=0; unann_agg_pos=0;
	gn_agg_neg=0; te_agg_neg=0; oth_agg_neg=0; unann_agg_neg=0;
	loc_perc_strand=list(); loc_perc_contents=list();
	
	pdf("piCluster-GenomeSummary.pdf")
	### FIXXX: set parameters! (even if it's just defaults)
	barplot(as.integer(piSumm_df[,10]), main="piCluster Size Comparison", ylab="Nucleotide Occupancy", names.arg=as.character(piSumm_df[,9])) # Barplot to compare piRNA cluster lengths	
	
	for(z in 1:length(piSumm_df[,1])) 
	{
		pos_loc=0; neg_loc=0; gn_loc=0; te_loc=0; oth_loc=0; unann_loc=0; loc_tot=0;
 
		##	Calculate Feature Occupancy of Positive (+) strand	##
		# Calculate Total Occupancy amongst piRNA clusters
		gn_agg_pos = gn_agg_pos + piSumm_df[z,1]; te_agg_pos = te_agg_pos + piSumm_df[z,3]; oth_agg_pos = oth_agg_pos + piSumm_df[z,5]; unann_agg_pos = unann_agg_pos + piSumm_df[z,7];
		##	Calculate Feature Occupancy of Negative (-) strand	##
		# Calculate Total Occupancy amongst piRNA clusters
		gn_agg_neg = gn_agg_neg + piSumm_df[z,2]; te_agg_neg = te_agg_neg + piSumm_df[z,4]; oth_agg_neg = oth_agg_neg + piSumm_df[z,6]; unann_agg_neg = unann_agg_neg + piSumm_df[z,8]; 
		
		# Calculate (or retrieve) percentages of each piRNA cluster occupancy (to be printed in stacked bar plot)
		pos_loc = piSumm_df[z,1] + piSumm_df[z,3] + piSumm_df[z,5] + piSumm_df[z,7]
		neg_loc = piSumm_df[z,2] + piSumm_df[z,4] + piSumm_df[z,6] + piSumm_df[z,8]

		gn_loc = piSumm_df[z,1] + piSumm_df[z,2]
		te_loc = piSumm_df[z,3] + piSumm_df[z,4]
		oth_loc = piSumm_df[z,5] + piSumm_df[z,6]
		unann_loc = piSumm_df[z,7] + piSumm_df[z,8]
		
		loc_tot = gn_loc + te_loc + oth_loc + unann_loc

		loc_perc_strand[[z]] = c(round((pos_loc/loc_tot)*100,0), round((neg_loc/loc_tot)*100,0))
		loc_perc_contents[[z]] = c(round((gn_loc/loc_tot)*100,0), round((te_loc/loc_tot)*100,0), round((oth_loc/loc_tot)*100,0), round((unann_loc/loc_tot)*100,0))
	}
	
	par(mfrow=c(1,1), mar=c(2, 3, 2, 7), oma=c(1,1,1,1))

	# Make stacked bar plot of features (going up) and piclust #'s (going right)	
	lbls = c("Gene", "TE", "Other", "Unannotated")
	barplot(as.matrix(as.data.frame(loc_perc_contents)), legend.text=lbls, args.legend = list(x = "bottomright", bty= "n", inset=c(-0.2, 0))) 

	# Barplot to show % positive and negative strand of feature calls	
	barplot(as.matrix(as.data.frame(loc_perc_strand)), legend.text=c("Sense", "Antisense"), args.legend = list(x = "bottomright", bty= "n", inset=c(-0.2, 0)))
	
	# Calculate Total Genomic Feature Occupancy
	gn_agg_tot = gn_agg_pos + gn_agg_neg
	te_agg_tot = te_agg_pos + te_agg_neg
	oth_agg_tot = oth_agg_pos + oth_agg_neg	
	unann_agg_tot = unann_agg_pos + unann_agg_neg
	tot = gn_agg_tot + te_agg_tot + oth_agg_tot + unann_agg_tot

	# Calculate Average Genomic Feature Occupancy (to compare with other species)
	gn_agg_perc = (gn_agg_tot/tot) * 100
	te_agg_perc = (te_agg_tot/tot) * 100
	oth_agg_perc = (oth_agg_tot/tot) * 100
	unann_agg_perc = (unann_agg_tot/tot) * 100
	ave_occ = c(gn_agg_perc, te_agg_perc, oth_agg_perc, unann_agg_perc)	

	# Stacked Barplot to show aggregate average Feature Occupancy of Top piRNA Clusters
	lbls2 = paste(lbls, " (", ave_occ, "%)", sep="")
	pie(ave_occ, labels=lbls2, main="Average Feature Occupancy amongst Top piRNA Clusters")
		
	# Print Text Table
	tab = read.table(fseq_xls)
	
	gr = GRanges(as.character(tab[,1]), IRanges(as.integer(tab[,2]),as.integer(tab[,3])))
	
	cluster_calls = length(tab[,1]) # Comparison of Number of piRNA Cluster Calls by F-seq
	wid = sum(width(gr))
	gen_perc = paste((wid/genome_size) * 100, "%", sep="")
	
	gen_df = data.frame(cluster_calls, wid, gen_perc) 
	colnames(gen_df) = c("Number of piRNA Cluster Calls", "piRNA Cluster Nucleotide Occupancy", "Percent of Genome")
	write.table(gen_df, paste(gid, "-AggregateSummary.xls"), quote=FALSE, col.names=TRUE, row.names=lbls, sep="\t") 
	
	# Line Plot of Percent Genome Occupancy
	####### Can append species that have been analyzed, when applicable
	plot(gen_perc, type="b", main="piRNA Cluster Percent Genome Occupancy")
	text(axTicks(1),labels = c("TEST"))  #### FIXXXXX: obviously

	dev.off()

	unlink("../*-TEMPORARY.xls")
	####### FIXXX: Can append species that have been analyzed, when applicable
	####### Comparison of Nucleotide Occupancy of piRNA Clusters between species
	####### Comparison of Percent Genome Occupancy of piRNA Clusters between species

	####### Aggregate comparison between species (ex: top 10 clusters)?
	return(cat("\n\n", "** piClusterBuster Run Complete.", "\n\n"))
}

print_help <-
function(error=NULL)
{
	system(paste("Rscript usage.R ", error, sep="")); 
}
