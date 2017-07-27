#!/usr/bin/Rscript
ags = commandArgs(TRUE)
xls=ags[1]; out=ags[2]; filt=ags[3]; wd=ags[4];

source(paste(wd, "bin/piClusterBuster_source.R", sep="")
library(GenomicRanges)

query<-NULL; subject<-NULL; perc_id<-NULL; len<-NULL; mis<-NULL; gap<-NULL; qstart<-NULL; qend<-NULL; sstart<-NULL; send<-NULL; sim<-NULL; bit<-NULL;

if( file.info(as.character(xls))$size==0 ) {  system(paste("cp", xls, out, sep=" ")); q(save="no")
} else {
	d <- read.table(xls);
		
	# A. aegypti insignificants
	d <- d[!(grepl("EF",as.character(d[,1]))),];
	d <- d[!(grepl("AC",as.character(d[,1]))),]
	d <- d[!(grepl("AF",as.character(d[,1]))),]
	d <- d[!(grepl("AJ",as.character(d[,1]))),]
	d <- d[!(grepl("AY",as.character(d[,1]))),]
	d <- d[!(grepl("EU",as.character(d[,1]))),]
	d <- d[!(grepl("XM",as.character(d[,1]))),]
	d <- d[!(grepl("CR",as.character(d[,1]))),]
	d <- d[!(grepl("FJ",as.character(d[,1]))),]
	d <- d[!(grepl("LTR",as.character(d[,1]))),]
		
	# D. melanogaster insignificants
	d <- d[!(grepl("CP",as.character(d[,1]))),]
	d <- d[!(grepl("CU",as.character(d[,1]))),]
	d <- d[!(grepl("AE",as.character(d[,1]))),]
	d <- d[!(grepl("V[0-9]",as.character(d[,1]))),]
	d <- d[!(grepl("M[0-9]",as.character(d[,1]))),]
	d <- d[!(grepl("AL",as.character(d[,1]))),]
	d <- d[!(grepl("DS",as.character(d[,1]))),]

	# M. musculus insignificants
	d <- d[!(grepl("CT",as.character(d[,1]))),]
	d <- d[!(grepl("AK",as.character(d[,1]))),]

	if(length(d[,1])==0) { system(paste("cp", xls, out, sep=" ")); q(save="no") }

	gr <- blast2GRanges(d)
	gr <- filterAnnotation(gr)

	for(ii in 1:length(gr)) 
	{
		query<-append(query, as.character(seqnames(gr)[ii])); subject<-append(subject,"picluster"); perc_id<-append(perc_id,as.character(filt)); len<-append(len,width(gr)[ii]); mis<-append(mis,"."); gap<-append(gap,"."); qstart<-append(qstart,"."); qend<-append(qend,"."); sim<-append(sim,elementMetadata(gr)[["sim"]][ii]); bit<-append(bit,"BLAST")
		if(as.logical(strand(gr)[ii] == "-")) { sstart<-append(sstart,end(gr)[ii]); send<-append(send,start(gr)[ii])
		} else { sstart<-append(sstart,start(gr)[ii]); send<-append(send,end(gr)[ii]) }
	}
	fin <- data.frame(query, subject, perc_id, len, mis, gap, qstart, qend, sstart, send, sim, bit)
}

write.table(fin, out, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
