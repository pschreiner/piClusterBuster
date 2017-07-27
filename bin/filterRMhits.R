#!/usr/bin/Rscript 
ags = commandArgs(TRUE)
xls<-ags[1]; out<-ags[2]; filt<-ags[3]; wd<-ags[4]
source(paste(wd, "bin/piClusterBuster_source.R", sep=""))
library(GenomicRanges)

query<-NULL; subject<-NULL; perc_id<-NULL; len<-NULL; mis<-NULL; gap<-NULL; qstart<-NULL; qend<-NULL; sstart<-NULL; send<-NULL; meta1<-NULL; bit<-NULL;

if( length(readLines(xls))>1 ) {  
	d <- read.table(xls, skip=2, row.names=NULL, fill=TRUE, header=FALSE);
	d <- d[!is.na(d[,6]),]

	gr <- rm2GRanges(d)
        gr <- filterAnnotation(gr)

        for(ii in 1:length(gr)) 
	{
		query<-append(query, as.character(elementMetadata(gr)[[2]][ii])); subject<-append(subject,"picluster"); perc_id<-append(perc_id,as.character(filt)); len<-append(len,width(gr)[ii]); mis<-append(mis,"."); gap<-append(gap,"."); qstart<-append(qstart,"."); qend<-append(qend,"."); meta1<-append(meta1,elementMetadata(gr)[[1]][ii]); bit<-append(bit,"RM") 
		if(as.logical(strand(gr)[ii] == "-")) { sstart<-append(sstart,end(gr)[ii]); send<-append(send,start(gr)[ii])
		} else { sstart<-append(sstart,start(gr)[ii]); send<-append(send,end(gr)[ii]) }
	}

        fin <- data.frame(query, subject, perc_id, len, mis, gap, qstart, qend, sstart, send, meta1, bit) 
	write.table(fin, out, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
} else { file.create(out) }
