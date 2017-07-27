#!/usr/bin/Rscript
ags = commandArgs(TRUE)
xls=ags[1]; out=ags[2];
source("../../bin/piClusterBuster_source.R")
library(GenomicRanges)

query<-NULL; subject<-NULL; perc_id<-NULL; len<-NULL; mis<-NULL; gap<-NULL; qstart<-NULL; qend<-NULL; sstart<-NULL; send<-NULL; sim<-NULL; bit<-NULL;

d <- read.table(xls);
gr <- blast2GRanges(d)
gr <- filterAnnotation(gr)

for(ii in 1:length(gr)) 
{
	query<-append(query, as.character(seqnames(gr)[ii])); subject<-append(subject,"picluster"); perc_id<-append(perc_id,as.character(elementMetadata(gr)[["feature"]][ii])); len<-append(len,width(gr)[ii]); mis<-append(mis,"."); gap<-append(gap,"."); qstart<-append(qstart,"."); qend<-append(qend,"."); sim<-append(sim,elementMetadata(gr)[["sim"]][ii]); bit<-append(bit,elementMetadata(gr)[["method"]][ii])
	if(as.logical(strand(gr)[ii] == "-")) { sstart<-append(sstart,end(gr)[ii]); send<-append(send,start(gr)[ii])
	} else { sstart<-append(sstart,start(gr)[ii]); send<-append(send,end(gr)[ii]) }
}
fin <- data.frame(query, subject, perc_id, len, mis, gap, qstart, qend, sstart, send, sim, bit)

write.table(fin, out, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
