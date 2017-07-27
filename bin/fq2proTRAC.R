#!/usr/bin/Rscript
ags <- commandArgs(TRUE)
input<-ags[1]; out<-ags[2]

library(ShortRead); library(seqinr)
fq <- readFastq(input)

tab <- table(sread(fq))
tab <- tab[order(-tab)]

write.fasta(as.list(rownames(tab)), as.list(as.integer(tab)), out)
