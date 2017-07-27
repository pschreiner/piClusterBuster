#!/bin/bash

outdir=${7}/${8}

cd $outdir
module load ncbi-blast

if [ ! -f  "$1.nsq" ];
then
 makeblastdb -dbtype 'nucl' -in $1 -out $1
fi

blastn -task $3 -db $1 -outfmt 6 -word_size 7 -evalue 1e-3 -query $2 -out $4 -parse_deflines -num_threads $9
Rscript "${10}bin/filterBLASThits.R" $4 $5 $6 $10 
