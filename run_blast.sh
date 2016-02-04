#!/bin/bash
#PBS -l nodes=1:ppn=8 -j oe

cd ~/scripts/piClusterBuster/
module load ncbi-blast

if [ ! -f  "$1.1.bt2" ];
then
 bowtie2-build $1 $1
fi

blastn -db $1 -outfmt 6 -evalue 1e-3 -query $2 -out $3
out2=${3/.BLAST/-FILTERED.BLAST} ## Manipulate original BLAST out filename to indicate it has been filtered
Rscript filterBLASThits.R $3 $out2 pwd
