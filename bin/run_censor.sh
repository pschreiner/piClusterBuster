#!/bin/bash

outdir=${7}/${8}
cd $outdir

module load censor
module load ncbi-blast/2.2.26

censor.ncbi $1 -lib $2 -bprg $3 -map $4
if [ -f $4 ];
	then Rscript "${9}bin/filterCENSORhits.R" $4 $5 $6 $9
	else touch $5
fi
