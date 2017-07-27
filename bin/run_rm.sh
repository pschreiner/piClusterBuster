#!/bin/bash
outdir=${3}/${4}

cd $outdir
module load RepeatMasker

RepeatMasker -pa $5 -lib $2 -dir . $1
cp "${4}.fa.out" $6

Rscript "${9}bin/filterRMhits.R" $6 $7 $8 $9
