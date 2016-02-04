#!/bin/bash
#PBS -l walltime=240:00:00,nodes=1:ppn=8 -j oe -N annotate_piCluster

cd ~/scripts/piClusterBuster/
Rscript annotate_piCluster.R $@
