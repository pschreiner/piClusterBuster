#!/bin/bash
if [ "$#" != 0  ]; then printf "
$@
"
fi

echo "
	piClusterBuster v1.0 - Arguments
		Required:
			Data input file (Only provide 1 of the 4 options):
			-fa		- character, input FASTA file containing piRNA cluster sequences
						+ FASTA file headers must be in the format -> chrom:start..end
							(e.g. X:21632661..21747236)
			-fq		- character, input FASTQ file (will override input FASTA file)
			-bed		- input BED file (will override fastq or BAM file)

			Databases:
			-x		- character, reference genome
			-gndb		- character, Organism-specific Gene Set
			-tedb		- character, Organism-specific Transposable Element Set

		Optional:
			-ncbidb		- character, indicate path to the NCBI nt database or \"NA\" if you would like to have piClusterBuster to retrieve the download
			-n		- integer, indicate the number of piRNA clusters to analyze
			-p		- integer, indicate number of processors to use
			--qsub		- indicate if you would like to submit jobs using Torque/Maui Resource Management (qsub)
			--srun		- indicate if you would like to submit jobs using Slurm Resource Management (srun)
			-d		- character, output directory
			-gid		- character, name used for the output of the run (Ex: Aaegypti, human, etc) - spaces in this name are not recommended
		
			--go		- option to run GO enrichment enalysis on identified gene fragments
						+ Species of interest must have GO terms annotated by the Gene Ontology Consortium
						+ gid must be specified as with the first letter of the genus and species name
							(e.g. "hsapiens", "dmelanogaster", "aaegypti", etc.)
			--all-srna	- option to include all sRNA, not just piRNAs, in analysis
			--verbose	- option to print the result of intermediate steps
			--help		- print the help menu
	Example:
	piClusterBuster -fq myfile.fastq -x reference_genome.fa -gndb Genes.fa -tedb TEs.fa -ncbidb nt -n 5 -p 6 -gid MySpecies --qsub
"
