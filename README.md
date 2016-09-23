#piClusterBuster
##A Program for Automated piRNA Cluster Characterization

        piClusterBuster v1.0
                Required:
                        Data input file (Only provide 1 of the 4 options):
                        -fa             - character, input FASTA file containing piRNA cluster sequences
                                                + FASTA file headers must be in the format -> chrom:start..end
                                                        (e.g. X:21632661..21747236)
                        -fq             - character, input FASTQ file (will override input FASTA file)
                        -bam            - character, input BAM file (will override input FASTQ file)
                                                + BAM file must be sorted and indexed
                        -bed            - input BED file (will override fastq or BAM file)

                        Databases:
                        -x              - character, reference genome
                        -gndb           - character, Organism-specific Gene Set
                        -tedb           - character, Organism-specific Transposable Element Set

                Optional:
                        -ncbidb         - character, indicate path to the NCBI nt database or \"NA\" if you would like to have piClusterBuster to retrieve the download
                        -n              - integer, indicate the number of piRNA clusters to analyze
                        -p              - integer, indicate number of processors to use
                        --qsub          - indicate if you would like to submit jobs via qsub
                        -d              - character, output directory
                        -gid            - character, name used for the output of the run (Ex: Aaegypti, human, etc) - spaces in this name are not recommended

                        --go            - option to run GO enrichment enalysis on identified gene fragments
                                                + Species of interest must have GO terms annotated by the Gene Ontology Consortium
                                                + gid must be specified as with the first letter of the genus and species name
                                                        (e.g. "hsapiens", "dmelanogaster", "aaegypti", etc.)
                        --all-srna      - option to include all sRNA, not just piRNAs, in analysis
                        --verbose       - option to print the result of intermediate steps
                        --help          - print the help menu
        Example:
        piClusterBuster -fq myfile.fastq -x reference_genome.fa -gndb Genes.fa -tedb TEs.fa -ncbidb nt -n 5 -p 6 -gid MySpecies --qsub
        
##Description

To be run on Linux... BLAH BLAH FINISH

##Required Files
1. FASTQ, BAM, or BED file
2. Reference Genome
3. Organism-specific Gene Set
	- Gene names must be in Entrez format for gene summary to be completed

###Optional Files
1. Gene-GO Association File
	- Made available by the [Gene Ontology Consortium] (http://geneontology.org/page/download-annotations)

##Required Software
1. Bowtie2
2. RepeatMasker
2. CENSOR
3. NCBI-BLAST 
	- versions 2.2.26 & 2.2.27
4. Bedtools

###Optional Software
5. F-seq
6. Samtools

##Installation

###Apache Ant Java
If you don't have the [Apache Ant Java] (http://ant.apache.org/) library (-bash: ant: command not found), you can download it by:

```
wget http://apache.mirrors.hoobly.com//ant/binaries/apache-ant-1.9.6-bin.tar.gz
tar -zxvf apache-ant-1.9.6-bin.tar.gz
```

####Now load Apache Ant into your path
You can find your path via

```echo $PATH```

then change to a directory in your path. FOR EXAMPLE (only):

```cd /usr/local/bin```

Create a symbolic link to the ant executable

```ln -s apache-ant-1.9.6/bin/ant```

###F-Seq
piClusterBuster uses [F-Seq] (http://fureylab.web.unc.edu/software/fseq/) software to make piRNA cluster calls based on sRNA density at a particular loci.  However, custom piRNA clusters can be analyzed by providing a BED file.

****** This is a modified excerpt from the F-seq README.txt file ******
```
~/piClusterBuster$ cd F-seq/
~/piClusterBuster/F-seq$ ant

This will build F-seq and package it in the dist~ folder. To then run F-seq:
~/piClusterBuster/F-seq$ cd dist~/
~/piClusterBuster/F-seq/dist~$ tar -xvf fseq.tgz
~/piClusterBuster/F-seq/dist~$ cd fseq/bin/
~/piClusterBuster/F-seq/dist~/fseq/bin$ ./fseq
```
Make sure 'bin/fseq' is executable:

```chmod 0755 bin/fseq```

#Load F-seq into your path, as done for ant above
Change to a directory in your path. FOR EXAMPLE (only):

```cd /usr/local/bin```

Create a symbolic link to the fseq executable

```ln -s ~/piClusterBuster/F-seq/dist~/fseq/bin/fseq```

****** END F-seq excerpt ******
