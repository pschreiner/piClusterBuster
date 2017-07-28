piClusterBusteR
===============
Description
-----------
piClusterBusteR is a series of R and bash scripts that interact together, along with other standalone programs, to perform piRNA cluster characterization and annotation.  piClusterBusteR is a software tool that accurately, automatically, and efficiently describes the contents of piRNA clusters in any biological system that utilizes the piRNA pathway.  piClusterBusteR performs nested annotation of piRNA cluster contents to ensure high-quality characterization, provide a quantitative representation of piRNA cluster composition by feature, and make available annotated and unannotated piRNA cluster sequence that can be utilized for downstream analysis. The data necessary to run piClusterBusteR and the skills necessary to execute this software on any species of interest are not overly burdensome for biological researchers.

piClusterBusteR has been utilized to compare the composition of top piRNA generating loci amongst Metazoan species.  Characterization and quantification of cluster composition allows for comparison within piRNA clusters of the same species and between piRNA clusters of different species. 

:Software for Automated Classification and Characterization of piRNA Cluster Loci:
 --------------------------------------------------------------------------------

        piClusterBuster v1.0
                Required:
                        Data input file (Only provide 1 of the 4 options):
                        -fa             - character, input FASTA file containing piRNA cluster sequences
                                                + FASTA file headers must be in the format -> chrom:start..end
                                                        (e.g. X:21632661..21747236)
                        -fq             - character, input FASTQ file (will override input FASTA file)
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
        
Required Files
--------------
1. Data file (FASTQ, BAM, FASTA, or BED file)
2. Reference Genome
3. Organism-specific Gene Set
	- Gene names must be in Entrez format for gene summary to be completed
4. Transposable Element Set

Optional Files:
1. NCBI nucleotide database (ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz)
2. Gene-GO Association File
	- Made available by the [Gene Ontology Consortium](http://geneontology.org/page/download-annotations)

Required Software
-----------------
1. [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
2. [CENSOR](http://www.girinst.org/downloads/software/censor/)
3. [R software](https://www.r-project.org/)
4. [RepeatMasker](http://www.repeatmasker.org/RMDownload.html)

Optional Software:

5. [proTRAC](http://www.smallrnagroup.uni-mainz.de/software.html)
	- necessary if piRNA cluster definitions aren't already made

R Packages Utilized
-------------------
1. [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
2. [doMC](https://cran.r-project.org/web/packages/doMC/index.html)
3. [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
4. [gProfileR](https://cran.r-project.org/web/packages/gProfileR/gProfileR.pdf)
5. [plyr](https://cran.r-project.org/web/packages/plyr/plyr.pdf)
6. [qcc](https://cran.r-project.org/web/packages/qcc/qcc.pdf)
7. [seqinr](http://seqinr.r-forge.r-project.org/)
8. [systemPipeR](https://bioconductor.org/packages/release/bioc/html/systemPipeR.html)

Installation
------------
1. Download zip file or git clone
2. From the main folder of piClusterBusteR, run "./piClusterBusteR" for the usage statement
