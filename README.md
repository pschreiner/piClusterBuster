#piClusterBuster
##Description
piClusterBusteR is a series of R and bash scripts that interact together, along with other standalone programs, to perform piRNA cluster characterization and annotation.  piClusterBusteR is a software tool that accurately, automatically, and efficiently describes the contents of piRNA clusters in any biological system that utilizes the piRNA pathway.  piClusterBusteR performs nested annotation of piRNA cluster contents to ensure high-quality characterization, provide a quantitative representation of piRNA cluster composition by feature, and make available annotated and unannotated piRNA cluster sequence that can be utilized for downstream analysis. The data necessary to run piClusterBusteR and the skills necessary to execute this software on any species of interest are not overly burdensome for biological researchers.
piClusterBusteR has been utilized to compare the composition of top piRNA generating loci amongst Metazoan species.  Characterization and quantification of cluster composition allows for comparison within piRNA clusters of the same species and between piRNA clusters of different species. 

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
        
##Required Files
1. Data file (FASTQ, BAM, FASTA, or BED file)
2. Reference Genome
3. Organism-specific Gene Set
	- Gene names must be in Entrez format for gene summary to be completed
4. Transposable Element Set

###Optional Files
1. NCBI nucleotide database (ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz)
2. Gene-GO Association File
	- Made available by the [Gene Ontology Consortium] (http://geneontology.org/page/download-annotations)

##Required Software
1. BLAST+ (versions 2.2.7 or higher)
2. CENSOR
3. R software
4. RepeatMasker
5. SAMtools

###Optional Software
6. Bowtie2 (necessary with FASTQ file input)
7. proTRAC (necessary if piRNA cluster definitions aren't already made)

## R Packages Utilized
1. Biostrings
2. doMC
3. GenomicRanges
4. gProfileR
5. Plyr
6. qcc
7. seqinr
8. systemPipeR

##Installation
1. Download zip file or git clone
2. From the main folder of piClusterBusteR, run "./piClusterBusteR" for the usage statement
