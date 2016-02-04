# piClusterBuster
A Program for Automated piRNA Cluster Characterization

piClusterBuster v1.0

Arguments:
  -i              - input fastq file
  -b              - input SAM/BAM file (will override input fastq file)
  -x              - reference genome
  -n              - indicate the number of piRNA clusters to analyze
  -ncbidb         - indicate path to the NCBI nt database or "NA" if you would like to have piClusterBuster to retrieve the download
  -gn             - Organism-specific Gene Set
  -te             - Organism-specific Transposable Element Set
  -p              - indicate number of processors to use
  --qsub          - option to submit jobs via qsub
  -d              - output directory
  -gid            - name used for the output of the run (Ex: Aaegypti, human, etc) - spaces in this name are not recommended
  --help          - print the help menu
  
Example:
  piClusterBuster -i myfile.fastq -x reference_genome.fa -n 15 -ncbidb ./nr -p 6 -gid MySpecies



***** Continue Full Instructions *****
