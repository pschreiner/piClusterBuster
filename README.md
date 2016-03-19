# piClusterBuster
A Program for Automated piRNA Cluster Characterization

	piClusterBuster v1.0
		Arguments:
		-fq		- character, input fastq file
		-bam		- character, input BAM file (will override input fastq file)
		-bed		- input BED file (will override fastq or BAM file)
		-x		- character, reference genome
		-n		- integer, indicate the number of piRNA clusters to analyze
		-ncbidb		- character, indicate path to the NCBI nt database or \"NA\" if you would like to have piClusterBuster to retrieve the download
		-gndb		- character, Organism-specific Gene Set
		-tedb		- character, Organism-specific Transposable Element Set
		-p		- integer, indicate number of processors to use
		--qsub		- indicate if you would like to submit jobs via qsub
		-d		- character, output directory
		-gid		- character, name used for the output of the run (Ex: Aaegypti, human, etc) - spaces in this name are not recommended
		--verbose	- option to print the result of intermediate steps
		--help		- print the help menu
	Example:
	piClusterBuster -fq myfile.fastq -x reference_genome.fa -gndb Genes.fa -tedb TEs.fa -ncbidb nt -n 5 -p 6 -gid MySpecies --qsub



***** Continue Full Instructions *****

If you don't have the Apache Ant Java library (-bash: ant: command not found), you can download it by:
	wget http://apache.mirrors.hoobly.com//ant/binaries/apache-ant-1.9.6-bin.tar.gz
	tar -zxvf apache-ant-1.9.6-bin.tar.gz

Now load Apache Ant into your path
You can find your path via
	echo $PATH
then change to a directory in your path. FOR EXAMPLE (only):
	cd /usr/local/bin
Create a symbolic link to the ant executable
	ln -s apache-ant-1.9.6/bin/ant

****** This is a modified excerpt from the F-seq README.txt file ******

~/piClusterBuster$ cd F-seq/
~/piClusterBuster/F-seq$ ant

This will build F-seq and package it in the dist~ folder. To then run F-seq:
~/piClusterBuster/F-seq$ cd dist~/
~/piClusterBuster/F-seq/dist~$ tar -xvf fseq.tgz
~/piClusterBuster/F-seq/dist~$ cd fseq/bin/
~/piClusterBuster/F-seq/dist~/fseq/bin$ ./fseq

#Make sure 'bin/fseq' is executable:
chmod 0755 bin/fseq

# Load F-seq into your path, as done for ant above
# Change to a directory in your path. FOR EXAMPLE (only):
cd /usr/local/bin
# Create a symbolic link to the fseq executable
ln -s ~/piClusterBuster/F-seq/dist~/fseq/bin/fseq

****** END F-seq excerpt ******
