Documentation

--------------------------
    Author Information
--------------------------
Patrick Schreiner
pschr001@ucr.edu
University of California, Riverside

---------------------
    Accessibility
---------------------
To download and execute piClusterBuster from GitHub:

git clone git clone https://github.com/piClusterBuster.git
cd piClusterBuster
./piClusterBuster	# Prints usage statement

------------------
    Parameters
------------------
Accepted input data file formats:
1. FASTQ
2. FASTA
3. BED (http://genome.ucsc.edu/FAQ/FAQformat#format1)

Required:
1. -fq, -fa, OR -bed (input data file)
2. -x (reference genome in FASTA format)
3. -tedb (transposable element database in FASTA format)
4. -gndb (organism-specific gene set in FASTA format)
5. -ncbidb (NCBI nt databasein FASTA format: ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz)
6. -n (number of piRNA clusters to analyze)

Optional:
1. -gid (provide a name for the output directory/ guide GO term analysis)
2. --all-srna (observe all sRNAs, not just piRNAs)
3. --verbose (retain intermediate files)

-------------------------------
    Performance Enhancement
-------------------------------
Optional Parameters are available for multithreading and parallel computing:
1. --qsub (job submission)
2. -p (multithreading)

----------------------------
    Interpreting Results
----------------------------
Upon completion, a new directory exists with your Genome Identification (gid) followed by "_results".  If you did not provide a gid, "genome" is used by default.

In this directory, there is the BED file used to define the piRNA clusters in the analysis, 3 genome-level summary files, and additional directories for each individual piRNA cluster under observation.

1. BED file: Contains the genomic coordinates (and proTRAC's normalized read count, if available) for each individual piRNA cluster

2. Genome-level Summary Files:
	a. piClusters-GenomeSummary.pdf:	
		Contains plots that describe the size of the piRNA clusters, percent of content attributed to TE, Genic, Other, or undefined origin, degree of sense/antisense orientation of the feature calls, average nucleotide occupancy of features across all piRNA clusters observed, and a comparison to the average feature occupancy of piRNA clusters in other organisms.
	b. piClusters-GenomeSummary.xls: 
		Data in text format that was utilized to create the plots in piClusters-GenomeSummary.pdf.
	c. "gid"-AggregateSummary.xls: 
		Data in text format describing the number of piRNA clusters analyzed, the nucleotide occupancy occupied by those piRNA clusters, and the percent of the genome that was occupied by the piRNA clusters under observation.

3. Picluster-level directories: 
		Each piRNA cluster under observation has an individual directory containing the details of its characterization.
	a. piclusterXX-Summary.pdf: 
		Contains plots that describe the number of features calls in the piRNA cluster, nucleotide occupancy of the feature calls, degree of sense/antisense orientation of the feature calls, the number of calls associated with a TE superfamily, and the number of specific TE calls for the top 5 TEs in the top 5 TE superfamilies that were characterized within the piRNA cluster.
	b. piclusterXX-Summary.xls: 
		Contains the data that was used the generate the pdf
	c. piclusterXX_FEATURElandscape-ALL.xls: 
		Contains the filtered data for each feature that was characterized in each piRNA cluster.
	d. piclusterXX-FinalAnnotation.xls: 
		File containing the final, filtered description of piRNA cluster contents
	e. piclusterXX-UnannotatedSeqs.fa: 
		FASTA file containing sequence within the piRNA cluster from which the origin could not be identified
	f. piclusterXX_Otherlandscape-FILTERED.BLAST: 
		File containing the filtered BLAST output, to account for redundant characterization, within the NCBI nt database
	g. piclusterXX-GnTE.xls: 
		All of the filtered Gene and TE calls within the piRNA cluster
	h. piclusterXX_FEATURElandscape-FILTERED.CENSOR: 
		File containing the filtered CENSOR output, to account for redundant characterization, for the described feature
	i. piclusterXX_FEATURElandscape-FILTERED.RM: 
		File containing the filtered RepeatMasker output, to account for redundant characterization, for the described feature
	j. piclusterXX_FEATURElandscape.CENSOR: 
		File containing the unfiltered CENSOR output for the described feature
	k. piclusterXX_FEATURElandscape.RM: 
		File containing the unfiltered RepeatMasker output for the described feature

* All of the above files are in the standard output format for RepeatMasker, CENSOR, and BLAST

---------------
    License
---------------
piClusterBusteR is free software that is licensed by the GNU General Public License 3.0 license (https://www.gnu.org/licenses/gpl-3.0.txt).  Redistribution and/or modification of piClusterBusteR is permitted under the terms of the GNU General Public License.

This software comes with absolutely no warranty.
