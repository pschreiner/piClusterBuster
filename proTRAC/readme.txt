## Credit to: 
## Rosenkranz D, Zischler H. proTRAC - a software for probabilistic piRNA cluster detection, visualization and analysis. BMC Bioinformatics 2012 13(1):5.
## https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-5
## https://sourceforge.net/projects/protrac/
## http://www.smallrnagroup.uni-mainz.de/software.html

latest version: 2.1

proTRAC 2.0 and later versions do no longer calculate probabilities
based on a hypotetical random distribution of mapped sequence reads.
Instead proTRAC calculates normalized hit counts, which makes it much
faster and thus suitable for really large datasets. Yet it still
considers all the relevant features of mapped sequence reads like
the amount of 1U/10A, size and strandiness to ensure a high specificity.

 - automatic search for a set of transcription factor binding sites.
 - include RepeatMasker annotation (www.repeatmasker.org).
 - include gene annotation (download .GTF files from www.ensembl.org).
 - starts direct from the command line.
 - OS independent (Requires Perl which is commonly preinstalled on Linux
   and Mac machines. You can find a free Perl distribution for Windows at
   strawberryperl.com).å
 - comes with a special piRNA maåpping tool (sRNAmapper.pl) that implies
   biological knowledge to find most probable alignments.


For more information and program options run proTRAC with the following
command:

 perl proTRAC -help


For older versions of proTRAC contact: rosenkranz@uni-mainz.de
#########################################################################
                       QUICK START INFORMATION:

You need...

 - a genome in FASTA format (available at ensembl.org or ncbi.nlm.nih.gov)
 - piRNA sequences in FASTA format. The FASTA headers must refer to the
   number of sequence reads for the sequence e.g:


   >23
   TAACGCGCGTTTCGATCGACTGCTAGTAC
   >122
   TTTCAGTATTATGGCGATTATACCC
   >2
   TGGGCACAAAAAGTCGGATAAAAAAAGGGC


   You can convert your FASTA file to the required format using the collapse
   tool from the NGS TOOLBOX: http://sourceforge.net/projects/ngs-toolbox/


Now that you have got your genome and and a set of piRNAs you must map the
piRNAs to the genome. proTRAC requires a map file that contains the (sorted)
coordinates of mapped piRNAs in the following format (ELAND3 tab-delimited
table):

Chr1  7400394  TTGCTACGTCAGATCGTGCGGGTAA  12  TTGCTACGTCAGATCGTGCGGGTCC  2  +

where the columns refer to chromosome, coordinate, target sequence, number
of reads, query sequence, number of mismatch, strand.


There are two simple ways to produce such a file:
1. Use the small RNA mapping tool sRNAmapper.pl that comes with the proTRAC
   software and start it from the command line with the following command:

   perl sRNAmapper.pl -i piRNAs.fasta -g genome.fas -o piRNAs.map

2. Use the SeqMap software from Jiang and Wong 2008 Bioinformatics, 24(20).
   SeqMap is available at: http://www-personal.umich.edu/~jianghui/seqmap/
   Start SeqMap from the command line with the following command:

   On Windows systems type
   seqmap 0 piRNAs.fasta genome.fas mapfile.eland3 /output_all_matches

   On Linux systems type
   ./seqmap 0 piRNAs.fasta genome.fas mapfile.eland3 /output_all_matches



Now you can start proTRAC from the command line. A typical command would look
like this:

perl proTRAC_2.1.pl -genome hsapiens.fas -map piRNAs.map
