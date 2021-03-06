Digital_RADs
============
Digital_RADs.py Version 1.03: Feb 2013

The sequencing of restriction site-associated DNA (RAD) markers enables the simultaneous targeting
and genotyping of thousands of homologous loci across a set of individual samples (Baird et al.
2008). By reducing the number of targeted loci in the genome more samples can be run in a single
reaction, therefore striking a balance between harnessing the power of next-gen throughput and
reducing sequencing costs. When designing a RAD-seq experiment it is important to consider the goals
of the study and estimate the number of necessary RAD markers.  Determining how to reach this number
is best achieved by conducting digital digests of the most appropriate reference genome with a panel
of restriction enzymes that cut at different frequencies (Davey et al. 2011). Additional reduction
can be achieved by using two enzymes, thus targeting only regions where the enzymes cut adjacent to
each other within a specified size window (Hohenlohe et al. 2012, Peterson et al. 2012, DaCosta and
Sorenson 2014).

This python script conducts digital digests of a reference genome (or any fasta file) with either
one or two enzymes, and can be used to optimize the experimental design of a RAD-seq study. The
script was designed and tested in the python 3 environment, but should also run in python 2.6 and
2.7.


RUNNING THE SCRIPT


The Digital_RADs.py script must be in the same directory as the fasta input file. Note that the
script will parse the fasta file into multiple files (see below), so you may want to have a
dedicated directory. Open a terminal window and navigate to this directory. Then run the command
line for either 1 or 2 enzymes (see below). The output file(s) will be written to the same
directory.


CONVERSION OF FASTA FILE


The first step of the script is to parse the fasta file into individual files for each
chromosome/contig/sequence, each of which will be a single line with a long string of the sequence.
The names of these files are taken from the names of the sequences in the fasta file (after the >
character). If you are running the script multiple times (e.g. testing different restriction
enzymes) it is recommended that you keep these files so future runs of the script will finish
faster.


COMMAND LINE WITH 1 ENZYME


python3 Digital_RADs.py infile outfile 1 enzyme_seq bases

infile		Name of infile. The infile must be in fasta format. There is no limit on the number of
			separate sequences (e.g. chromosomes, contigs, scaffolds) in the file.

outfile		Name of outfile. The outfile is a tab-delimited file that contains the following data
			on the RAD markers recovered:

			contig: name of sequence/chromosome/contig
			position: start position
			downstream: sequence downstream of enzyme, length determined by 'bases' parameter
			down_GC: GC content of the upstream sequence
			upstream: reversecomp sequence upstream of enzyme, length determined by 'bases' parameter

1			Number of enzymes

enzyme_seq	The sequence of the restriction enzyme (e.g. GAATTC for EcoRI)

bases		The length of the down/upstream sequence (including enzyme) reported in the outfile


COMMAND LINE WITH 2 ENZYMES


python3 Digital_RADs.py infile outfile 2 enzyme1_seq enzyme2_seq window_start window_end

infile		The infile must be in fasta format. There is no limit on the number of separate
			sequences (e.g. chromosomes, contigs, scaffolds) in the file.

outfile		The outfile is a tab-delimited file that contains the following data on the RAD
			markers recovered:

			contig: name of sequence/chromosome/contig
			position: start position of enzyme1 end of RAD marker
			length: length of RAD marker
			direction: direction of locus
						1 = enzyme2 downstream of enzyme1
						-1 = enzyme2 upstream of enzyme1
			sequence: sequence of RAD marker
			GC: GC content of the RAD marker

2			Number of enzymes

enzyme1_seq	The sequence of the primary restriction enzyme (e.g. GAATTC for EcoRI)

enzyme2_seq	The sequence of the secondary restriction enzyme (e.g. GAATTC for EcoRI)

window_start	The lower threshold for the size selection window*

window_end		The upper threshold for the size selection window*

*This size selection step is used in a double-digest protocol to further reduce the number of
targeted loci.  Note that the size window is for only the insert fragment (i.e., insert minus the
adapters). For example, if your protocol size selects for ligated fragments of 300-400 bp and your
adapters comprise 125 bp then you should input 175 and 275 for your window_start and window_end
values, respecitvely.


REFERENCES

Baird NA, Etter PD, Atwood TS, Currey MC, Shiver AL, Lewis ZA, Selker EU, Cresko WA, Johnson EA
(2008) Rapid SNP discovery and genetic mapping using sequenced RAD markers. PLoS One. 3, e3376.

DaCosta JM, Sorenson MD (2014). Amplification biases and consistent recovery of loci in a
double-digest RAD-seq protocol. PLoS One. 9, e106713.

Davey JW, Hohenlohe PA, Etter PD, Boone JQ, Catchen JM, Blaxter ML (2011) Genome-wide genetic marker
discovery and genotyping using next-generation sequencing. Nature Reviews Genetics. 12, 499-510.

Hohenlohe PA, Bassham S, Currey M, Cresko WA (2012) Extensive linkage disequilibrium and parallel
adaptive divergence across threespine stickleback genomes. Phil Trans Roy Soc B. 367, 395-408.

Peterson BK, Weber JN, Kay EH, Fisher HS, Hoekstra HE (2012). Double digest RADseq: an inexpensive
method for do novo SNP discovery and genotyping in model and non-model species. PLoS One. 7, e37135.

