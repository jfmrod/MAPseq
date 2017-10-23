
# MAPseq v1.2 (23 Oct 2017)
by Joao F. Matias Rodrigues, Thomas S.B. Schmidt, Janko Tackmann, and Christian von Mering  
Institute of Molecular Life Sciences, University of Zurich, Switzerland

Reference:
Matias Rodrigues JF, Schmidt TSB, Tackmann J, Mering C von. MAPseq: bringing speed, accuracy and consistency to metagenomic ribosomal RNA analysis. submitted.

---
## Table of contents

1. Installation  
2. MAPseq usage instructions  
 a. default reference  
 b. custom user-provided reference  
3. File output  
4. History  

---

MAPseq is a set of fast and accurate sequence read classification tools designed to assign taxonomy and
OTU classifications to ribosomal RNA sequences. This is done by using a
reference set of full-length ribosomal RNA sequences for which known taxonomies are known,
and for which a set of high quality OTU clusters has been previously generated.
For each read, the best guess and correspoding confidence in the assignment is shown at
each taxonomic and OTU level.

For bugs and more information contact: Joao F. Matias Rodrigues <jfmrod@gmail.com>


## 1. INSTALLATION

You can get the source code on github or binary packages at:

git clone https://github.com/jfmrod/MAPseq.git

http://meringlab.org/software/mapseq/


### i) Installing the binary package

To install the binary package simply unpack the contents of the mapseq tar.gz file, i.e.:

tar -xvzf mapseq-1.2-linux.tar.gz   # for the linux version  
or  
tar -xvzf mapseq-1.2-macosx.tar.gz  # for the MacOSX version  

The mapseq binary will be located in the created directory. You may move the whole directory to another location. Moving the binary elsewhere will break the installation though, as the data files are searched for in relation to the binary's path.

### ii) Installing from source

First make sure you have the gcc compiler / build tools installed. On Ubuntu you would need to run "apt-get install build-essential" with admin priviledges. You will also need the libncurses dev libraries installed: "apt-get install libncurses5-dev"

./setup.sh      # downloads eutils and the reference data files  
./bootstrap     # this step is only needed if you cloned the repository, you will also need to install autotools/autoconf packages  
./configure  
make  
make install  

In the directory where you unpacked the package contents.
Alternatively, if you want the program to be installed to another
location instead of the default system wide /usr/local/ directory,
you can change the ./configure command to:

./configure --prefix=$HOME/usr

This would install the program binaries to a directory usr/bin inside
your home directory (i.e.: $HOME/usr/bin/mapseq), after you type
the command "make install".


## 2. MAPseq usage instructions

### a) Default reference

MAPseq takes as input a fasta file with raw sequence data which should have been previously
demultiplexed and quality filtered usually from a fastq file produced by the sequencing machine.

If the input sequences can be found in the file "rawseqs.fa". Then to classify the reads
one simply has to run the following command:

mapseq rawseqs.fa > rawseqs.fa.mseq

This will classify all the sequences found in rawseqs against the standard reference dataset
provided with MAPseq.

You can change the number of threads that MAPseq uses with the -nthreads <no_threads> argument.

### b) Custom user-provided reference

You can use mapseq with your own fasta reference and taxonomy files with the following command:

mapseq rawseqs.fa <customref.fasta> <customref.tax> [customref.tax2 ...] > rawseqs.fa.mseq

Where customref.fasta is a nucleotide fasta file with your reference set and customref.tax, customref.tax2 are one or more taxonomic assignments for each sequence in the reference.

The taxonomy file should have a header (preceeded with the # character) with the identity cutoff parameters and description of the taxonomy followed by two tab-separated columns composed of the accession id and the taxonomy. For example:

#cutoff: 0.00:0.08 0.70:0.35 0.70:0.35 0.70:0.35 0.80:0.25 0.92:0.08 0.95:0.05  
#name: NCBI  
#levels: Kingdom Phylum Class Order Family Genus Species  
HE801216:78..1345       Bacteria;Proteobacteria;Gammaproteobacteria;Methylococcales;Methylococcaceae;Methylomonas;Methylomonas paludis  
HE802067:76740..77993   Bacteria;Actinobacteria;Actinobacteria;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium glutamicum  
HE804045:1012175..1013425       Bacteria;Actinobacteria;Actinobacteria;Pseudonocardiales;Pseudonocardiaceae;Saccharothrix;Saccharothrix espanaensis  


## 3. FILE OUTPUT

In the results output, each line indicates a classification of the read. Two output formats can be chosen ("simple" or "confidences") using the --outfmt option.

The "confidences" format outputs the confidence values for each of the taxonomic levels. For example:
SRR044946.347	GQ156763:1..1446	548	0.91985428	505	22	22	1	540	263	800	0.99072355	20		Bacteria	1	1	Firmicutes	0.55452305	1	Clostridia	0.55452305	1	Clostridiales	0.55452305	1	Ruminococcaceae	0.31190208	0.3119020760059357	Ruminococcus	0	0.2104288786649704	Ruminococcus gnavus	0	0.0604640431702137		Bacteria	0.58272612	1	F6159	0.22964814	1	G35588	0	1	S61033	0	0.7381679934055649	SS52094	0	0.2980887881680916

Each field is tab separated and indicates the following:  
Field  
1	Query sequence id  
2	Reference sequence id (highest alignment score)  
3	Alignment bitscore  
4	Pairwise identity  
5	Matches  
6	Mismatches  
7	Gaps  
8	Query start pos  
9	Query end pos  
10	Reference start pos  
11	Reference end pos  
12      Strand (+/-)  
13	[empty]  

After the first empty field the taxonomy classifications and confidences are shown, every taxonomy classification is separated by an empty field.
Although different fasta reference and taxonomy databases can be specified by the user, by default mapseq maps reads to the NCBI taxonomy and to OTU taxonomies

The combined confidence is computed based on a score confidence, used to control misclassification errors, and a identity cutoff confidence, used to ensure that the query isnt misclassified due to the inexistence of a sequence representative in the database of the true classification. The score confidence is
calculated by comparing the identity of the assigned taxonomy to the identity of the first sequence not matching the assigned taxonomy.
The identity cutoff confidence uses preoptimized cutoffs at each taxonomic level to calculate the confidence that the query is not too divergent from the assigned taxonomy.

We recommend using a combined confidence cutoff of 0.4, or 0.5 as this value yielded the highest F1/2-score for MAPseq in our benchmarks. Please see our article for further information.

The "simple" format gives the alignment information plus the taxonomy assignment for which the combined confidence at least 0.5. For example:  
query1	FJ560320:1..876	301	0.7369985	301	0	0	0	301	305	606	+		Archaea		Archaea;F94;G275  



## 4. HISTORY
1.2.1 (23 Oct 2017)
- Updated mapref to v2.2. Fixed several issues with v2.0.
- Dropped LTP taxonomy due to low coverage.


1.2 (16 July 2017)  
- Updated mapref to v2.0, now includes 1.5 million sequences.  
- Added assert checks.  

1.1 (24 April 2017)  
- Several improvements and bug fixes, updated to latest NCBI taxonomy.  

1.0 (14 October 2016)  
- First release of MAPseq.  


