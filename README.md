[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=mapseq)


# MAPseq v1.2.6 (24 Mar 2020)
by Joao F. Matias Rodrigues, Thomas S.B. Schmidt, Janko Tackmann, and Christian von Mering  
Institute of Molecular Life Sciences, University of Zurich, Switzerland

Reference:
Matias Rodrigues JF, Schmidt TSB, Tackmann J & von Mering C (2017) MAPseq: highly efficient k-mer search with confidence estimates, for rRNA sequence analysis. Bioinformatics. http://doi.org/10.1093/bioinformatics/btx517

---
## Table of contents

1. Installation  
2. MAPseq usage instructions
 a. Default reference  
 b. Custom user-provided reference  
 c. Single sample counts summar  
 d. OTU count table for multiple samples  
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

The binary packages can be found in the "Releases" page on GitHub:
https://github.com/jfmrod/MAPseq/releases


### i) Installing the binary package

To install the binary package simply unpack the contents of the mapseq tar.gz file, e.g.:

tar -xvzf mapseq-1.2-linux.tar.gz   # for the linux version  
or  
tar -xvzf mapseq-1.2-macosx.tar.gz  # for the MacOSX version  

The mapseq binary will be located in the created directory. You may move the whole directory to another location. Moving only the binary elsewhere will break the installation though, as the data files are searched for in relation to the binary's path.


### ii) Installing from source

To compile from the github source you will need:
- svn
- autotools/autoconf
- wget
- git
- libncurses5-dev
- libtool

On Ubuntu systems you can install these with the command:

sudo apt-get install build-essential wget subversion git libncurses5-dev libtool


You can then clone the mapseq repository with:

git clone https://github.com/jfmrod/MAPseq.git


Once you have cloned the repository you can type the following commands in the MAPseq directory


./setup.sh      # downloads eutils and the reference data files  
./bootstrap     # this step is only needed if you cloned the repository, you will also need to install autotools/autoconf packages  
./configure  
make  
make install  


If you want to install MAPseq to your home directory instead of the default system wide /usr/local/ directory,
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


### c) Single sample counts summary

mapseq -otucounts <sample1.mseq>
 
Provides you with summary counts of the MSEQ output files for each taxonomy and level. Example:

#sample.mseq	102301  
Taxonomy	TaxonomyLevel	Label	Counts  
0	0	Bacteria	102301  
0	1	Bacteria;Bacteroidetes	7586  
0	1	Bacteria;Firmicutes	2150  
0	1	Bacteria;Proteobacteria	112  
0	1	Bacteria;PHY_Coriobacteriia	7  
0	1	Bacteria;Actinobacteria	4  
0	1	Bacteria;Fusobacteria	1  
0	2	Bacteria;Bacteroidetes;Bacteroidia	7585  
0	2	Bacteria;Firmicutes;Clostridia	1330  
  


### d) OTU count table for multiple samples

mapseq -otutable <sample1.mseq> <sample2.mseq> ...  

Generates a tab separated value (tsv) file with the counts for each sample (column wise) and OTU or taxonomic labels (row wise).
Which taxomy (OTU or NCBI taxonomy) and levels in the taxonomy can be specified using the -ti and -tl, respectively.  


The generated table can be imported into R with the following R command:

myotutable <- read.table("map.otutable",sep="\t",header=TRUE)



## 3. MSEQ FILE OUTPUT

In the results output, each line indicates a classification of the read. Two output formats can be chosen ("simple" or "confidences") using the --outfmt option.

The "confidences" format outputs the confidence values for each of the taxonomic levels. For example:
SRR044946.347	GQ156763:1..1446	548	0.91985428	505	22	22	1	540	263	800	0.99072355	20		Bacteria	1	1	Firmicutes	0.55452305	1	Clostridia	0.55452305	1	Clostridiales	0.55452305	1	Ruminococcaceae	0.31190208	0.3119020760059357	Ruminococcus	0	0.2104288786649704	Ruminococcus gnavus	0	0.0604640431702137		Bacteria	0.58272612	1	F6159	0.22964814	1	G35588	0	1	S61033	0	0.7381679934055649	SS52094	0	0.2980887881680916

Each field is tab separated and indicates the following:  
Field  
1	Query sequence id  
2	Reference sequence id (highest alignment score)  
3	Alignment score  
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
1.2.6 (24 Mar 2020)
- Fixed "-otutable" option using only counts from first sample.

1.2.5 (12 Jul 2019)
- Added "-otucounts" and "-otutable" options to generate count summary for single mapseq (.mseq) files or an otu/taxa table for multiple .mseq files

1.2.4 (27 May 2019)
- Added "-ignoreEmptyTax" option. Default is off until thorough benchmarks are performed.
  Prevents second hits with missing taxonomic labels (uncertain annotation) from decreasing the confidence of the top hit assignment.

1.2.3 (2 Oct 2018)
- Fixed missing newline causing last sequence to be missed, added assert on empty sequences
- Fixed double hits reported when classifying long queries (>1200bp)

- Updated to mapref-2.2b:
  Removed low quality reference sequences that would cause issues when classifying low quality query sequences.

1.2.2 (30 Oct 2017)
- Fixed multithreaded race condition causing issues on some systems.

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

