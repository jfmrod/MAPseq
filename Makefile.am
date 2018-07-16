ACLOCAL_AMFLAGS=-I m4

SUBDIRS=$(subdirs)
DIST_SUBDIRS=$(subdirs)

AM_CXXFLAGS=-DDATAPATH='"$(datadir)/mapseq"'

if !USE_SYSEUTILS
  LDADD=eutils/libeutils.la
  AM_CXXFLAGS+= -I$(srcdir)/eutils
#  AM_CXXFLAGS+= -Ilibs
endif

#TESTS=tests/simple.sh

#mapref_files=data/mapref.fna.ltps119tax data/mapref.fna data/mapref.fna.mscluster data/mapref.fna.ncbitax data/mapref.fna.otutax # data/mapref.gold.fasta data/mapref.gold.fasta.tax
mapref_files=data/mapref-2.2.fna data/mapref-2.2.fna.mscluster data/mapref-2.2.fna.ncbitax data/mapref-2.2.fna.otutax # data/mapref.gold.fasta data/mapref.gold.fasta.tax

#noinst_LIBRARIES=libmapseqali.a
#libmapseqali_a_SOURCES=eseqali.h eseqali.cpp eseq.h eseq.cpp

dist_man_MANS=man/mapseq.1


bin_PROGRAMS=mapseq
#noinst_bin_PROGRAMS=kmerseqtables
#bin_SCRIPTS=tree-convert-json.esh

#LIBS=`eutils-config --libs || eutils/eutils-config --libs`
EXTRA_PROGRAMS=kmerseqtables
EXTRA_DIST=manual/mapseq-manual.tex manual/mapseq-manual.pdf data/make-transtable.sh data/prot_trans_table.h data/transl_table $(mapref_files)


maprefdir = $(datadir)/mapseq
mapref_DATA=$(mapref_files)

# data/mapseqref.fa data/mapseqref.ncbitax data/mapseqref.otutax 

mapseq_SOURCES=mapseq.cpp eseqali.h ekmerhashmap.h ekmerhashmap_dec.h kmerseqtables-data.h kmerseqtables-data.cpp eseqali.h eseqali.cpp eseq.h eseq.cpp
#mapseq_LDADD=libmapseqali.a eutils/libeutils.la
#mapseq_LDADD=eutils/libeutils.la
#mapseq_DEPENDENCIES=kmerseqtables-data.h

kmerseqtables_SOURCES=kmerseqtables.cpp
#kmerseqtables_LDADD=eutils/libeutils.la
#kmerseqtables_DEPENDENCIES=eutils/libeutils.la

# needed to force make to make this file before making kmerfilter
#BUILT_SOURCES=kmerseqtables-data.h

kmerseqtables-data.h : kmerseqtables$(EXE)
	./kmerseqtables

kmerseqtables-data.cpp : kmerseqtables$(EXE)
	./kmerseqtables

