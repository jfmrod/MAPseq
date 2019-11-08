#!/bin/bash

MSEQVER=1.2.5

./configure --prefix=`pwd`/mapseq-$MSEQVER-linux --bindir=`pwd`/mapseq-$MSEQVER-linux --enable-makestatic LIBS="-lgpm -ltermcap" CXXFLAGS="-O2"
make install

rm -rf mapseq-$MSEQVER-linux/{esh,eutils-config,include,lib}
rm -rf mapseq-$MSEQVER-linux/share/{aclocal,eutils,man}

tar -cvzf mapseq-$MSEQVER-linux.tar.gz mapseq-$MSEQVER-linux

#make install
#cp README mapseq-1.0-bin/
#rm -rf mapseq-1.0-bin/{esh,eutils-config,include,lib} 
#tar chof - mapseq-1.0-bin | GZIP=--best gzip -c >mapseq-1.0-bin.tar.gz
