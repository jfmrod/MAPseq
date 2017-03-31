#!/bin/bash

./configure --prefix=`pwd`/mapseq-1.0-bin --bindir=`pwd`/mapseq-1.0-bin --enable-static LIBS="-lgpm -ltermcap"
#make install
#cp README mapseq-1.0-bin/
#rm -rf mapseq-1.0-bin/{esh,eutils-config,include,lib} 
#tar chof - mapseq-1.0-bin | GZIP=--best gzip -c >mapseq-1.0-bin.tar.gz
