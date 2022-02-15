#!/bin/bash

VERSION=2.0.1alpha

./configure --prefix=`pwd`/mapseq-$VERSION-linux --bindir=`pwd`/mapseq-$VERSION-linux --enable-makestatic LIBS="-lgpm -ltermcap" CXXFLAGS="-O2"
#make install
#cp README mapseq-1.0-bin/
#rm -rf mapseq-1.0-bin/{esh,eutils-config,include,lib} 
#tar chof - mapseq-1.0-bin | GZIP=--best gzip -c >mapseq-1.0-bin.tar.gz
