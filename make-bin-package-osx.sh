#!/bin/bash

unset LD_LIBRARY_PATH
unset DYLD_LIBRARY_PATH
unset LIBRARY_PATH
unset CPLUS_INCLUDE_PATH
unset C_INCLUDE_PATH
unset OBJC_INCLUDE_PATH
unset OBJCPLUS_INCLUDE_PATH

MSEQVER=2.1

./configure --enable-makestatic --with-syseutils --prefix=`pwd`/mapseq-$MSEQVER-macosx --bindir=`pwd`/mapseq-$MSEQVER-macosx CXXFLAGS="-O2 `~/usr/static/bin/eutils-config --cxxflags`" LDFLAGS="`~/usr/static/bin/eutils-config --libs`" 
make clean
make -j4


#install_name_tool -change '`pwd`/mapseq-$MSEQVER-bin/lib/libeutils-1.0.dylib' '@loader_path/lib/libeutils-1.0.dylib' mapseq-1.0-bin/mapseq
