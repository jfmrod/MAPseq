#!/bin/bash

unset LD_LIBRARY_PATH
unset DYLD_LIBRARY_PATH
unset LIBRARY_PATH
unset CPLUS_INCLUDE_PATH
unset C_INCLUDE_PATH
unset OBJC_INCLUDE_PATH
unset OBJCPLUS_INCLUDE_PATH

./configure --enable-makestatic --prefix=`pwd`/mapseq-1.2.3-macosx --bindir=`pwd`/mapseq-1.2.3-macosx CXXFLAGS="-O2"


#install_name_tool -change '/Users/joao/work/libdev/mapseq/mapseq-1.0-bin/lib/libeutils-1.0.dylib' '@loader_path/lib/libeutils-1.0.dylib' mapseq-1.0-bin/mapseq
