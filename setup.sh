#!/bin/bash

wget http://www.microbeatlas.org/mapref/mapref-2.2b.tar.gz
tar -Cdata -xvzf mapref-2.2b.tar.gz && mv data/mapref-2.2b/* data/ && rmdir data/mapref-2.2b && touch data/mapref-2.2b.fna


if [ ! -d "libs/eutils" ]; then
  mkdir -p libs
  svn co -r 1194 https://www.konceptfx.com/svn/eutils libs/eutils
else
  pushd libs/eutils
  svn update
  popd
fi
