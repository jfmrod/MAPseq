#!/bin/bash

wget http://www.meringlab.org/software/mapseq/mapref-2.2.tar.gz
tar -Cdata -xvzf mapref-2.2.tar.gz && mv data/mapref-2.2/* data/ && rmdir data/mapref-2.2 && touch data/mapref.fna

svn co https://www.konceptfx.com/svn/eutils
pushd eutils
svn update


