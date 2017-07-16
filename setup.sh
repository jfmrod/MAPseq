#!/bin/bash

wget http://www.meringlab.org/software/mapseq/mapref-2.0.tar.gz
tar -Cdata -xvzf mapref-2.0.tar.gz && mv data/mapref-2.0/* data/ && rmdir data/mapref-2.0 && touch data/mapref.fna

svn co https://www.konceptfx.com/svn/eutils
pushd eutils
svn update


