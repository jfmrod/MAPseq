#ifndef ESEQDB_H
#define ESEQDB_H

#include "emseq_defs.h"

class etax
{
 public:
  estr name;
  earray<estr> levels;
  earray<earray<estr> > names;
  earray<estrhashof<int> > ind;
  ebasicarray<eseqtax*> seqs;
  efloatarray cutoff;
  efloatarray cutoffcoef;
};


class eseqdb
{
 public:
  estrarrayof<eseq> seqs;
  estrhashof<int> seqind;
  eintarray seqotu;
  earray<eintarray> otus;
  ebasicarray<ekmerarray> otukmers;
  ebasicarray<etax> taxa;

  void loadCluster(const estr& cfile);
  void makeCluster(const estr& cfile);

  void loadSequences(const estr& dbname);
  
  void loadTaxFormat1(efile& f,etax& tax);
  void loadTaxonomy(const estr& dbname);

  void init(const estr& dbname);
};




#endif

