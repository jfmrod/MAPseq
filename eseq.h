#ifndef ESEQ_H
#define ESEQ_H

#include <eutils/estr.h>
#include <eutils/eintarray.h>
#include <eutils/ebasicarray.h>

class eseqtaxlevel
{
 public:
  int tid;
  float cf;
  eseqtaxlevel(): tid(-1),cf(0.0) {}
  eseqtaxlevel(int _tid): tid(_tid),cf(0.0) {}
  eseqtaxlevel(int _tid,float _cf): tid(_tid),cf(_cf) {}
};

class eseqtax
{
 public:
//  int tid;
//  float cf;
//  float ncf;
  float bid;
//  int seqid;
  ebasicarray<eseqtaxlevel> tl;
  eseqtax(): bid(0.0) {}
//  eseqtax(): seqid(-1),bid(0.0) {}
//  eseqtax(int _seqid,float _bid): seqid(_seqid),bid(_bid) {}
  eseqtax(float _bid): bid(_bid) {}
};

class eseq
{
 public:
//  estr useq;
  estr seq;
  long seqlen;
  long seqstart;
  eintarray npos;

  eseq();
  eseq(const estr& seq);
  void setseq(const estr& seq);
  void setprot(estr& seq);
  void revcompl();
  void setrevcompl(const eseq& seq,long i,long e);
  estr print_seq(int i=0,int e=-1) const;
  estr print_pseq(int i=0,int e=-1) const;
  eseq subseq(int i,int l) const;
  bool operator<(const eseq& s);
  friend ostream& operator<<(ostream& stream,const eseq& seq);

  void serial(estr& sstr);
  long unserial(const estr& sstr,long i);
};

char cnuc2chr(uint32_t cc);
char cnuc2chru(uint32_t cc);

ostream& operator<<(ostream& stream,const eseq& seq);

#endif

