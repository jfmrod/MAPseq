#ifndef ESEQDB_H
#define ESEQDB_H

#include "emseq_defs.h"
#include "eseq.h"
#include "eseqali.h"

#include <eutils/ethread.h>
#include <eutils/ernd.h>

class eseqtax;

class ediag
{
 public:
  long i;
  long j;
  long i2;
  long j2;
  float V;
  ediag *bestseg;
  ediag(): i(-1),j(-1),i2(-1),j2(-1),V(-1),bestseg(0x00) {}
//  ediag(int _i,int _j,int _i2,int _j2): i(_i),j(_j),i2(_i2),j2(_j2),V(-1),bestseg(0x00) {}
  ediag(long _i,long _delta,long _len): i(_i),j(_i+_len),i2(_i+_delta),j2(_i+_delta+_len),V(-1),bestseg(0x00) {}
};



class epredinfotax {
 public:
  eintarray tax;
  efloatarray cf;
  efloatarray cutoff;
};

class epredinfo {
 public:
  int seqid;
  ealigndata tophit;
  ebasicarray<ealigndata> matchcounts;
//  earrayof<float,int> ptax;
  ebasicarray<epredinfotax> ptax;
  efloatarray cutoff;
  void clear() { matchcounts.clear(); tophit.seqid=-1; }
};

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

class eseqdb;

class esearchws
{
 public:
//  eintarray best;
  ernd rng;
  ealignws alignws;
  eintarray otukmerpos;
//  uint64_t *kmerbitmask;
//  uint64_t *bitmask;
  eintarray idcount;
  eintarray seqids;

  ebasicarray<uint32_t> idcount2;
  unsigned int maskid;
  euintarray kmermask;
  unsigned int offset;
  euintarray kmerpos;
  euintarray kmerposrev;
  unsigned int offset3;
  euintarray kmerpos3;
  euintarray kmerposrev3;
  unsigned int offset2;
  euintarray kmerpos2;
  ebasicarray<eintarray> taxcounts;

  esearchws();
  esearchws(const eseqdb& seqdb);
  void init(const eseqdb& seqdb);
  void initPaired(const eseqdb& seqdb);
  void initProt(const eseqdb& seqdb);
};

typedef estr (*outfmt_fd)(const etax&,const earrayof<double,int>&,const efloatarray&);

estr outfmt_simple(const etax& tax,const earrayof<double,int>& ptax,const efloatarray& mcfarr);
estr outfmt_confidences(const etax& tax,const earrayof<double,int>& ptax,const efloatarray& mcfarr);

class eseqdb
{
 public:
  int minscore;
  int tophits;
  int topotus;
  int otulim;

  estrarrayof<eseq> seqs;
  estrhashof<int> seqind;
  eintarray seqotu;
  earray<eintarray> otus;
  ebasicarray<ekmerarray> otukmers;
  ebasicarray<etax> taxa;

  unsigned int akmers[MAXSIZE]; // 6bases

  eseqdb();

  void loadCluster(const estr& cfile);
  void makeCluster(const estr& cfile);
  void makeClusterMT(ethreads& t); //,const estr& cfile)

  void loadSequencesBinary(const estr& filename);
  void loadSequences(const estr& dbname);
  void saveSequences(const estr& filename);
  
  void loadTaxFormat1(efile& f,etax& tax);
  void loadTaxonomy(const estr& dbname);

  void init(const estr& dbname,bool nocluster=false,outfmt_fd ofmt=outfmt_confidences);

  void seqsearch(const estr& str2id,eseq& s,earray<epredinfo>& pinfoarr,esearchws& sws);
  void seqsearch_global(const estr& str2id,eseq& s,earray<epredinfo>& pinfoarr,esearchws& sws);
  void seqsearchpair(const estr& id,eseq& s,eseq& s2,earray<epredinfo>& pinfoarr,esearchws& sws);
  void pseqsearch(const estr& str2id,eseq& s,earray<epredinfo>& pinfoarr,esearchws& sws);

  void seqalign_global(const estr& str2id,eseq& s,earray<epredinfo>& previnfoarr,earray<epredinfo>& pinfoarr,esearchws& sws);

  void processQueryFASTA(const estr& fname,void (*taskfunc)(),ethreads& t);
  void processQueryFASTQ(const estr& fname,void (*taskfunc)(),ethreads& t);
  void processQueryPairend(const estr& fname,const estr& fname2,void (*taskfunc)(),ethreads& t);

  void printSearchHeader();
};

struct emtdata {
  ebasicarray<estrarrayof<eseq>*> seqs;
  ebasicarray<estrarrayof<eseq>*> sbuffer;
  earray<estr> output;

  emutex m;
  econdsig seqsSignal;
  econdsig sbufferSignal;
  eseqdb *seqdb;
//  ebasicarray<etax> *taxa;

  int ipos;
  int ilen;
  int threadi;
  int dbotus;
  euintarray otukmersize;
  eintarray mapped;
  eintarray len_si;
  
  outfmt_fd outfmt;
  bool finished;
  bool print_align;
  bool print_hits;
  bool print_kmerhits;
  bool galign;

  emtdata(): outfmt(0x00),seqdb(0x00),galign(false),print_kmerhits(false) {}
};

extern emtdata mtdata;

extern float swmin;
extern float swmax;

extern int minid1;
extern int minid2;
extern int topotus;
extern int tophits;
//extern int minscore;
extern int otulim;
extern float sweight;
extern float sweightabs;

extern float cfthres;

extern int minlen;
extern float minqual;

void taxScoreE(earrayof<double,int>& ptax,efloatarray& mcfarr,ealigndata& adata,epredinfo& pinfo,edoublearray& taxscores,etax& tax,int slen);
void taxScoreSumE(edoublearray& taxscores,epredinfo& pinfo,etax& tax,ebasicarray<eintarray>& taxcounts,int slen);
void taxScore(earrayof<double,int>& ptax,efloatarray& mcfarr,ealigndata& adata,epredinfo& pinfo,edoublearray& taxscores,etax& tax,int slen);
void taxScoreSum(edoublearray& taxscores,epredinfo& pinfo,etax& tax,ebasicarray<eintarray>& taxcounts,int slen);

void otukmeradd(ebasicarray<ekmerarray>& otukmers,int i,eseq& s,eintarray& tmpkmers,int ti,unsigned int *akmers,unsigned long akmask);
void kmercount_single(ebasicarray<ekmerarray>& otukmers,eseq& s,ebasicarray<uint32_t>& idcount,eintarray& kmerpos);
void kmercount_single(ebasicarray<ekmerarray>& otukmers,eseq& s,ebasicarray<uint32_t>& idcount,eintarray& kmerpos,const euintarray& otukmersize);

void taskSearchPaired();
void taskSearch();
void taskSearchReturn();
void taskProtSearch();

#endif

