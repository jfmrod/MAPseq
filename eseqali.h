#ifndef ESEQALI_H
#define ESEQALI_H

#include <eutils/estr.h>
#include <eutils/logger.h>
#include "eseq.h"

extern uint8_t codon2prot[];
extern char aa[];

enum ealignelemtype {AT_NONE,AT_INS,AT_DEL,AT_MATCH,AT_MISS,AT_LEFT,AT_ID,AT_ALIGN,AT_RIGHT,AT_PAIR,AT_COMPAT};

class ealignelem {
 public:
  int count;
  ealignelemtype type;
  ealignelem();
  ealignelem(ealignelemtype type,int count);
  bool operator==(const ealignelem& e) const;
  bool operator!=(const ealignelem& e) const;
};

class ealignprofile {
 public:
  ebasicarray<ealignelem> elm;
  void add(ealignelemtype type,int count=1);
  void add(const ealignprofile& p);
  void addinv(const ealignprofile& p);
  void inv();
  bool operator==(const ealignprofile& p) const;
  estr str();
};

ostream& operator<<(ostream& stream,const ealignprofile& p);

class ealignscore {
 public:
  float dropoff;
  float match;
  float mismatch;
  float gapopen;
  float gapext;
  ealignscore(): match(0.0),mismatch(0.0),gapopen(0.0),gapext(0.0),dropoff(0.0) {}
  ealignscore(float _m,float _mi,float _go,float _ge,float _dropoff): match(_m),mismatch(_mi),gapopen(_go),gapext(_ge),dropoff(_dropoff) {}
};

class epalignscore {
 public:
  float dropoff;
  float pmatch;
  int8_t *smatrix;
  uint32_t smatrixwidth;
  float gapopen;
  float gapext;
  epalignscore(): smatrix(0x00),smatrixwidth(0),gapopen(0.0),gapext(0.0),dropoff(0.0),pmatch(0.0) {}
  epalignscore(int8_t *_m,int swidth,float _pmatch,float _go,float _ge,float _dropoff): smatrix(_m),smatrixwidth(swidth),pmatch(_pmatch),gapopen(_go),gapext(_ge),dropoff(_dropoff) {}
};

class ealigndata {
 public:
  double _eval;
  double _score;
  unsigned char* aln;
  int seqid;
  bool revcompl;
  int kmercount;
  int matches;
  int mismatches;
  int gaps;
  long s1;
  long e1;
  long s2;
  long e2;

  ealignprofile profile;
//  eintarray ainfo;
//  bool operator<(const ealigndata& d) const { return(score()<d.score()); }
  void operator+=(const ealigndata& d) { if (s1<=d.e1 && e1>=d.s1) return;  _score+=d._score; matches+=d.matches; mismatches+=d.mismatches; gaps+=d.gaps; profile.add(AT_PAIR); profile.add(d.profile); }
  bool operator==(const ealigndata& d) const { return(score()==d.score() && matches==d.matches && mismatches==d.mismatches && gaps==d.gaps); }
  bool operator<(const ealigndata& d) const { return(score()<d.score() || (score()==d.score() && identity()<d.identity())); }
//  bool operator<(const ealigndata& d) const { return(eval()>d.eval()); }
//  bool operator<(const ealigndata& d) const { return(identity()<d.identity() || identity()==d.identity() && score()<d.score()); }
  double eval() const { return(_eval); }
  double score() const { return(_score); }

  double identity() const { return(float(matches)/(matches+mismatches+gaps)); }
//  double identity() const;
//  double identity() const { return(float(matches)/(matches+mismatches)); }
  ealigndata(): _eval(1.0),_score(0.0),s1(-1),e1(-1),s2(-1),e2(-1),matches(0),mismatches(0),gaps(0),aln(0x00) {}
//  double score() const { return(matches*matchcost + mismatches*misscost + gaps*gapcost); }
  estr compress(const eseq& s1);
  estr align_str(const eseq& s1,const eseq& s2);
  estr palign_str(const eseq& s1,const eseq& s2);

  void globalalign(const ealignscore& as);
};

estr sali_decompress(const estr& sastr,const eseq& s2);

class ealignws
{
 public:
  double *mF;
  int size;
  ealignws();
  ~ealignws();
  void reserve(int _size);
};


float seqcalign(const eseq& a,int pa,int ea,const eseq& b,int pb,int eb,estr& as1,estr& as2);
float seqcalign_global(const eseq& a,int pa,int ea,const eseq& b,int pb,int eb,estr& as1,estr& as2,const ealignscore& as);
float seqcalign_global_noedgegap(const eseq& a,int pa,int ea,const eseq& b,int pb,int eb,estr& as1,estr& as2,const ealignscore& as);
float seqcalign_global_noleftedgegap(const eseq& a,int pa,int ea,const eseq& b,int pb,int eb,estr& as1,estr& as2,const ealignscore& as);
float seqcalign_global_norightedgegap(const eseq& a,int pa,int ea,const eseq& b,int pb,int eb,estr& as1,estr& as2,const ealignscore& as);

void seqcalign_global(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const ealignscore& as);
void seqcalign_global_noedgegap(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const ealignscore& as);
void seqcalign_global_noleftedgegap(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const ealignscore& as);
void seqcalign_global_norightedgegap(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const ealignscore& as);
void pseqcalign_global(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const epalignscore& as);

void pseqcalign_local_leftext(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const epalignscore& as);
void pseqcalign_local_rightext(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const epalignscore& as);


void seqcalign_local_leftext(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const ealignscore& as);
void seqcalign_local_rightext(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const ealignscore& as);

void seqcalign_nb_global_noedgegap(const eseq& a,int pa,int ea,const eseq& b,int pb,int eb,ealigndata& adata,unsigned char *tncount,ealignws& ws,const ealignscore& as);

void print_seqali(const estr& s1,const estr& s2);

void invertstr(estr& str);

inline int maxs(double a,double b,double c)
{
  if (a>b){
    if (c<a) return(0);
    return(2);
  }
  if (c<b) return(1);
  return(2);
}
inline double maxv(double a,double b,double c)
{
  if (a>b){
    if (c<a) return(a);
    return(c);
  }
  if (c<b) return(b);
  return(c);
}
inline void maxvg(double& retval,char& gapOpen,double gapx,double gapy,double align)
{
  if (gapx>gapy){
    if (align<gapx) { gapOpen=1; retval=gapx; return; }
    gapOpen=0;
    retval=align;
    return;
  }
  if (align<gapy) { gapOpen=2; retval=gapy; return; }
  gapOpen=0;
  retval=align;
}


template <class T,class K>
class eswalign
{
 public:
  int m_score[1u<<8u];


  int gap;
  int gapext;
  int match;
  int mismatch;

  int score;
  int aligned;
  int gaps;
  int matches;
  int mismatches;
  int maxsize;
  
  estr seq1;
  estr seq2;

  eswalign();
  void init();

  float align(const T& a,const T& b);
  void print();
  void print_info();
};

template <class T,class K>
eswalign<T,K>::eswalign(): match(1),mismatch(2),gap(5),gapext(2) {}

template <class T,class K>
void eswalign<T,K>::init()
{
  for (int i=0; i<(1u<<8u); ++i){
    m_score[i]=mismatch;
  }
  m_score[0x00u]=match;
}

template <class T,class K>
float eswalign<T,K>::align(const T& a,const T& b)
{
  score=0; aligned=0; gaps=0; matches=0; mismatches=0;
  maxsize=MAX(a.size(),b.size());
  if (a.size()==0 || b.size()==0) return(0.0);

  int w=a.size()+1;
  int h=b.size()+1;

  int i,j;
  double *mF=new double[w*h];
  char *gapOpen=new char[w*h];
//  char hgapOpen;
  
  mF[0]=0.0l;
  gapOpen[0]=0;

  // no gap penalty for beginning gaps
  for (i=1; i<a.len()+1; ++i)
    { mF[i*h]=0.0l; gapOpen[i*h]=0; }
  for (j=1; j<b.len()+1; ++j)
    { mF[j]=0.0l; gapOpen[j]=0; }

  double max_S=0.0;
  int max_i=-1;
  int max_j=-1;
  for (i=1; i<a.len(); ++i){
    for (j=1; j<b.len(); ++j) {
      maxvg(mF[i*h+j],gapOpen[i*h+j],mF[(i-1)*h+j]-(gapOpen[(i-1)*h+j]==1?gapext:gap),mF[i*h+j-1]-(gapOpen[i*h+j-1]==2?gapext:gap),mF[(i-1)*h+j-1]+m_score[a[i-1]^b[j-1]]);
      if (mF[i*h+j]>max_S){ max_S=mF[i*h+j]; max_i=i; max_j=j; }
    }
    // no gap penalty at end of seq
    maxvg(mF[i*h+j],gapOpen[i*h+j],mF[(i-1)*h+j],mF[i*h+j-1]-(gapOpen[i*h+j-1]==2?gapext:gap),mF[(i-1)*h+j-1]+m_score[a[i-1]^b[j-1]]);
    if (mF[i*h+j]>max_S){ max_S=mF[i*h+j]; max_i=i; max_j=j; }
  }

  // no gap penalties at end of seq
  for (j=1; j<b.len(); ++j) {
    maxvg(mF[i*h+j],gapOpen[i*h+j],mF[(i-1)*h+j]-(gapOpen[(i-1)*h+j]==1?gapext:gap),mF[i*h+j-1],mF[(i-1)*h+j-1]+m_score[a[i-1]^b[j-1]]);
    if (mF[i*h+j]>max_S){ max_S=mF[i*h+j]; max_i=i; max_j=j; }
  }
  maxvg(mF[i*h+j],gapOpen[i*h+j],mF[(i-1)*h+j],mF[i*h+j-1],mF[(i-1)*h+j-1]+m_score[a[i-1]^b[j-1]]);
  if (mF[i*h+j]>max_S){ max_S=mF[i*h+j]; max_i=i; max_j=j; }

/*
  for (i=0; i<=a.size(); ++i){
    for (j=0; j<=b.size(); ++j){
      printf("%7.1lf ",mF[i][j]);
    }
    printf("\n");
  }
*/

  ldieif(max_i==-1 || max_j==-1,"max_i or max_j are negative");
//  cout << "max_i: " << max_i << " max_j: " << max_j << " max_S: " << max_S << endl;

  aligned=0;
  matches=0;

  seq1.clear(); seq2.clear();

  double gapseqa,gapseqb;
//  i=a.len(); j=b.len(); // global alignment starting position
  i=max_i; j=max_j; // local alignment starting from maximum score
  while (i>0 && j>0){
    gapseqa=mF[(i-1)*h+j];
    if (i<a.len())
      gapseqb-=(gapOpen[i*h+j-1]==1?gapext:gap);
    gapseqb=mF[i*h+j-1];
    if (j<b.len())
      gapseqb-=(gapOpen[(i-1)*h+j]==2?gapext:gap);
    switch (maxs(gapseqa,gapseqb,mF[(i-1)*h+j-1]+m_score[a[i-1]^b[j-1]])){
      case 0: seq1 += a[i-1]; seq2 += '-'; --i; break;
      case 1: seq1 += '-'; seq2 += b[j-1]; --j; break;
      case 2: seq1 += a[i-1]; seq2 += b[j-1]; ++aligned; if (a[i-1]==b[j-1]) ++matches; --j; --i; break;
     default:
      ldie("unknown code");
    }
  }
  for (;i==0 && j>0;--j){
    seq1+='-'; seq2+=b[j-1];
  }
  for (;j==0 && i>0;--i){
    seq2+='-'; seq1+=a[i-1];
  }
  invertstr(seq1);
  invertstr(seq2);

/*
  i=a.size()-1; j=b.size()-1;
  while (i>0 && j>0){
    switch (maxs(mF[i-1][j],mF[i][j-1],mF[i-1][j-1])){
      case 0: --i; break;
      case 1: --j; break;
      case 2: ++aligned; if (a[i]==b[j]) ++matches; --j; --i; break;
    }
  }
  while (j>0 && i>=0){
    switch (maxs(-gap-(j-1)*gap,mF[i][j-1],-(j-1)*gap)){
      case 0: --i; break;
      case 1: --j; break;
      case 2: ++aligned; if (a[i]==b[j]) ++matches; --j; --i; break;
    }
  }
  while (j>=0 && i>0){
    switch (maxs(mF[i-1][j],-gap-(i-1)*gap,-(i-1)*gap)){
      case 0: --i; break;
      case 1: --j; break;
      case 2: ++aligned; if (a[i]==b[j]) ++matches; --j; --i; break;
    }
  }

  if (i==0 && j==0){
    ++aligned; if (a[i]==b[j]) ++matches;
  }
*/
  delete[] mF;
  delete[] gapOpen;

  mismatches=aligned-matches;
  gaps=maxsize-aligned;
  return((float)matches/aligned);
}

template <class T,class K>
void eswalign<T,K>::print_info() //const estr& seq,const estr& seq2)
{
  ldieif(seq1.len()!=seq2.len(),"seq mismatches"); 
  int gapscore=0;
  int gapopen=0;
  score=0; aligned=0; gaps=0; matches=0; mismatches=0;
  for (int i=0; i<seq1.len(); ++i){
    if (seq1[i]=='-' || seq2[i]=='-'){
      if (seq1[i]=='-'){
        if (gapopen==0 || gapopen==2) { gapopen=1; gapscore+=gap; }
        else gapscore+=gapext;
      }else{
        if (gapopen==0 || gapopen==1) { gapopen=2; gapscore+=gap; }
        else gapscore+=gapext;
      }
      ++gaps;
    }else if (seq1[i]==seq2[i])
      ++matches;
    else
      ++mismatches;
  }
  aligned=matches+mismatches;
  score=matches*match - mismatches*mismatch - gapscore;
  cout << "score: " << score << " identity: " << (float)matches/aligned << " aligned: " << aligned << " matches: " << matches << " mismatches: " << mismatches << " gaps: " << gaps << " aligned_seq_len: " << seq1.len() << endl; 
}

template <class T,class K>
void eswalign<T,K>::print() //const estr& seq,const estr& seq2)
{
  ldieif(seq1.len()!=seq2.len(),"seq1uences not aligned");
  if (seq2.len()==0) return;

  int i=0;
  estr tmpseq1;
  int pseq2=0;
  int pseq1=0;
  if (seq2[0]=='-')
    for (i=0; i<seq2.len() && seq2[i]=='-'; ++i,++pseq1);
  else if (seq1[0]=='-')
    for (i=0; i<seq1.len() && seq1[i]=='-'; ++i,++pseq2);

  for (; i<seq2.len(); i+=100){
    cout << seq2.substr(i,100) << endl;
    tmpseq1.clear();
    for (int j=i; j<i+100 && j<seq1.len(); ++j){
      if (seq2[j]!='-') ++pseq2; if (seq1[j]!='-') ++pseq1;
      if (seq2[j] == '-' || seq1[j] =='-') tmpseq1+=' ';
      else if (seq2[j] == seq1[j]) tmpseq1+='|';
      else tmpseq1+='.';
    }
    cout << tmpseq1 << endl;
    cout << seq1.substr(i,100) << endl;
    cout << endl;
  }
}





#endif /* ESEQALI_H */

