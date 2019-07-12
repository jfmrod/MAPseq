#include <eutils/emain.h>
#include <eutils/ebasichashmap.h>
#include <eutils/estrhashof.h>
#include <eutils/einthashof.h>
#include <eutils/efile.h>
#include <eutils/etimer.h>
#include <eutils/eheap.h>
#include <eutils/ernd.h>
#include <eutils/ethread.h>
#include <eutils/esystem.h>
#include <eutils/eoption.h>
#include <deque>
#include <unistd.h>

#define LDEBUG(key,cmd) LDEBUG2(key,cmd)
#define LDEBUG2(key,cmd) LDEBUG_##key(cmd)
#define LDEBUG_(cmd) {}
#define LDEBUG_TRUE(cmd) cmd

#define D_SEQALIGNMENT
#define D_SEQIDENT
#define D_SEQSEARCH
#define D_PROFILE

//#define D_SEQIDENT TRUE
//#define D_SEQALIGNMENT TRUE
//#define D_SEQSEARCH TRUE


#include <math.h>

#include <map>
using namespace std;

#include "mapseq-config.h"
#include "ekmerhashmap.h"
#include "eseqali.h"

#include "kmerseqtables-data.h"

#define KMERSIZE 10ul
#define KMERBITS (KMERSIZE*2ul)
#define MAXSIZE (1ul<<KMERBITS)
#define KMERMAX (1ul<<KMERBITS)
#define KMERMASK (MAXSIZE-1ul)

#define KMERSIZE2 8ul
#define KMERBITS2 (KMERSIZE2*2ul)
#define MAXSIZE2 (1ul<<KMERBITS2)
#define KMERMAX2 (1ul<<KMERBITS2)
#define KMERMASK2 (MAXSIZE2-1ul)



const uint32_t nuc[]={'a','t','g','c','u'};
const uint32_t compnuc[]={0x0u,0x1u,0x2u,0x3u,0x1u};
const uint32_t consnuc[]={0x0001u,0x0010u,0x0100u,0x1000u};
const unsigned long safe_shift[32]={0x0ul,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful};

unsigned int akmers[MAXSIZE]; // 6bases

int topotus=10;
int tophits=20;
int minscore=30;
int minid1=1;
int minid2=1;
int otulim=50;
float lambda=1.280;
float sweight=30.0;

/*
const float matchcost=2.0;
const float misscost=-1.0;
const float gapcost=-1.0;
const float gapopen=-5.0;
*/

//const ealignscore as(1.0,2.0,10.0,1.0);
//const ealignscore as(2.0,1.0,5.0,1.0,30.0);
const ealignscore as(1.0,1.0,5.0,1.0,20.0);
//const ealignscore as(2.0,1.0,5.0,1.0,30.0);

const int NCOUNT_MAXLEN=1600; // IMPORTANT: must be disivible by 8 to be 64bit aligned

float ti=0.0,ts=0.0,ti2=0.0,ts2=0.0,ta=0.0;
float tdp=0.0,tdp1=0.0,tdp2=0.0,tdpfl=0.0,tdpmd=0.0;
etimer t1,t2;

const uint32_t BMASK31=(1u<<31u)-1u;
const uint32_t BMASK16=(1u<<16u)-1u;
const uint32_t BMASK11=(1u<<11u)-1u;
const uint32_t BMASK10=(1u<<10u)-1u;


class etax
{
 public:
  estr name;
  earray<estr> levels;
  earray<estrarray> names;
  earray<estrhashof<int> > ind;
  ebasicarray<eseqtax*> seqs;
  efloatarray cutoff;
  efloatarray cutoffcoef;
};

typedef ebasicarray<unsigned int> euintarray;

class esearchws
{
 public:
//  eintarray best;
  ernd rng;
  ealignws alignws;
  eintarray otukmerpos;
  uint64_t *kmerbitmask;
  uint64_t *bitmask;
  eintarray idcount;
  ebasicarray<uint32_t> idcount2;
  unsigned int maskid;
  euintarray kmermask;
  unsigned int offset;
  euintarray kmerpos;
  euintarray kmerposrev;
  unsigned int offset2;
  euintarray kmerpos2;
  ebasicarray<eintarray> taxcounts;
};


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

/*
class ekmer
{
 public:
  int i;
//  eintarray pos;
  int pos;
  ekmer(): i(-1) {}
  ekmer(int _i,int _pos): i(_i),pos(_pos) {}
//  ekmer(int _i,eintarray& _pos): i(_i),pos(_pos) {}
};
*/

class emaxcount {
 public:
  int i;
  int maxcount;
  eintarray kwordaln;
  emaxcount(): i(-1), maxcount(0) {}
};

char oct2hex(unsigned char o)
{
  switch (o){
    case 0x00: return('0');
    case 0x01: return('1');
    case 0x02: return('2');
    case 0x03: return('3');
    case 0x04: return('4');
    case 0x05: return('5');
    case 0x06: return('6');
    case 0x07: return('7');
    case 0x08: return('8');
    case 0x09: return('9');
    case 0x0a: return('a');
    case 0x0b: return('b');
    case 0x0c: return('c');
    case 0x0d: return('d');
    case 0x0e: return('e');
    case 0x0f: return('f');
  }
  ldie("wrong octet value");
  return('-');
}


char* uint2hex(char *tmpstr,unsigned int c)
{
  tmpstr[0u]='0';
  tmpstr[1u]='x';
  for (unsigned int k=0; k<sizeof(c)*8u/4u; ++k)
    tmpstr[2u+k]=oct2hex((c>>(sizeof(c)*8u-4u-k*4u))&0xf);
  tmpstr[2u+sizeof(c)*8u/4u]='u';
  tmpstr[2u+sizeof(c)*8u/4u+1u]=0x00;
  return(tmpstr);
}


uint32_t chr2kmer(const char *str)
{
  uint32_t *tmp=(uint32_t*)str;
  return(seq_comp_table[tmp[0]&0xffffu]|(seq_comp_table[(tmp[0]>>16u)&0xffffu]<<4u)|(seq_comp_table[tmp[1]&0xffffu]<<8u)|(seq_comp_table[(tmp[1]>>16u)&0xffffu]<<12u));
}

#define IKMERSIZE 8ul
#define IKMERBITS (IKMERSIZE*2ul)
#define IKMERMAX (1ul<<IKMERBITS)
#define IKMERMASK (IKMERMAX-1ul)



estr long2str(unsigned long kmer)
{
  estr tmpstr;
  for (int i=0; i<64u; i+=2u,kmer>>=2u){
    switch (0x03u&kmer){
      case 0x00u: tmpstr+='A'; break;
      case 0x01u: tmpstr+='T'; break;
      case 0x02u: tmpstr+='G'; break;
      case 0x03u: tmpstr+='C'; break;
    }
  }
  return(tmpstr);
}
estr kmer2str(uint32_t kmer)
{
  estr tmpstr;
  for (int i=0; i<IKMERSIZE; ++i,kmer>>=2u){
    switch (0x03u&kmer){
      case 0x00u: tmpstr+='A'; break;
      case 0x01u: tmpstr+='T'; break;
      case 0x02u: tmpstr+='G'; break;
      case 0x03u: tmpstr+='C'; break;
    }
  }
  return(tmpstr);
}

class eseqdb
{
 public:
  estrarrayof<eseq> seqs;
  estrhashof<int> seqind;
  eintarray seqotu;
  earray<eintarray> otus;
  ebasicarray<deque<uint32_t> > otukmers;
  ebasicarray<etax> taxa;
};

unsigned long seqkmer(const eseq& s,long p1)
{
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s.seq._str);
  return(((pstr1[p1/32u]>>(2u*(p1%32u)))|((pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u]))&KMERMASK);
}

unsigned long seqrevkmer(const eseq& s,long p1)
{
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s.seq._str);
  p1=s.seqlen-p1;
  return(kmer_rev_lt[((pstr1[p1/32u]>>(2u*(p1%32u)))|((pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u]))&KMERMASK]);
}

void print_tncount(unsigned char* tncounts,int p,int e)
{
  uint64_t vid;

  for (int i=p; i<e; ++i){
    vid=tncounts[i];
    switch (vid&0x3u){ 
      case 0x0u: cout << "i: " << i << " -" << endl; break;
      case 0x1u: cout << "i: " << i << " m" << endl; break;
      case 0x2u: cout << "i: " << i << " x" << endl; break;
    }
  }
}

unsigned char d1_lt[1u<<16u];
unsigned char d2_lt[1u<<16u];
unsigned char id_lt[1u<<16u];


void initdlt()
{
  for (uint32_t i=0u; i<(1u<<16u); ++i){
    id_lt[i]=0u;
    d1_lt[i]=8u;
    d2_lt[i]=8u;
  }
  id_lt[0u]=1u;
  d2_lt[0u]=0u;
  for (uint32_t i=0u; i<8u; ++i)
    d1_lt[i]=i;
}

void seqident_seg_left(const eseq& s1,int p1,const eseq& s2,int p2,ealigndata& adata,const ealignscore& as){
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s1.seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(s2.seq._str);
  unsigned long v1,v2;
  unsigned char tmp[5];
  int k;
  int itmp;

  if (p1==0 || p2==0) return;

//  cout << "seg_left id: " << p1 << " :: " << p2 << " len1: " << p1 << " len2: " << p2 << endl;
  for (; p1>=32 && p2>=32; p1+=k-32u+IKMERSIZE,p2+=k-32u+IKMERSIZE){
    v1=pstr1[p1/32u-1u]>>(2u*(p1%32u));
    v2=pstr2[p2/32u-1u]>>(2u*(p2%32u));
    v1|=(pstr1[p1/32u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    v2|=(pstr2[p2/32u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
//    cout << p1 << " " << long2str(v1) << endl;
//    cout << p2 << " " << long2str(v2) << endl;
    for (k=32u-IKMERSIZE; k>=0; k-=IKMERSIZE,v1<<=IKMERSIZE*2u,v2<<=IKMERSIZE*2u){
//      cout << p1 << "+" << k << " " << kmer2str((v1>>((32u-IKMERSIZE)*2u))&IKMERMASK) << endl; 
//      cout << p2 << "+" << k << " " << kmer2str((v2>>((32u-IKMERSIZE)*2u))&IKMERMASK) << endl; 
      tmp[0]=seq_ident_table[((v1^v2)>>((32u-IKMERSIZE)*2u))&IKMERMASK];
      tmp[1]=seq_ident_table[((v1^(v2>>2ul))>>((32u-IKMERSIZE)*2u))&IKMERMASK];
      tmp[2]=seq_ident_table[((v1^(v2>>4ul))>>((32u-IKMERSIZE)*2u))&IKMERMASK];
      tmp[3]=seq_ident_table[((v2^(v1>>2ul))>>((32u-IKMERSIZE)*2u))&IKMERMASK];
      tmp[4]=seq_ident_table[((v2^(v1>>4ul))>>((32u-IKMERSIZE)*2u))&IKMERMASK];
      
      itmp=0;
      for (int i=1; i<5; ++i){
        if (tmp[i]>tmp[itmp]) itmp=i;
      }
      if (itmp!=0 && tmp[itmp]>tmp[0]+IKMERSIZE/4){
//        matches+=tmp[0]; mismatches+=IKMERSIZE-tmp[0];
        switch(itmp){
          case 1: p1-=1; adata.gaps+=1; break;
          case 2: p1-=2; adata.gaps+=2; break;
          case 3: p2-=1; adata.gaps+=1; break;
          case 4: p2-=2; adata.gaps+=2; break;
        }
        break;
      }else{
        adata.matches+=tmp[0]; adata.mismatches+=IKMERSIZE-tmp[0];
        adata._score+=tmp[0]*as.match-(IKMERSIZE-tmp[0])*as.mismatch;
      }
    }
  }
//  cout << "seg_left rest id: " << p1 << " :: " << p2 << " len1: " << p1 << " len2: " << p2 << endl;
  for (;p1>=int(IKMERSIZE) && p2>=int(IKMERSIZE); p1+=k+int(IKMERSIZE)-32,p2+=k+int(IKMERSIZE)-32){ // needed to handle gap shifts gracefully
//    cout << p1 << endl;
//    cout << p2 << endl;
    v1=(p1>=32?pstr1[p1/32u-1u]>>(2u*(p1%32u)):0ul);
    v2=(p2>=32?pstr2[p2/32u-1u]>>(2u*(p2%32u)):0ul);
    v1|=(pstr1[p1/32u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    v2|=(pstr2[p2/32u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
//    cout << p1 << " " << long2str(v1) << endl;
//    cout << p2 << " " << long2str(v2) << endl;
    for (k=32u-IKMERSIZE; k>=0 && p1+k>=int(IKMERSIZE) && p2+k>=int(IKMERSIZE); k-=IKMERSIZE,v1<<=IKMERSIZE*2u,v2<<=IKMERSIZE*2u){
//      cout << p1 << "+" << k << " " << kmer2str((v1>>((32u-IKMERSIZE)*2u))&IKMERMASK) << endl; 
//      cout << p2 << "+" << k << " " << kmer2str((v2>>((32u-IKMERSIZE)*2u))&IKMERMASK) << endl; 
      tmp[0]=seq_ident_table[((v1^v2)>>((32u-IKMERSIZE)*2u))&IKMERMASK];
      tmp[1]=seq_ident_table[((v1^(v2>>2ul))>>((32u-IKMERSIZE)*2u))&IKMERMASK];
      tmp[2]=seq_ident_table[((v1^(v2>>4ul))>>((32u-IKMERSIZE)*2u))&IKMERMASK];
      tmp[3]=seq_ident_table[((v2^(v1>>2ul))>>((32u-IKMERSIZE)*2u))&IKMERMASK];
      tmp[4]=seq_ident_table[((v2^(v1>>4ul))>>((32u-IKMERSIZE)*2u))&IKMERMASK];
      
      itmp=0;
      for (int i=1; i<5; ++i){
        if (tmp[i]>tmp[itmp]) itmp=i;
      }
      if (itmp!=0 && tmp[itmp]>tmp[0]+IKMERSIZE/4){
        switch(itmp){
          case 1: p1-=1; adata.gaps+=1; break;
          case 2: p1-=2; adata.gaps+=2; break;
          case 3: p2-=1; adata.gaps+=1; break;
          case 4: p2-=2; adata.gaps+=2; break;
        }
        break;
      }else{
        adata.matches+=tmp[0]; adata.mismatches+=IKMERSIZE-tmp[0];
        adata._score+=tmp[0]*as.match-(IKMERSIZE-tmp[0])*as.mismatch;
      }
    }
  }

  int l=(p1<p2?p1:p2);
//  cout << "missing: " << p1 << " l: " << l << endl;
  if (l<=0) return;
  p1-=l;
  p2-=l;
  v1=pstr1[p1/32u]>>(2u*(p1%32u));
  v2=pstr2[p2/32u]>>(2u*(p2%32u));
  v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
  v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
//  cout << p1 << " " << long2str(v1) << endl;
//  cout << p2 << " " << long2str(v2) << endl;

  tmp[0]=seq_ident_table[(v1^v2)&((0x1ul<<(l*2ul))-1ul)];
//  cout << long2str(v1&((0x1ul<<(l*2u))-1ul)) << " " << long2str(((0x1ul<<(l*2u))-1ul)) << endl;
//  cout << long2str(v2&((0x1ul<<(l*2u))-1ul)) << " " << long2str(((0x1ul<<(l*2u))-1ul)) << endl;
  adata.matches+=tmp[0]-(IKMERSIZE-l);
  adata.mismatches+=IKMERSIZE-tmp[0];
  adata._score+=(tmp[0]-(IKMERSIZE-l))*as.match-(IKMERSIZE-tmp[0])*as.mismatch;

/*
  tmp[0]=seq_ident_table[((v1^v2)>>((32u-IKMERSIZE)*2u))&~((0x1ul<<(IKMERSIZE*2u-l*2ul))-1ul)&IKMERMASK]; // put mask at end of 64bit int, i.e.: FF000000

  cout << long2str((v1>>((32u-IKMERSIZE)*2u))&(~((0x1ul<<(IKMERSIZE*2u-l*2ul))-1ul))&IKMERMASK) << " " << long2str((~((0x1ul<<(IKMERSIZE*2u-l*2ul))-1ul))&IKMERMASK) << endl;
  cout << long2str((v2>>((32u-IKMERSIZE)*2u))&~((0x1ul<<(IKMERSIZE*2u-l*2ul))-1ul)&IKMERMASK) << " " << long2str(~((0x1ul<<(IKMERSIZE*2u-l*2ul))-1ul)&IKMERMASK) << endl;
  matches+=tmp[0]-(IKMERSIZE-l);
  mismatches+=IKMERSIZE-tmp[0];
*/
}


bool galign=true;
bool ignoreEmptyTax=false;


inline void sumcounts(unsigned char* tncount,int p1,uint64_t tmpnuc)
{
  *reinterpret_cast<uint64_t*>(&tncount[p1])|=tmpnuc;
/*
  for (int i=0; i<8; ++i,tmpnuc>>=8u)
    tncount[p1+i]|=(tmpnuc&0xff);
*/
//  cout << uint2hex(tmpstr,*p) << " " << uint2hex(tmpstr2,tmpnuc) << " " << uint2hex(tmpstr3,tncount[p1]) << endl;
}

void seqncount(const eseq& s1,long p1,long e1,const eseq& s2,long p2,long e2,ealigndata& adata,const ealignscore& as)
{
  int k;
  adata.matches+=e1-p1;
  adata._score+=(e1-p1)*as.match;
  adata.profile.add(AT_MATCH,e1-p1);
  if (adata.aln==0x00) return;

/*
  uint64_t *pstr1=reinterpret_cast<uint64_t*>(s1.seq._str);
  uint64_t *pstr2=reinterpret_cast<uint64_t*>(s2.seq._str);
  uint64_t v1,v2;
*/

//  cout << "id: " << p1 << "," << e1 << " :: " << p2 << "," <<e2 << " len1: " << e1-p1 << " len2: " << e2-p2 << endl;
  for (; p2+IKMERSIZE<=e2; p2+=IKMERSIZE,p1+=IKMERSIZE){
//    cout << "sumcounts: " << p1 << endl;
//    v1=(pstr1[p1/32u]>>(2u*(p1%32u)))|((pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u]);
//    v2=(pstr2[p2/32u]>>(2u*(p2%32u)))|((pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u]);
    sumcounts(adata.aln,p1,seq_alignment_lt[0u]);
//    tmpas1+=kmer2str(v1);
//    tmpas2+=kmer2str(v2);
  }

  int l=e2-p2;
  if (l<=0) return;
//  if (l<=0) goto fend;
//  cout << "missing: p1: " << p1 << " p2: " << p2 << " l: " << l << endl;
  
  sumcounts(adata.aln,p1,seq_alignment_lt[0u]&((0x1ul<<(l*8u))-1ul));
//  v1=(pstr1[p1/32u]>>(2u*(p1%32u)))|((pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u]);
//  v2=(pstr2[p2/32u]>>(2u*(p2%32u)))|((pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u]);
//  tmpas1+=kmer2str(v1).substr(0,l);
//  tmpas2+=kmer2str(v2).substr(0,l);
// fend:
//  pas1=tmpas1+pas1;
//  pas2=tmpas2+pas2;
}



void seqident_seg_left_local(const eseq& s1,long p1,const eseq& s2,long p2,ealigndata& adata,const ealignscore& as){
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s1.seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(s2.seq._str);
  unsigned long v1,v2;
  unsigned char tmp[5];
  int k;
  int itmp;

  if (p1==0 || p2==0) {
    adata.s1=p1;
    adata.s2=p2;
    return;
  }

//  cout << "seg_left id: " << p1 << " :: " << p2 << " len1: " << p1 << " len2: " << p2 << endl;
  for (; p1>=32 && p2>=32; p1-=k,p2-=k){
    v1=pstr1[p1/32u-1u]>>(2u*(p1%32u));
    v2=pstr2[p2/32u-1u]>>(2u*(p2%32u));
    v1|=(pstr1[p1/32u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    v2|=(pstr2[p2/32u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
//    cout << p1 << " " << long2str(v1) << endl;
//    cout << p2 << " " << long2str(v2) << endl;
    for (k=IKMERSIZE; k<32; k+=IKMERSIZE,v1<<=IKMERSIZE*2u,v2<<=IKMERSIZE*2u){
//      cout << p1 << "+" << k << " " << kmer2str((v1>>((32u-IKMERSIZE)*2u))&IKMERMASK) << endl; 
//      cout << p2 << "+" << k << " " << kmer2str((v2>>((32u-IKMERSIZE)*2u))&IKMERMASK) << endl; 
      tmp[0]=seq_ident_table[((v1^v2)>>(64u-IKMERSIZE*2u))&IKMERMASK];
     
//      if (!galign && tmp[0]<IKMERSIZE/2){
      if (tmp[0]<IKMERSIZE/2){
        adata.s1=p1-k;
        adata.s2=p2-k;
        return;
      }
      adata.matches+=tmp[0]; adata.mismatches+=IKMERSIZE-tmp[0];
//      sumcounts(tncount,p1-k,seq_alignment_lt[((v1^v2)>>(64u-IKMERSIZE*2u))&IKMERMASK]);
    }
  }
//  cout << "seg_left rest id: " << p1 << " :: " << p2 << " len1: " << p1 << " len2: " << p2 << endl;
//  cout << "2nd last p1: " << p1 << " p2: " << p2 << endl;
  for (;p1>=int(IKMERSIZE) && p2>=int(IKMERSIZE); p1-=k,p2-=k){ // needed to handle gap shifts gracefully
//    cout << "- 2nd last p1: " << p1 << " p2: " << p2 << endl;
//    cout << p1 << endl;
//    cout << p2 << endl;
    v1=(p1>=32?pstr1[p1/32u-1u]>>(2u*(p1%32u)):0ul);
    v2=(p2>=32?pstr2[p2/32u-1u]>>(2u*(p2%32u)):0ul);
    v1|=(pstr1[p1/32u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    v2|=(pstr2[p2/32u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
//    cout << p1 << " " << long2str(v1) << endl;
//    cout << p2 << " " << long2str(v2) << endl;
    for (k=IKMERSIZE; p1>=k && p2>=k && k<=32; k+=IKMERSIZE,v1<<=IKMERSIZE*2u,v2<<=IKMERSIZE*2u){
//      cout << p1 << "+" << k << " " << kmer2str((v1>>((32u-IKMERSIZE)*2u))&IKMERMASK) << endl; 
//      cout << p2 << "+" << k << " " << kmer2str((v2>>((32u-IKMERSIZE)*2u))&IKMERMASK) << endl; 
      tmp[0]=seq_ident_table[((v1^v2)>>(64u-IKMERSIZE*2u))&IKMERMASK];
      
//      if (!galign && tmp[0]<=IKMERSIZE/2){
      if (tmp[0]<=IKMERSIZE/2){
        adata.s1=p1-k;
        adata.s2=p2-k;
//        cout << "break p1: " << adata.s1 << " p2: " << adata.s2 << endl;
        return;
      }
      adata.matches+=tmp[0]; adata.mismatches+=IKMERSIZE-tmp[0];
//      sumcounts(tncount,p1-k,seq_alignment_lt[((v1^v2)>>(64u-IKMERSIZE*2u))&IKMERMASK]);
    }
  }

  long l=(p1<p2?p1:p2);
//  ldieif(p1<0 || p2<0,"negative sequence positions, p1: " + estr(p1) +" p2: "+p2);
//  cout << "p1: " << p1 << " p2: " << p2 << " l: " << l << endl;
  if (l<=0) { 
    adata.s1=p1-l;
    adata.s2=p2-l;
    return;
  }

  p1+=IKMERSIZE-l; p2+=IKMERSIZE-l; // adjust p1 and p2, such that we get a full kmer which includes the left end, we will mask the right side afterwards

  v1=(p1>=32?pstr1[p1/32u-1u]>>(2u*(p1%32u)):0ul);
  v2=(p2>=32?pstr2[p2/32u-1u]>>(2u*(p2%32u)):0ul);
  v1|=(pstr1[p1/32u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
  v2|=(pstr2[p2/32u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
//  cout << p1 << " " << long2str(v1) << endl;
//  cout << p2 << " " << long2str(v2) << endl;

  tmp[0]=seq_ident_table[((v1^v2)>>(62u-IKMERSIZE*2u)) & ((0x1ul<<(l*2ul))-1ul) & IKMERMASK];
//  if (!galign && tmp[0]<l/2){
  if (tmp[0]<l/2){
    adata.s1=p1-IKMERSIZE+l;
    adata.s2=p2-IKMERSIZE+l;
    return;
  }

//  cout << long2str(v1&((0x1ul<<(l*2u))-1ul)) << " " << long2str(((0x1ul<<(l*2u))-1ul)) << endl;
//  cout << long2str(v2&((0x1ul<<(l*2u))-1ul)) << " " << long2str(((0x1ul<<(l*2u))-1ul)) << endl;
  adata.matches+=tmp[0]-(IKMERSIZE-l);
  adata.mismatches+=IKMERSIZE-tmp[0];
  adata.s1=p1-IKMERSIZE;
  adata.s2=p2-IKMERSIZE;
//  sumcounts(tncount,p1-IKMERSIZE,seq_alignment_lt[((v1^v2)>>(64u-IKMERSIZE*2u)) & IKMERMASK]&((0x1ul<<(l*8ul))-1ul));
//  cout << "final p1: " << p1 << " p2: " << p2 << " l: " << l << endl;

/*
  tmp[0]=seq_ident_table[((v1^v2)>>((32u-IKMERSIZE)*2u))&~((0x1ul<<(IKMERSIZE*2u-l*2ul))-1ul)&IKMERMASK]; // put mask at end of 64bit int, i.e.: FF000000

  cout << long2str((v1>>((32u-IKMERSIZE)*2u))&(~((0x1ul<<(IKMERSIZE*2u-l*2ul))-1ul))&IKMERMASK) << " " << long2str((~((0x1ul<<(IKMERSIZE*2u-l*2ul))-1ul))&IKMERMASK) << endl;
  cout << long2str((v2>>((32u-IKMERSIZE)*2u))&~((0x1ul<<(IKMERSIZE*2u-l*2ul))-1ul)&IKMERMASK) << " " << long2str(~((0x1ul<<(IKMERSIZE*2u-l*2ul))-1ul)&IKMERMASK) << endl;
  matches+=tmp[0]-(IKMERSIZE-l);
  mismatches+=IKMERSIZE-tmp[0];
*/
}

void seqident_seg_local_right(const eseq& s1,long p1,long e1,const eseq& s2,long p2,long e2,ealigndata& adata,ealignscore& as){
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s1.seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(s2.seq._str);
  unsigned long v1,v2;
  unsigned char tmp[5];
  int k;
  int itmp;

//  cout << "id: " << p1 << "," << e1 << " :: " << p2 << "," <<e2 << " len1: " << e1-p1 << " len2: " << e2-p2 << endl;
  for (; p1<e1-32 && p2<e2-32; p1+=k,p2+=k){
    v1=pstr1[p1/32u]>>(2u*(p1%32u));
    v2=pstr2[p2/32u]>>(2u*(p2%32u));
    v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
    for (k=0; k<32u-IKMERSIZE; k+=IKMERSIZE,v1>>=IKMERSIZE*2u,v2>>=IKMERSIZE*2u){
      tmp[0]=seq_ident_table[(v1^v2)&IKMERMASK];
//      if (!galign && tmp[0]<IKMERSIZE/2){
      if (tmp[0]<IKMERSIZE/2){
        adata.e1=p1+k;
        adata.e2=p2+k;
        return;
      }
      adata.matches+=tmp[0]; adata.mismatches+=IKMERSIZE-tmp[0];
//      sumcounts(tncount,p1+k,seq_alignment_lt[(v1^v2)&IKMERMASK]);
    }
  }

//  cout << "rest id: " << p1 << "," << e1 << " :: " << p2 << "," <<e2 << " len1: " << e1-p1 << " len2: " << e2-p2 << endl;
  for (; p1+IKMERSIZE<e1 && p2+IKMERSIZE<e2; p1+=k,p2+=k){
    v1=pstr1[p1/32u]>>(2u*(p1%32u));
    v2=pstr2[p2/32u]>>(2u*(p2%32u));
    v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
//    cout << p1 << " " << long2str(v1) << endl;
//    cout << p2 << " " << long2str(v2) << endl;
    for (k=0; p1+k+IKMERSIZE<e1 && p2+k+IKMERSIZE<e2; k+=IKMERSIZE,v1>>=IKMERSIZE*2u,v2>>=IKMERSIZE*2u){
//      cout << p1 << "+" << k << " " << kmer2str(v1) << endl;
//      cout << p2 << "+" << k << " " << kmer2str(v2) << endl;
      tmp[0]=seq_ident_table[(v1^v2)&IKMERMASK];
     
//      if (!galign && tmp[0]<IKMERSIZE/2){
      if (tmp[0]<IKMERSIZE/2){
        adata.e1=p1+k;
        adata.e2=p2+k;
        return;
      }
      adata.matches+=tmp[0]; adata.mismatches+=IKMERSIZE-tmp[0];
//      sumcounts(tncount,p1+k,seq_alignment_lt[(v1^v2)&IKMERMASK]);
    }
  }
  long l=(e1-p1<e2-p2?e1-p1:e2-p2);
//  cout << "missing: " << p1 << " l: " << l << endl;
//  ldieif(l<0,"error past end of sequence");
  if (l<=0) {
    adata.e1=p1+l;
    adata.e2=p2+l;
    return;
  }
  v1=pstr1[p1/32u]>>(2u*(p1%32u));
  v2=pstr2[p2/32u]>>(2u*(p2%32u));
//  cout << p1 << " " << long2str(v1) << endl;
//  cout << p2 << " " << long2str(v2) << endl;
  v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
  v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
//  cout << p1 << " " << long2str(v1) << endl;
//  cout << p2 << " " << long2str(v2) << endl;
  tmp[0]=seq_ident_table[(v1^v2)&((0x1ul<<(l*2ul))-1ul)];
//  cout << long2str(v1&((0x1ul<<(l*2u))-1ul)) << " " << long2str(((0x1ul<<(l*2u))-1ul)) << endl;
//  cout << long2str(v2&((0x1ul<<(l*2u))-1ul)) << " " << long2str(((0x1ul<<(l*2u))-1ul)) << endl;

//  if (!galign && tmp[0]<(IKMERSIZE-l)/2) {
  if (tmp[0]<(IKMERSIZE-l)/2) {
    adata.e1=p1;
    adata.e2=p2;
    return;
  }

  adata.matches+=tmp[0]-(IKMERSIZE-l);
  adata.mismatches+=IKMERSIZE-tmp[0];
//  sumcounts(tncount,p1,seq_alignment_lt[(v1^v2)&IKMERMASK]&((0x1ul<<(l*8u))-1ul));
  adata.e1=p1+l;
  adata.e2=p2+l;
}


void seqident_seg_nogaps(const eseq& s1,long p1,long e1,const eseq& s2,long p2,long e2,ealigndata& adata,const ealignscore& as){
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s1.seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(s2.seq._str);
  unsigned long v1,v2;
  unsigned char tmp;
  int k;
  int itmp;

  ealignprofile tmpprof;

//  cout << "id: " << p1 << "," << e1 << " :: " << p2 << "," <<e2 << " len1: " << e1-p1 << " len2: " << e2-p2 << endl;
  for (; p1<e1-32 && p2<e2-32; p1+=k,p2+=k){
    v1=pstr1[p1/32u]>>(2u*(p1%32u));
    v2=pstr2[p2/32u]>>(2u*(p2%32u));
    v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
    for (k=0; k<32u-IKMERSIZE; k+=IKMERSIZE,v1>>=IKMERSIZE*2u,v2>>=IKMERSIZE*2u){
      tmp=seq_ident_table[(v1^v2)&IKMERMASK];
      adata.matches+=tmp;
      adata.mismatches+=IKMERSIZE-tmp;
      tmpprof.add(AT_MATCH,tmp);
      tmpprof.add(AT_MISS,IKMERSIZE-tmp);
      adata._score+=tmp*as.match-(IKMERSIZE-tmp)*as.mismatch;
    }
  }

  for (; p1+IKMERSIZE<e1 && p2+IKMERSIZE<e2; p1+=k,p2+=k){
    v1=pstr1[p1/32u]>>(2u*(p1%32u));
    v2=pstr2[p2/32u]>>(2u*(p2%32u));
    v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
    for (k=0; p1+k+IKMERSIZE<e1 && p2+k+IKMERSIZE<e2; k+=IKMERSIZE,v1>>=IKMERSIZE*2u,v2>>=IKMERSIZE*2u){
      tmp=seq_ident_table[(v1^v2)&IKMERMASK];
      adata.matches+=tmp;
      adata.mismatches+=IKMERSIZE-tmp;
      tmpprof.add(AT_MATCH,tmp);
      tmpprof.add(AT_MISS,IKMERSIZE-tmp);
      adata._score+=tmp*as.match-(IKMERSIZE-tmp)*as.mismatch;
    }
  }

  long l=(e1-p1<e2-p2?e1-p1:e2-p2);
//  long ml=(e1-p1>e2-p2?e1-p1:e2-p2);
  if (l<=0) return;
  v1=pstr1[p1/32u]>>(2u*(p1%32u));
  v2=pstr2[p2/32u]>>(2u*(p2%32u));
  v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
  v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
  tmp=seq_ident_table[(v1^v2)&((0x1ul<<(l*2ul))-1ul)];
  adata.matches+=tmp-(IKMERSIZE-l);
  adata.mismatches+=IKMERSIZE-tmp;
  adata._score+=(tmp-(IKMERSIZE-l))*as.match-(IKMERSIZE-tmp)*as.mismatch;
  tmpprof.add(AT_MATCH,tmp-(IKMERSIZE-l));
  tmpprof.add(AT_MISS,IKMERSIZE-tmp);
  adata.profile.addinv(tmpprof);
}


void seqident_seg(const eseq& s1,long p1,long e1,const eseq& s2,long p2,long e2,ealigndata& adata,const ealignscore& as){
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s1.seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(s2.seq._str);
  unsigned long v1,v2;
  unsigned char tmp[5];
  int k;
  int itmp;

//  cout << "id: " << p1 << "," << e1 << " :: " << p2 << "," <<e2 << " len1: " << e1-p1 << " len2: " << e2-p2 << endl;
  for (; p1<e1-32 && p2<e2-32; p1+=k,p2+=k){
    v1=pstr1[p1/32u]>>(2u*(p1%32u));
    v2=pstr2[p2/32u]>>(2u*(p2%32u));
    v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
    for (k=0; k<32u-IKMERSIZE; k+=IKMERSIZE,v1>>=IKMERSIZE*2u,v2>>=IKMERSIZE*2u){
      tmp[0]=seq_ident_table[(v1^v2)&IKMERMASK];
      tmp[1]=seq_ident_table[(v1^(v2>>2ul))&IKMERMASK];
      tmp[2]=seq_ident_table[(v1^(v2>>4ul))&IKMERMASK];
      tmp[3]=seq_ident_table[(v2^(v1>>2ul))&IKMERMASK];
      tmp[4]=seq_ident_table[(v2^(v1>>4ul))&IKMERMASK];
      
      itmp=0;
      for (int i=1; i<5; ++i){
        if (tmp[i]>tmp[itmp]) itmp=i;
      }
      if (itmp!=0 && tmp[itmp]>tmp[0]+IKMERSIZE/4){
        switch(itmp){
          case 1: p2+=1; adata.gaps+=1; break;
          case 2: p2+=2; adata.gaps+=2; break;
          case 3: p1+=1; adata.gaps+=1; break;
          case 4: p1+=2; adata.gaps+=2; break;
        }
        break;
      }else{
        adata.matches+=tmp[0]; adata.mismatches+=IKMERSIZE-tmp[0];
        adata._score+=tmp[0]*as.match-(IKMERSIZE-tmp[0])*as.mismatch;
      }
    }
  }

//  cout << "rest id: " << p1 << "," << e1 << " :: " << p2 << "," <<e2 << " len1: " << e1-p1 << " len2: " << e2-p2 << endl;
  for (; p1+IKMERSIZE<e1 && p2+IKMERSIZE<e2; p1+=k,p2+=k){
    v1=pstr1[p1/32u]>>(2u*(p1%32u));
    v2=pstr2[p2/32u]>>(2u*(p2%32u));
    v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
//    cout << p1 << " " << long2str(v1) << endl;
//    cout << p2 << " " << long2str(v2) << endl;
    for (k=0; p1+k+IKMERSIZE<e1 && p2+k+IKMERSIZE<e2; k+=IKMERSIZE,v1>>=IKMERSIZE*2u,v2>>=IKMERSIZE*2u){
//      cout << p1 << "+" << k << " " << kmer2str(v1) << endl;
//      cout << p2 << "+" << k << " " << kmer2str(v2) << endl;
      tmp[0]=seq_ident_table[(v1^v2)&IKMERMASK];
      tmp[1]=seq_ident_table[(v1^(v2>>2ul))&IKMERMASK];
      tmp[2]=seq_ident_table[(v1^(v2>>4ul))&IKMERMASK];
      tmp[3]=seq_ident_table[(v2^(v1>>2ul))&IKMERMASK];
      tmp[4]=seq_ident_table[(v2^(v1>>4ul))&IKMERMASK];
      
      itmp=0;
      for (int i=1; i<5; ++i){
        if (tmp[i]>tmp[itmp]) itmp=i;
      }
      if (itmp!=0 && tmp[itmp]>tmp[0]+IKMERSIZE/4){
        switch(itmp){
          case 1: p2+=1; adata.gaps+=1; break;
          case 2: p2+=2; adata.gaps+=2; break;
          case 3: p1+=1; adata.gaps+=1; break;
          case 4: p1+=2; adata.gaps+=2; break;
        }
        break;
      }else{
        adata.matches+=tmp[0]; adata.mismatches+=IKMERSIZE-tmp[0];
        adata._score+=tmp[0]*as.match-(IKMERSIZE-tmp[0])*as.mismatch;
      }
    }
  }
  long l=(e1-p1<e2-p2?e1-p1:e2-p2);
  long ml=(e1-p1>e2-p2?e1-p1:e2-p2);
  if (e1!=s1.seqlen && e2!=s2.seqlen)
    adata.gaps+=ml-l;
//  cout << "missing: " << p1 << " l: " << l << endl;
  if (l<=0) return;
  v1=pstr1[p1/32u]>>(2u*(p1%32u));
  v2=pstr2[p2/32u]>>(2u*(p2%32u));
//  cout << p1 << " " << long2str(v1) << endl;
//  cout << p2 << " " << long2str(v2) << endl;
  v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
  v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
//  cout << p1 << " " << long2str(v1) << endl;
//  cout << p2 << " " << long2str(v2) << endl;
  tmp[0]=seq_ident_table[(v1^v2)&((0x1ul<<(l*2ul))-1ul)];
//  cout << long2str(v1&((0x1ul<<(l*2u))-1ul)) << " " << long2str(((0x1ul<<(l*2u))-1ul)) << endl;
//  cout << long2str(v2&((0x1ul<<(l*2u))-1ul)) << " " << long2str(((0x1ul<<(l*2u))-1ul)) << endl;
  adata.matches+=tmp[0]-(IKMERSIZE-l);
  adata.mismatches+=IKMERSIZE-tmp[0];
  adata._score+=(tmp[0]-(IKMERSIZE-l))*as.match-(IKMERSIZE-tmp[0])*as.mismatch;
}



void seqident_seg(const eseq& s1,long p1,long e1,const eseq& s2,long p2,long e2,ealigndata& adata,unsigned char *tncount){
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s1.seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(s2.seq._str);
  unsigned long v1,v2;
  unsigned char tmp[5];
  int k;
  int itmp;

//  cout << "id: " << p1 << "," << e1 << " :: " << p2 << "," <<e2 << " len1: " << e1-p1 << " len2: " << e2-p2 << endl;
  for (; p1<e1-32 && p2<e2-32; p1+=k,p2+=k){
    v1=pstr1[p1/32u]>>(2u*(p1%32u));
    v2=pstr2[p2/32u]>>(2u*(p2%32u));
    v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
    for (k=0; k<32u-IKMERSIZE; k+=IKMERSIZE,v1>>=IKMERSIZE*2u,v2>>=IKMERSIZE*2u){
      tmp[0]=seq_ident_table[(v1^v2)&IKMERMASK];
      tmp[1]=seq_ident_table[(v1^(v2>>2ul))&IKMERMASK];
      tmp[2]=seq_ident_table[(v1^(v2>>4ul))&IKMERMASK];
      tmp[3]=seq_ident_table[(v2^(v1>>2ul))&IKMERMASK];
      tmp[4]=seq_ident_table[(v2^(v1>>4ul))&IKMERMASK];
      
      itmp=0;
      for (int i=1; i<5; ++i){
        if (tmp[i]>tmp[itmp]) itmp=i;
      }
      if (itmp!=0 && tmp[itmp]>tmp[0]+IKMERSIZE/4){
        switch(itmp){
          case 1: p2+=1; adata.gaps+=1; break;
          case 2: p2+=2; adata.gaps+=2; break;
          case 3: p1+=1; adata.gaps+=1; break;
          case 4: p1+=2; adata.gaps+=2; break;
        }
        break;
      }else{
        adata.matches+=tmp[0]; adata.mismatches+=IKMERSIZE-tmp[0];
        sumcounts(tncount,p1+k,seq_alignment_lt[(v1^v2)&IKMERMASK]);
      }
    }
  }

//  cout << "rest id: " << p1 << "," << e1 << " :: " << p2 << "," <<e2 << " len1: " << e1-p1 << " len2: " << e2-p2 << endl;
  for (; p1+IKMERSIZE<e1 && p2+IKMERSIZE<e2; p1+=k,p2+=k){
    v1=pstr1[p1/32u]>>(2u*(p1%32u));
    v2=pstr2[p2/32u]>>(2u*(p2%32u));
    v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
//    cout << p1 << " " << long2str(v1) << endl;
//    cout << p2 << " " << long2str(v2) << endl;
    for (k=0; p1+k+IKMERSIZE<e1 && p2+k+IKMERSIZE<e2; k+=IKMERSIZE,v1>>=IKMERSIZE*2u,v2>>=IKMERSIZE*2u){
//      cout << p1 << "+" << k << " " << kmer2str(v1) << endl;
//      cout << p2 << "+" << k << " " << kmer2str(v2) << endl;
      tmp[0]=seq_ident_table[(v1^v2)&IKMERMASK];
      tmp[1]=seq_ident_table[(v1^(v2>>2ul))&IKMERMASK];
      tmp[2]=seq_ident_table[(v1^(v2>>4ul))&IKMERMASK];
      tmp[3]=seq_ident_table[(v2^(v1>>2ul))&IKMERMASK];
      tmp[4]=seq_ident_table[(v2^(v1>>4ul))&IKMERMASK];
      
      itmp=0;
      for (int i=1; i<5; ++i){
        if (tmp[i]>tmp[itmp]) itmp=i;
      }
      if (itmp!=0 && tmp[itmp]>tmp[0]+IKMERSIZE/4){
        switch(itmp){
          case 1: p2+=1; adata.gaps+=1; break;
          case 2: p2+=2; adata.gaps+=2; break;
          case 3: p1+=1; adata.gaps+=1; break;
          case 4: p1+=2; adata.gaps+=2; break;
        }
        break;
      }else{
        adata.matches+=tmp[0]; adata.mismatches+=IKMERSIZE-tmp[0];
        sumcounts(tncount,p1+k,seq_alignment_lt[(v1^v2)&IKMERMASK]);
      }
    }
  }
  long l=(e1-p1<e2-p2?e1-p1:e2-p2);
  long ml=(e1-p1>e2-p2?e1-p1:e2-p2);
  if (e1!=s1.seqlen && e2!=s2.seqlen)
    adata.gaps+=ml-l;
//  cout << "missing: " << p1 << " l: " << l << endl;
  if (l<=0) return;
  v1=pstr1[p1/32u]>>(2u*(p1%32u));
  v2=pstr2[p2/32u]>>(2u*(p2%32u));
//  cout << p1 << " " << long2str(v1) << endl;
//  cout << p2 << " " << long2str(v2) << endl;
  v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
  v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
//  cout << p1 << " " << long2str(v1) << endl;
//  cout << p2 << " " << long2str(v2) << endl;
  tmp[0]=seq_ident_table[(v1^v2)&((0x1ul<<(l*2ul))-1ul)];
//  cout << long2str(v1&((0x1ul<<(l*2u))-1ul)) << " " << long2str(((0x1ul<<(l*2u))-1ul)) << endl;
//  cout << long2str(v2&((0x1ul<<(l*2u))-1ul)) << " " << long2str(((0x1ul<<(l*2u))-1ul)) << endl;
  adata.matches+=tmp[0]-(IKMERSIZE-l);
  adata.mismatches+=IKMERSIZE-tmp[0];
  sumcounts(tncount,p1,seq_alignment_lt[(v1^v2)&IKMERMASK]&((0x1ul<<(l*8u))-1ul));
}

void seqident_rev_seg(const eseq& s1,long p1,long e1,const eseq& s2,long p2,long e2,ealigndata& adata){
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s1.seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(s2.seq._str);
  unsigned long v1,v2;
  unsigned char tmp[5];
  int k;
  int itmp;

//  cout << "id: " << p1 << "," << e1 << " :: " << p2 << "," <<e2 << " len1: " << e1-p1 << " len2: " << e2-p2 << endl;
  for (; p1<e1-32 && p2<e2-32; p1+=k,p2+=k){
    long rp2=s2.seqlen-p2-32;
    v1=pstr1[p1/32u]>>(2u*(p1%32u));
    v2=pstr2[rp2/32u]>>(2u*(rp2%32u));
    v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    v2|=(pstr2[rp2/32u+1u]<<(64u-2u*(rp2%32u)))&safe_shift[rp2%32u];
    for (k=0; k<32u-IKMERSIZE; k+=IKMERSIZE,v1>>=IKMERSIZE*2u,v2<<=IKMERSIZE*2u){
      tmp[0]=seq_ident_table[(v1^kmer_rev_lt[(v2>>((32u-IKMERSIZE)*2u))&IKMERMASK])&IKMERMASK];
      tmp[1]=seq_ident_table[(v1^kmer_rev_lt[(v2>>((32u-IKMERSIZE)*2u-2u))&IKMERMASK])&IKMERMASK];
      tmp[2]=seq_ident_table[(v1^kmer_rev_lt[(v2>>((32u-IKMERSIZE)*2u-4u))&IKMERMASK])&IKMERMASK];
      tmp[3]=seq_ident_table[(kmer_rev_lt[(v2>>((32u-IKMERSIZE)*2u))&IKMERMASK]^(v1>>2ul))&IKMERMASK];
      tmp[4]=seq_ident_table[(kmer_rev_lt[(v2>>((32u-IKMERSIZE)*2u))&IKMERMASK]^(v1>>4ul))&IKMERMASK];
      
      itmp=0;
      for (int i=1; i<5; ++i){
        if (tmp[i]>tmp[itmp]) itmp=i;
      }
      if (itmp!=0 && tmp[itmp]>tmp[0]+IKMERSIZE/4){
        switch(itmp){
          case 1: p2+=1; adata.gaps+=1; break;
          case 2: p2+=2; adata.gaps+=2; break;
          case 3: p1+=1; adata.gaps+=1; break;
          case 4: p1+=2; adata.gaps+=2; break;
        }
        break;
      }else{
        adata.matches+=tmp[0]; adata.mismatches+=IKMERSIZE-tmp[0];
      }
    }
  }

//  cout << "rest id: " << p1 << "," << e1 << " :: " << p2 << "," <<e2 << " len1: " << e1-p1 << " len2: " << e2-p2 << endl;
  for (; p1+IKMERSIZE<e1 && p2+IKMERSIZE<e2; p1+=k,p2+=k){
    v1=pstr1[p1/32u]>>(2u*(p1%32u));
    v2=pstr2[p2/32u]>>(2u*(p2%32u));
    v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
//    cout << p1 << " " << long2str(v1) << endl;
//    cout << p2 << " " << long2str(v2) << endl;
    for (k=0; p1+k+IKMERSIZE<e1 && p2+k+IKMERSIZE<e2; k+=IKMERSIZE,v1>>=IKMERSIZE*2u,v2>>=IKMERSIZE*2u){
//      cout << p1 << "+" << k << " " << kmer2str(v1) << endl;
//      cout << p2 << "+" << k << " " << kmer2str(v2) << endl;
      tmp[0]=seq_ident_table[(v1^v2)&IKMERMASK];
      tmp[1]=seq_ident_table[(v1^(v2>>2ul))&IKMERMASK];
      tmp[2]=seq_ident_table[(v1^(v2>>4ul))&IKMERMASK];
      tmp[3]=seq_ident_table[(v2^(v1>>2ul))&IKMERMASK];
      tmp[4]=seq_ident_table[(v2^(v1>>4ul))&IKMERMASK];
      
      itmp=0;
      for (int i=1; i<5; ++i){
        if (tmp[i]>tmp[itmp]) itmp=i;
      }
      if (itmp!=0 && tmp[itmp]>tmp[0]+IKMERSIZE/4){
        switch(itmp){
          case 1: p2+=1; adata.gaps+=1; break;
          case 2: p2+=2; adata.gaps+=2; break;
          case 3: p1+=1; adata.gaps+=1; break;
          case 4: p1+=2; adata.gaps+=2; break;
        }
        break;
      }else{
        adata.matches+=tmp[0]; adata.mismatches+=IKMERSIZE-tmp[0];
      }
    }
  }
  long l=(e1-p1<e2-p2?e1-p1:e2-p2);
//  cout << "missing: " << p1 << " l: " << l << endl;
  if (l<=0) return;
  v1=pstr1[p1/32u]>>(2u*(p1%32u));
  v2=pstr2[p2/32u]>>(2u*(p2%32u));
//  cout << p1 << " " << long2str(v1) << endl;
//  cout << p2 << " " << long2str(v2) << endl;
  v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
  v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
//  cout << p1 << " " << long2str(v1) << endl;
//  cout << p2 << " " << long2str(v2) << endl;
  tmp[0]=seq_ident_table[(v1^v2)&((0x1ul<<(l*2ul))-1ul)];
//  cout << long2str(v1&((0x1ul<<(l*2u))-1ul)) << " " << long2str(((0x1ul<<(l*2u))-1ul)) << endl;
//  cout << long2str(v2&((0x1ul<<(l*2u))-1ul)) << " " << long2str(((0x1ul<<(l*2u))-1ul)) << endl;
  adata.matches+=tmp[0]-(IKMERSIZE-l);
  adata.mismatches+=IKMERSIZE-tmp[0];
}



void seqident_rev(const eseq& s1,long offset1,eintarray& kmerpos1,const eseq& s2,long offset2,eintarray& kmerpos2,ealigndata& adata)
{
  long p1,p2;
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s1.seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(s2.seq._str);
  unsigned long v1,v2;
  long ndelta=-s1.seqlen-s2.seqlen;
  long lastkmerlen,lastkmerpos,lastndelta;

//  for (int i=0; i<11; ++i)
//    match[i]=0;
  adata.matches=0; adata.mismatches=0; adata.gaps=0;

  char tmp[5];
  char itmp;
  int k;

  ebasicarray<ediag> diags;

  if (s1.seqlen<s2.seqlen){
    lastkmerpos=-2*int(KMERSIZE);
    lastkmerlen=0;
    lastndelta=-s1.seqlen-s2.seqlen;
    for (p1=0; p1<s1.seqlen-32; p1+=k){
      v1=pstr1[p1/32u]>>(2u*(p1%32u));
      v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
      for (k=0; k<32u-KMERSIZE; ++k,v1>>=2u){
        int kpos2=kmerpos2[kmer_rev_lt[v1&KMERMASK]]; // -offset2;
        if (kpos2>offset2){
          kpos2=s2.seqlen-kpos2-offset2;
          int d=p1+k-lastkmerpos;
          if (d<=KMERSIZE && ((lastndelta==kpos2-p1-k) || (v1&KMERMASK)==seqrevkmer(s2,p1+k+lastndelta))){
            lastkmerlen+=d;
          }else{
            if (lastkmerlen){
//              cout << "p1: " << lastkmerpos-lastkmerlen+KMERSIZE<< "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
              diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE,lastndelta,lastkmerlen));
            }
            lastkmerlen=KMERSIZE;
            lastndelta=kpos2-p1-k;
          }
          lastkmerpos=p1+k;
        }
      }
    }
    if (lastkmerlen){
      diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE,lastndelta,lastkmerlen));
//      cout << "end p1: " << lastkmerpos-lastkmerlen+KMERSIZE << "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
    }
  }
/*
  else{
    lastkmerpos=-2*int(KMERSIZE);
    lastkmerlen=0;
    lastndelta=-s1.seqlen-s2.seqlen;
    for (p2=0; p2<s2.seqlen-32 && ndelta==-s1.seqlen-s2.seqlen; p2+=k){
      v2=pstr2[p2/32u]>>(2u*(p2%32u));
      v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
      for (k=0; k<32u-KMERSIZE; ++k,v2>>=2u){
        int kpos1=kmerpos1[kmer_rev_lt[v2&KMERMASK]]-offset1;
        if (kpos1>-1){
          int d=p2+k-lastkmerpos;
          if (d<2*KMERSIZE && ((lastndelta==kpos1-p2-k) || kmer_rev_lt[v2&KMERMASK]==seqkmer(s1,p2+k-lastndelta))){
            lastkmerlen+=d;
          }else{
            if (lastkmerlen){
//              cout << "p2: " << lastkmerpos-lastkmerlen+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
              diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE-lastndelta,lastndelta,lastkmerlen));
            }
            lastkmerlen=KMERSIZE;
            lastndelta=p2+k-kpos1;
          }
          lastkmerpos=p2+k;
        }
      }
    }
    if (lastkmerlen){
//      cout << "end p2: " << lastkmerpos-lastkmerlen+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
      diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE-lastndelta,lastndelta,lastkmerlen));
    }
  }
*/

  // sparse dynamic programming to find longest chain of segments
  multimap<int,ediag*> segments;
  for (int j=0; j<diags.size(); ++j){
    ediag& diag(diags.at(j));
    diag.V=diag.j-diag.i;
    segments.insert(pair<int,ediag*>(diag.i,&diag));
    segments.insert(pair<int,ediag*>(diag.j,&diag));
  }

  multimap<int,ediag*> Lsegments;
  multimap<int,ediag*>::iterator sit;
  multimap<int,ediag*>::iterator slb;
  ediag *bestseg=0x00;
  for (sit=segments.begin(); sit!=segments.end(); ++sit){
//      cout << "x: " << sit->first << " start: " << (sit->first==sit->second->i?"1":"0") << " i: " << sit->second->i << " i2: " << sit->second->i2 << " Lsegs: " << Lsegments.size() << endl;
    if (sit->first==sit->second->i) { // start of rectangle
      slb=Lsegments.upper_bound(sit->second->i2-1);
      if (Lsegments.size()==0 || slb==Lsegments.begin()) continue;
      --slb;
      sit->second->V=sit->second->j-sit->second->i+slb->second->V;
      sit->second->bestseg=slb->second;
//        cout << "linking to: " << slb->second->i << endl;
    } else {  // end of rectangle
      slb=Lsegments.lower_bound(sit->second->j2);
      if (slb==Lsegments.end() || sit->second->V > slb->second->V){
        while (slb!=Lsegments.end()){
          if (slb->second->j2>=sit->second->j2 && slb->second->V<=sit->second->V){
            Lsegments.erase(slb++);
            continue;
          }
          ++slb;
        }
        if (bestseg==0x00 || sit->second->V > bestseg->V)
          bestseg=sit->second;
        Lsegments.insert(pair<int,ediag*>(sit->second->j2,sit->second));
      }
    }
  }

//  seqident_seg_right(s1,s2,bestseg,matches,mismatches,gaps);
  if (bestseg==0x00) return;

//  cout << "s1: " << s1.seqlen << " s2: " << s2.seqlen << endl;

  seqident_rev_seg(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata);
//  cout << "rightid: " << bestseg->j << "," << s1.seqlen << " l: " << s1.seqlen-bestseg->j << " d: " << bestseg->i2-bestseg->i << " m: " << matches << " miss: " << mismatches << " gaps: " << gaps << endl;
  while (bestseg!=0x00 && bestseg->bestseg != 0x00){
    adata.matches+=bestseg->j-bestseg->i;
//    cout << "seg: "<< bestseg->i << "," << bestseg->j << " l: " << (bestseg->j-bestseg->i) << " m: " << matches << " miss: " << mismatches << " gaps: " << gaps << endl;
    seqident_rev_seg(s1,bestseg->bestseg->j,bestseg->i,s2,bestseg->bestseg->j2,bestseg->i2,adata);
//    cout << "segid: " << bestseg->i << "," << bestseg->j << " l: " << bestseg->j-bestseg->i << " d: " << bestseg->i2-bestseg->i << " m: " << matches << " miss: " << mismatches << " gaps: " << gaps << endl;
    bestseg=bestseg->bestseg;
  }
  adata.matches+=bestseg->j-bestseg->i;
//  cout << "seg: "<< bestseg->i << "," << bestseg->j << " l: " << (bestseg->j-bestseg->i) << " m: " << matches << " miss: " << mismatches << " gaps: " << gaps << endl;
//  seqident_rev_seg_left(s1,bestseg->i,s2,bestseg->i2,matches,mismatches,gaps);
//  cout << "leftid: "<< bestseg->i << " l: " << (bestseg->i<bestseg->i2?bestseg->i:bestseg->i2) << " d: " << bestseg->i2-bestseg->i << " m: " << matches << " miss: " << mismatches << " gaps: " << gaps << endl;

  return;
}

/*
void seqident(const eseq& s1,int offset1,eintarray& kmerpos1,const eseq& s2,int offset2,eintarray& kmerpos2,ealigndata& adata,unsigned char *tncounts)
{
  int p1,p2;
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s1.seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(s2.seq._str);
  unsigned long v1,v2;
  int ndelta=-s1.seqlen-s2.seqlen;
  int lastkmerlen,lastkmerpos,lastndelta;

//  for (int i=0; i<11; ++i)
//    match[i]=0;
  adata.matches=0; adata.mismatches=0; adata.gaps=0;
  adata.s1=-1; adata.e1=-1;
  adata.s2=-1; adata.e2=-1;

  char tmp[5];
  char itmp;
  int k;

  ebasicarray<ediag> diags;

  if (s1.seqlen<s2.seqlen){
    lastkmerpos=-2*int(KMERSIZE);
    lastkmerlen=0;
    lastndelta=-s1.seqlen-s2.seqlen;
    for (p1=0; p1<s1.seqlen-32; p1+=k){
      v1=pstr1[p1/32u]>>(2u*(p1%32u));
      v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
      for (k=0; k<32u-KMERSIZE; ++k,v1>>=2u){
        int kpos2=kmerpos2[v1&KMERMASK]-offset2;
        if (kpos2>-1){
          int d=p1+k-lastkmerpos;
          if (d<=KMERSIZE && ((lastndelta==kpos2-p1-k) || (v1&KMERMASK)==seqkmer(s2,p1+k+lastndelta))){
            lastkmerlen+=d;
          }else{
            if (lastkmerlen){
//              cout << "p1: " << lastkmerpos-lastkmerlen+KMERSIZE<< "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
              diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE,lastndelta,lastkmerlen));
            }
            lastkmerlen=KMERSIZE;
            lastndelta=kpos2-p1-k;
          }
          lastkmerpos=p1+k;
        }
      }
    }
    if (lastkmerlen){
      diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE,lastndelta,lastkmerlen));
//      cout << "end p1: " << lastkmerpos-lastkmerlen+KMERSIZE << "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
    }
  }else{
    lastkmerpos=-2*int(KMERSIZE);
    lastkmerlen=0;
    lastndelta=-s1.seqlen-s2.seqlen;
    for (p2=0; p2<s2.seqlen-32 && ndelta==-s1.seqlen-s2.seqlen; p2+=k){
      v2=pstr2[p2/32u]>>(2u*(p2%32u));
      v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
      for (k=0; k<32u-KMERSIZE; ++k,v2>>=2u){
        int kpos1=kmerpos1[v2&KMERMASK]-offset1;
        if (kpos1>-1){
          int d=p2+k-lastkmerpos;
          if (d<=KMERSIZE && ((lastndelta==kpos1-p2-k) || (v2&KMERMASK)==seqkmer(s1,p2+k-lastndelta))){
            lastkmerlen+=d;
          }else{
            if (lastkmerlen){
//              cout << "p2: " << lastkmerpos-lastkmerlen+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
              diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE-lastndelta,lastndelta,lastkmerlen));
            }
            lastkmerlen=KMERSIZE;
            lastndelta=p2+k-kpos1;
          }
          lastkmerpos=p2+k;
        }
      }
    }
    if (lastkmerlen){
//      cout << "end p2: " << lastkmerpos-lastkmerlen+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
      diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE-lastndelta,lastndelta,lastkmerlen));
    }
  }

  // sparse dynamic programming to find longest chain of segments
  multimap<int,ediag*> segments;
  for (int j=0; j<diags.size(); ++j){
    ediag& diag(diags.at(j));
    diag.V=diag.j-diag.i;
    segments.insert(pair<int,ediag*>(diag.i,&diag));
    segments.insert(pair<int,ediag*>(diag.j,&diag));
  }

  multimap<int,ediag*> Lsegments;
  multimap<int,ediag*>::iterator sit;
  multimap<int,ediag*>::iterator slb;
  ediag *bestseg=0x00;
  for (sit=segments.begin(); sit!=segments.end(); ++sit){
//      cout << "x: " << sit->first << " start: " << (sit->first==sit->second->i?"1":"0") << " i: " << sit->second->i << " i2: " << sit->second->i2 << " Lsegs: " << Lsegments.size() << endl;
    if (sit->first==sit->second->i) { // start of rectangle
      slb=Lsegments.upper_bound(sit->second->i2-1);
      if (Lsegments.size()==0 || slb==Lsegments.begin()) continue;
      --slb;
      sit->second->V=sit->second->j-sit->second->i+slb->second->V;
      sit->second->bestseg=slb->second;
//        cout << "linking to: " << slb->second->i << endl;
    } else {  // end of rectangle
      slb=Lsegments.lower_bound(sit->second->j2);
      if (slb==Lsegments.end() || sit->second->V > slb->second->V){
        while (slb!=Lsegments.end()){
          if (slb->second->j2>=sit->second->j2 && slb->second->V<=sit->second->V){
            Lsegments.erase(slb++);
            continue;
          }
          ++slb;
        }
        if (bestseg==0x00 || sit->second->V > bestseg->V)
          bestseg=sit->second;
        Lsegments.insert(pair<int,ediag*>(sit->second->j2,sit->second));
      }
    }
  }

//  seqident_seg_right(s1,s2,bestseg,matches,mismatches,gaps);
  if (bestseg==0x00) return;

//  cout << "s1: " << s1.seqlen << " s2: " << s2.seqlen << endl;

  int l=(s2.seqlen-bestseg->j2<s1.seqlen-bestseg->j?s2.seqlen-bestseg->j2:s1.seqlen-bestseg->j);
  adata.e1=bestseg->j+l;
  adata.e2=bestseg->j2+l;

  estr ps1,ps2;

  seqident_seg(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,tncounts);
//  cout << "rightid: " << bestseg->j << "," << s1.seqlen << " l: " << s1.seqlen-bestseg->j << " d: " << bestseg->i2-bestseg->i << " m: " << matches << " miss: " << mismatches << " gaps: " << gaps << endl;
  while (bestseg!=0x00 && bestseg->bestseg != 0x00){
    adata.matches+=bestseg->j-bestseg->i;
    seqncount(s1,bestseg->i,bestseg->j,s2,bestseg->i2,bestseg->j2,tncounts);
//    cout << "seg: "<< bestseg->i << "," << bestseg->j << " l: " << (bestseg->j-bestseg->i) << " m: " << matches << " miss: " << mismatches << " gaps: " << gaps << endl;
    seqident_seg(s1,bestseg->bestseg->j,bestseg->i,s2,bestseg->bestseg->j2,bestseg->i2,adata,tncounts);
//    cout << "segid: " << bestseg->i << "," << bestseg->j << " l: " << bestseg->j-bestseg->i << " d: " << bestseg->i2-bestseg->i << " m: " << matches << " miss: " << mismatches << " gaps: " << gaps << endl;
    bestseg=bestseg->bestseg;
  }
  adata.matches+=bestseg->j-bestseg->i;
  seqncount(s1,bestseg->i,bestseg->j,s2,bestseg->i2,bestseg->j2,tncounts);
//  cout << "seg: "<< bestseg->i << "," << bestseg->j << " l: " << (bestseg->j-bestseg->i) << " m: " << matches << " miss: " << mismatches << " gaps: " << gaps << endl;
  seqident_seg_left(s1,bestseg->i,s2,bestseg->i2,adata);
//  cout << "leftid: "<< bestseg->i << " l: " << (bestseg->i<bestseg->i2?bestseg->i:bestseg->i2) << " d: " << bestseg->i2-bestseg->i << " m: " << matches << " miss: " << mismatches << " gaps: " << gaps << endl;
  l=(bestseg->i2<bestseg->i?bestseg->i2:bestseg->i);
  adata.s1=bestseg->i-l;
  adata.s2=bestseg->i2-l;

  return;
}
*/

/*
void seqident_local_fast(const eseq& s1,eintarray& kmerpos1,const eseq& s2,ealigndata& adata,esearchws& sws,const ealignscore& as)
{
  long p1,p2;
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s1.seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(s2.seq._str);
  unsigned long v1,v2;
  long lastkmerlen,lastkmerpos,lastndelta;

//  for (int i=0; i<11; ++i)
//    match[i]=0;
  adata.matches=0; adata.mismatches=0; adata.gaps=0; adata._score=0.0;
  adata.s1=-1; adata.e1=-1;
  adata.s2=-1; adata.e2=-1;

  char tmp[5];
  char itmp;
  int k;

  ebasicarray<ediag> diags;
  ebasicarray<ediag> dheap;

  if (s1.seqlen<s2.seqlen){
    lastkmerpos=-2*long(KMERSIZE);
    lastkmerlen=0;
    lastndelta=-s1.seqlen-s2.seqlen;
    for (p1=0; p1<s1.seqlen-32; p1+=k){
      v1=pstr1[p1/32u]>>(2u*(p1%32u));
      v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
      for (k=0; k<32u-KMERSIZE; ++k,v1>>=2u){
        long kpos2=sws.kmerpos2[v1&KMERMASK]-sws.offset2;
        if (kpos2>-1){
          long d=p1+k-lastkmerpos;
          
          if (d<=KMERSIZE && ((lastndelta==kpos2-p1-k) || (p1+k+lastndelta+KMERSIZE<=s2.seqlen && (v1&KMERMASK)==seqkmer(s2,p1+k+lastndelta)))){
            lastkmerlen+=d;
          }else{
            if (lastkmerlen){
//              cout << "p1: " << lastkmerpos-lastkmerlen+KMERSIZE<< "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
              diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE,lastndelta,lastkmerlen));
            }
            lastkmerlen=KMERSIZE;
            lastndelta=kpos2-p1-k;
          }
          lastkmerpos=p1+k;
        }
      }
    }
    if (lastkmerlen){
      diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE,lastndelta,lastkmerlen));
//      cout << "end p1: " << lastkmerpos-lastkmerlen+KMERSIZE << "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
    }
  }else{
    lastkmerpos=-2*long(KMERSIZE);
    lastkmerlen=0;
    lastndelta=-s1.seqlen-s2.seqlen;
    for (p2=0; p2<s2.seqlen-32; p2+=k){
      v2=pstr2[p2/32u]>>(2u*(p2%32u));
      v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
      for (k=0; k<32u-KMERSIZE; ++k,v2>>=2u){
        long kpos1=kmerpos1[v2&KMERMASK]-sws.offset;
        if (kpos1>-1){
          long d=p2+k-lastkmerpos;
          if (d<=KMERSIZE && ((lastndelta==kpos1-p2-k) || (p2+k-lastndelta+KMERSIZE<=s1.seqlen && (v2&KMERMASK)==seqkmer(s1,p2+k-lastndelta)))){
            lastkmerlen+=d;
          }else{
            if (lastkmerlen){
//              cout << "p2: " << lastkmerpos-lastkmerlen+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
              diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE-lastndelta,lastndelta,lastkmerlen));
            }
            lastkmerlen=KMERSIZE;
            lastndelta=p2+k-kpos1;
          }
          lastkmerpos=p2+k;
        }
      }
    }
    if (lastkmerlen){
//      cout << "end p2: " << lastkmerpos-lastkmerlen+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
      diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE-lastndelta,lastndelta,lastkmerlen));
    }
  }

  // sparse dynamic programming to find longest chain of segments
  multimap<long,ediag*> segments;
  for (long j=0; j<diags.size(); ++j){
    ediag& diag(diags.at(j));
    diag.V=(diag.j-diag.i)*as.match;
    segments.insert(pair<long,ediag*>(diag.i,&diag));
    segments.insert(pair<long,ediag*>(diag.j,&diag));
  }

  multimap<long,ediag*> Lsegments;
  multimap<long,ediag*>::iterator sit;
  multimap<long,ediag*>::iterator slb;
  ediag *bestseg=0x00;
  for (sit=segments.begin(); sit!=segments.end(); ++sit){
//      cout << "x: " << sit->first << " start: " << (sit->first==sit->second->i?"1":"0") << " i: " << sit->second->i << " i2: " << sit->second->i2 << " Lsegs: " << Lsegments.size() << endl;
    if (sit->first==sit->second->i) { // start of rectangle
      slb=Lsegments.upper_bound(sit->second->i2-1);
      if (Lsegments.size()==0 || slb==Lsegments.begin()) continue;
      --slb;
//      sit->second->V=sit->second->j-sit->second->i+slb->second->V;
      int gapdiff=abs((sit->second->i-sit->second->i2)-(slb->second->i-slb->second->i2));
      sit->second->V=(sit->second->j-sit->second->i)*as.match+gapdiff*as.gapext+(gapdiff>0?as.gapopen:0)+ MIN(sit->second->i-slb->second->j,sit->second->i2-slb->second->j2)*(0.75*as.match+0.25*as.mismatch) +slb->second->V;
      sit->second->bestseg=slb->second;
//      cout << "linking to: " << sit->second->i << " " << slb->second->j << endl;
    } else {  // end of rectangle
      slb=Lsegments.lower_bound(sit->second->j2);
      if (slb==Lsegments.end() || sit->second->V > slb->second->V){
        while (slb!=Lsegments.end()){
          if (slb->second->j2>=sit->second->j2 && slb->second->V<=sit->second->V){
            Lsegments.erase(slb++);
            continue;
          }
          ++slb;
        }
        if (bestseg==0x00 || sit->second->V > bestseg->V)
          bestseg=sit->second;
        Lsegments.insert(pair<long,ediag*>(sit->second->j2,sit->second));
      }
    }
  }

  if (bestseg==0x00) return;

//  estr fas1,fas2;
//  cout << "full alignment" << endl;
//  cout << "ali: " << 0 <<","<<s2.seqlen<<" : " << 0 <<","<<s1.seqlen << endl; 
//  adata._score=seqcalign_global_noedgegap(s1,0,s1.seqlen,s2,0,s2.seqlen,fas1,fas2);
//  print_seqali(fas1,fas2);

//  estr pas1,pas2;
//  seqident_seg_local(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,tncounts);
//  seqident_seg_sw_right(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,tncounts);
//  estr as1,as2;
//  seqcalign_global_noedgegap(s1,0,s1.seqlen,s2,0,s2.seqlen,as1,as2);
//  print_seqali(as1,as2);

//  seqcalign_global_noedgegap(s1,0,s1.seqlen,s2,0,s2.seqlen,adata,alignws);
//  return;

  ldieif(s2.seqlen<bestseg->j2,"segment start larger than sequence length: "+estr(s2.seqlen)+" < "+estr(bestseg->j2));
  long minright=MIN(s1.seqlen-bestseg->j,s2.seqlen-bestseg->j2);
//  seqcalign_global_norightedgegap(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,alignws);

//  seqcalign_global_norightedgegap(s1,bestseg->j,MIN(s1.seqlen,bestseg->j+2*minright),s2,bestseg->j2,MIN(s2.seqlen,bestseg->j2+2*minright),adata,sws.alignws,as);
  seqident_seg(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,as);
  LDEBUG(D_SEQIDENT,cout << "rightid: " << bestseg->j << "," << s1.seqlen << " l: " << s1.seqlen-bestseg->j << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
  LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,0,s1.seqlen));
  while (bestseg!=0x00 && bestseg->bestseg != 0x00){
    seqncount(s1,bestseg->i,bestseg->j,s2,bestseg->i2,bestseg->j2,adata,as);
    LDEBUG(D_SEQIDENT,cout << "seg: "<< bestseg->i << "," << bestseg->j << " l: " << (bestseg->j-bestseg->i) << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
    LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,0,s1.seqlen));
    long gapdiff=abs((bestseg->bestseg->i-bestseg->bestseg->i2)-(bestseg->i-bestseg->i2));
    if (gapdiff==0)
      seqident_seg_nogaps(s1,bestseg->bestseg->j,bestseg->i,s2,bestseg->bestseg->j2,bestseg->i2,adata,as);
    else
      seqident_seg(s1,bestseg->bestseg->j,bestseg->i,s2,bestseg->bestseg->j2,bestseg->i2,adata,as);
//      seqcalign_global(s1,bestseg->bestseg->j,bestseg->i,s2,bestseg->bestseg->j2,bestseg->i2,adata,sws.alignws,as);
    LDEBUG(D_SEQIDENT,cout << "segid: " << bestseg->bestseg->j << "," << bestseg->i << " l: " << bestseg->i-bestseg->bestseg->j << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
    LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,0,s1.seqlen));
//    cout << "segid: " << bestseg->i << "," << bestseg->j << " l: " << bestseg->j-bestseg->i << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl;
    bestseg=bestseg->bestseg;
  }
  seqncount(s1,bestseg->i,bestseg->j,s2,bestseg->i2,bestseg->j2,adata,as);
  LDEBUG(D_SEQIDENT,cout << "seg: "<< bestseg->i << "," << bestseg->j << " l: " << (bestseg->j-bestseg->i) << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
  LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,0,s1.seqlen));
//  seqident_seg_left_local(s1,bestseg->i,s2,bestseg->i2,adata,tncounts);
  long minleft=MIN(bestseg->i,bestseg->i2);
  seqident_seg_left(s1,bestseg->i,s2,bestseg->i2,adata,as);
//  seqcalign_global_noleftedgegap(s1,MAX(0,bestseg->i-2*minleft),bestseg->i,s2,MAX(0,bestseg->i2-2*minleft),bestseg->i2,adata,sws.alignws,as);

//  estr as1,as2;
//  seqcalign_global_noleftedgegap(s1,0,bestseg->i,s2,0,bestseg->i2,as1,as2);
//  cout << endl;
//  cout << as1 << endl << as2 << endl;
  LDEBUG(D_SEQIDENT,cout << "leftid: "<< bestseg->i << " l: " << (bestseg->i<bestseg->i2?bestseg->i:bestseg->i2) << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
  LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,0,s1.seqlen));
//  cout << "partial alignment" << endl;
//  print_seqali(pas1,pas2);

  return;
}
*/

/*
void seqident_global(const eseq& s1,eintarray& kmerpos1,const eseq& s2,ealigndata& adata,esearchws& sws,const ealignscore& as)
{
  long p1,p2;
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s1.seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(s2.seq._str);
  unsigned long v1,v2;
  long lastkmerlen,lastkmerpos,lastndelta;

//  for (int i=0; i<11; ++i)
//    match[i]=0;
  adata.matches=0; adata.mismatches=0; adata.gaps=0; adata._score=0.0;
  adata.s1=-1; adata.e1=-1;
  adata.s2=-1; adata.e2=-1;

  char tmp[5];
  char itmp;
  int k;

  ebasicarray<ediag> diags;
  ebasicarray<ediag> dheap;

  if (s1.seqlen<s2.seqlen){
    lastkmerpos=-2*long(KMERSIZE);
    lastkmerlen=0;
    lastndelta=-s1.seqlen-s2.seqlen;
    for (p1=0; p1<s1.seqlen-32; p1+=k){
      v1=pstr1[p1/32u]>>(2u*(p1%32u));
      v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
      for (k=0; k<32u-KMERSIZE; ++k,v1>>=2u){
        long kpos2=sws.kmerpos2[v1&KMERMASK]-sws.offset2;
        if (kpos2>-1){
          long d=p1+k-lastkmerpos;
          
          if (d<=KMERSIZE && ((lastndelta==kpos2-p1-k) || (p1+k+lastndelta+KMERSIZE<=s2.seqlen && (v1&KMERMASK)==seqkmer(s2,p1+k+lastndelta)))){
            lastkmerlen+=d;
          }else{
            if (lastkmerlen){
//              cout << "p1: " << lastkmerpos-lastkmerlen+KMERSIZE<< "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
              diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE,lastndelta,lastkmerlen));
            }
            lastkmerlen=KMERSIZE;
            lastndelta=kpos2-p1-k;
          }
          lastkmerpos=p1+k;
        }
      }
    }
    v1=pstr1[p1/32u]>>(2u*(p1%32u));
    v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    for (k=0; p1+k+KMERSIZE<s1.seqlen; ++k,v1>>=2u){
      long kpos2=sws.kmerpos2[v1&KMERMASK]-sws.offset2;
      if (kpos2>-1){
        long d=p1+k-lastkmerpos;
        if (d<=KMERSIZE && ((lastndelta==kpos2-p1-k) || (p1+k+lastndelta+KMERSIZE<=s2.seqlen && (v1&KMERMASK)==seqkmer(s2,p1+k+lastndelta)))){
          lastkmerlen+=d;
        }else{
          if (lastkmerlen){
            diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE,lastndelta,lastkmerlen));
          }
          lastkmerlen=KMERSIZE;
          lastndelta=kpos2-p1-k;
        }
        lastkmerpos=p1+k;
      }
    }
    if (lastkmerlen){
      diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE,lastndelta,lastkmerlen));
//      cout << "end p1: " << lastkmerpos-lastkmerlen+KMERSIZE << "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
    }
  }else{
    lastkmerpos=-2*long(KMERSIZE);
    lastkmerlen=0;
    lastndelta=-s1.seqlen-s2.seqlen;
    for (p2=0; p2<s2.seqlen-32; p2+=k){
      v2=pstr2[p2/32u]>>(2u*(p2%32u));
      v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
      for (k=0; k<32u-KMERSIZE; ++k,v2>>=2u){
        long kpos1=kmerpos1[v2&KMERMASK]-sws.offset;
        if (kpos1>-1){
          long d=p2+k-lastkmerpos;
          if (d<=KMERSIZE && ((lastndelta==kpos1-p2-k) || (p2+k-lastndelta+KMERSIZE<=s1.seqlen && (v2&KMERMASK)==seqkmer(s1,p2+k-lastndelta)))){
            lastkmerlen+=d;
          }else{
            if (lastkmerlen){
//              cout << "p2: " << lastkmerpos-lastkmerlen+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
              diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE-lastndelta,lastndelta,lastkmerlen));
            }
            lastkmerlen=KMERSIZE;
            lastndelta=p2+k-kpos1;
          }
          lastkmerpos=p2+k;
        }
      }
    }
    v2=pstr2[p2/32u]>>(2u*(p2%32u));
    v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
    for (k=0; p2+k+KMERSIZE<s2.seqlen; ++k,v2>>=2u){
      long kpos1=kmerpos1[v2&KMERMASK]-sws.offset;
      if (kpos1>-1){
        long d=p2+k-lastkmerpos;
        if (d<=KMERSIZE && ((lastndelta==kpos1-p2-k) || (p2+k-lastndelta+KMERSIZE<=s1.seqlen && (v2&KMERMASK)==seqkmer(s1,p2+k-lastndelta)))){
          lastkmerlen+=d;
        }else{
          if (lastkmerlen){
            diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE-lastndelta,lastndelta,lastkmerlen));
          }
          lastkmerlen=KMERSIZE;
          lastndelta=p2+k-kpos1;
        }
        lastkmerpos=p2+k;
      }
    }
    if (lastkmerlen){
//      cout << "end p2: " << lastkmerpos-lastkmerlen+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
      diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE-lastndelta,lastndelta,lastkmerlen));
    }
  }

  if (diags.size()==0){
//    lderror("no shared segments found");
    return;
  }

  // sparse dynamic programming to find longest chain of segments
  ediag *bestseg=&diags[0];
  multimap<long,ediag*> segments;
  for (long j=0; j<diags.size(); ++j){
    ediag& diag(diags.at(j));
    diag.V=(diag.j-diag.i)*as.match;
    if (diag.V > bestseg->V) bestseg=&diag;
    segments.insert(pair<long,ediag*>(diag.i,&diag));
    segments.insert(pair<long,ediag*>(diag.j,&diag));
  }

  multimap<long,ediag*> Lsegments;
  multimap<long,ediag*>::iterator sit;
  multimap<long,ediag*>::iterator slb;
  for (sit=segments.begin(); sit!=segments.end(); ++sit){
//      cout << "x: " << sit->first << " start: " << (sit->first==sit->second->i?"1":"0") << " i: " << sit->second->i << " i2: " << sit->second->i2 << " Lsegs: " << Lsegments.size() << endl;
    if (sit->first==sit->second->i) { // start of rectangle
//      if (Lsegments.size()==0) continue;
//      slb=Lsegments.upper_bound(sit->second->i2-1);
//      if (slb==Lsegments.end()) continue;
      ediag *tmpbest=0x00;
      double tmpbestscore=(sit->second->j-sit->second->i)*as.match;
      for (slb=Lsegments.begin(); slb!=Lsegments.end(); ++slb){
        if (slb->first >= sit->second->i2) break;
//      sit->second->V=sit->second->j-sit->second->i+slb->second->V;
        long gapdiff=abs((sit->second->i-sit->second->i2)-(slb->second->i-slb->second->i2));
        double tmpscore=(sit->second->j-sit->second->i)*as.match - gapdiff*as.gapext-(gapdiff>0?as.gapopen:0) + slb->second->V;
        if (tmpscore>tmpbestscore) { tmpbestscore=tmpscore; tmpbest=slb->second; }
      }
      sit->second->bestseg=tmpbest;
      sit->second->V=tmpbestscore;
      
//      sit->second->V=(sit->second->j-sit->second->i)*as.match-gapdiff*as.gapext-(gapdiff>0?as.gapopen:0)+ MIN(sit->second->j-slb->second->i,sit->second->j2-slb->second->i2)*(9.0/10.0*as.match-1.0/10.0*as.mismatch) +slb->second->V;
*/
/*
      sit->second->bestseg=0x00;
      double newscore=-gapdiff*as.gapext-(gapdiff>0?as.gapopen:0)+ MIN(sit->second->j-slb->second->i,sit->second->j2-slb->second->i2)*(7.0/8.0*as.match-1.0/8.0*as.mismatch) +slb->second->V;
      if (sit->second->V+MIN(sit->second->i,sit->second->j)*(7.0/8.0*as.match-1.0/8.0*as.mismatch) < newscore){
        sit->second->V=newscore;
        sit->second->bestseg=slb->second;
      }else{
        sit->second->V+=MIN(sit->second->i,sit->second->j)*(7.0/8.0*as.match-1.0/8.0*as.mismatch);
      }
*/
/*
//        double newscore=(sit->second->j-sit->second->i)*as.match - gapdiff*as.gapext-(gapdiff>0?as.gapopen:0) + slb->second->V;
//      if (sit->second->V < newscore){
//        sit->second->V=newscore;
//        sit->second->bestseg=slb->second;
//      }
//      cout << "linking to: " << sit->second->i << " " << slb->second->j << endl;
    } else {  // end of rectangle
      slb=Lsegments.lower_bound(sit->second->j2);
      if (slb==Lsegments.end() || sit->second->V > slb->second->V){
        while (slb!=Lsegments.end()){
          if (slb->second->j2>=sit->second->j2 && slb->second->V<=sit->second->V){
            Lsegments.erase(slb++);
            continue;
          }
          ++slb;
        }
        if (bestseg==0x00 || sit->second->V > bestseg->V)
          bestseg=sit->second;
        Lsegments.insert(pair<long,ediag*>(sit->second->j2,sit->second));
      }
    }
  }

//  estr fas1,fas2;
//  cout << "full alignment" << endl;
//  cout << "ali: " << 0 <<","<<s2.seqlen<<" : " << 0 <<","<<s1.seqlen << endl; 
//  adata._score=seqcalign_global_noedgegap(s1,0,s1.seqlen,s2,0,s2.seqlen,fas1,fas2);
//  print_seqali(fas1,fas2);

//  estr pas1,pas2;
//  seqident_seg_local(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,tncounts);
//  seqident_seg_sw_right(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,tncounts);
//  estr as1,as2;
//  seqcalign_global_noedgegap(s1,0,s1.seqlen,s2,0,s2.seqlen,as1,as2);
//  print_seqali(as1,as2);
*/
/*  
  if (bestseg->bestseg==0x00){ // single segment hit, do full alignment
    seqcalign_global_noedgegap(s1,0,s1.seqlen,s2,0,s2.seqlen,adata,sws.alignws,as);
    return;
  }
*/
/*
  ldieif(s2.seqlen<bestseg->j2,"segment start larger than sequence length: "+estr(s2.seqlen)+" < "+estr(bestseg->j2));
  long minright=MIN(s1.seqlen-bestseg->j,s2.seqlen-bestseg->j2);
//  seqcalign_global_norightedgegap(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,alignws);
  seqcalign_global_norightedgegap(s1,bestseg->j,MIN(s1.seqlen,bestseg->j+2*minright),s2,bestseg->j2,MIN(s2.seqlen,bestseg->j2+2*minright),adata,sws.alignws,as);
  LDEBUG(D_SEQIDENT,cout << "rightid: " << bestseg->j << "," << s1.seqlen << " l: " << s1.seqlen-bestseg->j << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
  LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,0,s1.seqlen));
  while (bestseg!=0x00 && bestseg->bestseg != 0x00){
    seqncount(s1,bestseg->i,bestseg->j,s2,bestseg->i2,bestseg->j2,adata,as);
    LDEBUG(D_SEQIDENT,cout << "seg: "<< bestseg->i << "," << bestseg->j << " l: " << (bestseg->j-bestseg->i) << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
    LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,0,s1.seqlen));
    long gapdiff=abs((bestseg->bestseg->i-bestseg->bestseg->i2)-(bestseg->i-bestseg->i2));
    if (gapdiff==0)
      seqident_seg_nogaps(s1,bestseg->bestseg->j,bestseg->i,s2,bestseg->bestseg->j2,bestseg->i2,adata,as);
    else
      seqcalign_global(s1,bestseg->bestseg->j,bestseg->i,s2,bestseg->bestseg->j2,bestseg->i2,adata,sws.alignws,as);
    LDEBUG(D_SEQIDENT,cout << "segid: " << bestseg->bestseg->j << "," << bestseg->i << " l: " << bestseg->i-bestseg->bestseg->j << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << " gapdiff: " << gapdiff << endl);
    LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,0,s1.seqlen));
//    cout << "segid: " << bestseg->i << "," << bestseg->j << " l: " << bestseg->j-bestseg->i << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl;
    bestseg=bestseg->bestseg;
  }
  seqncount(s1,bestseg->i,bestseg->j,s2,bestseg->i2,bestseg->j2,adata,as);
  LDEBUG(D_SEQIDENT,cout << "seg: "<< bestseg->i << "," << bestseg->j << " l: " << (bestseg->j-bestseg->i) << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
  LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,0,s1.seqlen));
//  seqident_seg_left_local(s1,bestseg->i,s2,bestseg->i2,adata,tncounts);
  long minleft=MIN(bestseg->i,bestseg->i2);
  seqcalign_global_noleftedgegap(s1,MAX(0,bestseg->i-2*minleft),bestseg->i,s2,MAX(0,bestseg->i2-2*minleft),bestseg->i2,adata,sws.alignws,as);

//  estr as1,as2;
//  seqcalign_global_noleftedgegap(s1,0,bestseg->i,s2,0,bestseg->i2,as1,as2);
//  cout << endl;
//  cout << as1 << endl << as2 << endl;
  LDEBUG(D_SEQIDENT,cout << "leftid: "<< bestseg->i << " l: " << (bestseg->i<bestseg->i2?bestseg->i:bestseg->i2) << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
  LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,0,s1.seqlen));
//  cout << "partial alignment" << endl;
//  print_seqali(pas1,pas2);

  return;
}
*/

void seqident_local(const estr& s1id,const estr& s2id,const eseq& s1,euintarray& kmerpos1,const eseq& s2,ealigndata& adata,esearchws& sws,const ealignscore& as,long s1_start=0,long s1_end=0)
{
  long p1,p2;
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s1.seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(s2.seq._str);
  unsigned long v1,v2;
  long lastkmerlen,lastkmerpos,lastndelta;
  if (s1_end==0) s1_end=s1.seqlen;
  long s1_len=s1_end-s1_start;

//  for (int i=0; i<11; ++i)
//    match[i]=0;
  adata.matches=0; adata.mismatches=0; adata.gaps=0; adata._score=0.0;
  adata.s1=-1; adata.e1=-1;
  adata.s2=-1; adata.e2=-1;

  char tmp[5];
  char itmp;
  int k;

  ebasicarray<ediag> diags;
  ebasicarray<ediag> dheap;

  LDEBUG(D_PROFILE,t2.reset());
  if (s1_len<s2.seqlen){
    lastkmerpos=s1_start-2*long(KMERSIZE);
    lastkmerlen=0;
    lastndelta=-s1_end-s2.seqlen;
    for (p1=s1_start; p1<s1_end-32; p1+=k){
      v1=pstr1[p1/32u]>>(2u*(p1%32u));
      v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
      for (k=0; k<32u-KMERSIZE; ++k,v1>>=2u){
        long kpos2=sws.kmerpos2[v1&KMERMASK];
        if (kpos2>sws.offset2){
          kpos2-=sws.offset2;
          long d=p1+k-lastkmerpos;
          
          if (d<=KMERSIZE && ((lastndelta==kpos2-p1-k) || (p1+k+lastndelta+KMERSIZE<=s2.seqlen && (v1&KMERMASK)==seqkmer(s2,p1+k+lastndelta)))){
            lastkmerlen+=d;
          }else{
            if (lastkmerlen){
//              cout << "p1: " << lastkmerpos-lastkmerlen+KMERSIZE<< "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
              diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE,lastndelta,lastkmerlen));
            }
            lastkmerlen=KMERSIZE;
            lastndelta=kpos2-p1-k;
          }
          lastkmerpos=p1+k;
        }
      }
    }
    v1=pstr1[p1/32u]>>(2u*(p1%32u));
    v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    for (k=0; p1+k+KMERSIZE<s1_end; ++k,v1>>=2u){
      long kpos2=sws.kmerpos2[v1&KMERMASK];
      if (kpos2>sws.offset2){
        kpos2-=sws.offset2;
        long d=p1+k-lastkmerpos;
        if (d<=KMERSIZE && ((lastndelta==kpos2-p1-k) || (p1+k+lastndelta+KMERSIZE<=s2.seqlen && (v1&KMERMASK)==seqkmer(s2,p1+k+lastndelta)))){
          lastkmerlen+=d;
        }else{
          if (lastkmerlen){
            diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE,lastndelta,lastkmerlen));
          }
          lastkmerlen=KMERSIZE;
          lastndelta=kpos2-p1-k;
        }
        lastkmerpos=p1+k;
      }
    }
    if (lastkmerlen){
      diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE,lastndelta,lastkmerlen));
//      cout << "end p1: " << lastkmerpos-lastkmerlen+KMERSIZE << "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
    }
  }else{
    lastkmerpos=-2*long(KMERSIZE);
    lastkmerlen=0;
    lastndelta=-s1_end-s2.seqlen;
    for (p2=0; p2<s2.seqlen-32; p2+=k){
      v2=pstr2[p2/32u]>>(2u*(p2%32u));
      v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
      for (k=0; k<32u-KMERSIZE; ++k,v2>>=2u){
        long kpos1=kmerpos1[v2&KMERMASK];
        if (kpos1>sws.offset){
          kpos1+=s1_start-sws.offset;
          long d=p2+k-lastkmerpos;
          if (d<=KMERSIZE && ((lastndelta==kpos1-p2-k) || (p2+k-lastndelta+KMERSIZE<=s1_end && (v2&KMERMASK)==seqkmer(s1,p2+k-lastndelta)))){
            lastkmerlen+=d;
          }else{
            if (lastkmerlen){
//              cout << "p2: " << lastkmerpos-lastkmerlen+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
              diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE-lastndelta,lastndelta,lastkmerlen));
            }
            lastkmerlen=KMERSIZE;
            lastndelta=p2+k-kpos1;
          }
          lastkmerpos=p2+k;
        }
      }
    }
    v2=pstr2[p2/32u]>>(2u*(p2%32u));
    v2|=(pstr2[p2/32u+1u]<<(64u-2u*(p2%32u)))&safe_shift[p2%32u];
    for (k=0; p2+k+KMERSIZE<s2.seqlen; ++k,v2>>=2u){
      long kpos1=kmerpos1[v2&KMERMASK];
      if (kpos1>sws.offset){
        kpos1+=s1_start-sws.offset;
        long d=p2+k-lastkmerpos;
        if (d<=KMERSIZE && ((lastndelta==kpos1-p2-k) || (p2+k-lastndelta+KMERSIZE<=s1_end && (v2&KMERMASK)==seqkmer(s1,p2+k-lastndelta)))){
          lastkmerlen+=d;
        }else{
          if (lastkmerlen){
            diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE-lastndelta,lastndelta,lastkmerlen));
          }
          lastkmerlen=KMERSIZE;
          lastndelta=p2+k-kpos1;
        }
        lastkmerpos=p2+k;
      }
    }
    if (lastkmerlen){
//      cout << "end p2: " << lastkmerpos-lastkmerlen+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
      diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE-lastndelta,lastndelta,lastkmerlen));
    }
  }
  LDEBUG(D_PROFILE,tdp1=tdp1*0.99+t2.lap()*0.01);

  if (diags.size()==0){
//    lderror("no shared segments found");
    return;
  }

  // sparse dynamic programming to find longest chain of segments
  ediag *bestseg=&diags[0];
  multimap<long,ediag*> segments;
  for (long j=0; j<diags.size(); ++j){
    ediag& diag(diags.at(j));
    diag.V=(diag.j-diag.i)*as.match;
    if (diag.V > bestseg->V) bestseg=&diag;
    segments.insert(pair<long,ediag*>(diag.i,&diag));
    segments.insert(pair<long,ediag*>(diag.j,&diag));
  }

  multimap<long,ediag*> Lsegments;
  multimap<long,ediag*>::iterator sit;
  multimap<long,ediag*>::iterator slb;
  for (sit=segments.begin(); sit!=segments.end(); ++sit){
//      cout << "x: " << sit->first << " start: " << (sit->first==sit->second->i?"1":"0") << " i: " << sit->second->i << " i2: " << sit->second->i2 << " Lsegs: " << Lsegments.size() << endl;
    if (sit->first==sit->second->i) { // start of rectangle
//      if (Lsegments.size()==0) continue;
//      slb=Lsegments.upper_bound(sit->second->i2-1);
//      if (slb==Lsegments.end()) continue;
      ediag *tmpbest=0x00;
      double tmpbestscore=(sit->second->j-sit->second->i)*as.match;
      for (slb=Lsegments.begin(); slb!=Lsegments.end(); ++slb){
        if (slb->first >= sit->second->i2) break;
//      sit->second->V=sit->second->j-sit->second->i+slb->second->V;
        long gapdiff=abs((sit->second->i-sit->second->i2)-(slb->second->i-slb->second->i2));
        double tmpscore=(sit->second->j-sit->second->i)*as.match - gapdiff*as.gapext-(gapdiff>0?as.gapopen:0) + slb->second->V;
        if (tmpscore>tmpbestscore) { tmpbestscore=tmpscore; tmpbest=slb->second; }
      }
      sit->second->bestseg=tmpbest;
      sit->second->V=tmpbestscore;
      
//      sit->second->V=(sit->second->j-sit->second->i)*as.match-gapdiff*as.gapext-(gapdiff>0?as.gapopen:0)+ MIN(sit->second->j-slb->second->i,sit->second->j2-slb->second->i2)*(9.0/10.0*as.match-1.0/10.0*as.mismatch) +slb->second->V;
/*
      sit->second->bestseg=0x00;
      double newscore=-gapdiff*as.gapext-(gapdiff>0?as.gapopen:0)+ MIN(sit->second->j-slb->second->i,sit->second->j2-slb->second->i2)*(7.0/8.0*as.match-1.0/8.0*as.mismatch) +slb->second->V;
      if (sit->second->V+MIN(sit->second->i,sit->second->j)*(7.0/8.0*as.match-1.0/8.0*as.mismatch) < newscore){
        sit->second->V=newscore;
        sit->second->bestseg=slb->second;
      }else{
        sit->second->V+=MIN(sit->second->i,sit->second->j)*(7.0/8.0*as.match-1.0/8.0*as.mismatch);
      }
*/
//        double newscore=(sit->second->j-sit->second->i)*as.match - gapdiff*as.gapext-(gapdiff>0?as.gapopen:0) + slb->second->V;
//      if (sit->second->V < newscore){
//        sit->second->V=newscore;
//        sit->second->bestseg=slb->second;
//      }
//      cout << "linking to: " << sit->second->i << " " << slb->second->j << endl;
    } else {  // end of rectangle
      slb=Lsegments.lower_bound(sit->second->j2);
      if (slb==Lsegments.end() || sit->second->V > slb->second->V){
        while (slb!=Lsegments.end()){
          if (slb->second->j2>=sit->second->j2 && slb->second->V<=sit->second->V){
            Lsegments.erase(slb++);
            continue;
          }
          ++slb;
        }
        if (bestseg==0x00 || sit->second->V > bestseg->V)
          bestseg=sit->second;
        Lsegments.insert(pair<long,ediag*>(sit->second->j2,sit->second));
      }
    }
  }
  LDEBUG(D_PROFILE,tdp2=tdp2*0.99+t2.lap()*0.01);

//  estr fas1,fas2;
//  cout << "full alignment" << endl;
//  cout << "ali: " << 0 <<","<<s2.seqlen<<" : " << 0 <<","<<s1.seqlen << endl; 
//  adata._score=seqcalign_global_noedgegap(s1,0,s1.seqlen,s2,0,s2.seqlen,fas1,fas2);
//  print_seqali(fas1,fas2);

//  estr pas1,pas2;
//  seqident_seg_local(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,tncounts);
//  seqident_seg_sw_right(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,tncounts);
//  estr as1,as2;
//  seqcalign_global_noedgegap(s1,0,s1.seqlen,s2,0,s2.seqlen,as1,as2);
//  print_seqali(as1,as2);

/*  
  if (bestseg->bestseg==0x00){ // single segment hit, do full alignment
    seqcalign_global_noedgegap(s1,0,s1.seqlen,s2,0,s2.seqlen,adata,sws.alignws,as);
    return;
  }
*/


  float tmpflanks=0.0;
  ldieif(s2.seqlen<bestseg->j2,"segment start larger than sequence length: "+estr(s2.seqlen)+" < "+estr(bestseg->j2)+" s1 id: "+s1id+" s2 id: "+s2id);
  long minright=MIN(s1.seqlen-bestseg->j,s2.seqlen-bestseg->j2);
//  seqcalign_global_norightedgegap(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,alignws);
  seqcalign_local_rightext(s1,bestseg->j,MIN(s1.seqlen,bestseg->j+2*minright),s2,bestseg->j2,MIN(s2.seqlen,bestseg->j2+2*minright),adata,sws.alignws,as);
  LDEBUG(D_PROFILE,tmpflanks+=t2.lap());
  adata.profile.add(AT_RIGHT);
  LDEBUG(D_SEQIDENT,cout << "rightid: " << bestseg->j << "," << s1_end << " l: " << s1_end-bestseg->j << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
  LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,s1_start,s1_end));
  while (bestseg!=0x00 && bestseg->bestseg != 0x00){
    seqncount(s1,bestseg->i,bestseg->j,s2,bestseg->i2,bestseg->j2,adata,as);
    adata.profile.add(AT_ID);
    LDEBUG(D_SEQIDENT,cout << "seg: "<< bestseg->i << "," << bestseg->j << " l: " << (bestseg->j-bestseg->i) << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
    LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,s1_start,s1_end));
    long gapdiff=abs((bestseg->bestseg->i-bestseg->bestseg->i2)-(bestseg->i-bestseg->i2));
//    if (gapdiff==0) {/* fast alignment no gaps */
//      seqident_seg_nogaps(s1,bestseg->bestseg->j,bestseg->i,s2,bestseg->bestseg->j2,bestseg->i2,adata,as);
//      adata.profile.add(AT_ALIGNF);
//    } else { /* full sw alignment */
      seqcalign_global(s1,bestseg->bestseg->j,bestseg->i,s2,bestseg->bestseg->j2,bestseg->i2,adata,sws.alignws,as);
      adata.profile.add(AT_ALIGN);
//    }
    LDEBUG(D_SEQIDENT,cout << "segid: " << bestseg->bestseg->j << "," << bestseg->i << " l: " << bestseg->i-bestseg->bestseg->j << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << " gapdiff: " << gapdiff << endl);
    LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,s1_start,s1_end));
//    cout << "segid: " << bestseg->i << "," << bestseg->j << " l: " << bestseg->j-bestseg->i << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl;
    bestseg=bestseg->bestseg;
  }
  seqncount(s1,bestseg->i,bestseg->j,s2,bestseg->i2,bestseg->j2,adata,as);
  adata.profile.add(AT_ID);
  LDEBUG(D_SEQIDENT,cout << "seg: "<< bestseg->i << "," << bestseg->j << " l: " << (bestseg->j-bestseg->i) << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
  LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,s1_start,s1_end));


  
  LDEBUG(D_PROFILE,tdpmd=tdpmd*0.99+t2.lap()*0.01);
//  seqident_seg_left_local(s1,bestseg->i,s2,bestseg->i2,adata,tncounts);
  long minleft=MIN(bestseg->i,bestseg->i2);
  adata.s1=bestseg->i; adata.s2=bestseg->i2;
  seqcalign_local_leftext(s1,MAX(0,bestseg->i-2*minleft),bestseg->i,s2,MAX(0,bestseg->i2-2*minleft),bestseg->i2,adata,sws.alignws,as);
  adata.profile.add(AT_LEFT);
  LDEBUG(D_PROFILE,tmpflanks+=t2.lap());
  tdpfl=tdpfl*0.99+tmpflanks*0.01;

//  estr as1,as2;
//  seqcalign_global_noleftedgegap(s1,0,bestseg->i,s2,0,bestseg->i2,as1,as2);
//  cout << endl;
//  cout << as1 << endl << as2 << endl;
  LDEBUG(D_SEQIDENT,cout << "leftid: "<< bestseg->i << " l: " << (bestseg->i<bestseg->i2?bestseg->i:bestseg->i2) << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
  LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,s1_start,s1_end));
//  cout << "partial alignment" << endl;
//  print_seqali(pas1,pas2);

  return;
}



void kmercount_mask(estrarrayof<eseq>& seqs,eintarray& seqids,uint64_t *kmermask,unsigned int maskid,eintarray& idcount,unsigned int *akmers,unsigned long akmask)
{
  int k;
  long p;
  long kmerpos;
  for (int i=0; i<seqids.size(); ++i){
    if (seqids[i]<seqs.size()){
      eseq &s(seqs.values(seqids[i]));
      unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
      unsigned long v;
      kmerpos=-long(KMERSIZE2);
      for (p=0; p+32<s.seqlen; p+=32-KMERSIZE2){
        v=pstr[p/32u]>>(2u*(p%32u));
        v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
        for (k=0; k+KMERSIZE2<32; ++k,v>>=2u){
          if ((kmermask[(v&KMERMASK2)/64u]>>((v&KMERMASK2)%64u))&0x01ul){
            unsigned long d=p+k-kmerpos;
            if (d>KMERSIZE2) d=KMERSIZE2;
            idcount[i]+=d;
//            idcount[i]+=d1_lt[(p+k-kmerpos)&BMASK16];
            kmerpos=p+k;
          }
        }
      }
      v=pstr[p/32u]>>(2u*(p%32u));
      v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
      for (k=0; p+k+KMERSIZE2<s.seqlen; ++k,v>>=2u){
        if ((kmermask[(v&KMERMASK2)/64u]>>((v&KMERMASK2)%64u))&0x01ul){
          unsigned long d=p+k-kmerpos;
          if (d>KMERSIZE2) d=KMERSIZE2;
          idcount[i]+=d;
//          idcount[i]+=d1_lt[(p+k-kmerpos)&BMASK16];
          kmerpos=p+k;
        }
      }
    }else{
      eseq &s(seqs.values(seqids[i]-seqs.size()));
      unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
      unsigned long v;
      kmerpos=-long(KMERSIZE2);
      for (p=0; p+32<s.seqlen; p+=32-KMERSIZE2){
        v=pstr[p/32u]>>(2u*(p%32u));
        v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
        for (k=0; k+KMERSIZE2<32; ++k,v>>=2u){
          uint32_t tmpv=kmer_rev_lt2[v&KMERMASK2];
          if ((kmermask[tmpv/64u]>>(tmpv%64u))&0x01ul){
            unsigned long d=p+k-kmerpos;
            if (d>KMERSIZE2) d=KMERSIZE2;
            idcount[i]+=d;
//            idcount[i]+=d1_lt[(p+k-kmerpos)&BMASK16];
            kmerpos=p+k;
          }
        }
      }
      v=pstr[p/32u]>>(2u*(p%32u));
      v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
      for (k=0; p+k+KMERSIZE2<s.seqlen; ++k,v>>=2u){
        uint32_t tmpv=kmer_rev_lt[v&KMERMASK2];
        if ((kmermask[tmpv/64u]>>(tmpv%64u))&0x01ul){
          unsigned long d=p+k-kmerpos;
          if (d>KMERSIZE2) d=KMERSIZE2;
          idcount[i]+=d;
//          idcount[i]+=d1_lt[(p+k-kmerpos)&BMASK16];
          kmerpos=p+k;
        }
      }
    }
  }
}

void kmercount_mask(estrarrayof<eseq>& seqs,eintarray& seqids,euintarray& kmermask,unsigned int maskid,eintarray& idcount)
{
  int k;
  long p;
  long kmerpos;
  for (int i=0; i<seqids.size(); ++i){
    if (seqids[i]<seqs.size()){
      eseq &s(seqs.values(seqids[i]));
      unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
      unsigned long v;
      kmerpos=-long(KMERSIZE2);
      for (p=0; p+32<s.seqlen; p+=32-KMERSIZE2){
        v=pstr[p/32u]>>(2u*(p%32u));
        v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
        for (k=0; k+KMERSIZE2<32; ++k,v>>=2u){
          if (kmermask[v&KMERMASK2]==maskid){
//          if ((kmermask[(v&KMERMASK)/64u]>>((v&KMERMASK)%64u))&0x01ul){
            unsigned long d=p+k-kmerpos;
            if (d>KMERSIZE2) d=KMERSIZE2;
            idcount[i]+=d;
            kmerpos=p+k;
          }
        }
      }
      v=pstr[p/32u]>>(2u*(p%32u));
      v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
      for (k=0; p+k+KMERSIZE2<s.seqlen; ++k,v>>=2u){
        if (kmermask[v&KMERMASK2]==maskid){
//        if ((kmermask[(v&KMERMASK)/64u]>>((v&KMERMASK)%64u))&0x01ul){
          unsigned long d=p+k-kmerpos;
          if (d>KMERSIZE2) d=KMERSIZE2;
          idcount[i]+=d;
          kmerpos=p+k;
        }
      }
    }else{
      eseq &s(seqs.values(seqids[i]-seqs.size()));
      unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
      unsigned long v;
      kmerpos=-long(KMERSIZE2);
      for (p=0; p+32<s.seqlen; p+=32-KMERSIZE2){
        v=pstr[p/32u]>>(2u*(p%32u));
        v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
        for (k=0; k+KMERSIZE2<32; ++k,v>>=2u){
//          uint32_t tmpv=kmer_rev_lt[v&KMERMASK];
//          if ((kmermask[tmpv/64u]>>(tmpv%64u))&0x01ul){
          if (kmermask[kmer_rev_lt2[v&KMERMASK2]]==maskid){
            unsigned long d=p+k-kmerpos;
            if (d>KMERSIZE2) d=KMERSIZE2;
            idcount[i]+=d;
            kmerpos=p+k;
          }
        }
      }
      v=pstr[p/32u]>>(2u*(p%32u));
      v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
      for (k=0; p+k+KMERSIZE2<s.seqlen; ++k,v>>=2u){
//        uint32_t tmpv=kmer_rev_lt[v&KMERMASK];
//        if ((kmermask[tmpv/64u]>>(tmpv%64u))&0x01ul){
        if (kmermask[kmer_rev_lt2[v&KMERMASK2]]==maskid){
          unsigned long d=p+k-kmerpos;
          if (d>KMERSIZE2) d=KMERSIZE2;
          idcount[i]+=d;
          kmerpos=p+k;
        }
      }
    }
  }
}

void setkmermask(uint64_t *kmermask,eseq& s,unsigned int *akmers,unsigned long akmask)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
//  unsigned int akmask=0xFul;
  int k;
  long p;
  for (p=0; p<s.seqlen-32; p+=32-KMERSIZE2){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k<32-KMERSIZE2; ++k,v>>=2u){
      if (akmers[v&akmask]==0x0u) continue;
      uint32_t tmpv=v&KMERMASK2;
      kmermask[tmpv/64u]|=0x1u<<(tmpv%64u);
    }
  }
  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k<s.seqlen-long(KMERSIZE2); ++k,v>>=2u){
    if (akmers[v&akmask]==0x0u) continue;
    uint32_t tmpv=v&KMERMASK2;
    kmermask[tmpv/64u]|=0x1u<<(tmpv%64u);
  }
}

void setkmermask(eintarray& kmermask,eseq& s,int maskid)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
//  unsigned int akmask=0xFul;
  int k;
  long p;
  for (p=0; p<s.seqlen-32; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k<32-KMERSIZE; ++k,v>>=2u)
      kmermask[v&KMERMASK]=maskid;
  }
  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k<s.seqlen-long(KMERSIZE); ++k,v>>=2u)
    kmermask[v&KMERMASK]=maskid;
}

void setkmermask(euintarray& kmermask,eseq& s,unsigned int maskid,unsigned int *akmers,unsigned long akmask,unsigned long start,unsigned long end)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
//  unsigned int akmask=0xFul;
  int k;
  long p;
  for (p=start; p+32<end; p+=32-KMERSIZE2){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k+KMERSIZE2<32; ++k,v>>=2u){
      if (akmers[v&akmask]==0x0u) continue;
      kmermask[v&KMERMASK2]=maskid;
    }
  }
  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k+KMERSIZE2<end; ++k,v>>=2u){
    if (akmers[v&akmask]==0x0u) continue;
    kmermask[v&KMERMASK2]=maskid;
  }
}

void setkmermask(euintarray& kmermask,eseq& s,unsigned int maskid,unsigned int *akmers,unsigned long akmask)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
//  unsigned int akmask=0xFul;
  int k;
  long p;
  for (p=0; p+32<s.seqlen; p+=32-KMERSIZE2){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k<32-KMERSIZE2; ++k,v>>=2u){
      if (akmers[v&akmask]==0x0u) continue;
      kmermask[v&KMERMASK2]=maskid;
    }
  }
  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k<s.seqlen-int(KMERSIZE2); ++k,v>>=2u){
    if (akmers[v&akmask]==0x0u) continue;
    kmermask[v&KMERMASK2]=maskid;
  }
}

void setkmerpos(euintarray& kmerpos,eseq& s,unsigned int offset,unsigned long start,unsigned long end)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  int k;
  long p;
  for (p=start; p+32<end; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k+KMERSIZE<32; ++k,v>>=2u)
      kmerpos[v&KMERMASK]=offset+(p+k-start);
  }
  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k+KMERSIZE<end; ++k,v>>=2u)
    kmerpos[v&KMERMASK]=offset+(p+k-start);
}

void setkmerpos(euintarray& kmerpos,eseq& s,unsigned int offset)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  int k;
  long p;
  for (p=0; p+32<s.seqlen; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k+KMERSIZE<32; ++k,v>>=2u)
      kmerpos[v&KMERMASK]=offset+p+k;
  }
  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k+KMERSIZE<s.seqlen; ++k,v>>=2u)
    kmerpos[v&KMERMASK]=offset+p+k;
}



void kmercount_single3(ebasicarray<deque<int>*>& otukmers,eseq& s,uint32_t sp,uint32_t rp,ebasicarray<uint32_t>& idcount,eintarray& kmerpos)
//void kmercount_single2(ebasicarray<deque<int>*>& otukmers,eseq& s,uint32_t rp,eintarray& idcount,eintarray& kmerpos)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  int k;
  int p=0;
  sp<<=10u;
  for (; p<long(s.seqlen)-32; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k<32-KMERSIZE; ++k,v>>=2u){
      if (otukmers[v&KMERMASK]){
        deque<int> &otuind(*otukmers[v&KMERMASK]);
        for (unsigned int l=0; l<otuind.size(); ++l){
          unsigned int tl=otuind[l];
          unsigned int d=rp+p+k-kmerpos[tl];
          unsigned int &idc(idcount[tl]);
//          d=d1_lt[d&BMASK16]|d2_lt[(d>>16u)&BMASK16];
          if (d>KMERSIZE) d=KMERSIZE;
          uint32_t tmp=(id_lt[((idc^sp)>>10u)&BMASK11]*id_lt[((idc^sp)>>21u)&BMASK11]*idc&BMASK10);
          tmp=(((tmp+d)&BMASK10)<tmp?1024u:(tmp+d)&BMASK10);
//          if (d==1 || d==KMERSIZE)
          idcount[tl]=sp|tmp;
//          idcount[tl]+=d;
          kmerpos[tl]=rp+p+k;
        }
      }
    }
  }

//  cout << "p: " << p << " seqlen: " << s.seqlen << endl;

  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k<long(s.seqlen)-long(KMERSIZE); ++k,v>>=2u){
    if (otukmers[v&KMERMASK]){
      deque<int> &otuind(*otukmers[v&KMERMASK]);
      for (unsigned int l=0; l<otuind.size(); ++l){
        unsigned int tl=otuind[l];
        unsigned int d=rp+p+k-kmerpos[tl];
        unsigned int &idc(idcount[tl]);
        if (d>KMERSIZE) d=KMERSIZE;
//        d=d1_lt[d&BMASK16]|d2_lt[(d>>16u)&BMASK16];
//        idcount[tl]+=d;
        idcount[tl]=sp|((id_lt[((idc^sp)>>10u)&BMASK11]*id_lt[((idc^sp)>>21u)&BMASK11]*idc)&BMASK10)+d;
        kmerpos[tl]=rp+p+k;
      }
    }
  }
}

void kmercount_single2(ebasicarray<deque<int>*>& otukmers,eseq& s,uint32_t rp,ebasicarray<uint32_t>& idcount,eintarray& kmerpos)
//void kmercount_single2(ebasicarray<deque<int>*>& otukmers,eseq& s,uint32_t rp,eintarray& idcount,eintarray& kmerpos)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  int k;
  int p=0;
  for (; p<long(s.seqlen)-32; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k<32-KMERSIZE; ++k,v>>=2u){
      if (otukmers[v&KMERMASK]){
        deque<int> &otuind(*otukmers[v&KMERMASK]);
        for (unsigned int l=0; l<otuind.size(); ++l){
          unsigned int tl=otuind[l];
          unsigned int d=rp+p+k-kmerpos[tl];
//          d=d1_lt[d&BMASK16]|d2_lt[(d>>16u)&BMASK16];
          if (d>KMERSIZE) d=KMERSIZE;
//          if (d==1 || d==KMERSIZE)
//          idcount[tl]=sp|((id_lt[((idc^sp)>>10u)&BMASK11]*id_lt[((idc^sp)>>21u)&BMASK11]*idc)&BMASK10)+d;
          idcount[tl]+=d;
          kmerpos[tl]=rp+p+k;
        }
      }
    }
  }

//  cout << "p: " << p << " seqlen: " << s.seqlen << endl;

  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k<long(s.seqlen)-long(KMERSIZE); ++k,v>>=2u){
    if (otukmers[v&KMERMASK]){
      deque<int> &otuind(*otukmers[v&KMERMASK]);
      for (unsigned int l=0; l<otuind.size(); ++l){
        unsigned int tl=otuind[l];
        unsigned int d=rp+p+k-kmerpos[tl];
        unsigned int &idc(idcount[tl]);
        if (d>KMERSIZE) d=KMERSIZE;
//        d=d1_lt[d&BMASK16]|d2_lt[(d>>16u)&BMASK16];
        idcount[tl]+=d;
        kmerpos[tl]=rp+p+k;
      }
    }
  }
}

void kmercount_single(ebasicarray<deque<uint32_t> >& otukmers,eseq& s,ebasicarray<uint32_t>& idcount,eintarray& kmerpos)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  int k;
  long p=0;
  for (; p<long(s.seqlen)-32; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k<32-KMERSIZE; ++k,v>>=2u){
//      if (otukmers[v&KMERMASK]){
        deque<uint32_t> &otuind(otukmers[v&KMERMASK]);
        for (unsigned int l=0; l<otuind.size(); ++l){
          unsigned int tl=otuind[l];
          unsigned long d=p+k-kmerpos[tl];
          if (d>KMERSIZE) d=KMERSIZE;
//          if (d==1 || d==KMERSIZE)
            idcount[tl]+=d;
          kmerpos[tl]=p+k;
        }
//      }
    }
  }

//  cout << "p: " << p << " seqlen: " << s.seqlen << endl;

  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k<long(s.seqlen)-long(KMERSIZE); ++k,v>>=2u){
//    if (otukmers[v&KMERMASK]){
      deque<uint32_t> &otuind(otukmers[v&KMERMASK]);
      for (unsigned int l=0; l<otuind.size(); ++l){
        unsigned int tl=otuind[l];
        unsigned long d=p+k-kmerpos[tl];
        if (d>KMERSIZE) d=KMERSIZE;
//        if (d==1 || d==KMERSIZE)
          idcount[tl]+=d;
        kmerpos[tl]=p+k;
      }
//    }
  }
}

uint32_t bmask_lt[2]={0x0u,0xffffffffu};
int16_t idc_lt[2]={1,-1};

void kmercount_both_nopos2_skip(int scount,ebasicarray<deque<uint32_t> >& otukmers,eseq& s,unsigned long start,unsigned long end,eintarray& idcount,eintarray& kmerpos,unsigned int *akmers,unsigned long akmask,int skipn)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  int k;
  long p=start;
  for (; p+32<end; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k<32-KMERSIZE; ++k,v>>=2u){
      if ((p+k)%skipn) continue;
//      cout << kmer2str(v&KMERMASK) << " " << kmer2str(kmer_rev_lt[v&KMERMASK]) << endl;
//      if (otukmers[v&KMERMASK]){
      if (akmers[v&akmask]){
        deque<uint32_t> &otuind(otukmers[v&KMERMASK]);
        for (unsigned int l=0; l<otuind.size(); ++l)
          ++idcount[(otuind[l]&BMASK31)];
//          idcount[(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
      }
      if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
//      if (otukmers[kmer_rev_lt[v&KMERMASK]]){
        deque<uint32_t> &otuind(otukmers[kmer_rev_lt[v&KMERMASK]]);
        for (unsigned int l=0; l<otuind.size(); ++l)
//          idcount[scount+(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
          ++idcount[scount+(otuind[l]&BMASK31)];
      }
    }
  }

//  cout << "p: " << p << " seqlen: " << s.seqlen << endl;

  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k+KMERSIZE<end; ++k,v>>=2u){
    if ((p+k)%skipn) continue;
//    if (otukmers[v&KMERMASK]){
    if (akmers[v&akmask]){
      deque<uint32_t> &otuind(otukmers[v&KMERMASK]);
      for (unsigned int l=0; l<otuind.size(); ++l)
        ++idcount[(otuind[l]&BMASK31)];
//        idcount[(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
    }
//    if (otukmers[kmer_rev_lt[v&KMERMASK]]){
    if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
      deque<uint32_t> &otuind(otukmers[kmer_rev_lt[v&KMERMASK]]);
      for (unsigned int l=0; l<otuind.size(); ++l)
        ++idcount[scount+(otuind[l]&BMASK31)];
//        idcount[scount+(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
    }
  }
}

void kmercount_both_nopos2(int scount,ebasicarray<deque<uint32_t> >& otukmers,eseq& s,eintarray& idcount,eintarray& kmerpos,unsigned int *akmers,unsigned long akmask)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  int k;
  long p=0;
  for (; p+32<s.seqlen; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k<32-KMERSIZE; ++k,v>>=2u){
//      cout << kmer2str(v&KMERMASK) << " " << kmer2str(kmer_rev_lt[v&KMERMASK]) << endl;
//      if (otukmers[v&KMERMASK]){
      if (akmers[v&akmask]){
        deque<uint32_t> &otuind(otukmers[v&KMERMASK]);
        for (unsigned int l=0; l<otuind.size(); ++l)
          idcount[(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
      }
      if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
//      if (otukmers[kmer_rev_lt[v&KMERMASK]]){
        deque<uint32_t> &otuind(otukmers[kmer_rev_lt[v&KMERMASK]]);
        for (unsigned int l=0; l<otuind.size(); ++l)
          idcount[scount+(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
      }
    }
  }

//  cout << "p: " << p << " seqlen: " << s.seqlen << endl;

  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k+KMERSIZE<s.seqlen; ++k,v>>=2u){
//    if (otukmers[v&KMERMASK]){
    if (akmers[v&akmask]){
      deque<uint32_t> &otuind(otukmers[v&KMERMASK]);
      for (unsigned int l=0; l<otuind.size(); ++l)
        idcount[(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
    }
//    if (otukmers[kmer_rev_lt[v&KMERMASK]]){
    if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
      deque<uint32_t> &otuind(otukmers[kmer_rev_lt[v&KMERMASK]]);
      for (unsigned int l=0; l<otuind.size(); ++l)
        idcount[scount+(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
    }
  }
}

void kmercount_both_nopos(int scount,ebasicarray<deque<int> >& otukmers,eseq& s,eintarray& idcount,eintarray& kmerpos,unsigned int *akmers,unsigned long akmask)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  int k;
  long p=0;
  for (; p+32<s.seqlen; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k<32-KMERSIZE; ++k,v>>=2u){
//      cout << kmer2str(v&KMERMASK) << " " << kmer2str(kmer_rev_lt[v&KMERMASK]) << endl;
//      if (otukmers[v&KMERMASK]){
      if (akmers[v&akmask]){
        deque<int> &otuind(otukmers[v&KMERMASK]);
        for (unsigned int l=0; l<otuind.size(); ++l)
          ++idcount[otuind[l]];
      }
      if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
//      if (otukmers[kmer_rev_lt[v&KMERMASK]]){
        deque<int> &otuind(otukmers[kmer_rev_lt[v&KMERMASK]]);
        for (unsigned int l=0; l<otuind.size(); ++l)
          ++idcount[scount+otuind[l]];
      }
    }
  }

//  cout << "p: " << p << " seqlen: " << s.seqlen << endl;

  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k+KMERSIZE<s.seqlen; ++k,v>>=2u){
//    if (otukmers[v&KMERMASK]){
    if (akmers[v&akmask]){
      deque<int> &otuind(otukmers[v&KMERMASK]);
      for (unsigned int l=0; l<otuind.size(); ++l)
        ++idcount[otuind[l]];
    }
//    if (otukmers[kmer_rev_lt[v&KMERMASK]]){
    if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
      deque<int> &otuind(otukmers[kmer_rev_lt[v&KMERMASK]]);
      for (unsigned int l=0; l<otuind.size(); ++l)
        ++idcount[scount+otuind[l]];
    }
  }
}

void kmercount_both2(int scount,ebasicarray<deque<int> >& otukmers,eseq& s,eintarray& idcount,eintarray& kmerpos,unsigned int *akmers,unsigned long akmask)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  int k;
  long p=0;
  for (; p+32<s.seqlen; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k<32-KMERSIZE; ++k,v>>=2u){
//      cout << kmer2str(v&KMERMASK) << " " << kmer2str(kmer_rev_lt[v&KMERMASK]) << endl;
//      if (otukmers[v&KMERMASK]){
      if (akmers[v&akmask]){
        deque<int> &otuind(otukmers[v&KMERMASK]);
        for (unsigned int l=0; l<otuind.size(); ++l){
          unsigned int tl=otuind[l];
          unsigned int d=p+k+KMERSIZE-kmerpos[tl];
          if (d>KMERSIZE) d=KMERSIZE;
//          idcount[tl]+=d1_lt[(p+k+KMERSIZE-kmerpos[tl])&BMASK16];
          idcount[tl]+=d;
          kmerpos[tl]=p+k+KMERSIZE;
        }
      }
      if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
//      if (otukmers[kmer_rev_lt[v&KMERMASK]]){
        deque<int> &otuind2(otukmers[kmer_rev_lt[v&KMERMASK]]);
        for (unsigned int l=0; l<otuind2.size(); ++l){
          uint32_t tl=scount+otuind2[l];
          unsigned int d=p+k+KMERSIZE-kmerpos[tl];
          if (d>KMERSIZE) d=KMERSIZE;
          idcount[tl]+=d;
          kmerpos[tl]=p+k+KMERSIZE;
        }
      }
    }
  }

//  cout << "p: " << p << " seqlen: " << s.seqlen << endl;

  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k+KMERSIZE<s.seqlen; ++k,v>>=2u){
//    if (otukmers[v&KMERMASK]){
    if (akmers[v&akmask]){
      deque<int> &otuind(otukmers[v&KMERMASK]);
      for (unsigned int l=0; l<otuind.size(); ++l){
        uint32_t tl=otuind[l];
        unsigned int d=p+k+KMERSIZE-kmerpos[tl];
        if (d>KMERSIZE) d=KMERSIZE;
        idcount[tl]+=d;
        kmerpos[tl]=p+k+KMERSIZE;
      }
    }
//    if (otukmers[kmer_rev_lt[v&KMERMASK]]){
    if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
      deque<int> &otuind2(otukmers[kmer_rev_lt[v&KMERMASK]]);
      for (unsigned int l=0; l<otuind2.size(); ++l){
        uint32_t tl=scount+otuind2[l];
        unsigned int d=p+k+KMERSIZE-kmerpos[tl];
        if (d>KMERSIZE) d=KMERSIZE;
        idcount[tl]+=d;
        kmerpos[tl]=p+k+KMERSIZE;
      }
    }
  }
}

void kmercount_both(int scount,ebasicarray<deque<int> >& otukmers,eseq& s,ebasicarray<uint32_t>& idcount,uint64_t *bitmask,unsigned int *akmers,unsigned long akmask)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  int k;
  long p=0;
  for (; p+32<s.seqlen; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k<32-KMERSIZE; ++k,v>>=2u){
//      cout << kmer2str(v&KMERMASK) << " " << kmer2str(kmer_rev_lt[v&KMERMASK]) << endl;
      if (akmers[v&akmask]){
        deque<int> &otuind(otukmers[v&KMERMASK]);
        for (unsigned int l=0; l<otuind.size(); ++l){
          unsigned int tl=otuind[l];
          uint16_t bt=(bitmask[tl/64u]>>(tl%64u))&0x1ul;
          bitmask[tl/64u]|=(0x1ul<<(tl%64u));
//          unsigned int d=p+k+KMERSIZE-(kmerpos[tl]&bmask_lt[bt]);
//          uint32_t idc=idcount[tl]&bmask_lt[bt];
          uint32_t idc=idcount[tl]&bmask_lt[bt];
          uint32_t d=d1_lt[p+k+KMERSIZE-(idc>>16u)];
//          if (d>KMERSIZE) d=KMERSIZE;
          idcount[tl]=((idc+d)&BMASK16)|((p+k+KMERSIZE)<<16u);
//          if (d>KMERSIZE) d=KMERSIZE;
//          idcount[tl]=((idc+d)&BMASK16)|((p+k+KMERSIZE)<<16u);
//          kmerpos[tl]=p+k+KMERSIZE;
//          idcount[tl]+=d;
//          kmerpos[tl]=p+k;
        }
      }
      if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
        deque<int> &otuind(otukmers[kmer_rev_lt[v&KMERMASK]]);
        for (unsigned int l=0; l<otuind.size(); ++l){
          uint32_t tl=scount+otuind[l];
          uint16_t bt=(bitmask[tl/64u]>>(tl%64u))&0x1ul;
          bitmask[tl/64u]|=(0x1ul<<(tl%64u));
          uint32_t idc=idcount[tl]&bmask_lt[bt];
          uint32_t d=d1_lt[p+k+KMERSIZE-(idc>>16u)];
//          uint32_t d=p+k+KMERSIZE-(idc>>16u);
//          if (d>KMERSIZE) d=KMERSIZE;
          idcount[tl]=((idc+d)&BMASK16)|((p+k+KMERSIZE)<<16u);
//          uint32_t idc=idcount[tl]&bmask_lt[bt];
//          uint32_t d=p+k+KMERSIZE-(idc>>16u);
//          if (d>KMERSIZE) d=KMERSIZE;
//          idcount[tl]=((idc+d)&BMASK16)|((p+k+KMERSIZE)<<16u);
        }
      }
    }
  }

//  cout << "p: " << p << " seqlen: " << s.seqlen << endl;

  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k+KMERSIZE<s.seqlen; ++k,v>>=2u){
    if (akmers[v&akmask]){
      deque<int> &otuind(otukmers[v&KMERMASK]);
      for (unsigned int l=0; l<otuind.size(); ++l){
        uint32_t tl=otuind[l];
        uint16_t bt=(bitmask[tl/64u]>>(tl%64u))&0x1ul;
        bitmask[tl/64u]|=(0x1ul<<(tl%64u));
        uint32_t idc=idcount[tl]&bmask_lt[bt];
        uint32_t d=d1_lt[p+k+KMERSIZE-(idc>>16u)];
//        uint32_t d=p+k+KMERSIZE-(idc>>16u);
//        if (d>KMERSIZE) d=KMERSIZE;
        idcount[tl]=((idc+d)&BMASK16)|((p+k+KMERSIZE)<<16u);
//        uint32_t idc=idcount[tl]&bmask_lt[bt];
//        uint32_t d=p+k+KMERSIZE-(idc>>16u);
//        idcount[tl]=((idc+d)&BMASK16)|((p+k+KMERSIZE)<<16u);
      }
    }
    if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
      deque<int> &otuind(otukmers[kmer_rev_lt[v&KMERMASK]]);
      for (unsigned int l=0; l<otuind.size(); ++l){
        uint32_t tl=scount+otuind[l];
        uint16_t bt=(bitmask[tl/64u]>>(tl%64u))&0x1ul;
        bitmask[tl/64u]|=(0x1ul<<(tl%64u));
        uint32_t idc=idcount[tl]&bmask_lt[bt];
        uint32_t d=d1_lt[p+k+KMERSIZE-(idc>>16u)];
//        uint32_t d=p+k+KMERSIZE-(idc>>16u);
//        if (d>KMERSIZE) d=KMERSIZE;
        idcount[tl]=((idc+d)&BMASK16)|((p+k+KMERSIZE)<<16u);
//        uint32_t d=p+k+KMERSIZE-(idc>>16u);
      }
    }
  }
}

void kmercount(earray<eintarray>& otus,ebasicarray<eintarray*>& otukmers,int seqcount,eseq& s,eintarray& idcount,eintarray& kmerpos,eintarray& tmpkmers,eintarray& seqkmers,eintarray& revtmpkmers,eintarray& revseqkmers)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  int k;
  int p=0;
  for (; p<int(s.seqlen)-32; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k<32-KMERSIZE; ++k,v>>=2u){
      if (otukmers[v&KMERMASK]){
        eintarray &otuind(*otukmers[v&KMERMASK]);
        for (unsigned int l=0; l<otuind.size(); ++l){
          unsigned int tl=otuind[l];
          unsigned int d=p+k-kmerpos[tl];
          if (d>KMERSIZE) d=KMERSIZE;
          idcount[tl]+=d;
          kmerpos[tl]=p+k;
        }
      }
      if (tmpkmers[v&KMERMASK]!=seqcount){
        tmpkmers[v&KMERMASK]=seqcount;
        seqkmers[v&KMERMASK]=p+k;
      }
      if (otukmers[kmer_rev_lt[v&KMERMASK]]){
        eintarray &otuind(*otukmers[kmer_rev_lt[v&KMERMASK]]);
        for (unsigned int l=0; l<otuind.size(); ++l){
          unsigned int tl=otus.size()+otuind[l];
          unsigned int d=p+k-kmerpos[tl];
          if (d>KMERSIZE) d=KMERSIZE;
          idcount[tl]+=d;
          kmerpos[tl]=p+k;
        }
      }
      if (revtmpkmers[kmer_rev_lt[v&KMERMASK]]!=seqcount){
        revtmpkmers[kmer_rev_lt[v&KMERMASK]]=seqcount;
        revseqkmers[kmer_rev_lt[v&KMERMASK]]=s.seqlen-p-k-1-KMERSIZE;
      }
    }
  }

//  cout << "p: " << p << " seqlen: " << s.seqlen << endl;

  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k<int(s.seqlen)-int(KMERSIZE); ++k,v>>=2u){
    if (otukmers[v&KMERMASK]){
      eintarray &otuind(*otukmers[v&KMERMASK]);
      for (unsigned int l=0; l<otuind.size(); ++l){
        unsigned int tl=otuind[l];
        unsigned int d=p+k-kmerpos[tl];
        if (d>KMERSIZE) d=KMERSIZE;
        idcount[tl]+=d;
        kmerpos[tl]=p+k;
      }
    }
    if (tmpkmers[v&KMERMASK]!=seqcount){
      tmpkmers[v&KMERMASK]=seqcount;
      seqkmers[v&KMERMASK]=p+k;
    }
    if (otukmers[kmer_rev_lt[v&KMERMASK]]){
      eintarray &otuind(*otukmers[kmer_rev_lt[v&KMERMASK]]);
      for (unsigned int l=0; l<otuind.size(); ++l){
        unsigned int tl=otus.size()+otuind[l];
        unsigned int d=p+k-kmerpos[tl];
        if (d>KMERSIZE) d=KMERSIZE;
        idcount[tl]+=d;
        kmerpos[tl]=p+k;
      }
    }
    if (revtmpkmers[kmer_rev_lt[v&KMERMASK]]!=seqcount){
      revtmpkmers[kmer_rev_lt[v&KMERMASK]]=seqcount;
      revseqkmers[kmer_rev_lt[v&KMERMASK]]=s.seqlen-p-k-1-KMERSIZE;
    }
  }
}


void kmercount_single_nopos(ebasicarray<deque<uint32_t> >& otukmers,eseq& s,uint64_t *bitmask,ebasicarray<uint32_t>& idcount,eintarray& tmpkmers,int ti,int& kmercount,uint32_t& bid,short& bcount,unsigned int *akmers,unsigned long akmask)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  int k;
  long p=0;
  kmercount=0;
  bcount=0; bid=0;
  for (; p<long(s.seqlen)-32; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k<32-KMERSIZE; ++k,v>>=2u){
      unsigned long kmer=v&KMERMASK;
      if (akmers[v&akmask]==0x0u || tmpkmers[kmer]==ti) continue;
      ++kmercount;
      tmpkmers[kmer]=ti;
      deque<uint32_t> &otuind(otukmers[kmer]);
      for (unsigned int l=0; l<otuind.size(); ++l){
        uint32_t oi=otuind[l]&BMASK31;
        uint16_t bt=(bitmask[oi/64u]>>(oi%64u))&0x1ul;
        bitmask[oi/64u]|=(0x1ul<<(oi%64u));
        idcount[oi]=(idcount[oi]&bmask_lt[bt])+idc_lt[otuind[l]>>31u];
        if (bcount<idcount[oi]){
          bcount=idcount[oi];
          bid=oi;
        }
      }
    }
  }

//  cout << "p: " << p << " seqlen: " << s.seqlen << endl;

  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k<int(s.seqlen)-int(KMERSIZE); ++k,v>>=2u){
    unsigned long kmer=v&KMERMASK;
    if (akmers[v&akmask]==0x0u || tmpkmers[kmer]==ti) continue;
    ++kmercount;
    tmpkmers[kmer]=ti;
    deque<uint32_t> &otuind(otukmers[kmer]);
    for (unsigned int l=0; l<otuind.size(); ++l){
      uint32_t oi=otuind[l]&BMASK31;
      uint16_t bt=(bitmask[oi/64u]>>(oi%64u))&0x1ul;
      bitmask[oi/64u]|=(0x1ul<<(oi%64u));
      idcount[oi]=(idcount[oi]&bmask_lt[bt])+idc_lt[otuind[l]>>31u];
      if (bcount<idcount[oi]){
        bcount=idcount[oi];
        bid=oi;
      }
    }
  }
}



void otuaddkmerdiff(ebasicarray<deque<uint32_t> >& otukmers,eseq& s1,eseq& s2,eintarray& tmpkmers,int ti,int pi,unsigned int *akmers,unsigned long akmask)
{
  unsigned long v;
  long p;
  int k;
  unsigned long *pstr;
  
  pstr=reinterpret_cast<unsigned long*>(s1.seq._str);
  for (p=0; p<s1.seqlen-32; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (int k=0; k<32-KMERSIZE; ++k,v>>=2u){
      if (akmers[v&akmask]==0x0u || tmpkmers[v&KMERMASK]==-ti) continue; // only add allowed kmers 
      if (tmpkmers[v&KMERMASK]!=ti) // kmer not in new sequence
        otukmers[v&KMERMASK].push_back(pi|(1u<<31u));
      tmpkmers[v&KMERMASK]=-ti;
    }
  }
  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (int k=0; p+k<long(s1.seqlen)-long(KMERSIZE); ++k,v>>=2u){
    if (akmers[v&akmask]==0x0u || tmpkmers[v&KMERMASK]==-ti) continue; // only add allowed kmers 
    if (tmpkmers[v&KMERMASK]!=ti) // kmer not in new sequence
      otukmers[v&KMERMASK].push_back(pi|(1u<<31u));
    tmpkmers[v&KMERMASK]=-ti;
  }

  // add kmer increments in new leaf
  pstr=reinterpret_cast<unsigned long*>(s2.seq._str);
  for (p=0; p<s2.seqlen-32; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (int k=0; k<32-KMERSIZE; ++k,v>>=2u){
      if (akmers[v&akmask]==0x0u || tmpkmers[v&KMERMASK]==-ti) continue; // only add allowed kmers 
      tmpkmers[v&KMERMASK]=-ti;
      otukmers[v&KMERMASK].push_back(pi);
    }
  }
  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (int k=0; p+k<long(s2.seqlen)-long(KMERSIZE); ++k,v>>=2u){
    if (akmers[v&akmask]==0x0u || tmpkmers[v&KMERMASK]==-ti) continue; // only add allowed kmers 
    tmpkmers[v&KMERMASK]=-ti;
    otukmers[v&KMERMASK].push_back(pi);
  }
}




void otukmeradd(ebasicarray<deque<uint32_t> >& otukmers,int i,eseq& s,eintarray& tmpkmers,int ti,unsigned int *akmers,unsigned long akmask)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  long p;
  int k;
  for (p=0; p<s.seqlen-32; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (int k=0; k<32-KMERSIZE; ++k,v>>=2u){
      if (akmers[v&akmask]==0x0u) continue; // only add allowed kmers 
      if (tmpkmers[v&KMERMASK]==ti) continue;
      tmpkmers[v&KMERMASK]=ti;
//      if (otukmers[v&KMERMASK]==0x00){
//        otukmers[v&KMERMASK]=new deque<int>();
//        otukmers[v&KMERMASK]->reserve(1000);
//      }
      otukmers[v&KMERMASK].push_back(i);
    }
  }
  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (int k=0; p+k<long(s.seqlen)-long(KMERSIZE); ++k,v>>=2u){
    if (akmers[v&akmask]==0x0u) continue; // only add allowed kmers 
    if (tmpkmers[v&KMERMASK]==ti) continue;
    tmpkmers[v&KMERMASK]=ti;
//    if (otukmers[v&KMERMASK]==0x00){
//      otukmers[v&KMERMASK]=new deque<int>();
//      otukmers[v&KMERMASK]->reserve(1000);
//    }
    otukmers[v&KMERMASK].push_back(i);
  }
}

void otukmeradd(ebasicarray<deque<int> >& otukmers,int i,eseq& s,eintarray& tmpkmers,int ti,unsigned int *akmers,unsigned long akmask)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  long p;
  int k;
  for (p=0; p<s.seqlen-32; p+=32-KMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (int k=0; k<32-KMERSIZE; ++k,v>>=2u){
      if (akmers[v&akmask]==0x0u) continue; // only add allowed kmers 
      if (tmpkmers[v&KMERMASK]==ti) continue;
      tmpkmers[v&KMERMASK]=ti;
//      if (otukmers[v&KMERMASK]==0x00){
//        otukmers[v&KMERMASK]=new deque<int>();
//        otukmers[v&KMERMASK]->reserve(1000);
//      }
      otukmers[v&KMERMASK].push_back(i);
    }
  }
  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (int k=0; p+k<int(s.seqlen)-int(KMERSIZE); ++k,v>>=2u){
    if (akmers[v&akmask]==0x0u) continue; // only add allowed kmers 
    if (tmpkmers[v&KMERMASK]==ti) continue;
    tmpkmers[v&KMERMASK]=ti;
//    if (otukmers[v&KMERMASK]==0x00){
//      otukmers[v&KMERMASK]=new deque<int>();
//      otukmers[v&KMERMASK]->reserve(1000);
//    }
    otukmers[v&KMERMASK].push_back(i);
  }
}

void randomize(ernd& prnd,eintarray& seqids)
{
  for (int i=seqids.size()-1; i>0; --i){
    int r=int(prnd.uniform()*(i+1));
    if (r!=i) seqids.swap(r,i);
  }
}

/*
void seqsearch_benchmark(const estr& str2id,estrarrayof<eseq>& seqs,earray<eintarray>& otus,ebasicarray<eintarray*>& otukmers,eseq& s,eintarray& otukmerpos,eintarray& idcount,int& maskid,eintarray& kmermask,int& offset,eintarray& kmerpos,eintarray& kmerposrev,int& offset2,eintarray& kmerpos2,ebasicarray<ealigndata>& matchcounts,eintarray& truetax)
{
  eseq srev(s);
  srev.revcompl();

  matchcounts.clear();
  otukmerpos.init(otus.size()*2,-int(KMERSIZE));
  idcount.init(otus.size()*2,0);

  kmercount_both(otus.size(),otukmers,s,idcount,otukmerpos);

//  eintarray best;
//  ebasicarray<ealigndata> matchcounts;

  eintarray best;
//  best.clear();
  best.add(0);
  for (int l=1; l<idcount.size(); ++l){
//    if (idcount[l]<idcount[best[best.size()-1]] && (idcount[l]<0.8*idcount[ibest] && best.size()>=20 || best.size()==100)) continue; // worst id than the bottom of the list
    if (idcount[l]<idcount[best[best.size()-1]] && best.size()==10) continue; // worst id than the bottom of the list

    int j;
    for (j=0; j<best.size() && idcount[l]<idcount[best[j]]; ++j);

    if (best.size()<10)
      best.add(l);
    for (int tj=best.size()-1; tj>j; --tj)
      best[tj]=best[tj-1];
    best[j]=l;
  }

  eintarray seqids;
  for (int l=0; l<best.size(); ++l){
    int ibest=best[l];
    if (ibest<otus.size()){
      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest][0]) << "\t" << idcount[ibest] << endl;
      for (int l2=0; l2<otus[ibest].size(); ++l2)
        seqids.add(otus[ibest][l2]);
    }else{
      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest-otus.size()][0]) << "\t" << idcount[ibest] << " reversed" << endl;
      ibest-=otus.size();
      for (int l2=0; l2<otus[ibest].size(); ++l2)
        seqids.add(seqs.size()+otus[ibest][l2]);
    }
  }
//  randomize(seqids);
  idcount.init(seqids.size(),0);

  if (maskid+1<maskid){
    maskid=-1;
    kmermask.init(KMERSIZE,-1);
  }
  ++maskid;
  setkmermask(kmermask,s,maskid,akmers,0xFul);
  kmercount_mask(seqs,seqids,kmermask,maskid,idcount);

  best.clear();
  for (int l=0; l<idcount.size(); ++l){
    eseq& tmps(seqs.values(seqids[l]<seqs.size()?seqids[l]:seqids[l]-seqs.size()));
    ldieif(tmps.tax!=0x00 && tmps.tax->size()==0,"error in tax");
    if ((best.size()>0 && idcount[l]<idcount[best[best.size()-1]] && best.size()==20) || (tmps.tax!=0x00 && (*tmps.tax)[tmps.tax->size()-1]==truetax[truetax.size()-1])) continue; // worst id than the bottom of the list, ignore same species

    int j;
    for (j=0; j<best.size() && idcount[l]<idcount[best[j]]; ++j);

    if (best.size()<20)
      best.add(l);
    for (int tj=best.size()-1; tj>j; --tj)
      best[tj]=best[tj-1];
    best[j]=l;
  }
  if (idcount[best[0]]<=0) return;  // no hits in db found, skip query

  if (offset+s.seqlen<offset){
    offset=0;
    kmerpos.init(MAXSIZE,-1);
    kmerposrev.init(MAXSIZE,-1);
  }
  setkmerpos(kmerpos,s,offset);
  setkmerpos(kmerposrev,srev,offset);
  
  for (int l=0; l<best.size(); ++l){
    int sbest=(seqids[best[l]]<seqs.size()?seqids[best[l]]:seqids[best[l]]-seqs.size());
    if (offset2+seqs.values(sbest).seqlen<offset2){
      offset2=0;
      kmerpos2.init(MAXSIZE,-1);
    }
    setkmerpos(kmerpos2,seqs.values(sbest),offset2);

    ealigndata adata;
    adata.seqid=sbest;
    adata.revcompl=(seqids[best[l]]>=seqs.size());
    adata.kmercount=idcount[best[l]];
    if (seqids[best[l]]<seqs.size())
      seqident_local(s,offset,kmerpos,seqs.values(sbest),offset2,kmerpos2,adata); //matches,mismatches,gaps);
    else
      seqident_local(srev,offset,kmerposrev,seqs.values(sbest),offset2,kmerpos2,adata); //matches,mismatches,gaps);
    matchcounts.add(adata);
    offset2+=seqs.values(sbest).seqlen;
  }
  offset+=s.seqlen;

  heapsort(matchcounts);
}
*/

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


/*
void calc_score_bs(unsigned char* tncounts,int p,int e,edoublearray &scores,unsigned char *bscoefs,int bscount)
{
//  uint64_t *pcount=reinterpret_cast<uint64_t*>(tncounts);
  uint64_t vid;
  scores.init(scores.size(),0.0l);
  double tmp;
  char tmpstr[32];

  for (int i=p; i<e; ++i){
    vid=tncounts[i];
    switch (vid&0x3u){ 
      case 0x0u: tmp=gapcost; break;
      case 0x1u: tmp=matchcost; break;
      case 0x2u: tmp=misscost; break;
    }
    for (int j=0; j<scores.size(); ++j) scores[j]+=double(bscoefs[j+i*bscount])*tmp;
  }
}
*/

const int SEQSEGSIZE=1500;

void seqsearch(const estr& str2id,eseqdb& db,eseq& s,earray<epredinfo>& pinfoarr,esearchws& sws)
{
  t1.reset();
  pinfoarr.clear();
//  memset(sws.bitmask,0,((db.otus.size()*2)/64+1)*sizeof(uint64_t));
//  kmercount_both(db.otus.size(),db.otukmers,s,sws.idcount2,sws.bitmask,akmers,0x0Fu);
//  kmercount_both2(db.otus.size(),db.otukmers,s,sws.idcount,sws.otukmerpos,akmers,0x0Fu);
//  kmercount_both_nopos2(db.otus.size(),db.otukmers,s,sws.idcount,sws.otukmerpos,akmers,0x0Fu);
  for (long segi=0; segi<s.seqlen/SEQSEGSIZE+1; ++segi){
    epredinfo pinfo;
    long s_start=segi*SEQSEGSIZE;
    long s_end=(s_start+SEQSEGSIZE<s.seqlen?s_start+SEQSEGSIZE:s.seqlen);
    long s_len=s_end-s_start;
//    int step=1;
    int step=s_len/150+1;
//    if (s_len<150) step=
    eseq srev;
    srev.setrevcompl(s,s_start,s_end);
//    cout << "# " << str2id << " start: " << s_start << " end: " << s_end << " len: " << s_len << " seqlen: " << s.seqlen << " srev: " << srev.seqlen << endl;
  
    sws.otukmerpos.init(db.otus.size()*2,0);
    sws.idcount.init(db.otus.size()*2,0);

//    cout << "# counting" << endl;
    kmercount_both_nopos2_skip(db.otus.size(),db.otukmers,s,s_start,s_end,sws.idcount,sws.otukmerpos,akmers,0x0Fu,step);
    ti=ti*0.99+t1.lap()*0.01;
  
    eintarray best;
  //  ebasicarray<ealigndata> matchcounts;
  
    // choosing sequences for kmercounting step
  //  eintarray& best(sws.best);
  //  eintarray bestcount;
  //  best.clear();
    int l;
    for (l=0; l<sws.idcount.size(); ++l){
  //    uint16_t bt=(sws.bitmask[l/64u]>>(l%64u))&0x1ul;
      if (sws.idcount[l]>=minid1) { best.add(l); ++l; break; }
    }
    for (; l<sws.idcount.size(); ++l){
  //    uint16_t bt=(sws.bitmask[l/64u]>>(l%64u))&0x1ul;
      if (sws.idcount[l]<minid1 || (sws.idcount[l]<sws.idcount[best[best.size()-1]] && best.size()==3*topotus)) continue; // worst id than the bottom of the list
  
      int j;
      for (j=0; j<best.size() && sws.idcount[l]<sws.idcount[best[j]]; ++j);
  
      if (best.size()<3*topotus)
        best.add(l);
      for (int tj=best.size()-1; tj>j; --tj)
        best[tj]=best[tj-1];
      best[j]=l;
    }
    if (best.size()==0) { // no overlap found to any seq group representatives
      continue;
//      pinfo.seqid=-2;
//      return;
    }
    ts=ts*0.99+t1.lap()*0.01;

    eintarray seqids;
    for (int l=0; l<best.size() && l<topotus; ++l){
      int ibest=best[l];
      if (ibest<db.otus.size()){
  //      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest][0]) << "\t" << idcount[ibest] << endl;
        for (int l2=0; l2<db.otus[ibest].size() && (otulim==0 || l2<otulim); ++l2)
          seqids.add(db.otus[ibest][l2]);
      }else{
  //      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest-otus.size()][0]) << "\t" << idcount[ibest] << " reversed" << endl;
        ibest-=db.otus.size();
        for (int l2=0; l2<db.otus[ibest].size() && (otulim==0 || l2<otulim); ++l2)
          seqids.add(db.seqs.size()+db.otus[ibest][l2]);
      }
    }
  
    // add representatives from all non-chosen OTUS (may improve confidence and novel otu estimation)
  /*
    for (int l=topotus; l<best.size(); ++l){
      int ibest=best[l];
      if (ibest<db.otus.size()){
        if (db.otus[ibest].size()==0 || sws.idcount[ibest]==0) continue;
        seqids.add(db.otus[ibest][0]);
      } else {
  //      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest-otus.size()][0]) << "\t" << idcount[ibest] << " reversed" << endl;
        ibest-=db.otus.size();
        if (db.otus[ibest].size()==0 || sws.idcount[ibest]==0) continue;
        seqids.add(db.seqs.size()+db.otus[ibest][0]);
      }
    }
  */
  
//    cout << "# 2nd counting" << endl;

    randomize(sws.rng,seqids);
    sws.idcount.init(seqids.size(),0);
  
    if (sws.maskid+1u<sws.maskid){
      sws.maskid=1u;
      sws.kmermask.init(KMERSIZE,0u);
    }
    ++sws.maskid;
  
    ti2=ti2*0.99+t1.lap()*0.01;
  //  memset(sws.kmerbitmask,0,(MAXSIZE2/64+1)*sizeof(uint64_t));
  //  setkmermask(sws.kmerbitmask,s,akmers,0xFul);
  //  kmercount_mask(db.seqs,seqids,sws.kmerbitmask,sws.maskid,sws.idcount,akmers,0xFul);
//    cout << "# 2nd counting -- setkmermask" << endl;
    setkmermask(sws.kmermask,s,sws.maskid,akmers,0xFul,s_start,s_end);
//    cout << "# 2nd counting -- kmercount_mask" << endl;
    kmercount_mask(db.seqs,seqids,sws.kmermask,sws.maskid,sws.idcount);
    ts2=ts2*0.99+t1.lap()*0.01;
  
  /*
    int ibest=0;
    for (int l=1; l<sws.idcount.size(); ++l)
      if (sws.idcount[l]>sws.idcount[ibest]) ibest=l;
  
    if (sws.idcount[ibest]<1) {
      pinfo.seqid=-2;
      return;  // no hits in db found, skip query
    }
  */
  
    best.clear();
    for (l=0; l<sws.idcount.size(); ++l){
      if (sws.idcount[l]>=minid2){ best.add(l); ++l; break; }
    }
    for (; l<sws.idcount.size(); ++l){
  //    if (idcount[l]<idcount[best[best.size()-1]] && (idcount[l]<0.8*idcount[ibest] && best.size()>=20 || best.size()==100)) continue; // worse id than the bottom of the list
      if (sws.idcount[l]<minid2 || (sws.idcount[l]<sws.idcount[best[best.size()-1]] && best.size()>=tophits)) continue; // zero counts or worse id than the bottom of the list
  
      int j;
      for (j=0; j<best.size() && sws.idcount[l]<sws.idcount[best[j]]; ++j);
  
  //    if (best.size()<20 || idcount[l]<0.8*idcount[ibest] && best.size()<100)
      if (best.size()<tophits)
        best.add(l);
      for (int tj=best.size()-1; tj>j; --tj)
        best[tj]=best[tj-1];
      best[j]=l;
    }
    if (best.size()==0){
//      pinfo.seqid=-2;
//      return;  // no hits in db found, skip query
      continue;
    }

/*  
    for (l=0; l<sws.idcount.size(); ++l){
      if (db.seqs.keys(seqids[l]%db.seqs.size())==str2id){
        cout << "# topmatch2: " << sws.idcount[best[0]] << endl;
        cout << "#" << str2id << " selfmatch2: " << sws.idcount[l] << endl;
        break;
      }
    }
*/

    if (sws.offset+(unsigned int)(s_len)<sws.offset){ // need an int here otherwise the comparison is made in long and the offset is not correctly reset
      sws.offset=1u;
      sws.kmerpos.init(MAXSIZE,0u);
      sws.kmerposrev.init(MAXSIZE,0u);
    }
//    cout << "# 2nd counting -- setkmerpos" << endl;
    setkmerpos(sws.kmerpos,s,sws.offset,s_start,s_end);
//    cout << "# 2nd counting -- setkmerposrev" << endl;
    setkmerpos(sws.kmerposrev,srev,sws.offset);
   
  /*
    // add last added seqs (worst kmercounts) to improve confidence estimation by keeping lower scoring seqs for full alignment
    // TODO: might not need full alignment but just estimated alignment score from shared kmers
    eintarray chosenflag;
    chosenflag.init(seqids.size(),0);
  
    for (int l=0; l<best.size(); ++l) chosenflag[best[l]]=1;
    for (int l=0,c=seqids.size()-1; c>=0 && l<10; --c){
      if (chosenflag[c]==1 || sws.idcount[c]==0) continue; // no shared kmer or already in list
      chosenflag[c]=1;
      best.add(c);
      ++l;
    }
  */
  
  //  eintarray seqboth;
  //  seqboth.init(db.seqs.size(),-1);
//    cout << "# aligning" << endl;
  
    for (int l=0; l<best.size(); ++l){
      int sbest=seqids[best[l]]%db.seqs.size();
      eseq &sdb(db.seqs.values(sbest));
      if (sws.offset2+(unsigned int)(db.seqs.values(sbest).seqlen)<sws.offset2){  // need an unsigned int here otherwise the comparison is made in long and the offset is not correctly reset, signed int overflows are undefined so this cannot be done with signed ints either
        sws.offset2=1u;
        sws.kmerpos2.init(MAXSIZE,0);
      }
      setkmerpos(sws.kmerpos2,db.seqs.values(sbest),sws.offset2);
  
      ealigndata adata;
      adata.seqid=sbest;
      adata.revcompl=(seqids[best[l]]>=db.seqs.size());
      adata.kmercount=sws.idcount[best[l]];
    
      if (seqids[best[l]]<db.seqs.size()){
//        cout << "# forward align" << endl;
        seqident_local(str2id,db.seqs.keys(sbest),s,sws.kmerpos,sdb,adata,sws,as,s_start,s_end);
  //      seqident_global(s,sws.kmerpos,sdb,adata,sws,as);
  //      seqcalign_global_noedgegap(s,0,s.seqlen,sdb,0,sdb.seqlen,adata,sws.alignws,as);
  /*
        seqcalign_global_noedgegap(s,0,s.seqlen,sdb,0,sdb.seqlen,adata,sws.alignws,as);
        seqident_local(s,sws.kmerpos,sdb,adata2,sws,as);
  
        cout << str2id << endl;
        cout << "# " << " M:"<<adata.matches << " X:" << adata.mismatches << " G:"<<adata.gaps << " S:"<< adata._score << " "<< adata.profile << endl;
        if (!(adata==adata2))
          cout << "#!" << " M:"<<adata2.matches << " X:" << adata2.mismatches << " G:"<<adata2.gaps << " S:"<<adata2._score << " K:" << adata2.kmercount << " " << adata2.profile << endl;
  */
      }else{
//        cout << "# reverse align" << endl;
        seqident_local(str2id,db.seqs.keys(sbest),srev,sws.kmerposrev,sdb,adata,sws,as);
        // flip 
        int tmp=srev.seqlen-adata.s1+s_start; adata.s1=srev.seqlen-adata.e1+s_start; adata.e1=tmp; 
        adata.revcompl=true;
  //      seqident_global(srev,sws.kmerposrev,sdb,adata,sws,as);
  //      seqcalign_global_noedgegap(srev,0,srev.seqlen,sdb,0,sdb.seqlen,adata,sws.alignws,as);
  /*
        seqcalign_global_noedgegap(srev,0,srev.seqlen,sdb,0,sdb.seqlen,adata,sws.alignws,as);
        seqident_local(srev,sws.kmerposrev,sdb,adata2,sws,as);
  
        cout << str2id << endl;
        cout << "# " << " M:"<<adata.matches << " X:" << adata.mismatches << " G:"<<adata.gaps << " S:" <<adata._score << " " << adata.profile << endl;
        if (!(adata==adata2))
          cout << "#!" << " M:"<<adata2.matches << " X:" << adata2.mismatches << " G:"<<adata2.gaps << " S:" << adata2._score << " K:" << adata2.kmercount << " " << adata2.profile << endl;
  */
  
      }
  //    LDEBUG(D_SEQALIGNMENT,print_tncount(&tncounts[l*NCOUNT_MAXLEN],0,s.seqlen));
      if (adata.matches+adata.mismatches>0 && adata.score()>=minscore){
        adata._eval=sdb.seqlen*exp(-lambda*adata.score()); // for K-A stats we need (*s.seqlen) but this is constant
        pinfo.matchcounts.add(adata);
      }
  
      sws.offset2+=db.seqs.values(sbest).seqlen;
    }
    sws.offset+=s_len;
 
    if (pinfo.matchcounts.size()>0){
      heapsort(pinfo.matchcounts);
    //  for (int i=0; i<pinfo.matchcounts.size(); ++i)
    //    cout << "#best: " << i << " " << seqs.keys(pinfo.matchcounts[i].seqid) << " " << pinfo.matchcounts[i].score() << endl;
      pinfo.tophit=pinfo.matchcounts[pinfo.matchcounts.size()-1];
      pinfo.seqid=pinfo.tophit.seqid;
      pinfoarr.add(pinfo);
    }
  }
  ta=ta*0.99+t1.lap()*0.01;
}

void findtaxcutoff(efloatarray& taxcutoff,epredinfo& pinfo2,earrayof<double,int>& predtax,etax& tax)
{
  ealigndata refhit=pinfo2.matchcounts[pinfo2.matchcounts.size()-1];
//      cout << "# ref: " << str2id << "\t" << seqs.keys(refhit.seqid) << "\t" << refhit.score() << endl;
  for (int l=pinfo2.matchcounts.size()-1; l>=0; --l){
    int sbest=pinfo2.matchcounts[l].seqid;
    if (tax.seqs[sbest]==0x00) continue;
    eseqtax &taxhit(*tax.seqs[sbest]);
    for (int k=0; k<taxhit.tl.size(); ++k){
      if (taxhit.tl[k].tid!=predtax.keys(k) && taxcutoff[k]==0.0) taxcutoff[k]=pinfo2.matchcounts[l].identity();
    }
  }
  // if computed thresholds are lower than what was predicted for the silver hit, then increase them to silver hit values
//        for (int k=0; k<tmptaxhit.tl.size(); ++k){
//          if (taxcutoff[k]>tmptaxhit.tl[k].ncf && tmptaxhit.tl[k].ncf>0.0) taxcutoff[k]=tmptaxhit.tl[k].ncf;
//        }
}

void taxscore(earrayof<double,int>& ptax,epredinfo& pinfo,etax& tax,ebasicarray<eintarray>& taxcounts)
{
  edoublearray taxscores;
  taxscores.init(tax.names.size(),0.0);
  for (int i=0; i<tax.names.size(); ++i)
    ptax.add(-1,0.0);

  int tophitl=-1;
  for (int l=pinfo.matchcounts.size()-1; l>=0; --l){
    int sbest=pinfo.matchcounts[l].seqid;
    if (tax.seqs[sbest]!=0x00) { tophitl=l; break; }
  }
  if (tophitl==-1) return;

  ealigndata &tophit(pinfo.matchcounts[tophitl]);
  pinfo.tophit=tophit; 
  ldieif(tax.seqs[tophit.seqid]==0x00,"top hit does not have taxonomy: "+estr(tophit.seqid));
  eseqtax &toptax(*tax.seqs[tophit.seqid]);

  double topscore=tophit.score();

  int taxid=0;
  taxcounts.init(tax.names.size());
  for (int i=0; i<taxcounts.size(); ++i)
    taxcounts[i].init(tax.names[i].size(),-1);
  for (int l=pinfo.matchcounts.size()-1; l>=0; --l){
    int sbest=pinfo.matchcounts[l].seqid;
    if (tax.seqs[sbest]==0x00) continue;
//    if (pinfo.matchcounts[l].score()<=0.0) break; // do not use alignments with less than or zero score
    eseqtax &taxhit(*tax.seqs[sbest]);
    for (int k=0; k<taxhit.tl.size(); ++k){
      if (ignoreEmptyTax && tax.names[k][taxhit.tl[k].tid].len()==0) continue;
      if (taxcounts[k][taxhit.tl[k].tid]==taxid) continue;
      taxcounts[k][taxhit.tl[k].tid]=taxid;
      taxscores[k]+=exp((1.0l-topscore/pinfo.matchcounts[l].score())*sweight);
//      taxscores[k]+=exp(log(pinfo.matchcounts[l].identity()*100.0)*30.0);
    }
/*
      if (seqids[best[imc[l]]]<seqs.size())
        cout << str2id << "\t" << seqs.keys(sbest) << "\t" << matchcounts[imc[l]] << "\t" << s.seqlen << endl;
      else
        cout << str2id << "\t" << seqs.keys(sbest) << "\t" << matchcounts[imc[l]] << "\t" << s.seqlen << "\treversed" << endl;
*/
  }

  
  for (int l=0; l<toptax.tl.size(); ++l){
    ldieif(toptax.tl[l].tid>=tax.names[l].size(),estr("key out of tax: ")+tophit.seqid+" "+estr(l)+" "+toptax.tl[l].tid+" "+tax.names[l].size());
    ptax.keys(l)=toptax.tl[l].tid;
    ptax.values(l)=1.0l/taxscores[l];
//    cout << "# i: " << tophit.identity() << " exp: " << exp(log(tophit.identity()*100.0l)*30.0l) << " sum: " << taxscores[l] << " = " << pinfo.predtax2.values(l) << endl;
//    pinfo.predtax2.add(toptax[l],exp(log(tophit.score())*30.0l)/taxscores[l]);
//    cout << "\t" << tax[l].values(toptax[l]) << "\t" << exp(log(tophit.score())*30.0l)/taxscores[l] << "\t" << (noveltaxa?estr(taxcutoff[l]/hitlen):estr("-"));
  }
}

void taxscore2(earrayof<double,int>& ptax,epredinfo& pinfo,etax& tax,ebasicarray<eintarray>& taxcounts)
{
  edoublearray taxscores;
  taxscores.init(tax.names.size(),0.0);
  for (int i=0; i<tax.names.size(); ++i)
    ptax.add(-1,0.0);

  int tophitl=-1;
  for (int l=pinfo.matchcounts.size()-1; l>=0; --l){
    int sbest=pinfo.matchcounts[l].seqid;
    if (tax.seqs[sbest]!=0x00) { tophitl=l; break; }
  }
  if (tophitl==-1) return;

  ealigndata &tophit(pinfo.matchcounts[tophitl]);
  pinfo.tophit=tophit; 
  ldieif(tax.seqs[tophit.seqid]==0x00,"top hit does not have taxonomy: "+estr(tophit.seqid));
  eseqtax &toptax(*tax.seqs[tophit.seqid]);

  int taxid=0;
  taxcounts.init(tax.names.size());
  for (int i=0; i<taxcounts.size(); ++i)
    taxcounts[i].init(tax.names[i].size(),-1);
  for (int l=pinfo.matchcounts.size()-1; l>=0; --l){
    int sbest=pinfo.matchcounts[l].seqid;
    if (tax.seqs[sbest]==0x00) continue;
//    if (pinfo.matchcounts[l].score()<=0.0) break; // do not use alignments with less than or zero score
    eseqtax &taxhit(*tax.seqs[sbest]);
    for (int k=0; k<taxhit.tl.size(); ++k){
      if (taxcounts[k][taxhit.tl[k].tid]==taxid) continue;
      taxcounts[k][taxhit.tl[k].tid]=taxid;
//      taxscores[k]+=exp(log(pinfo.matchcounts[l].score())*30.0);
      taxscores[k]+=exp(log(pinfo.matchcounts[l].identity()*100.0)*30.0);
    }
/*
      if (seqids[best[imc[l]]]<seqs.size())
        cout << str2id << "\t" << seqs.keys(sbest) << "\t" << matchcounts[imc[l]] << "\t" << s.seqlen << endl;
      else
        cout << str2id << "\t" << seqs.keys(sbest) << "\t" << matchcounts[imc[l]] << "\t" << s.seqlen << "\treversed" << endl;
*/
  }

  
  for (int l=0; l<toptax.tl.size(); ++l){
    ldieif(toptax.tl[l].tid>=tax.names[l].size(),estr("key out of tax: ")+tophit.seqid+" "+estr(l)+" "+toptax.tl[l].tid+" "+tax.names[l].size());
    ptax.keys(l)=toptax.tl[l].tid;
    ptax.values(l)=exp(log(tophit.identity()*100.0l)*30.0l)/taxscores[l];
//    cout << "# i: " << tophit.identity() << " exp: " << exp(log(tophit.identity()*100.0l)*30.0l) << " sum: " << taxscores[l] << " = " << pinfo.predtax2.values(l) << endl;
//    pinfo.predtax2.add(toptax[l],exp(log(tophit.score())*30.0l)/taxscores[l]);
//    cout << "\t" << tax[l].values(toptax[l]) << "\t" << exp(log(tophit.score())*30.0l)/taxscores[l] << "\t" << (noveltaxa?estr(taxcutoff[l]/hitlen):estr("-"));
  }
}

void load_taxa(const estr& taxfile,eseqdb& db)
{
  efile f;
  estr line;
  estrarray parts,parts2;
  int notfound=0;
  f.open(taxfile,"r");
  etax& tax(db.taxa.add(etax()));
  tax.seqs.init(db.seqs.size(),0x00);
  int taxind=0;
  while (!f.eof()){
    f.readln(line);
    if (line.len()==0) continue;
    if (line[0]=='#'){
      parts=line.explode(" ");
      if (parts[0]=="#cutoff:"){
        ldieif(tax.cutoff.size()>0,"duplicate #cutoff lines!!");
        for (int i=1; i<parts.size(); ++i){
          parts2=parts[i].explode(":");
          ldieif(parts2.size()<2,"not enough parts on #cutoff line, i.e.: 0.97:0.02");
          tax.cutoff.add(parts2[0].f());
          tax.cutoffcoef.add(parts2[1].f());
        }
      }else if (parts[0]=="#name:" && parts.size()>1){
        tax.name=parts[1];
      }else if (parts[0]=="#levels:" && parts.size()>1){
        for (int i=1; i<parts.size(); ++i)
          tax.levels.add(parts[i]);
      }
      continue; 
    }
    parts=line.explode("\t");
    ldieif(parts.size()<2,"loading taxonomy, not enough fields in line: "+line);
    if (db.seqind.exists(parts[0])){
      if (parts.size()==2){ // simple taxonomy file
        parts2=parts[1].explode(";");
        eseqtax *newtax=new eseqtax();
        newtax->tl.reserve(parts2.size());
        for (int i=0; i<parts2.size(); ++i){
          if (i>=tax.ind.size() || i>=tax.names.size()){
            tax.names.add(estrarray());
            tax.ind.add(estrhashof<int>());
          }
          if (!tax.ind[i].exists(parts2[i])){
            tax.ind[i].add(parts2[i],tax.names[i].size());
            tax.names[i].add(parts2[i]);
          }
          newtax->tl.add(eseqtaxlevel(tax.ind[i][parts2[i]]));
        }
        tax.seqs[db.seqind[parts[0]]]=newtax;
      }else if (parts.size()>11) { // indirect taxonomy file with confidences
        
        if (taxind==0){
          for (int i=1; i<parts.size(); ++i)
            if (parts[i].len()==0)
              { taxind=i+1; break; }
        }
        int dbhit=(db.seqind.exists(parts[1])?db.seqind[parts[1]]:-1);
//        lerrorif(dbhit==-1,"gold hit not found: "+parts[1]);
        eseqtax *newtax=new eseqtax(parts[3].f());
        int taxfields=(parts.size()-taxind)/3;
        newtax->tl.reserve(taxfields);
//        cout << parts[0] << " " << parts[3] << " " << taxfields;
        for (int i=0; i<taxfields; ++i){
          int field=i*3+taxind;
          if (i>=tax.ind.size() || i>=tax.names.size()){
            tax.names.add(estrarray());
            tax.ind.add(estrhashof<int>());
          }
          if (!tax.ind[i].exists(parts[field])){
            tax.ind[i].add(parts[field],tax.names[i].size());
            tax.names[i].add(parts[field]);
          }
          newtax->tl.add(eseqtaxlevel(tax.ind[i][parts[field]],parts[field+2].f()));
//          cout << parts[field] << "(" << taxind[i][parts[field]] << ") " << parts[field+1] << " " << parts[field+2] << endl;
        }
//        cout << endl;
        tax.seqs[db.seqind[parts[0]]]=newtax;
      }
    }else{
      ++notfound;
    }
  }
  if (tax.name.len()==0)
    tax.name=::basename(taxfile);
  lwarnif(notfound>0,"loading taxonomy, "+estr(notfound)+" sequences not found in sequence database");
  f.close();
}

struct emtdata {
  ebasicarray<estrarrayof<eseq>*> seqs;
  ebasicarray<estrarrayof<eseq>*> sbuffer;
  emutex m;
  econdsig seqsSignal;
  econdsig sbufferSignal;
  eseqdb *seqdb;
//  ebasicarray<etax> *taxa;
  bool noveltaxa;
  bool finished;
  bool print_align;
  bool print_hits;
} mtdata;

class epredtax
{
 public:
  earrayof<double,int> tax;
  efloatarray cutoff;
  
};

estr tax2str(emtdata& mtdata,long seqid)
{
  estr tmpstr;
  for (int t=0; t<mtdata.seqdb->taxa.size(); ++t){
    etax& tax(mtdata.seqdb->taxa.at(t));
    if (tax.seqs[seqid]==0x00)
      tmpstr+="-";
    else
      for (int l=0; l<tax.seqs[seqid]->tl.size(); ++l)
        tmpstr+=estr(";")+(tax.seqs[seqid]->tl[l].tid==-1?estr("NA"):tax.names[l].values(tax.seqs[seqid]->tl[l].tid));
    tmpstr+="\t";
  }
  tmpstr.del(0,1);
  return(tmpstr);
}

typedef estr (*outfmt_fd)(const etax&,const earrayof<double,int>&,const efloatarray&);

eoption<outfmt_fd> outfmt;

float cf=0.5;

estr outfmt_simple(const etax& tax,const earrayof<double,int>& ptax,const efloatarray& mcfarr)
{
  estr res;
  for (int l=0; l<ptax.size() && l<tax.names.size(); ++l){
    if (mcfarr[l]<0.5) break;
    res+=(ptax.keys(l)==-1?estr("NA"):tax.names[l].values(ptax.keys(l)))+";";
  }
  res.del(-1);
  return(res);
}

estr outfmt_confidences(const etax& tax,const earrayof<double,int>& ptax,const efloatarray& mcfarr)
{
  estr res;
  for (int l=0; l<ptax.size() && l<tax.names.size(); ++l)
    res+=(ptax.keys(l)==-1?estr("NA"):tax.names[l].values(ptax.keys(l)))+"\t"+mcfarr[l]+"\t"+ptax.values(l)+"\t";
  res.del(-1);
  return(res);
}



void taskSeqsearch()
{
  estr outstr;
  esearchws searchws;
//  searchws.seqkmers.init(MAXSIZE,-1);
//  searchws.revseqkmers.init(MAXSIZE,-1);
  searchws.kmerbitmask=new uint64_t[MAXSIZE/64+1];
  searchws.bitmask=new uint64_t[(mtdata.seqdb->otus.size()*2)/64+1];
  searchws.otukmerpos.init(mtdata.seqdb->otus.size()*2,0);
  searchws.idcount.init(mtdata.seqdb->otus.size()*2,0);
  searchws.idcount2.init(mtdata.seqdb->otus.size()*2,0);
  searchws.kmerpos.init(MAXSIZE,0u);
  searchws.kmerpos2.init(MAXSIZE,0u);
  searchws.kmerposrev.init(MAXSIZE,0u);
  searchws.offset=1u;
  searchws.offset2=1u;
  searchws.kmermask.init(MAXSIZE,0u);
  searchws.maskid=1u;
  while(1){
    mtdata.m.lock();
    while (mtdata.seqs.size()==0 && !mtdata.finished) mtdata.seqsSignal.wait(mtdata.m);
    if (mtdata.seqs.size()==0 && mtdata.finished) {mtdata.m.unlock(); return;}

    estrarrayof<eseq> *pbuf=mtdata.seqs[0];
    mtdata.seqs.erase(0);
    mtdata.m.unlock();
    outstr.clear();

    for (int i=0; i<pbuf->size(); ++i){
      eseq& s(pbuf->values(i));
      earray<epredinfo> pinfoarr;
  //    ebasicarray<ealigndata> matchcounts;
      seqsearch(pbuf->keys(i),*mtdata.seqdb,s,pinfoarr,searchws);
      if (pinfoarr.size()==0) continue;
//      if (pinfo.tophit.seqid<0) continue;
  
      epredinfo *topinfo=&pinfoarr[0];
      for (int pi=1; pi<pinfoarr.size(); ++pi){
        if (pinfoarr[pi].tophit.score()>topinfo->tophit.score()) topinfo=&pinfoarr[pi];
      }

      epredinfo& pinfo(*topinfo);
      float taxcutoffmin=pinfo.matchcounts[0].identity();
      float bid=pinfo.tophit.identity();
      if (mtdata.seqdb->taxa.size() && mtdata.seqdb->taxa.at(0).seqs[pinfo.tophit.seqid]!=0x00){
        eseqtax &tmptaxhit(*mtdata.seqdb->taxa.at(0).seqs[pinfo.tophit.seqid]);
        // adjust id to closest gold hit
//        if (tmptaxhit.bid>0.0) bid=tmptaxhit.bid*bid;
        if (tmptaxhit.bid>0.0 && bid>tmptaxhit.bid) bid=tmptaxhit.bid;
      }
      
//        outstr+=pbuf->keys(i)+"\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+"\t"+pinfo.tophit.score()+"\t"+bid+"\t"+pinfo.tophit.matches+"\t"+pinfo.tophit.mismatches+"\t"+pinfo.tophit.gaps+"\t"+pinfo.tophit.s2+"\t"+pinfo.tophit.e2+"\t"+taxcutoffmin+"\t"+pinfo.matchcounts.size()+"\t";
  
//        outstr+=pbuf->keys(i)+"\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+"\t"+pinfo.tophit.score()+"\t"+bid+"\t"+pinfo.tophit.matches+"\t"+pinfo.tophit.mismatches+"\t"+pinfo.tophit.gaps+"\t"+(s.seqstart+pinfo.tophit.s1)+"\t"+(s.seqstart+pinfo.tophit.e1)+"\t"+pinfo.tophit.s2+"\t"+pinfo.tophit.e2+"\t"+taxcutoffmin+"\t"+pinfo.matchcounts.size()+"\t";
      outstr+=pbuf->keys(i)+"\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+"\t"+pinfo.tophit.score()+"\t"+pinfo.tophit.identity()+"\t"+pinfo.tophit.matches+"\t"+pinfo.tophit.mismatches+"\t"+pinfo.tophit.gaps+"\t"+(s.seqstart+pinfo.tophit.s1)+"\t"+(s.seqstart+pinfo.tophit.e1)+"\t"+pinfo.tophit.s2+"\t"+pinfo.tophit.e2+"\t"+(pinfo.tophit.revcompl?"-":"+")+"\t";
       
      if (mtdata.seqdb->taxa.size()==0){
        double topscore=pinfo.tophit.score();
        double tscore=0.0;
        for (int t=0; t<pinfo.matchcounts.size(); ++t)
          tscore+=exp((1.0l-topscore/pinfo.matchcounts[t].score())*sweight);
         outstr+="\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+"\t"+estr(1.0/tscore);
      }

      for (int t=0; t<mtdata.seqdb->taxa.size(); ++t){
        etax& tax(mtdata.seqdb->taxa.at(t));
    
        earrayof<double,int> ptax;
        taxscore(ptax,pinfo,tax,searchws.taxcounts);
    
        bid=pinfo.tophit.identity();
        if (mtdata.seqdb->taxa.size() && mtdata.seqdb->taxa.at(t).seqs[pinfo.tophit.seqid]!=0x00){
          eseqtax &tmptaxhit(*mtdata.seqdb->taxa.at(t).seqs[pinfo.tophit.seqid]);
          // adjust id to closest gold hit
//          if (tmptaxhit.bid>0.0 && bid>tmptaxhit.bid) bid=tmptaxhit.bid;
//          if (tmptaxhit.bid>0.0) bid=0.5*(bid<tmptaxhit.bid?bid:tmptaxhit.bid)+0.5*bid*tmptaxhit.bid; 
          if (tmptaxhit.bid>0.0) {
            bid=bid*tmptaxhit.bid;
            for (int l=0; l<ptax.size(); ++l)
              ptax.values(l)=ptax.values(l)*tmptaxhit.tl[l].cf;
//            if (ptax.values(l)>tmptaxhit.tl[l].cf) ptax.values(l)=tmptaxhit.tl[l].cf;
//          ptax.values(l)=ptax.values(l)*tmptaxhit.tl[l].cf;
          }
        }

/*
        efloatarray taxcutoff;
        if (mtdata.noveltaxa){
          taxcutoff.init(tax.names.size(),0.0);
          if (pinfo2.matchcounts.size()>0)
            findtaxcutoff(taxcutoff,pinfo2,ptax,tax);
        }
*/
    
//        int hitlen=pinfo.tophit.e2-pinfo.tophit.s2;
   
        efloatarray mcfarr;
        mcfarr.init(ptax.size());
        float lastmcf=0.0;
        //  adjust computed confidences
        for (int l=MIN(ptax.size(),tax.names.size())-1; l>=0; --l){
          ldieif(ptax.keys(l)!=-1 && ptax.keys(l)>=tax.names[l].size(),"key out of tax: "+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+" "+estr(l)+" "+ptax.keys(l)+" "+tax.names[l].size());
          mcfarr[l]=ptax.values(l);

          // if both fixed thresholds and precalculated identity thresholds exist, mix both in 3:1 ratio
          if (tax.cutoff.size()>0){ // if only fixed id threshold exists
            float ncf=(bid-tax.cutoff[l]+0.02)/tax.cutoffcoef[l];
            mcfarr[l]=ncf<ptax.values(l)?ncf:ptax.values(l);
          }
//            mcf=ptax.values(l)+0.5*(bid-adjcutoff+0.02)/0.2;
//            mcf=0.5*ptax.values(l)+0.5*(bid-adjcutoff+0.02)/0.05;
//            mcf=0.5*ptax.values(l)+0.5*(bid-adjcutoff+0.02)/0.1;
//            mcf=0.5*ptax.values(l)+0.5*(bid-adjcutoff+0.02)/0.2;

          if (mcfarr[l]>1.0) mcfarr[l]=1.0; else if (mcfarr[l]<0.0) mcfarr[l]=0.0;
          if (mcfarr[l]<lastmcf) mcfarr[l]=lastmcf;  // do not let confidences get smaller
          lastmcf=mcfarr[l];
        }

/*      
          for (int l=0; l<ptax.size() && l<tax.names.size(); ++l){
  //          outstr+=estr("\t")+(ptax.keys(l)==-1?estr("NA"):tax.names[l].values(ptax.keys(l)))+"\t"+ptax.values(l)+"\t"+(taxcutoff.size()?taxcutoff[l]:estr("-"));
            outstr+=estr("\t")+(ptax.keys(l)==-1?estr("NA"):tax.names[l].values(ptax.keys(l)))+"\t"+mcfarr[l]+"\t"+ptax.values(l);
          }
*/
        outstr+="\t";
        outstr+=outfmt.value()(tax,ptax,mcfarr);
        outstr+="\t";
      }
      outstr+="\n";
      if (mtdata.print_hits){
        outstr+="#\t";
        etax& tax(mtdata.seqdb->taxa.at(0));
        for (int l=pinfo.matchcounts.size()-1; l>0; --l){
          ealigndata& adata(pinfo.matchcounts[l]);
//          outstr+="# "+mtdata.seqdb->seqs.keys(adata.seqid)+"\t"+adata.score()+"\t"+adata.identity()+"\t"+adata.matches+"\t"+adata.mismatches+"\t"+adata.gaps+"\t"+adata.s2+"\t"+adata.e2+"\t"+tax2str(mtdata,adata.seqid)+"\n";
          outstr+=mtdata.seqdb->seqs.keys(adata.seqid)+"\t"+adata.score()+"\t"+adata.identity()+"\t"+adata.kmercount+"\t";
        }
        outstr+="\n";
      }
//      outstr+=pinfo.tophit.profile.str() + "\n";
      if (mtdata.print_align){
        pinfo.tophit.profile.inv();
        outstr+=pinfo.tophit.profile.str() + "\n";
        outstr+=pinfo.tophit.align_str(s,mtdata.seqdb->seqs.values(pinfo.tophit.seqid));
        outstr+="\n";
      }
      mtdata.m.lock();
      cout << outstr; outstr.clear();
      mtdata.m.unlock();

//      if (pinfo.tophit.pair)
//        outstr+=estr("#pair: ")+pinfo.tophit.pair->score()+"\t"+pinfo.tophit.pair->identity()+"\t"+pinfo.tophit.pair->matches+"\t"+pinfo.tophit.pair->mismatches+"\t"+pinfo.tophit.pair->gaps+"\t"+pinfo.tophit.pair->s2+"\t"+pinfo.tophit.pair->e2+"\n";
//      outstr+="# "+pinfo.tophit.profile.str()+"\n";

/*
      for (int j=pinfo.matchcounts.size()-1; j>=0 && j>=pinfo.matchcounts.size()-10; --j){
        ealigndata &ad(pinfo.matchcounts[j]);
        outstr+="# "+mtdata.seqdb->seqs.keys(ad.seqid)+" "+ad.identity()+" "+ad.score()+" ";
        ldieif(mtdata.seqdb->taxa.size()==0,"no taxa");
        ldieif(ad.seqid>mtdata.seqdb->taxa[0].seqs.size(),"seqid larger than taxa: "+estr(ad.seqid)+" > "+mtdata.seqdb->taxa[0].seqs.size());
        if (mtdata.seqdb->taxa[0].seqs[ad.seqid]==0x00) outstr+="-";
        else{ 
          eseqtax& st(*mtdata.seqdb->taxa[0].seqs[ad.seqid]);
          for (int k=0; k<st.tl.size(); ++k)
            outstr+=";"+mtdata.seqdb->taxa[0].names[k][st.tl[k].tid];
        }
        outstr+="\n";
      }
*/
    }

    mtdata.m.lock();
    cout << outstr;
    mtdata.sbuffer.add(pbuf);
    mtdata.sbufferSignal.signal();
    mtdata.m.unlock();
  }
  delete searchws.bitmask;
}

template <class T>
void heapsortr(T& arr){
  heapsort(arr);
  for (int i=0; i<arr.size()/2; ++i)
    arr.swap(i,arr.size()-1-i);
}

void actionOTUTable()
{
  ldieif(getParser().args.size()<2,"syntax: sample1.mseq [sample2.mseq ...]");

  int ti=0;
  int tl=3;
  epregisterI(ti,"<integer> choose which taxonomy to make table for, usually OTU (0) or NCBI (1)");
  epregisterI(tl,"<integer> choose which taxonomic level to make table for, usually OTU 97% (3)");
  eparseArgs();

  int taxind=-1;
  egzfile f;

  estr tmptax;
  int sind,tl1,tli;
  float cf;
  earray<earray<estrhashof<eintarray> > > tax;
  tax.add(earray<estrhashof<eintarray> >());

  eintarray mseqlines;
  earray<estr> samples;

  for (int l=1; l<getParser().args.size(); ++l){
    estr run=getParser().args[l];
    estrarray arr=run.explode("/");
    samples.add(arr[arr.size()-1]);

    mseqlines.add(0);
    f.open(getParser().args[1],"r");
    while (!f.eof() && f.readarr(arr,"\t")){
      if (arr.size()==0 || arr[0].len()==0 || arr[0][0]=='#') continue;
      if (taxind==-1){
        for (taxind=0; taxind<arr.size(); ++taxind){
          if (arr[taxind].len()==0){
            taxind=taxind+1;
  //          cout << "taxind: " << taxind << endl;
            break;
          }
        }
      }
      ++mseqlines[l-1];
      sind=taxind;
      tmptax.clear();
      tli=0;
      for (int i=taxind; i<arr.size(); i+=3){
        if (arr[i].len()==0){
          ++i;
          if (i>=arr.size()) break;
          sind=i;
          ++tli;
          if (tax.size()<=tli) tax.add(earray<estrhashof<eintarray> >());
          tmptax.clear();
        }
        tl1=(i-sind)/3;
        cf=arr[i+1].f();
  //      cout << arr[i] << " " << cf << endl;
        if (cf<0.5) continue;
        tmptax+=";"+arr[i];
  //      cout << tli << " " << tl << " " << tmptax << " " << cf << endl;
  
        earray<estrhashof<eintarray> > &taxonc(tax[tli]);
  
        if (taxonc.size()<=tl1) taxonc.add(estrhashof<eintarray>());
        estrhashof<eintarray> &taxc(taxonc[tl1]);
        if (!(taxc.exists(tmptax)))
          taxc.add(tmptax,eintarray());
        eintarray& taxciarr(taxc[tmptax]);
        while (taxciarr.size()<l) taxciarr.add(0);
        ++taxciarr[l-1];
      }
    }
  }

  estrarrayof<int> tmpabs;
//  cout << ">" << run << "." << readsample << "\t" << mseqlines << endl;


  cout << "#TotalCounts:";
  for (int l=0; l<mseqlines.size(); ++l)
    cout << "\t" << mseqlines[l];

  cout << endl;
  for (int l=0; l<samples.size(); ++l)
    cout << "\t" << samples[l];
  cout << endl;


  ldieif(ti>=tax.size(),"chosen taxonomy number does not exist, please choose a number smaller than "+estr(tax.size()));
  earray<estrhashof<eintarray> > &taxonc(tax[ti]);
  ldieif(tl>=taxonc.size(),"chosen taxonomic level does not exist, please choose a number smaller than "+estr(taxonc.size()));
  estrhashof<eintarray> &taxc(taxonc[tl]);
  for (int t=0; t<taxc.size(); ++t){
    cout << taxc.keys(t).substr(1);
    eintarray& taxciarr(taxc.values(t));
    for (int l=0; l<mseqlines.size(); ++l)
      cout << "\t" << (l<taxciarr.size()?taxciarr[l]:0);
    cout << endl;
  }

  exit(0);
}

void actionOTUCounts()
{
  ldieif(getParser().args.size()<2,"syntax: sample1.mseq");

  int taxind=-1;
  egzfile f;

  estr run=getParser().args[1];
  estrarray arr=run.explode("/");
  run=arr[arr.size()-1];

  estr tmptax;
  int sind,tl,tli;
  float cf;
  earray<earray<estrhashof<int> > > tax;
  tax.add(earray<estrhashof<int> >());

  int mseqlines=0;
  f.open(getParser().args[1],"r");
  while (!f.eof() && f.readarr(arr,"\t")){
    if (arr.size()==0 || arr[0].len()==0 || arr[0][0]=='#') continue;
    if (taxind==-1){
      for (taxind=0; taxind<arr.size(); ++taxind){
        if (arr[taxind].len()==0){
          taxind=taxind+1;
//          cout << "taxind: " << taxind << endl;
          break;
        }
      }
    }
    ++mseqlines;
    sind=taxind;
    tmptax.clear();
    tli=0;
    for (int i=taxind; i<arr.size(); i+=3){
      if (arr[i].len()==0){
        ++i;
        if (i>=arr.size()) break;
        sind=i;
        ++tli;
        if (tax.size()<=tli) tax.add(earray<estrhashof<int> >());
        tmptax.clear();
      }
      tl=(i-sind)/3;
      cf=arr[i+1].f();
//      cout << arr[i] << " " << cf << endl;
      if (cf<0.5) continue;
      tmptax+=";"+arr[i];
//      cout << tli << " " << tl << " " << tmptax << " " << cf << endl;

      earray<estrhashof<int> > &taxonc(tax[tli]);

      if (taxonc.size()<=tl) taxonc.add(estrhashof<int>());
      estrhashof<int> &taxc(taxonc[tl]);
      if (!(taxc.exists(tmptax)))
        taxc.add(tmptax,0);
      ++taxc[tmptax];
    }
  }

  estrarrayof<int> tmpabs;
//  cout << ">" << run << "." << readsample << "\t" << mseqlines << endl;
  cout << "#" << run << "\t" << mseqlines << endl;
  cout << "Taxonomy\tTaxonomyLevel\tLabel\tCounts" << endl;
  for (int i=0; i<tax.size(); ++i){
    earray<estrhashof<int> > &taxonc(tax[i]);
    for (int j=0; j<taxonc.size(); ++j){
      estrhashof<int> &taxc(taxonc[j]);
      tmpabs.clear();
      for (int k=0; k<taxc.size(); ++k)
        tmpabs.add(taxc.keys(k).substr(1),taxc.values(k));
      heapsortr(tmpabs);
      for (int k=0; k<tmpabs.size(); ++k)
        cout << i << "\t" << j << "\t" << tmpabs.keys(k) << "\t" << tmpabs.values(k) << endl;
    }
  }
  exit(0);
}


void help()
{
  printf("MAPseq v%s\n",MAPSEQ_PACKAGE_VERSION);
  printf("by Joao F. Matias Rodrigues, Thomas S. B. Schmidt, Janko Tackmann, and Christian von Mering\n");
  printf("Institute of Molecular Life Sciences, University of Zurich, Switzerland\n");
  printf("Matias Rodrigues JF, Schmidt TSB, Tackmann J, and von Mering C (2017) MAPseq: highly efficient k-mer search with confidence estimates, for rRNA sequence analysis. Bioinformatics. http://doi.org/10.1093/bioinformatics/btx517\n");
  printf("\n");
  printf("Usage:\n");
  printf("    %s input.fa [<db> <tax1> <tax2> ...]\n",efile(getParser().args[0]).basename()._str);
  printf("\n");
  printf("Classify a fasta file containing sequence reads to the default NCBI taxonomy and OTU classifications.\n");
  printf("Example: mapseq -nthreads 4 rawreads.fa\n"); 
  printf("\n"); 
  printf("Optional arguments:\n");
  printf("%10s    %s\n","-nthreads","number of threads to use [default: 4]");
  printf("\n");

//  printf("After clustering:\n");
//  printf("%10s    %s\n","-makeotus <alignment> <mergelog> <threshold>","generate an OTU file at a given threshold");
//  printf("%10s    %s\n","-makeotus_mothur <alignment> <mergelog> <threshold>","generate a MOTHUR compatible OTU file at a given threshold");
//  printf("%10s    %s\n","-makereps <alignment> <otu>","generate a fasta file of OTU representatives. Sequences chosen have the minimum average distance to other sequences in the OTU.");
//  printf("\n");

  printf("Report bugs to: joao.rodrigues@imls.uzh.ch\n");
  printf("http://meringlab.org/software/mapseq/\n");

  exit(0);
}

int emain()
{
  getParser().onHelp=help;
  cout << "# mapseq v"<< MAPSEQ_PACKAGE_VERSION << " (" << __DATE__ << ")" << endl;
  initdlt();
  bool print_hits=false;
  bool print_align=false;
  bool nocluster=false;
  estr cutoffs;
  epregister(sweight);
  epregister(nocluster);
  epregister(otulim);
  epregister(lambda);
  epregister(print_hits);
  epregister(print_align);
  epregister(minscore);
  epregister(minid1);
  epregister(minid2);
  epregister(cf);
  epregister2(as.match,"match");
  epregister2(as.mismatch,"mismatch");
  epregister2(as.gapext,"gapext");
  epregister2(as.gapopen,"gapopen");
  epregister2(as.dropoff,"dropoff");
  epregister(cutoffs);
  estr dbfilter;
  epregister(dbfilter);
  outfmt.choice=0;
  outfmt.add("confidences",outfmt_confidences);
  outfmt.add("simple",outfmt_simple);

  epregisterClassInheritance(eoption<outfmt_fd>,ebaseoption);
//  epregisterClassMethod4(eoption<outfmt_fd>,operator=,int,(const estr& val),"=");

  epregister(outfmt);

  epregisterFunc(help);
  
  epregister(ignoreEmptyTax);
//  epregister(galign);
  bool benchmark=false;
//  epregister(benchmark);

  bool noveltaxa=false;
//  epregister(noveltaxa);

  int nthreads=4;
  epregister(nthreads);

  epregister(tophits);
  epregister(topotus);
  int step=0;


  getParser().actions.add("otucounts",actionOTUCounts);
  getParser().actions.add("otutable",actionOTUTable);

//  epregister(step);

  estr ignfile;
//  epregister(ignfile);

//  epregister(kmer);
  eparseArgs();
  mtdata.print_align=print_align;
  mtdata.print_hits=print_hits;

  if(getParser().args.size()<2){
    cout << "syntax: mapseq <query> [db] [tax] [tax2] ..." << endl;
    return(0);
  }
  estr cfile;

  for (unsigned int i=0; i<MAXSIZE; ++i)
    akmers[i]=0x00u;
  akmers[0x08ul]=1u; // AG
  akmers[0x0Cul]=1u; // AC
  akmers[0x09ul]=1u; // TG
  akmers[0x0Dul]=1u; // TC
  
  for (unsigned int i=0u; i<(1u<<4u); ++i){
    if (akmers[i]==1u) continue;
    for (unsigned int j=0u; j<(1u<<4u); ++j){
      if (akmers[j]!=1u) continue;
      akmers[(i<<4u)|j]=1u;
//      cout << kmer2str((i<<4u)|j) << endl;
    }
  }

//  cout << "4mer shift" << endl;

  for (unsigned int i=0u; i<(1u<<8u); ++i){
    if (akmers[i]==1u) continue;
    for (unsigned int j=0u; j<(1u<<8u); ++j){
      if (akmers[j]!=1u) continue;
      akmers[(i<<8u)|j]=1u;
//      cout << kmer2str((i<<8u)|j) << endl;
    }
  }


  
  efile f;
  estr line;
  estrarray args;

  eseqdb db;

//  estrarrayof<eseq> seqs;
//  ebasicarray<etax> taxa;
//  ebasicarray<eintarray> taxcounts;


  estrhash ignseqs;
  if (ignfile.len()>0){
    f.open(ignfile,"r");
    while (!f.eof()){
      f.readln(line);
      ignseqs.add(line,estr());
    }
  }

//  estrarrayof<eseq> query;
/*
  while(f.readln(line)){
    if (line.len()==0 || line[0]=='#') continue;
    args=line.explode(" ");
    ldieif(args.size()<2,"id and sequence not found in line: "+line);
    seqs.add(args[0],args[1]);
  }

  f.open(argv[2],"r");
  while(f.readln(line)){
    if (line.len()==0 || line[0]=='#') continue;
    args=line.explode(" ");
    ldieif(args.size()<2,"id and sequence not found in line: "+line);
    query.add(args[0],args[1]);
  }
*/

  estr dbfile=estr(DATAPATH)+"/mapref-2.2b.fna";
  if (getParser().args.size()>2)
    dbfile=getParser().args[2];
  else if (!efile(dbfile).exists())
    dbfile=dirname(getSystem().getExecutablePath())+"/share/mapseq/mapref-2.2b.fna";
  if (!efile(dbfile).exists()) ldie("fasta db not found: "+dbfile);

//  estrhashof<int> seqind;
  estr str2id,str2seq;
  int maxseqlen=0;
  f.open(dbfile,"r");
  f.readln(line);
  while (!f.eof()){
 //f.readln(str2id)){
    str2id=line;
    ldieif(str2id.len()==0 || str2id[0]!='>',"Unexpected line: "+str2id);

    str2seq.clear();
    while (f.readln(line) && line.len() && line[0]!='>') str2seq+=line;
    str2id.del(0,1);

    int i=str2id.findchr(" \t");
    if (i!=-1l) str2id.del(i); // only keep id up to first white space
    if (dbfilter.len() && str2id!=dbfilter) continue;

    ldieif(str2seq.len()==0,"Empty sequence in database: "+str2id);
    db.seqind.add(str2id,db.seqs.size());
    db.seqs.add(str2id,eseq(str2seq));

    if (str2seq.len()>maxseqlen) maxseqlen=str2seq.len();
//    cout << ">" << str2id << endl;
//    cout << str2seq << endl;
//    cout << seqs.values(seqs.size()-1) << endl;
  }
  f.close();

  cerr << "# loaded " << db.seqs.size() << " sequences" << endl;
  ldieif(db.seqs.size()==0,"empty database");

  etimer t1;
  t1.reset();

//  earray<eintarray> otus;
//  ebasicarray<eintarray*> otukmers;

  db.otukmers.init(MAXSIZE);

  eintarray tmpkmers;
//  eintarray kmercounts;
  int i;

  if (getParser().args.size()>3) {
    for (int i=3; i<getParser().args.size(); ++i){
      cerr << "# loading taxonomy file: " << getParser().args[i] << endl;
      load_taxa(getParser().args[i],db);
    }
  }else if (getParser().args.size()<=2){ // do not automatically load taxonomy if database is specified
    if (efile(estr(DATAPATH)+"/mapref-2.2b.fna.ncbitax").exists()){
      load_taxa(estr(DATAPATH)+"/mapref-2.2b.fna.ncbitax",db);
      load_taxa(estr(DATAPATH)+"/mapref-2.2b.fna.otutax",db);
//      load_taxa(estr(DATAPATH)+"/mapref.fna.ltps119tax",db);
    }else if (efile(dirname(getSystem().getExecutablePath())+"/share/mapseq/mapref-2.2b.fna.ncbitax").exists()){
      load_taxa(dirname(getSystem().getExecutablePath())+"/share/mapseq/mapref-2.2b.fna.ncbitax",db);
      load_taxa(dirname(getSystem().getExecutablePath())+"/share/mapseq/mapref-2.2b.fna.otutax",db);
//      load_taxa(dirname(getSystem().getExecutablePath())+"/share/mapseq/mapref.fna.ltps119tax",db);
    }
  }
  if (cutoffs.size()>0 && db.taxa.size()>0){
    ldieif(db.taxa.size()>1,"more than one taxonomy provided, cant use -cutoffs to specify cutoffs");
    etax& tax(db.taxa[0]);
    estrarray parts2;
    estrarray parts=cutoffs.explode(" ");
    ldieif(parts.size()<tax.cutoff.size(),"not enough cutoff levels: "+estr(tax.cutoff.size()));
    tax.cutoff.clear(); tax.cutoffcoef.clear();
    for (int i=0; i<parts.size(); ++i){
      parts2=parts[i].explode(":");
      ldieif(parts2.size()<2,"not enough parts on cutoff argument, i.e.: 0.97:0.02");
      tax.cutoff.add(parts2[0].f());
      tax.cutoffcoef.add(parts2[1].f());
    }
  }


  
  cerr << "# taxonomies: " << db.taxa.size() << endl;
  for (int i=0; i<db.taxa.size(); ++i){
    cerr << "# tax levels: " << db.taxa[i].names.size() << endl;
//    taxcounts.init(taxa[i].names.size());
    for (int j=0; j<db.taxa[i].names.size(); ++j){
      cerr << "# tax: " << i << " level: " << j << " ("<< db.taxa[i].names[j].size() << ")" << endl;
//      taxcounts[i][j].init(taxa[i].names[j].size(),-1);
    }
  }



  eintarray kmerpos;
  ebasicarray<uint32_t> idcount;

//  eintarray seqotu;
  db.seqotu.init(db.seqs.size(),-1);

  tmpkmers.init(MAXSIZE,-1);



  const int MAXSEQS=1000000;
  uint64_t bitmask[MAXSEQS/64+1];


  cfile=dbfile+".mscluster";
  if (nocluster){
    db.otus.reserve(db.seqs.size());
    db.otus.add(eintarray(0));

    tmpkmers.init(MAXSIZE,-1);
    otukmeradd(db.otukmers,0,db.seqs.values(db.otus[0][0]),tmpkmers,0,akmers,0xFul);
    for (i=1; i<db.seqs.size(); ++i){
      eseq& s(db.seqs.values(i));
      db.otus.add(eintarray(i));
      int ibest=db.otus.size()-1;
      db.seqotu[i]=ibest;
      otukmeradd(db.otukmers,ibest,s,tmpkmers,ibest,akmers,0xFul);
    }
  } else if (dbfilter.len()==0 && efile(cfile).exists()){ // load clustering
    efile f;
    f.open(cfile,"r");
    estr line;
    tmpkmers.init(MAXSIZE,0);
    while (!f.eof() && f.readln(line)){
      if (line.len()==0) continue;
      estrarray parts(line.explode(" "));

      int tmpi,repid=-1;
      for (i=1; i<parts.size(); ++i){
        tmpi=parts[i].i();
//        if (ignseqs.size() && !ignseqs.exists(seqs.keys(tmpi)) || seqs.values(tmpi).tax!=0x00) {repid=tmpi; break;}
        if (ignseqs.size() && ignseqs.exists(db.seqs.keys(tmpi))) continue;
        repid=tmpi;
        break;
      }
      if (repid==-1) continue; // no sequence found with taxonomic annotation

      db.otus.add(eintarray());
      eintarray &tmpo(db.otus[db.otus.size()-1]);
      tmpo.reserve(parts.size()-1);
      tmpi=parts[1].i();
//      otukmeradd(otukmers,otus.size()-1,seqs.values(tmpi),tmpkmers,tmpi,akmers,0xFFFFul);
//      if (db.otus.size()==0)
//        otukmeradd(db.otukmers,db.otus.size()-1,db.seqs.values(tmpi),tmpkmers,tmpi,akmers,0xFul);
//      else{
//      memset(bitmask,0x0u,(db.seqs.size()/64+1)*sizeof(uint64_t));
      uint32_t ibest=0u;
      short bcount=0;
      int kmerlen=0;
      eseq &s(db.seqs.values(tmpi));
//      kmercount_single_nopos(db.otukmers,s,bitmask,idcount,tmpkmers,db.otus.size(),kmerlen,ibest,bcount,akmers,0xFul);
//      if (idcount.size()==0){
      otukmeradd(db.otukmers,db.otus.size()-1,s,tmpkmers,db.otus.size(),akmers,0xFul);
//        tmpkmers.init(MAXSIZE,0);
//        kmerhash(db.otukmers,db.kmerlast,s,tmpkmers,ti++,db.seqs.size(),akmers,0xFul);
//      }else
//        otuaddkmerdiff(db.otukmers,db.seqs.values(db.otus[ibest][0]),s,tmpkmers,i,db.otus.size()-1,akmers,0xFul);
//      idcount.add(0);

      

      for (i=1; i<parts.size(); ++i){
        tmpi=parts[i].i();
//        if (ignseqs.size() && !ignseqs.exists(seqs.keys(tmpi)) || seqs.values(tmpi).tax==0x00) continue; // do not add sequence if there is no taxonomic annotation
        if (ignseqs.size() && ignseqs.exists(db.seqs.keys(tmpi))) continue;
        ldieif(tmpi>=db.seqs.size(),"cluster file has more sequence ids than original file, please remove cluster file: "+cfile);
        tmpo.add(tmpi);
        db.seqotu[tmpi]=db.otus.size()-1;
      }
    }
    ldieif(db.otus.size()==0,"no clusters in database");
  }else{ // perform clustering

    fprintf(stderr,"# no clustering file found, performing clustering\n");

    // phase 1: find all cluster seeds, but do not add non-seeds in first phase
//    eintarray dbotus;
    eintarray len_si(iheapsort(db.seqs));
    db.otus.reserve(db.seqs.size()); // worst case scenario
    db.otus.add(eintarray(len_si[0]));
    db.seqotu[len_si[0]]=0;
//    dbotus.add(len_si[0]);

//    fprintf(stderr,"# adding kmers from %i\n",db.otus[0][0]);

    tmpkmers.init(MAXSIZE,-1);
    otukmeradd(db.otukmers,0,db.seqs.values(db.otus[0][0]),tmpkmers,0,akmers,0xFul);
    idcount.reserve(db.seqs.size());
    kmerpos.reserve(db.seqs.size());
    idcount.init(db.otus.size(),0u);
    kmerpos.init(db.otus.size(),0u);
//    uint32_t sp=0u;
//    uint32_t kp=8u;
  
    for (i=1; i<db.seqs.size(); ++i){
      eseq& s(db.seqs.values(len_si[i]));
      if (i%100==0) 
        fprintf(stderr,"\rseq: %i otus: %li",i,db.otus.size());
//      if (sp+1<sp){
        idcount.init(db.otus.size(),0);
//        sp=0;
//      }
//      if (kp+s.seqlen<kp){
        kmerpos.init(db.otus.size(),0u);
//        kp=8u;
//      }

//      kmercount_single3(db.otukmers,s,sp,kp,idcount,kmerpos);
//      kmercount_single2(db.otukmers,s,kp,idcount,kmerpos);
      kmercount_single(db.otukmers,s,idcount,kmerpos);
//      kp+=s.seqlen;

      int ibest=0;
      for (int l=1; l<idcount.size(); ++l){
        if (idcount[l]>idcount[ibest]) ibest=l;
      }
/*
      uint32_t bestcount=0;
      int l;
      for (l=0; l<idcount.size(); ++l){
        uint32_t ti=idcount[l]>>10u;
        uint32_t tc=idcount[l]&BMASK10;
        if (ti==sp) {
          ibest=l;
          bestcount=tc;
          break;
        }
      }
        
      for (; l<idcount.size(); ++l){
        uint32_t ti=idcount[l]>>10u;
        uint32_t tc=idcount[l]&BMASK10;
        if (ti==sp && tc>bestcount) { ibest=l; bestcount=tc; }
//        if (idcount[l]>idcount[ibest]) ibest=l;
      }
*/
      if (float(idcount[ibest])/s.seqlen >= 0.80) {
//        cout << "# " << float(idcount[ibest])/seqs.values(len_si[i]).seqlen << endl;
        // then merge sequence to otu
//        otus[ibest].add(len_si[i]);
//        seqotu[len_si[i]]=ibest;
      } else {
//        cout << "+ " << float(idcount[ibest])/seqs.values(len_si[i]).seqlen << endl;
        // create new otu
//        idcount.add(0);
//        kmerpos.add(0);
        db.otus.add(eintarray(len_si[i]));
//        otukmercount.add(kmercount);
        ibest=db.otus.size()-1;
        db.seqotu[len_si[i]]=ibest;
        otukmeradd(db.otukmers,ibest,s,tmpkmers,ibest,akmers,0xFul);
      }
//      ++sp;
    }
    fprintf(stderr,"\rseq: %i clusters: %li\n",i,db.otus.size());

    fprintf(stderr,"phase 2:\n");
    // phase 2: after finding all seeds, add all non-seeds to the most similar clusters
    for (i=0; i<db.seqotu.size(); ++i){
      eseq& s(db.seqs.values(i));
      if (i%100==0)
        fprintf(stderr,"\rseq: %i",i);
      if (db.seqotu[i]!=-1) continue;

//      if (sp+1<sp){
        idcount.init(db.otus.size(),0);
//        sp=0;
//      }
//      if (kp+s.seqlen<kp){
        kmerpos.init(db.otus.size(),0u);
//        kp=8u;
//      }

//      kmercount_single3(db.otukmers,s,sp,kp,idcount,kmerpos);
//      kmercount_single2(db.otukmers,s,kp,idcount,kmerpos);
      kmercount_single(db.otukmers,s,idcount,kmerpos);
//      kp+=s.seqlen;

      int ibest=0;
/*
      uint32_t bestcount=0;
      int l;
      for (l=0; l<idcount.size(); ++l){
        uint32_t ti=idcount[l]>>10u;
        uint32_t tc=idcount[l]&BMASK10;
        if (ti==sp) {
          ibest=l;
          bestcount=tc;
          break;
        }
      }
        
      for (; l<idcount.size(); ++l){
        uint32_t ti=idcount[l]>>10u;
        uint32_t tc=idcount[l]&BMASK10;
        if (ti==sp && tc>bestcount) { ibest=l; bestcount=tc; }
//        if (idcount[l]>idcount[ibest]) ibest=l;
      }
*/
      for (int l=1; l<idcount.size(); ++l){
        if (idcount[l]>idcount[ibest]) ibest=l;
      }
      db.otus[ibest].add(i);
      db.seqotu[i]=ibest;
//      ++sp;
    }
//    return(0);
//    db.otus.init(dbotus.size());
//    for (int i=0; i<db.seqotu.size(); ++i)
//      db.otus[db.seqotu[i]].add(i);

    if (dbfilter.len()==0){ // do not save cluster file for a filtered db
      fprintf(stderr,"\rseq: %li clusters: %li\n",db.seqs.size(),(long)db.otus.size());
      if (cfile.len()){
        efile f(cfile,"w");
        for (int i=0; i<db.otus.size(); ++i){
          f.write(estr(i));
          for (int j=0; j<db.otus[i].size(); ++j)
            f.write(" "+estr(db.otus[i][j]));
          f.write("\n");
        }
        f.close();
      }
    }
  }


  int fcount=0;
  int okmercount=0;
  for (i=0; i<db.otukmers.size(); ++i){
    if (db.otukmers[i].size()>0)
      ++okmercount;
    if (db.otukmers[i].size()>1000 && db.otukmers[i].size()>db.otus.size()*0.50){
      db.otukmers[i].clear();
      ++fcount;
    }
  }
  cerr << "# fcount: " << fcount << " otukmercount: " << okmercount << endl;


//  earray<estrarray> tax;
//  earray<estrhashof<int> > taxind;

/*
  int alen;
  int match[11];
  int j,k;
*/

  t1.reset();

/*
  eintarray seqkmers;
  eintarray revseqkmers,revtmpkmers;
  eintarray tmpkmers2;
*/

  ethreads t;
  emutex m;
 
  if (estr(getParser().args[1])=="-")
    f.open(stdin);
  else
    f.open(getParser().args[1],"r");

  mtdata.seqdb=&db;
  mtdata.noveltaxa=noveltaxa;
  mtdata.finished=false;


  t1.reset();
//  taskSeqsearchRead(f);
//  cerr << "# done processing " << " seqs (" << t1.lap()*0.001 << "s)" << endl;
  
  t.run(taskSeqsearch,evararray(),nthreads);

  mtdata.sbuffer.reserve(nthreads*2);
  for (int i=0; i<nthreads*2; ++i){
    estrarrayof<eseq> *sarr=new estrarrayof<eseq>;
    for (int j=0; j<100; ++j) sarr->add(estr(),eseq());
    mtdata.sbuffer.add(sarr);
  }

  cout << "#query\tdbhit\tbitscore\tidentity\tmatches\tmismatches\tgaps\tquery_start\tquery_end\tdbhit_start\tdbhit_end\tstrand\t";
  if (outfmt.key()=="simple"){
    for (int i=0; i<db.taxa.size(); ++i)
      cout << "\t" << db.taxa[i].name << "\t";
  }else if (outfmt.key()=="confidences"){
    for (int i=0; i<db.taxa.size(); ++i){
      etax& tax(db.taxa[i]);
      if (tax.levels.size()==0){
        if (tax.names.size())
          cout << "\t"<<tax.name<<":taxlevel0\tcombined_cf\tscore_cf";
        for (int j=1; j<tax.names.size(); ++j)
          cout << "\ttaxlevel" << j << "\tcombined_cf\tscore_cf";
      }else{
        cout << "\t"<< tax.name <<":" << tax.levels[0] << "\tcombined_cf\tscore_cf";
        for (int j=1; j<tax.levels.size(); ++j)
          cout << "\t" << tax.levels[j] << "\tcombined_cf\tscore_cf";
      }
      cout << "\t";
    }
  }
  cout << endl;


//  cerr << "# processing input... ";
  long seqcount=0l;
  fprintf(stderr,"# processing input... ");
  fflush(stderr);
  lerrorif(!f.readln(line),"unable to read query file");
  estrarrayof<eseq> *cbuf=0x00;
  int cbufind=0;
  const int MAXSEQLEN=100000;
  long seqstart;
  str2seq.reserve(MAXSEQLEN);
  while (!f.eof()){
    if (cbuf==0x00){
      mtdata.m.lock();
      while (mtdata.sbuffer.size()==0) mtdata.sbufferSignal.wait(mtdata.m);
      cbuf=mtdata.sbuffer[mtdata.sbuffer.size()-1];
      mtdata.sbuffer.erase(mtdata.sbuffer.size()-1);
      mtdata.m.unlock();
      cbufind=0;
    }

    seqstart=0l;
    str2id=line;
    ldieif(str2id.len()==0 || str2id[0]!='>',"Unexpected line: "+str2id);
    str2id.del(0,1);
    int i=str2id.findchr(" \t");
    if (i!=-1l) str2id.del(i); // only keep id up to first white space

    str2seq.clear();
    //TODO: read a limited number of nucleotides each time, so the code does not depend on line breaks. For every nucleotide "chunk", compress the nucleotides and add to sequence
    while (f.readln(line) && line.len() && line[0]!='>'){
      if (str2seq.len()+line.len()>MAXSEQLEN){
        cbuf->keys(cbufind)=str2id;
        eseq& s(cbuf->values(cbufind));
        s.setseq(str2seq);
        s.seqstart=seqstart;
        ++cbufind;
        if (cbufind==cbuf->size()){
          mtdata.m.lock();
          mtdata.seqs.add(cbuf);
          mtdata.seqsSignal.signal();
          while (mtdata.sbuffer.size()==0) mtdata.sbufferSignal.wait(mtdata.m);
          cbuf=mtdata.sbuffer[mtdata.sbuffer.size()-1];
          mtdata.sbuffer.erase(mtdata.sbuffer.size()-1);
          mtdata.m.unlock();
          cbufind=0;
        }
        seqstart+=s.seqlen;
        str2seq.clear();
      }
      str2seq+=line;
    }

    cbuf->keys(cbufind)=str2id;
    eseq& s(cbuf->values(cbufind));
    s.setseq(str2seq);
    s.seqstart=seqstart;
    ++cbufind;
    if (cbufind==cbuf->size()){
      mtdata.m.lock();
      mtdata.seqs.add(cbuf);
      mtdata.seqsSignal.signal();
      mtdata.m.unlock();
      cbuf=0x00;
    }
    ++seqcount;
    if (seqcount%100==0 && isatty(2))
      fprintf(stderr,"\r# processing input... %li",seqcount);
//      fprintf(stderr,"\r# processing input... %li  ti: %f ts: %f ti2: %f ts2: %f ta: %f tdp1: %f tdp2: %f tdpfl: %f tdpmd: %f",seqcount,ti,ts,ti2,ts2,ta,tdp1,tdp2,tdpfl,tdpmd);
 }
  mtdata.m.lock();
  if (cbuf!=0x00){
    for (int i=cbuf->size()-1; i>=cbufind; --i) cbuf->erase(i);
    mtdata.seqs.add(cbuf);
    cbuf=0x00;
  }
  mtdata.finished=true;
  mtdata.seqsSignal.broadcast();
  mtdata.m.unlock();
  
  fprintf(stderr,"\r# processing input... %li\n",seqcount);
  f.close();
  t.wait();
  cerr << "# done processing " << seqcount << " seqs (" << t1.lap()*0.001 << "s)" << endl;
  exit(0);

  return(0);
}
