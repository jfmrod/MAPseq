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
#include <eutils/edaemon.h>
#include <eutils/earrayof.h>
#include <eutils/eblockarray.h>
#include <deque>
#include <unistd.h>



#include <math.h>

#include <map>
using namespace std;

#include "mapseq-config.h"
#include "emseq_defs.h"
#include "ekmerhashmap.h"
#include "eseqdb.h"


int nthreads=4;
bool nocluster=false;
estr cutoffs;
float lambda=1.280;

bool galign=false;

eoption<outfmt_fd> outfmt;


/*
const float matchcost=2.0;
const float misscost=-1.0;
const float gapcost=-1.0;
const float gapopen=-5.0;
*/


float ti=0.0,ts=0.0,ti2=0.0,ts2=0.0,ta=0.0;
float tdp=0.0,tdp1=0.0,tdp2=0.0,tdpfl=0.0,tdpmd=0.0;
etimer t1,t2;


#include "emseq_defs.h"



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
estr ckmer2str(uint32_t kmer)
{
  estr tmp;
  tmp+="("+kmer2str(kmer)+")";
  kmer=kmer_prot_lt[kmer&PKMERMASK];
  for (int i=0; i<PKMERSIZE; i+=3,kmer>>=6u){
    tmp+=aa[kmer&((1u<<6u)-1u)];
  }
  return(tmp);
}
estr pkmer2str(uint32_t kmer)
{
  estr tmp;
  for (int i=0; i<PKMERSIZE; i+=3,kmer>>=6u){
    tmp+=aa[kmer&((1u<<6u)-1u)];
  }
  return(tmp);
}

class eseqcluster
{
 public:
  eintarray seqotu;
  ebasicarray<deque<uint32_t> > otukmers;
  earray<eintarray> otus;
};

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


eseqdb db;

/*
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

}
*/



/*
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
*/
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

void kmercount_mask_single(eseq& s,euintarray& kmermask,unsigned int maskpos,int idcount[2])
{
  int k;
  long p;
  long kmerpos,kmerposrev;
  idcount[0]=0;
  idcount[1]=0;
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  kmerpos=-long(KMERSIZE2);
  kmerposrev=-long(KMERSIZE2);
  for (p=0; p+32<s.seqlen; p+=32-KMERSIZE2){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k+KMERSIZE2<32; ++k,v>>=2u){
      if (kmermask[v&KMERMASK2]>maskpos){
        unsigned long d=p+k-kmerpos;
        if (d>KMERSIZE2) d=KMERSIZE2;
        idcount[0]+=d;
        kmerpos=p+k;
      }
      if (kmermask[kmer_rev_lt2[v&KMERMASK2]]>maskpos){
        unsigned long d=p+k-kmerposrev;
        if (d>KMERSIZE2) d=KMERSIZE2;
        idcount[1]+=d;
        kmerposrev=p+k;
      }
    }
  }
  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k+KMERSIZE2<s.seqlen; ++k,v>>=2u){
    if (kmermask[v&KMERMASK2]>maskpos){
      unsigned long d=p+k-kmerpos;
      if (d>KMERSIZE2) d=KMERSIZE2;
      idcount[0]+=d;
      kmerpos=p+k;
    }
    if (kmermask[kmer_rev_lt2[v&KMERMASK2]]>maskpos){
      unsigned long d=p+k-kmerposrev;
      if (d>KMERSIZE2) d=KMERSIZE2;
      idcount[1]+=d;
      kmerposrev=p+k;
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


uint32_t bmask_lt[2]={0x0u,0xffffffffu};
int16_t idc_lt[2]={1,-1};

// 6 frame translation kmer search

void kmercount_both_nopos2(int scount,ebasicarray<ekmerarray>& otukmers,eseq& s,eintarray& idcount,eintarray& kmerpos,unsigned int *akmers,unsigned long akmask)
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
        ekmerarray &otuind(otukmers[v&KMERMASK]);
        for (unsigned int l=0; l<otuind.size(); ++l)
          idcount[(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
      }
      if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
//      if (otukmers[kmer_rev_lt[v&KMERMASK]]){
        ekmerarray &otuind(otukmers[kmer_rev_lt[v&KMERMASK]]);
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
      ekmerarray &otuind(otukmers[v&KMERMASK]);
      for (unsigned int l=0; l<otuind.size(); ++l)
        idcount[(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
    }
//    if (otukmers[kmer_rev_lt[v&KMERMASK]]){
    if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
      ekmerarray &otuind(otukmers[kmer_rev_lt[v&KMERMASK]]);
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


void kmercount_single_nopos(ebasicarray<ekmerarray>& otukmers,eseq& s,uint64_t *bitmask,ebasicarray<uint32_t>& idcount,eintarray& tmpkmers,int ti,int& kmercount,uint32_t& bid,short& bcount,unsigned int *akmers,unsigned long akmask)
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
      ekmerarray &otuind(otukmers[kmer]);
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
    ekmerarray &otuind(otukmers[kmer]);
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



void otuaddkmerdiff(ebasicarray<ekmerarray>& otukmers,eseq& s1,eseq& s2,eintarray& tmpkmers,int ti,int pi,unsigned int *akmers,unsigned long akmask)
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




void otukmerprotadd(ebasicarray<ekmerarray>& otukmers,int i,eseq& s,eintarray& tmpkmers,int ti,unsigned int *akmers,unsigned long akmask)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  long p;
  int k;
  for (p=0; p<s.seqlen-32; p+=32-PKMERSIZE){
//    p=(p/3)*3;
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (int k=0; k<32-PKMERSIZE; k+=3,v>>=6u){
//      if (akmers[v&akmask]==0x0u) continue; // only add allowed kmers 
      if (tmpkmers[v&PKMERMASK]==ti) continue;
      tmpkmers[v&PKMERMASK]=ti;
//      if (otukmers[v&KMERMASK]==0x00){
//        otukmers[v&KMERMASK]=new deque<int>();
//        otukmers[v&KMERMASK]->reserve(1000);
//      }
      otukmers[v&PKMERMASK].push_back(i);
    }
  }
//  p=(p/3)*3;
  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (int k=0; p+k<long(s.seqlen)-long(PKMERSIZE); k+=3,v>>=6u){
//    if (akmers[v&akmask]==0x0u) continue; // only add allowed kmers 
    if (tmpkmers[v&PKMERMASK]==ti) continue;
    tmpkmers[v&PKMERMASK]=ti;
//    if (otukmers[v&KMERMASK]==0x00){
//      otukmers[v&KMERMASK]=new deque<int>();
//      otukmers[v&KMERMASK]->reserve(1000);
//    }
    otukmers[v&PKMERMASK].push_back(i);
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

/*
void seqalign(eseq& s,eseq& s2,epredinfo& pinfo,esearchws& sws)
{
  long s_start=0;
  long s_end=s.seqlen;
  long s_len=s_end-s_start;
  eseq srev;
  srev.setrevcompl(s,s_start,s_end);
  
  if (sws.offset+(unsigned int)(s_len)<sws.offset){ // need an int here otherwise the comparison is made in long and the offset is not correctly reset
    sws.offset=1u;
    sws.kmerpos.init(MAXSIZE,0u);
    sws.kmerposrev.init(MAXSIZE,0u);
  }
  setkmerpos(sws.kmerpos,s,sws.offset,s_start,s_end);
  setkmerpos(sws.kmerposrev,srev,sws.offset);
 
  ealigndata adata;
  adata.seqid=-1;
//  adata.revcompl=(seqids[best[l]]>=db.seqs.size());
  adata.revcompl=false;
  adata.kmercount=0;
//  adata.kmercount=sws.idcount[best[l]];
    
//  if (seqids[best[l]]<db.seqs.size()){
    seqident_local(estr(),estr(),s,sws.kmerpos,s2,adata,sws,as,s_start,s_end);
//  }else{
//    seqident_local(str2id,db.seqs.keys(sbest),srev,sws.kmerposrev,s2,adata,sws,as);
      // flip 
//    int tmp=srev.seqlen-adata.s1+s_start; adata.s1=srev.seqlen-adata.e1+s_start; adata.e1=tmp; 
//    adata.revcompl=true;
//  }
  if (adata.matches+adata.mismatches>0 && adata.score()>=minscore){
    adata._eval=s2.seqlen*exp(-lambda*adata.score()); // for K-A stats we need (*s.seqlen) but this is constant
    pinfo.matchcounts.add(adata);
  }
  sws.offset+=s_len;
}
*/


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
      if (taxhit.tl[k].tid==-1 || taxcounts[k][taxhit.tl[k].tid]==taxid) continue;
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
    ldieif(toptax.tl[l].tid>=long(tax.names[l].size()),estr("key out of tax: ")+tophit.seqid+" "+estr(l)+" "+toptax.tl[l].tid+" "+tax.names[l].size());
    ptax.keys(l)=toptax.tl[l].tid;
    ptax.values(l)=exp(log(tophit.identity()*100.0l)*30.0l)/taxscores[l];
//    cout << "# i: " << tophit.identity() << " exp: " << exp(log(tophit.identity()*100.0l)*30.0l) << " sum: " << taxscores[l] << " = " << pinfo.predtax2.values(l) << endl;
//    pinfo.predtax2.add(toptax[l],exp(log(tophit.score())*30.0l)/taxscores[l]);
//    cout << "\t" << tax[l].values(toptax[l]) << "\t" << exp(log(tophit.score())*30.0l)/taxscores[l] << "\t" << (noveltaxa?estr(taxcutoff[l]/hitlen):estr("-"));
  }
}


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
        tmpstr+=estr(";")+(tax.seqs[seqid]->tl[l].tid==-1?estr("N/A"):tax.names[l].at(tax.seqs[seqid]->tl[l].tid));
    tmpstr+="\t";
  }
  tmpstr.del(0,1);
  return(tmpstr);
}

void taskClusterCompress()
{
  estr outstr;
  eseqcluster scluster;
  estrarrayof<eseq> tmpdbseqs;
  eintarray tmpkmers;
  int iotucount;
  tmpkmers.init(MAXSIZE,-1);
  scluster.otukmers.init(MAXSIZE);

  esearchws searchws(*mtdata.seqdb);
/*
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
*/

  while(1){
    mtdata.m.lock();
    while (mtdata.seqs.size()==0 && !mtdata.finished) mtdata.seqsSignal.wait(mtdata.m);
    if (mtdata.seqs.size()==0 && mtdata.finished) {mtdata.m.unlock(); return;}

    estrarrayof<eseq> *pbuf=mtdata.seqs[0];
    mtdata.seqs.erase(0);
    mtdata.m.unlock();
    outstr.clear();
    tmpdbseqs.clear();
    estr cstr;

    for (int i=0; i<pbuf->size(); ++i){
      eseq& s(pbuf->values(i));
      earray<epredinfo> pinfoarr;

      mtdata.seqdb->seqsearch(pbuf->keys(i),s,pinfoarr,searchws);
      if (pinfoarr.size()>0 && pinfoarr[0].tophit.score()>=s.seqlen*0.8){
        outstr+=">"+pbuf->keys(i)+"\n";
        epredinfo& pinfo(pinfoarr[0]);
        pinfo.tophit.profile.inv();
        
        cstr=estr(pinfo.tophit.seqid)+"\t"+(pinfo.tophit.revcompl?"-":"+")+"\t"+pinfo.tophit.compress(s);
        for (int j=0; j<s.npos.size(); ++j)
          cstr+="\t"+estr(s.npos[j]);
        if (cstr.len()<s.seqlen*0.7)
          outstr+=cstr;
        else
          outstr+=s.print_seq();
        outstr+="\n";

        int hitotu=mtdata.seqdb->seqotu[pinfo.tophit.seqid];
        mtdata.seqdb->seqotu.add(hitotu);
        mtdata.seqdb->otus[hitotu].add(mtdata.seqdb->seqs.size());
        mtdata.seqdb->seqs.add(pbuf->keys(i),pbuf->values(i));
        continue;
      }

      mtdata.seqdb->seqotu.add(mtdata.seqdb->otus.size());
      mtdata.seqdb->otus.add(eintarray(mtdata.seqdb->seqs.size()));
      otukmeradd(mtdata.seqdb->otukmers,mtdata.seqdb->otus.size()-1,s,tmpkmers,mtdata.seqdb->otus.size(),mtdata.seqdb->akmers,0xFul);
      mtdata.seqdb->seqs.add(pbuf->keys(i),pbuf->values(i));

      outstr+=">"+pbuf->keys(i)+"\n";

      if (pinfoarr.size()>0){
        epredinfo& pinfo(pinfoarr[0]);
        pinfo.tophit.profile.inv();
        cstr=estr(pinfo.tophit.seqid)+"\t"+(pinfo.tophit.revcompl?"-":"+")+"\t"+pinfo.tophit.compress(s);
        for (int j=0; j<s.npos.size(); ++j)
          cstr+="\t"+estr(s.npos[j]);
        if (cstr.len()<s.seqlen*0.7)
          outstr+=cstr;
        else
          outstr+=s.print_seq();
        outstr+="\n";
      }else
        outstr+=s.print_seq()+"\n";
    }

    mtdata.m.lock();
    cout << outstr;
    mtdata.sbuffer.add(pbuf);
    mtdata.sbufferSignal.signal();
    mtdata.m.unlock();
  }
//  delete searchws.bitmask;
}


void taskCluster()
{
  estr outstr;
  eseqcluster scluster;
  estrarrayof<eseq> tmpdbseqs;
  eintarray tmpkmers;
  int iotucount;
  tmpkmers.init(MAXSIZE,-1);
  scluster.otukmers.init(MAXSIZE);

  esearchws searchws(*mtdata.seqdb);
//  searchws.seqkmers.init(MAXSIZE,-1);
//  searchws.revseqkmers.init(MAXSIZE,-1);
/*
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
*/
  while(1){
    mtdata.m.lock();
    while (mtdata.seqs.size()==0 && !mtdata.finished) mtdata.seqsSignal.wait(mtdata.m);
    if (mtdata.seqs.size()==0 && mtdata.finished) {mtdata.m.unlock(); return;}

/*
    iotucount=mtdata.seqdb->otus.size();
    if (mtdata.seqdb->otus.size()>scluster.otus.size()){
      scluster.seqotu=mtdata.seqdb->seqotu;
      scluster.otukmers=mtdata.seqdb->otukmers;
      scluster.otus=mtdata.seqdb->otus;
    }
*/
    estrarrayof<eseq> *pbuf=mtdata.seqs[0];
    mtdata.seqs.erase(0);
    mtdata.m.unlock();
    outstr.clear();
    tmpdbseqs.clear();

    for (int i=0; i<pbuf->size(); ++i){
      eseq& s(pbuf->values(i));
      earray<epredinfo> pinfoarr;

      mtdata.seqdb->seqsearch(pbuf->keys(i),s,pinfoarr,searchws);
      if (pinfoarr.size()>0 && pinfoarr[0].tophit.identity()>=0.98) {
        epredinfo& pinfo(pinfoarr[0]);
        cout << pbuf->keys(i) << "\t" << mtdata.seqdb->seqs.keys(pinfo.tophit.seqid) << "\t" << pinfo.tophit.score() << "\t" << pinfo.tophit.identity() << "\t" << pinfo.tophit.matches << "\t" << pinfo.tophit.mismatches << "\t" << pinfo.tophit.gaps << "\t" << (s.seqstart+pinfo.tophit.s1) << "\t" << (s.seqstart+pinfo.tophit.e1) << "\t" << pinfo.tophit.s2 << "\t" << pinfo.tophit.e2 << "\t" << (pinfo.tophit.revcompl?"-":"+") << endl;
        continue;
      }
      

      mtdata.seqdb->seqotu.add(mtdata.seqdb->otus.size());
      mtdata.seqdb->otus.add(eintarray(mtdata.seqdb->otus.size()));
      otukmeradd(mtdata.seqdb->otukmers,mtdata.seqdb->otus.size()-1,s,tmpkmers,mtdata.seqdb->otus.size(),mtdata.seqdb->akmers,0xFul);
      mtdata.seqdb->seqs.add(pbuf->keys(i),pbuf->values(i));

/*
      scluster.seqotu.add(scluster.otus.size());
      scluster.otus.add(eintarray(scluster.seqotu.size()-1));
      otukmeradd(scluster.otukmers,scluster.otus.size()-1,s,tmpkmers,scluster.otus.size(),akmers,0xFul);
      tmpdbseqs.add(pbuf->keys(i),pbuf->values(i));
*/
      cout << pbuf->keys(i) << "\t-\t" << pinfoarr.size() << "\t" << (pinfoarr.size()>0?mtdata.seqdb->seqs.keys(pinfoarr[0].tophit.seqid):estr()) << "\t" << endl;
 
/* 
      outstr+=">"+pbuf->keys(i)+"\n";
      epredinfo& pinfo(pinfoarr[0]);
      pinfo.tophit.profile.inv();
      outstr+=estr(pinfo.tophit.seqid)+"\t"+(pinfo.tophit.revcompl?"-":"+")+"\t"+pinfo.tophit.compress(s);
      for (int j=0; j<s.npos.size(); ++j)
        outstr+="\t"+estr(s.npos[j]);
      outstr+="\n";
*/
    }

    mtdata.m.lock();
/*
    if (iotucount<scluster.otus.size()){
      if (mtdata.seqdb->otus.size()>iotucount){
        ldie("mapseq not ready for multithreaded clustering");

        for (int i=iotucount; i<scluster.otus.size(); ++i){
          for (int j=iotucount; j<mtdata.seqdb->otus.size(); ++j){
            // align tmpdbseqs[i] and db.seqs[db.otus[j][0]] and check if reference should be discarded or added
//            if (seqid<0.98){
              //add to reference
//              break;
//            }
          }
        }
      }else{
        for (int i=0; i<tmpdbseqs.size(); ++i)
          mtdata.seqdb->seqs.add(tmpdbseqs.keys(i),tmpdbseqs.values(i));
        mtdata.seqdb->seqotu=scluster.seqotu;
        mtdata.seqdb->otukmers=scluster.otukmers;
        mtdata.seqdb->otus=scluster.otus;
      }
    }
*/
    cout << "# seqdb: " << mtdata.seqdb->seqs.size() << " otus: " << mtdata.seqdb->otus.size() << endl;
    cout << outstr;
    mtdata.sbuffer.add(pbuf);
    mtdata.sbufferSignal.signal();
    mtdata.m.unlock();
  }
//  delete searchws.bitmask;
}



void taskCompress()
{
  estr outstr;
  esearchws searchws(*mtdata.seqdb);
/*
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
*/
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

      mtdata.seqdb->seqsearch(pbuf->keys(i),s,pinfoarr,searchws);
      if (pinfoarr.size()==0){
        outstr+=">"+pbuf->keys(i)+"\n";
        outstr+=pbuf->values(i).print_seq()+"\n";
        continue;
      }
  
      outstr+=">"+pbuf->keys(i)+"\n";
      epredinfo& pinfo(pinfoarr[0]);
      pinfo.tophit.profile.inv();
      outstr+=estr(pinfo.tophit.seqid)+"\t"+(pinfo.tophit.revcompl?"-":"+")+"\t"+pinfo.tophit.compress(s);
      for (int j=0; j<s.npos.size(); ++j)
        outstr+="\t"+estr(s.npos[j]);
      outstr+="\n";
    }

    mtdata.m.lock();
    cout << outstr;
    mtdata.sbuffer.add(pbuf);
    mtdata.sbufferSignal.signal();
    mtdata.m.unlock();
  }
//  delete searchws.bitmask;
}


/*

void taskProtsearch()
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
  searchws.kmerpos2.reserve(MAXSIZE);
  searchws.kmerpos.init(MAXSIZE,-1);
  searchws.kmerposrev.init(MAXSIZE,-1);
  searchws.offset=0;
  searchws.kmerpos2.init(PMAXSIZE,-1);
  searchws.offset2=0;
  searchws.kmermask.init(MAXSIZE,-1);
  searchws.maskid=-1;
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
      pseqsearch(pbuf->keys(i),*mtdata.seqdb,s,pinfoarr,searchws);
      if (pinfoarr.size()==0) continue;
//      if (pinfo.tophit.seqid<0) continue;
  
      for (int pi=0; pi<pinfoarr.size(); ++pi){
        epredinfo& pinfo(pinfoarr[pi]);
        float taxcutoffmin=pinfo.matchcounts[0].identity();
        float bid=pinfo.tophit.identity();
        if (mtdata.seqdb->taxa.size() && mtdata.seqdb->taxa.at(0).seqs[pinfo.tophit.seqid]!=0x00){
          eseqtax &tmptaxhit(*mtdata.seqdb->taxa.at(0).seqs[pinfo.tophit.seqid]);
          // adjust id to closest gold hit
  //        if (tmptaxhit.bid>0.0) bid=tmptaxhit.bid*bid;
          if (tmptaxhit.bid>0.0 && bid>tmptaxhit.bid) bid=tmptaxhit.bid;
        }
      
//        outstr+=pbuf->keys(i)+"\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+"\t"+pinfo.tophit.score()+"\t"+bid+"\t"+pinfo.tophit.matches+"\t"+pinfo.tophit.mismatches+"\t"+pinfo.tophit.gaps+"\t"+pinfo.tophit.s2+"\t"+pinfo.tophit.e2+"\t"+taxcutoffmin+"\t"+pinfo.matchcounts.size()+"\t";
  
        outstr+=pbuf->keys(i)+"\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+"\t"+pinfo.tophit.score()+"\t"+bid+"\t"+pinfo.tophit.matches+"\t"+pinfo.tophit.mismatches+"\t"+pinfo.tophit.gaps+"\t"+(s.seqstart+pinfo.tophit.s1)+"\t"+(s.seqstart+pinfo.tophit.e1)+"\t"+pinfo.tophit.s2+"\t"+pinfo.tophit.e2+"\t"+taxcutoffmin+"\t"+pinfo.matchcounts.size()+"\t";
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
      
          for (int l=0; l<ptax.size() && l<tax.names.size(); ++l){
  //          outstr+=estr("\t")+(ptax.keys(l)==-1?estr("NA"):tax.names[l].values(ptax.keys(l)))+"\t"+ptax.values(l)+"\t"+(taxcutoff.size()?taxcutoff[l]:estr("-"));
            outstr+=estr("\t")+(ptax.keys(l)==-1?estr("NA"):tax.names[l].values(ptax.keys(l)))+"\t"+mcfarr[l]+"\t"+ptax.values(l);
          }
          outstr+="\t";
        }
        outstr+="\n";
        mtdata.m.lock();
        cout << outstr; outstr.clear();
        mtdata.m.unlock();

        if (mtdata.print_hits){
          outstr+="#\t";
          etax& tax(mtdata.seqdb->taxa.at(0));
          for (int l=pinfo.matchcounts.size()-1; l>0; --l){
            ealigndata& adata(pinfo.matchcounts[l]);
  //          outstr+="# "+mtdata.seqdb->seqs.keys(adata.seqid)+"\t"+adata.score()+"\t"+adata.identity()+"\t"+adata.matches+"\t"+adata.mismatches+"\t"+adata.gaps+"\t"+adata.s2+"\t"+adata.e2+"\t"+tax2str(mtdata,adata.seqid)+"\n";
            outstr+=mtdata.seqdb->seqs.keys(adata.seqid)+"\t"+adata.score()+"\t"+adata.identity()+"\t";
          }
          outstr+="\n";
        }
  //      outstr+=pinfo.tophit.profile.str() + "\n";
        if (mtdata.print_align){
          pinfo.tophit.profile.inv();
          outstr+=pinfo.tophit.profile.str() + "\n";
          outstr+=pinfo.tophit.palign_str(s,mtdata.seqdb->seqs.values(pinfo.tophit.seqid));
          outstr+="\n";
        }
//      if (pinfo.tophit.pair)
//        outstr+=estr("#pair: ")+pinfo.tophit.pair->score()+"\t"+pinfo.tophit.pair->identity()+"\t"+pinfo.tophit.pair->matches+"\t"+pinfo.tophit.pair->mismatches+"\t"+pinfo.tophit.pair->gaps+"\t"+pinfo.tophit.pair->s2+"\t"+pinfo.tophit.pair->e2+"\n";
//      outstr+="# "+pinfo.tophit.profile.str()+"\n";

      }
    }

    mtdata.m.lock();
    cout << outstr;
    mtdata.sbuffer.add(pbuf);
    mtdata.sbufferSignal.signal();
    mtdata.m.unlock();
  }
  delete searchws.bitmask;
}
*/

template <class T>
void heapsortr(T& arr){
  heapsort(arr);
  for (int i=0; i<arr.size()/2; ++i)
    arr.swap(i,arr.size()-1-i);
}

void initDB(eseqdb& db,int argi=2);
void loadSequences(eseqdb& db,int argi=2);
void loadTaxonomy(eseqdb& db,int argi=3);


void actionASVOTUTable()
{
  ldieif(getParser().args.size()<2,"syntax: sample1.asv.csv [sample2.asv.csv ...]");

  eintarray tl;
  epregisterI(tl,"[<integer>,<integer>,...] choose which level to make table for, usually 6,5 for species and 97% OTU");
  eparseArgs();
  cout << tl << endl;

  int taxind=-1;
  egzfile f;

  estr tmptax;
  int sind,tl1,tli;
  float cf;
  earray<earray<estrhashof<eintarray> > > tax;
  tax.add(earray<estrhashof<eintarray> >());

  loadSequences(db,2);
  cerr << "# loaded " << db.seqs.size() << " sequences" << endl;
  ldieif(db.seqs.size()==0,"empty database");
  loadTaxonomy(db,3);
  initDB(db,2);

  esearchws searchws(db);
  mtdata.seqdb=&db;

  eintarray mseqlines;
  earray<estr> samples;
  estr line;
  estrarray arr;
  eseq s;

  earray<estrhashof<eintarray> > otusamples;
  otusamples.init(mtdata.seqdb->taxa.size());


  efile asvotuf;

  asvotuf.open("asv.otumap","w");

  estrhash otuasv;

  estr samplef=getParser().args[1];
  estrarray samplesf=samplef.explode(",");
  for (int i=0; i<samplesf.size(); ++i){
    f.open(samplesf[i],"r");
    if (f.eof() || !f.readarr(line,arr," ")) { lerror("empty first line"); exit(0); }
    int si=samples.size();
    for (int k=1; k<arr.size(); ++k)
      samples.add(arr[k]);
    for (int k=0; k<otusamples.size(); ++k){
      for (int l=0; l<otusamples[k].size(); ++l){
        while (otusamples[k].values(l).size() < samples.size())
          otusamples[k].values(l).add(0);
      }
    }
    while (!f.eof() && f.readarr(line,arr," ")){
      earray<epredinfo> pinfoarr;
      s.setseq(arr[0]);
      db.seqsearch(arr[0],s,pinfoarr,searchws);

      if (pinfoarr.size()==0) continue;

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
        if (tmptaxhit.bid>0.0 && bid>tmptaxhit.bid) bid=tmptaxhit.bid;
      }
 
      estr taxstr;
      earray<edoublearray> taxscores;
      taxscores.init(mtdata.seqdb->taxa.size());
      for (int t=0; t<mtdata.seqdb->taxa.size(); ++t){
        etax& tax(mtdata.seqdb->taxa.at(t));
    
        efloatarray mcfarr;
        earrayof<double,int> ptax;
        efloatarray tmpmcfarr;
        earrayof<double,int> tmptax;
        taxScoreSumE(taxscores[t],pinfo,tax,searchws.taxcounts,s.seqlen);
        earrayof<double,int> mixedtax; // best mixed prediction over all taxonomies
        efloatarray mixedmcfarr;

        taxScoreE(ptax,mcfarr,pinfo.tophit,pinfo,taxscores[t],tax,s.seqlen);
        mixedtax=ptax;
        mixedmcfarr=mcfarr;
        for (int l=pinfo.matchcounts.size()-2; l>=0; --l){
          ealigndata& adata(pinfo.matchcounts[l]);
          taxScoreE(tmptax,tmpmcfarr,adata,pinfo,taxscores[t],tax,s.seqlen);
          for (int tl=0; tl<tmptax.size(); ++tl){ // choose top hit in list which is not always best aligned (when including evidence or indirect taxonomy)
            if (tmpmcfarr[tl]>mixedmcfarr[tl]){ // should check per level or just for top?
              ldieif(mixedtax.size()!=tmptax.size(),"size mismatch: "+estr(mixedtax.size())+" "+tmptax.size());
              mixedtax.values(tl)=tmptax.values(tl); 
              mixedtax.keys(tl)=tmptax.keys(tl); 
              mixedmcfarr[tl]=tmpmcfarr[tl];
            }
            if (tmpmcfarr[tl]>mcfarr[tl]){ // should check per level or just for top?
              ptax=tmptax; 
              mcfarr=tmpmcfarr;
            }
          }
        }

        estr tmpstr,tmpstrfull;
        int k=0;
        for (k=0; k<mixedtax.size(); ++k){
          if (mixedmcfarr[k]<0.5) break;
          tmpstrfull+=";";
          tmpstrfull+=tax.names[k].at(mixedtax.keys(k));
        }
        for (k=0; k<mixedtax.size(); ++k){
          if (mixedmcfarr[k]<0.5 || t<tl.size() && k>=tl[t]) break;
          tmpstr+=";";
          tmpstr+=tax.names[k].at(mixedtax.keys(k));
        }
        tmpstr.del(0,1);
        if (t>=tl.size() || k>=tl[t]) {
          if (otuasv.exists(tmpstr))
            otuasv[tmpstr]+=","+arr[0];
          else
            otuasv.add(tmpstr,arr[0]);
          if (!(otusamples[t].exists(tmpstr)))
            otusamples[t].add(tmpstr,eintarray()).init(samples.size(),0);
          for (int k=1; k<arr.size(); ++k)
            otusamples[t][tmpstr][si+k-1]+=arr[k].i();
        }else{
          if (!(otusamples[t].exists("unmapped")))
            otusamples[t].add("unmapped",eintarray()).init(samples.size(),0);
          for (int k=1; k<arr.size(); ++k)
            otusamples[t]["unmapped"][si+k-1]+=arr[k].i();
        }
        taxstr+="\t";
        taxstr+=tmpstrfull.substr(1);
      }
      asvotuf.write(arr[0] + taxstr + "\n");
      // mapseq arr[0]
      // for each samplecounts arr[1..N] add to sample OTU counts table to mapped OTU
    }
    f.close();
  }

  efile otuasvf;
  otuasvf.open("otu.asvmap","w");
  for (int i=0; i<otuasv.size(); ++i)
    otuasvf.write(otuasv.keys(i) + "\t" + otuasv.values(i) + "\n");
  otuasvf.close();

  for (int i=0; i<otusamples.size(); ++i){
    cout << "### " << i << endl;
    for (int k=0; k<samples.size(); ++k)
      cout << "\t" << samples[k];
    cout << endl;
    for (int j=0; j<otusamples[i].size(); ++j){
      cout << otusamples[i].keys(j);
      for (int m=0; m<otusamples[i].values(j).size(); ++m)
        cout << "\t" << otusamples[i].values(j)[m];
      cout << endl;
    }
  }
  exit(0);
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
    f.open(getParser().args[l],"r");
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
    f.close();
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
  int ti=-1;
  int tl=-1;
  epregister(ti);
  epregister(tl);
  eparseArgs();
  ldieif(ti==-1 && tl>=0 || ti>=0 && tl==-1,"only one tl or ti option provided, but both are needed");

  egzfile f;

  for (int tfi=1; tfi<getParser().args.size(); ++tfi){
    int taxind=-1;
    estr run=getParser().args[tfi];
    estrarray arr=run.explode("/");
    run=arr[arr.size()-1];
  
    estr tmptax;
    int sind;
    float cf;
    earray<earray<estrhashof<int> > > tax;
    tax.add(earray<estrhashof<int> >());
  
    int mseqlines=0;
    f.open(getParser().args[tfi],"r");
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
      int tli=0;
      for (int i=taxind; i<arr.size(); i+=3){
        if (arr[i].len()==0){
          ++i;
          if (i>=arr.size()) break;
          sind=i;
          ++tli;
          if (tax.size()<=tli) tax.add(earray<estrhashof<int> >());
          tmptax.clear();
        }
        int tl=(i-sind)/3;
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
    if (ti!=-1 && tl!=-1){
      cout << ">" << run << "\t" << mseqlines << endl;
      earray<estrhashof<int> > &taxonc(tax[ti]);
      estrhashof<int> &taxc(taxonc[tl]);
      tmpabs.clear();
      for (int k=0; k<taxc.size(); ++k)
        tmpabs.add(taxc.keys(k).substr(1),taxc.values(k));
      heapsortr(tmpabs);
      for (int k=0; k<tmpabs.size(); ++k)
        cout << tmpabs.keys(k) << "\t" << tmpabs.values(k) << endl;
    }else{
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
    }
    f.close();
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
  printf("\n  Paired end sequence search:\n");
  printf("    %s -paired input1.fa input2.fa [<db> <tax1> <tax2> ...]\n",efile(getParser().args[0]).basename()._str);
  printf("\n");
  printf("Classify a fasta file containing sequence reads to the default NCBI taxonomy and OTU classifications.\n");
  printf("Example: mapseq -nthreads 4 rawreads.fa\n"); 
  printf("\n"); 
  printf("Optional arguments:\n");
  printf("%20s   %5s  %s\n","-nthreads","<int>","number of threads to use [default: 4]");
  printf("\n"); 
  printf("Performance/sensitivity:\n");
  printf("%20s   %5s  %s\n","-tophits","<int>","number of reference sequences to include in alignment phase [default: 20]");
  printf("%20s   %5s  %s\n","-topotus","<int>","number of internal reference otus to include in alignment phase [default: 10]");
  printf("\n"); 
  printf("Search parameters:\n"); 
  printf("%20s   %5s  %s\n","-minscore","<int>","minimum score cutoff to consider for a classification, should be reduced when searching very small sequences, i.e.: primer search [default: 30]");
  printf("%20s   %5s  %s\n","-minid1","<int>","minimum number of shared kmers to consider hit in second phase kmer search [default: 1]");
  printf("%20s   %5s  %s\n","-minid2","<int>","minimum number of shared kmers to consider hit in alignment phase [default: 1]");
  printf("%20s   %5s  %s\n","-otulim","<int>","number of sequences per internal cluster to include in alignment phase [default: 50]");
  printf("\n"); 
  printf("Extra information:\n"); 
  printf("%20s   %5s  %s\n","-print_hits","","outputs list of top hits for each input sequence");
  printf("%20s   %5s  %s\n","-print_kmerhits","","outputs list of top kmer hits for each input sequence");
  printf("%20s   %5s  %s\n","-print_align","","outputs alignments");
  printf("\n"); 
  printf("Generating count summaries from mapseq output:\n"); 
  printf("%20s   %s\n","-otucounts","<sample1.mseq>");
  printf("%20s   %5s  %s\n","","","computes summary of classification counts from the classification output file");
  printf("%20s   %s\n","-otutable","<sample1.mseq> [sample2.mseq [...]]");
  printf("%20s   %5s  %s\n","","","generates a tsv file with taxonomic labels as rows and samples as columns from classification output files");
  printf("\n");

  printf("Report bugs to: joao.rodrigues@imls.uzh.ch\n");
  printf("http://meringlab.org/software/mapseq/\n");
  printf("http://github.org/jfmrod/MAPseq/\n");

  exit(0);
}


void loadProtSequences(eseqdb& db,int argi=2)
{
  estr dbfile=estr(DATAPATH)+"/mapref-2.2b.fna";
  if (getParser().args.size()>argi)
    dbfile=getParser().args[argi];
  else if (!efile(dbfile).exists())
    dbfile=dirname(getSystem().getExecutablePath())+"/share/mapseq/mapref-2.2b.fna";
  if (!efile(dbfile).exists()) ldie("fasta db not found: "+dbfile);

  estr str2id,str2seq,line;
  egzfile f;
  f.open(dbfile,"r");
  f.readln(line);
  while (!f.eof()){
    str2id=line;
    ldieif(str2id.len()==0 || str2id[0]!='>',"Unexpected line: "+str2id);

    str2seq.clear();
    while (f.readln(line) && line.len() && line[0]!='>') str2seq+=line;
    str2id.del(0,1);

    int i=str2id.findchr(" \t");
    if (i!=-1l) str2id.del(i); // only keep id up to first white space
//    if (dbfilter.len() && str2id!=dbfilter) continue;

    eseq s;
    s.setprot(str2seq);
    db.seqind.add(str2id,db.seqs.size());
    db.seqs.add(str2id,s);
  }
  f.close();

  db.seqotu.init(db.seqs.size(),-1);
  db.otukmers.init(PMAXSIZE);
//  estr cfile=dbfile+".mscluster";
//  if (nocluster){
    eintarray tmpkmers;
    db.otus.reserve(db.seqs.size());
    db.otus.add(eintarray(0));
    tmpkmers.init(PMAXSIZE,-1);
    otukmerprotadd(db.otukmers,0,db.seqs.values(db.otus[0][0]),tmpkmers,0,db.akmers,0xFul);
    for (int i=1; i<db.seqs.size(); ++i){
      eseq& s(db.seqs.values(i));
      db.otus.add(eintarray(i));
      int ibest=db.otus.size()-1;
      db.seqotu[i]=ibest;
      otukmerprotadd(db.otukmers,ibest,s,tmpkmers,ibest,db.akmers,0xFul);
    }
/*
//  } else if (dbfilter.len()==0 && efile(cfile).exists()){ // load clustering
  } else if (efile(cfile).exists()){ // load clustering
    loadCluster(db,cfile);
    ldieif(db.otus.size()==0,"no clusters in database");
  }else{ // perform clustering
    fprintf(stderr,"# no clustering file found, performing clustering\n");
    makeCluster(db,cfile);
  }
*/
  int fcount=0;
  int okmercount=0;
  for (int i=0; i<db.otukmers.size(); ++i){
    if (db.otukmers[i].size()>0)
      ++okmercount;
    if (db.otukmers[i].size()>1000 && db.otukmers[i].size()>db.otus.size()*0.50){
      db.otukmers[i].clear();
      ++fcount;
    }
  }
//  cerr << "# fcount: " << fcount << " otukmercount: " << okmercount << endl;
}


void loadSequences(eseqdb& db,int argi)
{
  estr dbfile=estr(DATAPATH)+"/mapref-2.2b.fna";
  if (getParser().args.size()>argi)
    dbfile=getParser().args[argi];
  else if (!efile(dbfile).exists())
    dbfile=dirname(getSystem().getExecutablePath())+"/share/mapseq/mapref-2.2b.fna";
  if (!efile(dbfile).exists()) ldie("fasta db not found: "+dbfile);

  db.loadSequences(dbfile);
}

void initDB(eseqdb& db,int argi){
  estr dbfile=estr(DATAPATH)+"/mapref-2.2b.fna";
  if (getParser().args.size()>argi)
    dbfile=getParser().args[argi];
  else if (!efile(dbfile).exists())
    dbfile=dirname(getSystem().getExecutablePath())+"/share/mapseq/mapref-2.2b.fna";
  if (!efile(dbfile).exists()) ldie("fasta db not found: "+dbfile);

  db.init(dbfile,nocluster,outfmt.value());

//  cerr << "# fcount: " << fcount << " otukmercount: " << okmercount << endl;
}


/*
int minlen=75;
float minqual=0.05;

char qlt[256];

unsigned char lt[256];
unsigned char lt16[1u<<16u];

void fastq_filter(const estr& id,const estr& qual,estr& seq){
  ldieif(qual.len()!=seq.len(),"qual and seq string length mismatch: "+id+" "+seq.len()+" "+qual.len());
  if (id.len()==0 || qual.len()==0) return;
  int lastbad=0;
  int bad=0;
  int e=seq.len();
  for (int j=0; j<qual.len(); ++j){
    ldieif(qlt[qual[j]]==-1,"unknown quality char: "+id+" '"+estr(qual[j])+"'");
    if (qlt[qual[j]]<=10 || seq[j]=='N'){
      if (lastbad==1){
        e=j-1;
        --bad;
        break;
      }
      ++bad;
      lastbad=1;
    }else
      lastbad=0;
  }
 
  if (float(bad)/e<=minqual && e>=minlen)
//    cout << '>' << id << endl;
//    cout << (e==seq.len()?seq:seq.substr(0,e)) << endl;
    seq.del(e);
  else
    seq.clear();
}


void processQueryFASTQ(eseqdb& db,const estr& fname,void (*taskfunc)()){
  estr qualchrs="!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
  for (int i=0; i<256; ++i){
    lt[i]=1;
    qlt[i]=-1;
  }
  for (int i=0; i<qualchrs.len(); ++i)
    qlt[(unsigned char)qualchrs[i]]=i;

  lt['N']=0;
  lt['n']=0;
  lt['A']=0;
  lt['T']=0;
  lt['G']=0;
  lt['C']=0;
  lt['a']=0;
  lt['t']=0;
  lt['g']=0;
  lt['c']=0;

  for (unsigned int i=0u; i<(1u<<16u); ++i)
    lt16[i]=lt[i&0xff]+lt[(i>>8u)&0xff];

  etimer t1;
  t1.reset();

  egzfile f;
  estr line,line2;
  estrarray args;

  t1.reset();
  ethreads t;
  emutex m;
  
  ldieif(fname=="-","reading from stdin not supported with paired end data");
 
  f.open(fname,"r");

  mtdata.seqdb=&db;
  mtdata.finished=false;

  t1.reset();
  t.run(taskfunc,evararray(),nthreads);

  mtdata.sbuffer.reserve(nthreads*2);
  for (int i=0; i<nthreads*2; ++i){
    estrarrayof<eseq> *sarr=new estrarrayof<eseq>;
    for (int j=0; j<100; ++j) sarr->add(estr(),eseq());
    mtdata.sbuffer.add(sarr);
  }

  estr str2id,str2seq,strqual;

  long seqcount=0l;
  fprintf(stderr,"# processing FASTQ input... ");
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
    ldieif(str2id.len()==0 || str2id[0]!='@',"Unexpected line: "+str2id);
    str2id.del(0,1);
    int i=str2id.findchr(" \t");
    if (i!=-1l) str2id.del(i); // only keep id up to first white space

    str2seq.clear();
    strqual.clear();
    //TODO: read a limited number of nucleotides each time, so the code does not depend on line breaks. For every nucleotide "chunk", compress the nucleotides and add to sequence
    while (f.readln(line) && line.len()){
      ldieif(lt[line[0]]!=0,"expected nucleotide in: "+line);
      str2seq+=line;
    }

    f.readln(line);
    ldieif(line[0]!='+',"expected + in: "+line);

    while (f.readln(line) && line.len()){
      ldieif(qlt[line[0]]==-1,"expected quality char in: "+line);
      strqual+=line;
    }

    fastq_filter(str2id,strqual,str2seq);
    if (str2seq.len()==0) { ++seqcount; continue; }

    cbuf->keys(cbufind)=str2id;
    eseq& s(cbuf->values(cbufind));
    s.setseq(str2seq);
    s.seqstart=seqstart;
    ++cbufind;
    if (cbufind==cbuf->size()){
      mtdata.m.lock();
      mtdata.seqs.add(cbuf);
      mtdata.seqsSignal.signal();
      cbuf=0x00;
      cbufind=0;
    }
    ++seqcount;
    if (seqcount%100==0 && isatty(2))
      fprintf(stderr,"\r# processing FASTQ input... %li",seqcount);
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
  
  fprintf(stderr,"\r# processing FASTQ input... %li\n",seqcount);
  f.close();
}


void processQueryFASTA(eseqdb& db,const estr& fname,void (*taskfunc)())
{
  etimer t1;
  t1.reset();

  egzfile f;
  estr line;
  estrarray args;

  t1.reset();
  ethreads t;
  emutex m;
 
  if (fname=="-")
    f.open(stdin,"r");
  else
    f.open(fname,"r");

  mtdata.seqdb=&db;
  mtdata.finished=false;

  t1.reset();
  t.run(taskfunc,evararray(),nthreads);

  mtdata.sbuffer.reserve(nthreads*2);
  for (int i=0; i<nthreads*2; ++i){
    estrarrayof<eseq> *sarr=new estrarrayof<eseq>;
    for (int j=0; j<100; ++j) sarr->add(estr(),eseq());
    mtdata.sbuffer.add(sarr);
  }


//  cerr << "# processing input... ";
  long seqcount=0l;
  fprintf(stderr,"# processing input... ");
  fflush(stderr);
  lerrorif(!f.readln(line),"unable to read query file");
  estrarrayof<eseq> *cbuf=0x00;
  int cbufind=0;
  const int MAXSEQLEN=100000;
  long seqstart;
  estr str2seq,str2id;
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
}
*/

void saveTaxonomy(eseqdb& db,etax& tax,const estr& fname)
{
  efile f;
  f.open(fname,"w");
  estr tmpstr;
  tax.name.serial(tmpstr);
  serialint(tax.levels.size(),tmpstr);
  for (int j=0; j<tax.levels.size(); ++j)
    tax.levels[j].serial(tmpstr);
  serialint(tax.names.size(),tmpstr);
  for (int i=0; i<tax.names.size(); ++i){
    serialint(tax.names[i].size(),tmpstr);
    for (int j=0; j<tax.names[i].size(); ++j)
      tax.names[i][j].serial(tmpstr);
  }
  f.write(tmpstr);
  tmpstr.clear();
/*
  for (int i=0; i<tax.names.size(); ++i){
    serialint(tax.names[i].size(),tmpstr);
    for (int j=0; j<tax.names[i].size(); ++j)
      tax.names[i][j].serial(tmpstr);
  }
  f.write(tmpstr);
  tmpstr.clear();
*/
  int staxcount=0;
  for (int i=0; i<tax.seqs.size(); ++i){
    if (tax.seqs[i]==0x00) continue;
    ++staxcount;
  }
  serialint(staxcount,tmpstr);
  f.write(tmpstr);
  tmpstr.clear();
  for (int i=0; i<tax.seqs.size(); ++i){
    if (tax.seqs[i]==0x00) continue;
    serialint(i,tmpstr);
    eseqtax &stax(*tax.seqs[i]);
    serialfloat(stax.bid,tmpstr);
    for (int j=0; j<stax.tl.size(); ++j){
      serialint(stax.tl[j].tid,tmpstr);
      serialfloat(stax.tl[j].cf,tmpstr);
    }
    f.write(tmpstr);
    tmpstr.clear();
  }
  f.close();
}

void loadTaxonomyBinary(eseqdb& db,const estr& fname)
{
  efile f;
  f.open(fname,"r");
  estr tmpstr;
  tmpstr.reserve(100000000);
  long si=0,tsi;
  f.read(tmpstr,100000000);

  etax tax;
  tsi=tax.name.unserial(tmpstr,si);
  ldieif(tsi==-1,"error unserializing name");
  si=tsi;
  int len,len2;
  si=unserialint(len,tmpstr,si);
  ldieif(si==-1,"error unserializing size");
  cout << "reading levels: " << len << endl;
  tax.levels.init(len);
  for (int i=0; i<len; ++i){
    si=tax.levels[i].unserial(tmpstr,si);
    ldieif(si==-1,"error unserializing float");
  }
  si=unserialint(len,tmpstr,si);
  ldieif(si==-1,"error unserializing size");
  cout << "reading names: " << len << endl;
  tax.names.init(len);
  for (int i=0; i<len; ++i){
    si=unserialint(len2,tmpstr,si);
    ldieif(si==-1,"error unserializing size");
    tax.names[i].init(len2);
    for (int j=0; j<len2; ++j){
      si=tax.names[i][j].unserial(tmpstr,si);
      ldieif(si==-1,"error unserializing name items");
    }
  }
  si=unserialint(len,tmpstr,si);
  ldieif(si==-1,"error unserializing size");
  tax.seqs.init(db.seqs.size(),0x00);
  int seqi,tid;
  float bid,cf;
  eseqtax *stax;
  cout << "reading seqtax: " << len << endl;
  for (int i=0; i<len; ++i){
    si=unserialint(seqi,tmpstr,si);
    si=unserialfloat(bid,tmpstr,si);
    ldieif(si==-1 || seqi<0 || seqi>=db.seqs.size(),"error loading taxonomy in binary format, si: "+estr(si)+" seqi: "+estr(seqi));
    stax=new eseqtax();
    stax->bid=bid;
    for (int j=0; j<tax.levels.size(); ++j){
      si=unserialint(seqi,tmpstr,si);
      si=unserialfloat(cf,tmpstr,si);
      ldieif(si==-1,"error loading taxonomy in binary format, si: "+estr(si));
      stax->tl.add(eseqtaxlevel(tid,cf));
    }
    tax.seqs[seqi]=stax;
  }
  db.taxa.add(tax);
  f.close();
}

void loadTaxonomy(eseqdb& db,int argi)
{
  if (getParser().args.size()>argi) {
    for (int i=argi; i<getParser().args.size(); ++i){
      cerr << "# loading taxonomy file: " << getParser().args[i] << endl;
      db.loadTaxonomy(getParser().args[i]);
    }
  }else if (getParser().args.size()<argi){ // do not automatically load taxonomy if database is specified
    if (efile(estr(DATAPATH)+"/mapref-2.2b.fna.ncbitax").exists()){
      db.loadTaxonomy(estr(DATAPATH)+"/mapref-2.2b.fna.ncbitax");
      db.loadTaxonomy(estr(DATAPATH)+"/mapref-2.2b.fna.otutax");
//      load_taxa(estr(DATAPATH)+"/mapref.fna.ltps119tax",db);
    }else if (efile(dirname(getSystem().getExecutablePath())+"/share/mapseq/mapref-2.2b.fna.ncbitax").exists()){
      db.loadTaxonomy(dirname(getSystem().getExecutablePath())+"/share/mapseq/mapref-2.2b.fna.ncbitax");
      db.loadTaxonomy(dirname(getSystem().getExecutablePath())+"/share/mapseq/mapref-2.2b.fna.otutax");
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
}

void actionCluster()
{
  cerr << "Clustering: loading database" << endl;

  db.otukmers.init(MAXSIZE);
  mtdata.seqdb=&db;

  nthreads=1;
  ethreads t;
  t.setThreads(nthreads);

//  loadSequences(db);
//  cerr << "# loaded " << db.seqs.size() << " sequences" << endl;
//  ldieif(db.seqs.size()==0,"empty database");
  if (db.processQueryFASTA(getParser().args[1],taskCluster,t)!=0)
    ldie("failed processing data");
  exit(0);
}

void actionClusterDB()
{
  cerr << "ClusterDB: loading database" << endl;

  db.otukmers.init(MAXSIZE);
  mtdata.seqdb=&db;
  mtdata.finished=false;

//  nthreads=1;

  loadSequences(db,1);
  cerr << "# loaded " << db.seqs.size() << " sequences" << endl;
  ldieif(db.seqs.size()==0,"empty database");

  ethreads t;
  t.setThreads(nthreads);
  db.makeClusterMT(t);

//  processQuery(db,getParser().args[1],taskCluster);
  exit(0);
}




void actionClusterCompress()
{
  cerr << "Compressing: denovo reference" << endl;

  db.otukmers.init(MAXSIZE);
  mtdata.seqdb=&db;
  nthreads=1;

  loadSequences(db);
  initDB(db);
  cerr << "# loaded " << db.seqs.size() << " sequences" << endl;
//  ldieif(db.seqs.size()==0,"empty database");
  ethreads t;
  t.setThreads(nthreads);

  if (db.processQueryFASTA(getParser().args[1],taskClusterCompress,t)!=0)
    ldie("process query");
  exit(0);
}


void actionPairend()
{
  cerr << "# mapseq v"<< MAPSEQ_PACKAGE_VERSION << " (" << __DATE__ << ")" << endl;
  ldieif(getParser().args.size()<3,"not enough arguments, syntax: mapseq -paired <paired1.fna> <paired2.fna> [db] [[tax1] [tax2] ...]");
  cerr << "# threads: "<< nthreads << endl;

  loadSequences(db,3);
  cerr << "# loaded " << db.seqs.size() << " sequences" << endl;
  ldieif(db.seqs.size()==0,"empty database");
  loadTaxonomy(db,4);

  initDB(db,3);

  db.printSearchHeader();

  ethreads t;
  t.setThreads(nthreads);

  if (db.processQueryPairend(getParser().args[1],getParser().args[2],taskSearchPaired,t)!=0)
    ldie("process query");

  exit(0);
}

void actionChimera()
{
  cerr << "# mapseq v"<< MAPSEQ_PACKAGE_VERSION << " (" << __DATE__ << ")" << endl;
  ldieif(getParser().args.size()<2,"not enough arguments, syntax: mapseq -chimera <seqs.fna>");

  loadSequences(db,1);
  initDB(db,1);

  ethreads t;
  t.setThreads(nthreads);

  esearchws searchws(db);
  earray<epredinfo> pinfoarr,pinfoarrfull,pinfoarr2,pinfoarrpart;
  ldieif(db.seqs.size()==0,"empty db");

  int si=0;
  if (getParser().args.size()>2){
    if (!db.seqind.exists(getParser().args[2]))
      ldie("sequence not found: "+getParser().args[2]);
    si=db.seqind[getParser().args[2]];
  }

/*
  otulim=0;
  tophits=100;
  topotus=100;
*/

  for (; si<db.seqs.size(); ++si){
    eseq s(db.seqs.values(si).subseq(100,300));
    db.seqsearch_global(db.seqs.keys(si),s,pinfoarr,searchws);
    db.seqsearch_global(db.seqs.keys(si),db.seqs.values(si),pinfoarr2,searchws);
    db.seqalign_global(db.seqs.keys(si),db.seqs.values(si),pinfoarr,pinfoarrfull,searchws);
    db.seqalign_global(db.seqs.keys(si),s,pinfoarr2,pinfoarrpart,searchws);
  //  db.seqsearch(db.seqs.keys(0),db.seqs.values(0),pinfoarr,searchws);
  
  
    ldieif(pinfoarr.size()==0,"no matches");
    ldieif(pinfoarr.size()!=pinfoarrfull.size(),"pinfoarr.size() != pinfoarrfull.size()");
    ldieif(pinfoarr[0].matchcounts.size()!=pinfoarrfull[0].matchcounts.size(),"pinfoarr.matchcounts.size() != pinfoarrfull.matchcounts.size() : "+estr(pinfoarrfull[0].matchcounts.size()));
   
    epredinfo& pinfo2(pinfoarr2[0]);
    epredinfo& pinfopart(pinfoarrpart[0]);
    for (int l=pinfo2.matchcounts.size()-1; l>=0 && l>=pinfo2.matchcounts.size()-50; --l){
      ealigndata& adata2(pinfo2.matchcounts[l]);
      ealigndata& adatap(pinfopart.matchcounts[l]);
      cout << db.seqs.keys(si)+"\tglobal "+(pinfo2.matchcounts.size()-l-1)+"\t"+db.seqs.keys(adata2.seqid)+"\t"+adata2.score()+"\t"+adata2.identity()+"\t"+adatap.score()+"\t"+adatap.identity()+"\t"+adata2.matches+"\t"+adata2.mismatches+"\t"+adata2.gaps+"\t"+(s.seqstart+adata2.s1)+"\t"+(s.seqstart+adata2.e1)+"\t"+adata2.s2+"\t"+adata2.e2+"\t"+(adata2.revcompl?"-":"+") << endl;
    }
    epredinfo& pinfo(pinfoarr[0]);
    epredinfo& pinfofull(pinfoarrfull[0]);
    for (int l=pinfo.matchcounts.size()-1; l>=0 && l>=pinfo.matchcounts.size()-50; --l){
      ealigndata& adata(pinfo.matchcounts[l]);
      ealigndata& adatafull(pinfofull.matchcounts[l]);
      cout << db.seqs.keys(si)+"\t"+(pinfo.matchcounts.size()-l-1)+"\t"+db.seqs.keys(adata.seqid)+"\t"+adata.score()+"\t"+adata.identity()+"\t"+adatafull.score()+"\t"+adatafull.identity()+"\t"+adatafull.matches+"\t"+adatafull.mismatches+"\t"+adatafull.gaps+"\t"+(s.seqstart+adatafull.s1)+"\t"+(s.seqstart+adatafull.e1)+"\t"+adatafull.s2+"\t"+adatafull.e2+"\t"+(adatafull.revcompl?"-":"+")+"\t"+adata.matches+"\t"+adata.mismatches+"\t"+adata.gaps+"\t"+(s.seqstart+adata.s1)+"\t"+(s.seqstart+adata.e1)+"\t"+adata.s2+"\t"+adata.e2+"\t"+(adata.revcompl?"-":"+") << endl;
  
  
  /*
      int tmpmatch,tmpmismatch,tmpgaps;
      double tmpid,tmpid2,tmpid3,tmpscore,tmpscore2,tmpscore3;
      for (int i=100; i+100<s.seqlen; i+=100){
        adata.partscore_global(0,i,tmpscore,tmpid,tmpmatch,tmpmismatch,tmpgaps,as);
        adata.partscore_global(i,s.seqlen,tmpscore2,tmpid2,tmpmatch,tmpmismatch,tmpgaps,as);
        adata.partscore_global(i-100,i+100,tmpscore3,tmpid3,tmpmatch,tmpmismatch,tmpgaps,as);
  
        cout << i << "\t" << tmpid << "\t" << tmpid2 << "\t" << tmpid3 << "\t" << tmpid/tmpid2 << endl;
      }
  */
    }
  }
  exit(0);
}

estr mapseq(const estr& fastr,eseqdb& db,estrhash& sids)
{
  if (!fastr.len()) return("{\"error\":\"empty query\"}"); 

  estr sid,seqstr;
  estr line;
  estrarrayof<eseq> seqs;
  eseq tmps;
  if (fastr[0]=='>'){
    estr tmpfastr(fastr);
    tmpfastr.getline(line);
    while (tmpfastr.len()){
      sid=line;
      sid.del(0,1);
      seqstr.clear();
      while (tmpfastr.getline(line) && line.len() && line[0]!='>')
        seqstr+=line;
      tmps.setseq(seqstr);
      if (seqstr.len()==0) return("{\"error\":\"missing sequence: "+sid+"\"}");
      tmps.setseq(seqstr);
      seqs.add(sid,tmps);
    }
  }else if (sids.exists(fastr) && db.seqs.findkey(sids[fastr])>=0){ // accession id to sequence in 
    seqs.add(fastr,db.seqs[sids[fastr]]);
  }else{
    estr tmpfastr(fastr);
    tmpfastr.replace("\n","");
    tmps.setseq(tmpfastr);
    seqs.add("query",tmps);
  }
  if (seqs.size()==0) return("[]");
  
  estr outstr;
  esearchws searchws(db);

  cerr << "# query: " << seqs.size() << endl;
  
  outstr.clear();
  outstr+="[";
  for (int i=0; i<seqs.size(); ++i){
    outstr+="{\"queryid\":\""+seqs.keys(i)+"\",\"hits\":[";
    eseq& s(seqs.values(i));
    earray<epredinfo> pinfoarr;
    db.seqsearch(seqs.keys(i),s,pinfoarr,searchws);
    if (pinfoarr.size()==0) { outstr+="]},"; continue; } // return("{\"error\":\"No hits found for your search\"}");
  
    epredinfo *topinfo=&pinfoarr[0];
    for (int pi=1; pi<pinfoarr.size(); ++pi){
      if (pinfoarr[pi].tophit.score()>topinfo->tophit.score()) topinfo=&pinfoarr[pi];
    }
  
    epredinfo& pinfo(*topinfo);
    float taxcutoffmin=pinfo.matchcounts[0].identity();
    float bid=pinfo.tophit.identity();
    if (db.taxa.size() && db.taxa.at(0).seqs[pinfo.tophit.seqid]!=0x00){
      eseqtax &tmptaxhit(*db.taxa.at(0).seqs[pinfo.tophit.seqid]);
      if (tmptaxhit.bid>0.0 && bid>tmptaxhit.bid) bid=tmptaxhit.bid;
    }
      
    outstr+="{\"hitid\":\""+db.seqs.keys(pinfo.tophit.seqid)+"\",\"score\":"+pinfo.tophit.score()+",\"identity\":"+pinfo.tophit.identity()+",\"matches\":"+pinfo.tophit.matches+",\"mismatches\":"+pinfo.tophit.mismatches+",\"gaps\":"+pinfo.tophit.gaps+",\"qstart\":"+(s.seqstart+pinfo.tophit.s1)+",\"qend\":"+(s.seqstart+pinfo.tophit.e1)+",\"dstart\":"+pinfo.tophit.s2+",\"dend\":"+pinfo.tophit.e2+",\"strand\":"+(pinfo.tophit.revcompl?"1":"0");
  //    outstr+="query\t"+db.seqs.keys(pinfo.tophit.seqid)+"\t"+pinfo.tophit.score()+"\t"+pinfo.tophit.identity()+"\t"+pinfo.tophit.matches+"\t"+pinfo.tophit.mismatches+"\t"+pinfo.tophit.gaps+"\t"+(s.seqstart+pinfo.tophit.s1)+"\t"+(s.seqstart+pinfo.tophit.e1)+"\t"+pinfo.tophit.s2+"\t"+pinfo.tophit.e2+"\t"+(pinfo.tophit.revcompl?"-":"+")+"\t";
       
    pinfo.tophit.profile.inv();
    outstr+=",\"sumalign\":\""+pinfo.tophit.profile.str()+"\"";
    outstr+=",\"alignment\":\""+pinfo.tophit.align_str(s,db.seqs.values(pinfo.tophit.seqid)).replace("\n","\\n")+"\"";
  
    if (db.taxa.size()==0){
      double topscore=pinfo.tophit.score();
      double tscore=0.0;
      for (int t=0; t<pinfo.matchcounts.size(); ++t)
        tscore+=exp((pinfo.matchcounts[t].score()-topscore)*sweightabs);
  //      tscore+=exp((1.0l-topscore/pinfo.matchcounts[t].score())*sweight);
      outstr+=",\"taxa\":\""+db.seqs.keys(pinfo.tophit.seqid)+"\",\"confidence\":"+(1.0/tscore);
  //      outstr+="\t"+db.seqs.keys(pinfo.tophit.seqid)+"\t"+estr(1.0/tscore);
    }
  
    edoublearray taxscores;
    for (int t=0; t<db.taxa.size(); ++t){
      etax& tax(db.taxa.at(t));
      
      earrayof<double,int> ptax;
      edoublearray taxscores;
      efloatarray mcfarr;
      taxScoreSum(taxscores,pinfo,tax,searchws.taxcounts,s.seqlen);
      taxScore(ptax,mcfarr,pinfo.tophit,pinfo,taxscores,tax,s.seqlen);
     
      if (ptax.size()>3){
        outstr+=",\"taxa\":\""+tax.names[0].at(ptax.keys(0))+";"+tax.names[1].at(ptax.keys(1))+";"+tax.names[2].at(ptax.keys(2))+";"+tax.names[3].at(ptax.keys(3))+"\",\"confidence\":"+mcfarr[3];
      }else{
        outstr+=",\"error\":\"unexpected number of taxa levels\"";
      }
      outstr+="}";
  
      efloatarray tmpmcfarr;
      earrayof<double,int> tmptax;
      for (int l=pinfo.matchcounts.size()-2; l>=0; --l){
        ealigndata& adata(pinfo.matchcounts[l]);
        taxScore(tmptax,tmpmcfarr,adata,pinfo,taxscores,tax,s.seqlen);
        outstr+=",{\"hitid\":\""+db.seqs.keys(adata.seqid)+"\",\"score\":"+adata.score()+",\"identity\":"+adata.identity()+",\"matches\":"+adata.matches+",\"mismatches\":"+adata.mismatches+",\"gaps\":"+adata.gaps+",\"qstart\":"+(s.seqstart+adata.s1)+",\"qend\":"+(s.seqstart+adata.e1)+",\"dstart\":"+adata.s2+",\"dend\":"+adata.e2+",\"strand\":"+(adata.revcompl?"1":"0");
        adata.profile.inv();
        outstr+=",\"sumalign\":\""+adata.profile.str()+"\"";
        outstr+=",\"alignment\":\""+adata.align_str(s,db.seqs.values(adata.seqid)).replace("\n","\\n")+"\"";
        outstr+=",\"taxa\":\""+tax.names[0].at(tmptax.keys(0))+";"+tax.names[1].at(tmptax.keys(1))+";"+tax.names[2].at(tmptax.keys(2))+";"+tax.names[3].at(tmptax.keys(3))+"\",\"confidence\":"+tmptax.values(3);
  //         outstr+=",\"taxa\":\""+(stax.tl.size()>3?tax.names[0].at(stax.tl[0].tid)+";"+tax.names[1].at(stax.tl[1].tid)+";"+tax.names[2].at(stax.tl[2].tid)+";"+tax.names[3].at(stax.tl[3].tid):"")+"\",\"confidence\":"+tmptax.values(3)+"}";
        outstr+="}";
      }
    }
    outstr+="]},";
  }
  outstr.del(-1);
  outstr+="]";

  return(outstr);
}

void actionDaemon()
{
  estrhash sids;
  epregisterFuncD(mapseq,evararray(evarRef(db),evarRef(sids)));
  eparseArgs();

  loadSequences(db,1);
  cerr << "# loaded " << db.seqs.size() << " sequences" << endl;
  ldieif(db.seqs.size()==0,"empty database");
  loadTaxonomy(db,2);
  initDB(db,1);

  for (int i=0; i<db.seqs.size(); ++i)
    sids.add(db.seqs.keys(i).explode(":")[0],db.seqs.keys(i));
  
  cerr << "# waiting for queries..." << endl;
  
  getSystem().run();
}

void actionProtSearch()
{
  ldieif(getParser().args.size()<3,"syntax: <query.fna> <db.faa> [db.tax]");
  loadProtSequences(db);
  cerr << "# loaded " << db.seqs.size() << " sequences" << endl;
  ldieif(db.seqs.size()==0,"empty database");

//  loadTaxonomy(db);
  db.printSearchHeader();

  ethreads t;
  t.setThreads(nthreads);

  if (db.processQueryFASTA(getParser().args[1],taskProtSearch,t)!=0)
    ldie("process query");

  exit(0);
}

void actionLoadTaxBinary()
{
  cerr << "Binary: loading taxonomy" << endl;
  ldieif(getParser().args.size()<2,"syntax: <db.mstax>");

//  nocluster=true;

  db.loadSequencesBinary(getParser().args[1]);
  cerr << "# seqs laoded: " << db.seqs.size() << endl;
  loadTaxonomyBinary(db,getParser().args[2]);
  cerr << "# taxonomies: " << db.taxa.size() << endl;
  for (int i=0; i<db.taxa.size(); ++i){
    cerr << "# tax levels: " << db.taxa[i].names.size() << endl;
    for (int j=0; j<db.taxa[i].names.size(); ++j){
      cerr << "# tax: " << i << " level: " << j << " ("<< db.taxa[i].names[j].size() << ")" << endl;
    }
  }
  initDB(db,1);
  cerr << "# init db done" << endl;

  exit(0);
}

void actionLoadBinary()
{
  cerr << "Convert: loading database" << endl;
  ldieif(getParser().args.size()<2,"syntax: <db.msfna>");

  db.loadSequencesBinary(getParser().args[1]);
  cerr << "# loaded " << db.seqs.size() << " sequences" << endl;
  ldieif(db.seqs.size()==0,"empty database");
  exit(0);
}


void actionConvertTax()
{
  cerr << "ConvertTax: loading database" << endl;
  ldieif(getParser().args.size()<2,"syntax: <db.fna> <db.tax>");

  nocluster=true;

  loadSequences(db,1);
  loadTaxonomy(db,2);
  cerr << "# loaded " << db.seqs.size() << " sequences" << endl;
  ldieif(db.seqs.size()==0,"empty database");
  initDB(db,1);
  saveTaxonomy(db,db.taxa[0],getParser().args[1]+".mstax");
  exit(0);
}


void actionConvert()
{
  cerr << "Convert: loading database" << endl;
  ldieif(getParser().args.size()<2,"syntax: <db.fna>");

  nocluster=true;

  loadSequences(db,1);
  initDB(db,1);

  cerr << "# loaded " << db.seqs.size() << " sequences" << endl;
  ldieif(db.seqs.size()==0,"empty database");
  db.saveSequences(getParser().args[1]+".msfna");
  exit(0);
}


void actionCompress()
{
  cerr << "Compressing: loading database" << endl;

  loadSequences(db);
  initDB(db);
  cerr << "# loaded " << db.seqs.size() << " sequences" << endl;
  ldieif(db.seqs.size()==0,"empty database");
  ethreads t;
  t.setThreads(nthreads);

  if (db.processQueryFASTA(getParser().args[1],taskCompress,t)!=0)
    ldie("process query");
  exit(0);
}

void actionDecompress()
{
  cerr << "loading database" << endl;

  loadSequences(db);
  initDB(db);
  cerr << "# loaded " << db.seqs.size() << " sequences" << endl;
  ldieif(db.seqs.size()==0,"empty database");

  egzfile f;
  estr line;
  estr seqstr;
  eseq seq;
  estrarray arr;
  f.open(getParser().args[1],"r");
  while (!f.eof()){
    f.readln(line);
    if (line.len()==0 || line[0]=='#') continue;
    if (line[0]=='>') { 
      cout << line << endl;
      continue;
    }
    int i=line.find("\t");
    if (i==-1) {
      cout << line << endl;
      continue;
    }
    arr=line.explode("\t");
    int sid=arr[0].i();
    ldieif(sid<0 || sid>=db.seqs.size(),"reference sequence not found, make sure to use the same database used in compression");
    seqstr=sali_decompress(arr[2],db.seqs.values(sid));
    for (int j=3; j<arr.size(); ++j){
      ldieif(arr[j].i()<0 || arr[j].i()>seqstr.len(),"N pos out of bounds");
      seqstr[arr[j].i()]='N';
    }
    if (arr[1]=="-"){
      seq.setseq(seqstr);
      seq.revcompl();
      cout << seq.print_seq() << endl;
    }else
      cout << seqstr << endl;
  }
  exit(0);
}

int emain()
{
//  epregisterClassInheritance(eoption<outfmt_fd>,ebaseoption);

  getParser().onHelp=help;
  initdlt();
  bool fastq=false;
  bool sens=false;

  epregister(sens);

  epregister(minqual);
  epregister(minlen);
  epregister(fastq);
  epregister(galign);
  epregister(sweight);
  epregister(sweightabs);
  epregister(swmin);
  epregister(swmax);
  epregister(nocluster);
  epregister2(db.otulim,"otulim");
  epregister(lambda);
  epregister2(mtdata.print_hits,"print_hits");
  epregister2(mtdata.print_kmerhits,"print_kmerhits");
  epregister2(mtdata.print_align,"print_align");
  epregister2(db.minscore,"minscore");
  epregister(minid1);
  epregister(minid2);
  epregister(cfthres);
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
  epregister(outfmt);

  daemonArgs(actionDaemon);

  epregisterAction("paired",actionPairend);
  getParser().actions.add("clusterdb",actionClusterDB);
  getParser().actions.add("cluster",actionCluster);
  getParser().actions.add("clustercompress",actionClusterCompress);
  getParser().actions.add("compress",actionCompress);
  getParser().actions.add("convert",actionConvert);
  getParser().actions.add("converttax",actionConvertTax);
  getParser().actions.add("loadbinary",actionLoadBinary);
  getParser().actions.add("loadtaxbinary",actionLoadTaxBinary);
  getParser().actions.add("decompress",actionDecompress);
//  getParser().actions.add("paired",actionPairend);
  getParser().actions.add("psearch",actionProtSearch);
//  getParser().actions.add("daemon",actionDaemon);

  getParser().actions.add("otucounts",actionOTUCounts);
  getParser().actions.add("otutable",actionOTUTable);
  getParser().actions.add("asvotutable",actionASVOTUTable);
  getParser().actions.add("chimera",actionChimera);



//  epregisterClassMethod4(eoption<outfmt_fd>,operator=,int,(const estr& val),"=");


  epregisterFunc(help);
  
//  epregister(galign);
  bool benchmark=false;
//  epregister(benchmark);

  epregister(nthreads);

  epregister2(db.maxhits,"maxhits");
  epregister2(db.tophits,"tophits");
  epregister2(db.topotus,"topotus");
  int step=0;

//  epregister(step);

  estr ignfile;
//  epregister(ignfile);

//  epregister(kmer);
  eparseArgs();
  if (sens){
    db.otulim=10000;
    db.tophits=10000;
    db.topotus=60;
  }
  cerr << "# mapseq v"<< MAPSEQ_PACKAGE_VERSION << " (" << __DATE__ << ")" << endl;

  if(getParser().args.size()<2){
    cout << "syntax: mapseq <query> [db] [tax] [tax2] ..." << endl;
    cout << "mapseq -h for more help" << endl;
    return(0);
  }
  cerr << "# threads: "<< nthreads << endl;

/*
  estrhash ignseqs;
  if (ignfile.len()>0){
    f.open(ignfile,"r");
    while (!f.eof()){
      f.readln(line);
      ignseqs.add(line,estr());
    }
  }
*/

  loadSequences(db);
  cerr << "# loaded " << db.seqs.size() << " sequences" << endl;
  ldieif(db.seqs.size()==0,"empty database");

  loadTaxonomy(db);
  initDB(db);

//  eintarray kmerpos;
//  ebasicarray<uint32_t> idcount;

//  const int MAXSEQS=1000000;
//  uint64_t bitmask[MAXSEQS/64+1];

  db.printSearchHeader();
  ethreads t;
  t.setThreads(nthreads);

  int procret=0;
  if (!fastq)
    procret=db.processQueryFASTA(getParser().args[1],taskSearch,t);
  else
    procret=db.processQueryFASTQ(getParser().args[1],taskSearch,t);

  ldieif(procret!=0,"process query");

  exit(0);

  return(0);
}
