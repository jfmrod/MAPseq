#include "eseqali.h"

void invertstr(estr& str)
{
  int i;
  char tmp;
  for (i=0; i<str.len()/2; ++i){
    tmp=str[i];
    str[i]=str[str.len()-i-1];
    str[str.len()-i-1]=tmp;
  }
}

char nuccomplement(char c)
{
  switch(c){
    case 't': return('a');
    case 'a': return('t');
    case 'g': return('c');
    case 'c': return('g');
   default:
    return(c);
  }
}

estr seqcomplement(const estr& str)
{
  int i;
  estr tmp;
  for (i=str.size()-1; i>=0; --i)
    tmp+=nuccomplement(str[i]);
  return(tmp);
}


float match=1;
float penalty=1;
float gap=5;
float gapext=1;

inline double S(const char &a,const char &b)
{
  if (a==b) return(match);
  return(-penalty);
}


float salign(const estr& a,const estr& b)
{
  int i,j;
  double F[a.len()][b.len()];
  
  F[0][0]=maxv(-gap,-gap,S(a[0],b[0]));
  for (i=1; i<a.len(); ++i)
    F[i][0]=maxv(F[i-1][0]-gap,-gap-i*gap,-i*gap);
  for (j=1; j<b.len(); ++j)
    F[0][j]=maxv(-gap-j*gap,F[0][j-1]-gap,-j*gap);
  for (i=1; i<a.len(); ++i){
    for (j=1; j<b.len(); ++j)
      F[i][j]=maxv(F[i-1][j]-gap,F[i][j-1]-gap,F[i-1][j-1]+S(a[i],b[j]));
  }

/*
  for (i=0; i<a.len(); ++i){
    for (j=0; j<b.len(); ++j){
      printf("%5.1lf ",F[i][j]);
    }
    printf("\n");
  }
*/

  estr seq1,seq2;
  seq1.clear(); seq2.clear();

  int alncount=0;
  int idcount=0;

  i=a.len()-1; j=b.len()-1;
  while (i>0 && j>0){
    switch (maxs(F[i-1][j],F[i][j-1],F[i-1][j-1])){
      case 0: seq1 += a[i]; seq2 += '-'; --i; break;
      case 1: seq1 += '-'; seq2 += b[j]; --j; break;
      case 2: seq1 += a[i]; seq2 += b[j]; ++alncount; if (a[i]==b[j]) ++idcount; --j; --i; break;
    }
  }
  while (j>0 && i>=0){
    switch (maxs(-gap-(j-1)*gap,F[i][j-1],-(j-1)*gap)){
      case 0: seq1 += a[i]; seq2 += '-'; --i; break;
      case 1: seq1 += '-'; seq2 += b[j]; --j; break;
      case 2: seq1 += a[i]; seq2 += b[j]; ++alncount; if (a[i]==b[j]) ++idcount; --j; --i; break;
    }
  }
  while (j>=0 && i>0){
    switch (maxs(F[i-1][j],-gap-(i-1)*gap,-(i-1)*gap)){
      case 0: seq1 += a[i]; seq2 += '-'; --i; break;
      case 1: seq1 += '-'; seq2 += b[j]; --j; break;
      case 2: seq1 += a[i]; seq2 += b[j]; ++alncount; if (a[i]==b[j]) ++idcount; --j; --i; break;
    }
  }

  for (;i==-1 && j>=0;--j){
    seq1+='-'; seq2+=b[j];
  }
  for (;j==-1 && i>=0;--i){
    seq2+='-'; seq1+=a[i];
  }
  if (i==0 && j==0){
    seq1+=a[i]; seq2+=b[j];
    ++alncount; if (a[i]==b[j]) ++idcount;
  }


/*
  for (i=al_a.len()-1; i>=0; --i)
    cout << al_a[i];
  cout << endl;

  for (i=al_b.len()-1; i>=0; --i)
    cout << al_b[i];
  cout << endl;

  cout << " Identity: " << (double)idcount/(double)alncount<< " (" << idcount << "/"<< alncount << ")"<<endl;
*/
  return((float)idcount/alncount);
}

void seq_print_ref(const estr& seq,const estr& ref)
{
  ldieif(seq.len()!=ref.len(),"sequences not aligned");
  int i;
  estr tmpseq;
  for (i=0; i<ref.len(); ++i)
    if (ref[i]!='-')
      tmpseq+=seq[i];
  cout << tmpseq << endl;
}

inline char nuc2chr(unsigned char t)
{
  switch(t){
    case 0x0u: return('a');
    case 0x1u: return('t');
    case 0x2u: return('g');
    case 0x3u: return('c');
  }
  return('x');
}

inline int maxv3(double v0,double v1,double v2)
{
  if (v0>v1) { if (v0>v2) return(0); else return(2); }
  if (v2>v1) return(2);
  return(1);
}

#define MAX3(a,b,c) (a>b?(a>c?a:c):(b>c?b:c))

inline void updateF(int i,int j,int h,double *mF,float ms,const ealignscore& as)
{
  mF[i*h*3+j*3+1]=MAX(mF[(i-1)*h*3+j*3]-as.gapopen,mF[(i-1)*h*3+j*3+1])-as.gapext;
  mF[i*h*3+j*3+2]=MAX(mF[i*h*3+(j-1)*3]-as.gapopen,mF[i*h*3+(j-1)*3+2])-as.gapext;
  mF[i*h*3+j*3]=MAX3(mF[(i-1)*h*3+(j-1)*3],mF[(i-1)*h*3+(j-1)*3+1],mF[(i-1)*h*3+(j-1)*3+2])+ms;
}

inline void updateF(int i,int j,int h,double *mF,double *mFi,double *mFj,float ms,const ealignscore& as)
{
  mFi[i*h+j]=MAX(mF[(i-1)*h+j]-as.gapopen-as.gapext,mFi[(i-1)*h+j]-as.gapext);
  mFj[i*h+j]=MAX(mF[i*h+j-1]-as.gapopen-as.gapext,mFj[i*h+j-1]-as.gapext);
  mF[i*h+j]=MAX3(mF[(i-1)*h+j-1],mFi[(i-1)*h+(j-1)],mFj[(i-1)*h+(j-1)])+ms;
}

// No Edge Gap
inline double updateFNEG(int i,int j,int w,int h,double *mF,float ms,const ealignscore& as)
{
  if (j==h-1)
    mF[i*h*3+j*3+1]=MAX(mF[(i-1)*h*3+j*3],mF[(i-1)*h*3+j*3+1]);
  else
    mF[i*h*3+j*3+1]=MAX(mF[(i-1)*h*3+j*3]-as.gapopen,mF[(i-1)*h*3+j*3+1])-as.gapext;
  if (i==w-1)
    mF[i*h*3+j*3+2]=MAX(mF[i*h*3+(j-1)*3],mF[i*h*3+(j-1)*3+2]);
  else
    mF[i*h*3+j*3+2]=MAX(mF[i*h*3+(j-1)*3]-as.gapopen,mF[i*h*3+(j-1)*3+2])-as.gapext;
  mF[i*h*3+j*3]=MAX3(mF[(i-1)*h*3+(j-1)*3],mF[(i-1)*h*3+(j-1)*3+1],mF[(i-1)*h*3+(j-1)*3+2])+ms;
  return(MAX3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2]));
}
inline double updateFNEG2(int i,int j,int w,int h,double *mF,float ms,const ealignscore& as)
{
  mF[i*h*3+j*3+1]=MAX(mF[(i-1)*h*3+j*3]-as.gapopen,mF[(i-1)*h*3+j*3+1])-as.gapext;
  mF[i*h*3+j*3+2]=MAX(mF[i*h*3+(j-1)*3]-as.gapopen,mF[i*h*3+(j-1)*3+2])-as.gapext;
  mF[i*h*3+j*3]=MAX3(mF[(i-1)*h*3+(j-1)*3],mF[(i-1)*h*3+(j-1)*3+1],mF[(i-1)*h*3+(j-1)*3+2])+ms;
  return(MAX3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2]));
}
inline void updateFNEGh(int i,int j,int w,int h,double *mF,float ms,const ealignscore& as)
{
  mF[i*h*3+j*3+1]=MAX(mF[(i-1)*h*3+j*3],mF[(i-1)*h*3+j*3+1]);
  mF[i*h*3+j*3+2]=MAX(mF[i*h*3+(j-1)*3]-as.gapopen,mF[i*h*3+(j-1)*3+2])-as.gapext;
  mF[i*h*3+j*3]=MAX3(mF[(i-1)*h*3+(j-1)*3],mF[(i-1)*h*3+(j-1)*3+1],mF[(i-1)*h*3+(j-1)*3+2])+ms;
}
inline void updateFNEGw(int i,int j,int w,int h,double *mF,float ms,const ealignscore& as)
{
  mF[i*h*3+j*3+1]=MAX(mF[(i-1)*h*3+j*3]-as.gapopen,mF[(i-1)*h*3+j*3+1])-as.gapext;
  mF[i*h*3+j*3+2]=MAX(mF[i*h*3+(j-1)*3],mF[i*h*3+(j-1)*3+2]);
  mF[i*h*3+j*3]=MAX3(mF[(i-1)*h*3+(j-1)*3],mF[(i-1)*h*3+(j-1)*3+1],mF[(i-1)*h*3+(j-1)*3+2])+ms;
}
inline void updateFNEGwh(int i,int j,int w,int h,double *mF,float ms,const ealignscore& as)
{
  mF[i*h*3+j*3+1]=MAX(mF[(i-1)*h*3+j*3],mF[(i-1)*h*3+j*3+1]);
  mF[i*h*3+j*3+2]=MAX(mF[i*h*3+(j-1)*3],mF[i*h*3+(j-1)*3+2]);
  mF[i*h*3+j*3]=MAX3(mF[(i-1)*h*3+(j-1)*3],mF[(i-1)*h*3+(j-1)*3+1],mF[(i-1)*h*3+(j-1)*3+2])+ms;
}

// No Edge Gap
inline void updateFNEG(int i,int j,int w,int h,double *mF,double *mFi,double *mFj,float ms,const ealignscore& as)
{
  if (j==h-1)
    mFi[i*h+j]=MAX(mF[(i-1)*h+j],mFi[(i-1)*h+j]);
  else
    mFi[i*h+j]=MAX(mF[(i-1)*h+j]-as.gapopen-as.gapext,mFi[(i-1)*h+j]-as.gapext);
  if (i==w-1)
    mFj[i*h+j]=MAX(mF[i*h+j-1],mFj[i*h+j-1]);
  else
    mFj[i*h+j]=MAX(mF[i*h+j-1]-as.gapopen-as.gapext,mFj[i*h+j-1]-as.gapext);
  mF[i*h+j]=MAX3(mF[(i-1)*h+j-1],mFi[(i-1)*h+(j-1)],mFj[(i-1)*h+(j-1)])+ms;
}

float seqcalign_global_noleftedgegap(const eseq& a,int pa,int ea,const eseq& b,int pb,int eb,estr& as1,estr& as2,const ealignscore& as)
{
  uint64_t *psa=reinterpret_cast<uint64_t*>(a.seq._str),*psb=reinterpret_cast<uint64_t*>(b.seq._str);
  float mc_score[4]={as.match,-as.mismatch,-as.mismatch,-as.mismatch};
  float score=0;
  int aligned=0,gaps=0,matches=0,mismatches=0;
  int maxsize=MAX(ea-pa,eb-pb);
  if (maxsize==0) return(0.0);

  int w=(ea-pa)+1;
  int h=(eb-pb)+1;

  int i,j;
  double *mF=new double[w*h];
  double *mFi=new double[w*h];
  double *mFj=new double[w*h];
//  char *gapOpen=new char[w*h];
//  char hgapOpen;
  

  // no gap penalty for beginning gaps
  for (i=0; i<w; ++i){
    mFi[i*h]=0.0l;
    mFj[i*h]=-1.0e100l;
    mF[i*h]=-1.0e100l;
  }
  for (j=0; j<h; ++j){
    mFi[j]=-1.0e100l;
    mFj[j]=0.0l;
    mF[j]=-1.0e100l;
  }
  mF[0]=0.0l;
  mFi[0]=-as.gapopen;
  mFj[0]=-as.gapopen;

  uint64_t ca,cb;
  for (i=1; i<w; ++i){
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    for (j=1; j<h; ++j) {
      cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      updateF(i,j,h,mF,mFi,mFj,mc_score[ca^cb],as);
    }
  }

  aligned=0;
  matches=0;

  as1.clear(); as2.clear();
  as1.reserve(a.seqlen+b.seqlen);
  as2.reserve(a.seqlen+b.seqlen);

  i=w-1; j=h-1;
  while (i>0 && j>0){
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    switch (maxv3(mF[i*h+j],mFi[i*h+j],mFj[i*h+j])){
      case 0: as1 += nuc2chr(ca); as2 += nuc2chr(cb); ++aligned; if (ca==cb) ++matches; --j; --i; break;
      case 1: as1 += nuc2chr(ca); as2 += '-'; --i; break;
      case 2: as1 += '-'; as2 += nuc2chr(cb); --j; break;
     default:
      ldie("unknown code");
    }
  }
  for (;i==0 && j>0;--j){
    cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    as1+='-'; as2+=nuc2chr(cb);
  }
  for (;j==0 && i>0;--i){
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    as2+='-'; as1+=nuc2chr(ca);
  }
  invertstr(as1);
  invertstr(as2);
//  cout << as1 << endl;
//  cout << as2 << endl;

  delete[] mF;
  delete[] mFi;
  delete[] mFj;

  mismatches=aligned-matches;
  gaps=maxsize-aligned;
  return((float)matches/aligned);
}

float seqcalign_global_norightedgegap(const eseq& a,int pa,int ea,const eseq& b,int pb,int eb,estr& as1,estr& as2,const ealignscore& as)
{
  uint64_t *psa=reinterpret_cast<uint64_t*>(a.seq._str),*psb=reinterpret_cast<uint64_t*>(b.seq._str);
  float mc_score[4]={as.match,-as.mismatch,-as.mismatch,-as.mismatch};
  float score=0;
  int aligned=0,gaps=0,matches=0,mismatches=0;
  int maxsize=MAX(ea-pa,eb-pb);
  if (maxsize==0) return(0.0);

  int w=(ea-pa)+1;
  int h=(eb-pb)+1;

  int i,j;
  double *mF=new double[w*h];
  double *mFi=new double[w*h];
  double *mFj=new double[w*h];
//  char *gapOpen=new char[w*h];
//  char hgapOpen;
  

  // no gap penalty for beginning gaps
  for (i=0; i<w; ++i){
    mFi[i*h]=-as.gapopen-i*as.gapext;
    mFj[i*h]=-1.0e100l;
    mF[i*h]=-1.0e100l;
  }
  for (j=0; j<h; ++j){
    mFi[j]=-1.0e100l;
    mFj[j]=-as.gapopen-j*as.gapext;
    mF[j]=-1.0e100l;
  }
  mF[0]=0.0l;
  mFi[0]=-as.gapopen;
  mFj[0]=-as.gapopen;

  uint64_t ca,cb;
  for (i=1; i<w; ++i){
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    for (j=1; j<h; ++j) {
      cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      updateFNEG(i,j,w,h,mF,mFi,mFj,mc_score[ca^cb],as);
    }
  }

  aligned=0;
  matches=0;

  as1.clear(); as2.clear();
  as1.reserve(a.seqlen+b.seqlen);
  as2.reserve(a.seqlen+b.seqlen);

  i=w-1; j=h-1;
  while (i>0 && j>0){
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    switch (maxv3(mF[i*h+j],mFi[i*h+j],mFj[i*h+j])){
      case 0: as1 += nuc2chr(ca); as2 += nuc2chr(cb); ++aligned; if (ca==cb) ++matches; --j; --i; break;
      case 1: as1 += nuc2chr(ca); as2 += '-'; --i; break;
      case 2: as1 += '-'; as2 += nuc2chr(cb); --j; break;
     default:
      ldie("unknown code");
    }
  }
  for (;i==0 && j>0;--j){
    cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    as1+='-'; as2+=nuc2chr(cb);
  }
  for (;j==0 && i>0;--i){
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    as2+='-'; as1+=nuc2chr(ca);
  }
  invertstr(as1);
  invertstr(as2);
//  cout << as1 << endl;
//  cout << as2 << endl;

  delete[] mF;
  delete[] mFi;
  delete[] mFj;

  mismatches=aligned-matches;
  gaps=maxsize-aligned;
  return((float)matches/aligned);
}

float seqcalign_global_noedgegap(const eseq& a,int pa,int ea,const eseq& b,int pb,int eb,estr& as1,estr& as2,const ealignscore& as)
{
  uint64_t *psa=reinterpret_cast<uint64_t*>(a.seq._str),*psb=reinterpret_cast<uint64_t*>(b.seq._str);
  float mc_score[4]={as.match,-as.mismatch,-as.mismatch,-as.mismatch};
  float score=0;
  int aligned=0,gaps=0,matches=0,mismatches=0;
  int maxsize=MAX(ea-pa,eb-pb);
  if (maxsize==0) return(0.0);

  int w=(ea-pa)+1;
  int h=(eb-pb)+1;

  int i,j;
  double *mF=new double[w*h];
  double *mFi=new double[w*h];
  double *mFj=new double[w*h];
//  char *gapOpen=new char[w*h];
//  char hgapOpen;
  

  // no gap penalty for beginning gaps
  for (i=0; i<w; ++i){
    mFi[i*h]=0.0l;
    mFj[i*h]=-1.0e100l;
    mF[i*h]=-1.0e100l;
  }
  for (j=0; j<h; ++j){
    mFi[j]=-1.0e100l;
    mFj[j]=0.0l;
    mF[j]=-1.0e100l;
  }
  mF[0]=0.0l;
  mFi[0]=-as.gapopen;
  mFj[0]=-as.gapopen;

  uint64_t ca,cb;
  for (i=1; i<w; ++i){
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    for (j=1; j<h; ++j) {
      cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      updateFNEG(i,j,w,h,mF,mFi,mFj,mc_score[ca^cb],as);
    }
  }

  aligned=0;
  matches=0;

  as1.clear(); as2.clear();
  as1.reserve(a.seqlen+b.seqlen);
  as2.reserve(a.seqlen+b.seqlen);

  i=w-1; j=h-1;
  while (i>0 && j>0){
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    switch (maxv3(mF[i*h+j],mFi[i*h+j],mFj[i*h+j])){
      case 0: as1 += nuc2chr(ca); as2 += nuc2chr(cb); ++aligned; if (ca==cb) ++matches; --j; --i; break;
      case 1: as1 += nuc2chr(ca); as2 += '-'; --i; break;
      case 2: as1 += '-'; as2 += nuc2chr(cb); --j; break;
     default:
      ldie("unknown code");
    }
  }
  for (;i==0 && j>0;--j){
    cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    as1+='-'; as2+=nuc2chr(cb);
  }
  for (;j==0 && i>0;--i){
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    as2+='-'; as1+=nuc2chr(ca);
  }
  invertstr(as1);
  invertstr(as2);
//  cout << as1 << endl;
//  cout << as2 << endl;

  delete[] mF;
  delete[] mFi;
  delete[] mFj;

  mismatches=aligned-matches;
  gaps=maxsize-aligned;
  return((float)matches/aligned);
}

float seqcalign_global(const eseq& a,int pa,int ea,const eseq& b,int pb,int eb,estr& as1,estr& as2,const ealignscore& as)
{
  uint64_t *psa=reinterpret_cast<uint64_t*>(a.seq._str),*psb=reinterpret_cast<uint64_t*>(b.seq._str);
  float mc_score[4]={as.match,-as.mismatch,-as.mismatch,-as.mismatch};
  float score=0;
  int aligned=0,gaps=0,matches=0,mismatches=0;
  int maxsize=MAX(ea-pa,eb-pb);
  if (maxsize==0) return(0.0);

  int w=(ea-pa)+1;
  int h=(eb-pb)+1;

  int i,j;
  double *mF=new double[w*h];
  double *mFi=new double[w*h];
  double *mFj=new double[w*h];
//  char *gapOpen=new char[w*h];
//  char hgapOpen;
  

  // no gap penalty for beginning gaps
  for (i=0; i<w; ++i){
    mFi[i*h]=-as.gapopen-i*as.gapext;
    mFj[i*h]=-1.0e100l;
    mF[i*h]=-1.0e100l;
  }
  for (j=0; j<h; ++j){
    mFi[j]=-1.0e100l;
    mFj[j]=-as.gapopen-j*as.gapext;
    mF[j]=-1.0e100l;
  }
  mF[0]=0.0l;
  mFi[0]=-as.gapopen;
  mFj[0]=-as.gapopen;

  uint64_t ca,cb;
  for (i=1; i<w; ++i){
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    for (j=1; j<h; ++j) {
      cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      updateF(i,j,h,mF,mFi,mFj,mc_score[ca^cb],as);
    }
  }

  aligned=0;
  matches=0;

  as1.clear(); as2.clear();
  as1.reserve(a.seqlen+b.seqlen);
  as2.reserve(a.seqlen+b.seqlen);

  i=w-1; j=h-1;
  double maxS=MAX3(mF[i*h+j],mFi[i*h+j],mFj[i*h+j]);
  while (i>0 && j>0){
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    switch (maxv3(mF[i*h+j],mFi[i*h+j],mFj[i*h+j])){
      case 0: as1 += nuc2chr(ca); as2 += nuc2chr(cb); ++aligned; if (ca==cb) ++matches; --j; --i; break;
      case 1: as1 += nuc2chr(ca); as2 += '-'; --i; break;
      case 2: as1 += '-'; as2 += nuc2chr(cb); --j; break;
     default:
      ldie("unknown code");
    }
  }
  for (;i==0 && j>0;--j){
    cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    as1+='-'; as2+=nuc2chr(cb);
  }
  for (;j==0 && i>0;--i){
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    as2+='-'; as1+=nuc2chr(ca);
  }
  invertstr(as1);
  invertstr(as2);
//  cout << as1 << endl;
//  cout << as2 << endl;

  delete[] mF;
  delete[] mFi;
  delete[] mFj;

  mismatches=aligned-matches;
  gaps=maxsize-aligned;
//  return((float)matches/aligned);
  return(maxS);
}

float seqcalign(const eseq& a,int pa,int ea,const eseq& b,int pb,int eb,estr& as1,estr& as2)
{
  uint64_t *psa=reinterpret_cast<uint64_t*>(a.seq._str),*psb=reinterpret_cast<uint64_t*>(b.seq._str);
  float mc_score[4]={match,-penalty,-penalty,-penalty};
  float score=0;
  int aligned=0,gaps=0,matches=0,mismatches=0;
  int maxsize=MAX(ea-pa,eb-pb);
  if (maxsize==0) return(0.0);

  int w=(ea-pa)+1;
  int h=(eb-pb)+1;

  int i,j;
  double *mF=new double[w*h];
  char *gapOpen=new char[w*h];
//  char hgapOpen;
  
  mF[0]=0.0l;
  gapOpen[0]=0;

  // no gap penalty for beginning gaps
  for (i=1; i<w; ++i)
    { mF[i*h]=0.0l; gapOpen[i*h]=0; }
  for (j=1; j<h; ++j)
    { mF[j]=0.0l; gapOpen[j]=0; }

  uint64_t ca,cb;
  double max_S=-100.0;
  int max_i=-1;
  int max_j=-1;
  for (i=1; i<w-1; ++i){
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    for (j=1; j<h-1; ++j) {
      cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      maxvg(mF[i*h+j],gapOpen[i*h+j],mF[(i-1)*h+j]-(gapOpen[(i-1)*h+j]==1?gapext:gap),mF[i*h+j-1]-(gapOpen[i*h+j-1]==2?gapext:gap),mF[(i-1)*h+j-1]+mc_score[ca^cb]);
      if (mF[i*h+j]>max_S){ max_S=mF[i*h+j]; max_i=i; max_j=j; }
    }
    // no gap penalty at end of seq
    cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    maxvg(mF[i*h+j],gapOpen[i*h+j],mF[(i-1)*h+j],mF[i*h+j-1]-(gapOpen[i*h+j-1]==2?gapext:gap),mF[(i-1)*h+j-1]+mc_score[ca^cb]);
    if (mF[i*h+j]>max_S){ max_S=mF[i*h+j]; max_i=i; max_j=j; }
  }

  // no gap penalties at end of seq
  ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
  for (j=1; j<h-1; ++j) {
    cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    maxvg(mF[i*h+j],gapOpen[i*h+j],mF[(i-1)*h+j]-(gapOpen[(i-1)*h+j]==1?gapext:gap),mF[i*h+j-1],mF[(i-1)*h+j-1]+mc_score[ca^cb]);
    if (mF[i*h+j]>max_S){ max_S=mF[i*h+j]; max_i=i; max_j=j; }
  }
  cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
  maxvg(mF[i*h+j],gapOpen[i*h+j],mF[(i-1)*h+j],mF[i*h+j-1],mF[(i-1)*h+j-1]+mc_score[ca^cb]);
  if (mF[i*h+j]>max_S){ max_S=mF[i*h+j]; max_i=i; max_j=j; }

  ldieif(max_i==-1 || max_j==-1,"max_i or max_j are negative");
//  cout << "max_i: " << max_i << " max_j: " << max_j << " max_S: " << max_S << endl;

  aligned=0;
  matches=0;

  as1.clear(); as2.clear();
  as1.reserve(a.seqlen+b.seqlen);
  as2.reserve(a.seqlen+b.seqlen);

  double gapseqa,gapseqb;
//  i=a.len(); j=b.len(); // global alignment starting position
//  i=max_i; j=max_j; // local alignment starting from maximum score
  i=w-1; j=h-1; // local alignment starting from maximum score
  while (i>0 && j>0){
    gapseqa=mF[(i-1)*h+j];
    if (i<w-1)
      gapseqb-=(gapOpen[i*h+j-1]==1?gapext:gap);
    gapseqb=mF[i*h+j-1];
    if (j<h-1)
      gapseqb-=(gapOpen[(i-1)*h+j]==2?gapext:gap);
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    switch (maxs(gapseqa,gapseqb,mF[(i-1)*h+j-1]+mc_score[ca^cb])){
      case 0: as1 += nuc2chr(ca); as2 += '-'; --i; break;
      case 1: as1 += '-'; as2 += nuc2chr(cb); --j; break;
      case 2: as1 += nuc2chr(ca); as2 += nuc2chr(cb); ++aligned; if (ca==cb) ++matches; --j; --i; break;
     default:
      ldie("unknown code");
    }
  }
  for (;i==0 && j>0;--j){
    cb=((j-1+pb)%32u==0u?psb[(j-1+pb)/32u]:psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    as1+='-'; as2+=nuc2chr(cb);
  }
  for (;j==0 && i>0;--i){
    ca=((i-1+pa)%32u==0u?psa[(i-1+pa)/32u]:psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    as2+='-'; as1+=nuc2chr(ca);
  }
  invertstr(as1);
  invertstr(as2);
//  cout << as1 << endl;
//  cout << as2 << endl;

  delete[] mF;
  delete[] gapOpen;

  mismatches=aligned-matches;
  gaps=maxsize-aligned;
  return((float)matches/aligned);
}


void print_seqali(const estr& seq1,const estr& seq2)
{
  ldieif(seq1.len()!=seq2.len(),"seq mismatches"); 
  if (seq2.len()==0) return;

  int i=0,e=seq1.len()-1;
  int si1=0,si2=0;

  if (seq2[0]=='-')
    for (; i<seq2.len() && seq2[i]=='-'; ++i,++si1);
  else if (seq1[0]=='-')
    for (; i<seq1.len() && seq1[i]=='-'; ++i,++si2);
  
  if (seq1[seq1.len()-1]=='-')
    for (; e>=i && seq1[e]=='-'; --e);
  else if (seq2[seq2.len()-1]=='-')
    for (; e>=i && seq2[e]=='-'; --e);

  int gapscore=0;
  int gapopen=0;
  int score=0,aligned=0,gaps=0,matches=0,mismatches=0;
  for (int ti=i; ti<=e; ++ti){
    if (seq1[ti]=='-' || seq2[ti]=='-'){
      if (seq1[ti]=='-'){
        if (gapopen==0 || gapopen==2) { gapopen=1; gapscore+=gap+gapext; }
        else gapscore+=gapext;
      }else{
        if (gapopen==0 || gapopen==1) { gapopen=2; gapscore+=gap+gapext; }
        else gapscore+=gapext;
      }
      ++gaps;
    }else if (seq1[ti]==seq2[ti])
      ++matches;
    else
      ++mismatches;
  }
  aligned=matches+mismatches;
  score=matches*match - mismatches*penalty - gapscore;
  cout << "score: " << score << " identity: " << (float)matches/aligned << " aligned: " << aligned << " matches: " << matches << " mismatches: " << mismatches << " gaps: " << gaps << " aligned_seq_len: " << seq1.len() << endl; 

  estr tmpseq1;


  for (; i<=e; i+=100){
    printf("%4.i ",si2);
    for (int ti=0; ti<100 && i+ti<=e; ++ti)
      if (seq2[ti]!='-') ++si2;
    for (int k=0; k<100 && i+k<=e; k+=10)
      cout << seq2.substr(i+k,10) << " ";
    cout << endl;
    tmpseq1.clear();
    for (int j=i; j<i+100 && j<=e; ++j){
      if (seq2[j] == '-' || seq1[j] =='-') tmpseq1+=' ';
      else if (seq2[j] == seq1[j]) tmpseq1+='|';
      else tmpseq1+='.';
    }
    printf("     ");
    for (int k=0; k<100 && k<=e; k+=10)
      cout << tmpseq1.substr(k,10) << " ";
    cout << endl;
    printf("%4.i ",si1);
    for (int ti=0; ti<100 && i+ti<=e; ++ti)
      if (seq1[ti]!='-') ++si1;
    for (int k=0; k<100 && i+k<=e; k+=10)
      cout << seq1.substr(i+k,10) << " ";
    cout << endl;
    cout << endl;
  }
}

void seqcalign_global(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const ealignscore& as)
{
  uint64_t *psa=reinterpret_cast<uint64_t*>(a.seq._str),*psb=reinterpret_cast<uint64_t*>(b.seq._str);
  float mc_score[4]={as.match,-as.mismatch,-as.mismatch,-as.mismatch};
  float score=0;
  int aligned=0,gaps=0,matches=0,mismatches=0;
  long maxsize=MAX(ea-pa,eb-pb);
  long minsize=MIN(ea-pa,eb-pb);
  if (maxsize==0) return;
  if (minsize==0) {
    if (adata.e1==-1 || adata.e1<pa) adata.e1=pa;
    if (adata.e2==-1 || adata.e2<pb) adata.e2=pb;
    adata._score+=-as.gapopen-maxsize*as.gapext;
    adata.gaps+=maxsize;
    if (ea==pa)
      adata.profile.add(AT_DEL,maxsize);
    else
      adata.profile.add(AT_INS,maxsize);
    return;
  }

  int w=(ea-pa)+1;
  int h=(eb-pb)+1;

  int i,j;
  ws.reserve(w*h*3);
  double *&mF(ws.mF);

  // no gap penalty for beginning gaps
  for (i=0; i<w; ++i){
    mF[i*h*3+1]=-as.gapopen-i*as.gapext;
    mF[i*h*3+2]=-1.0e100l;
    mF[i*h*3]=-1.0e100l;
  }
  for (j=0; j<h; ++j){
    mF[j*3+1]=-1.0e100l;
    mF[j*3+2]=-as.gapopen-j*as.gapext;
    mF[j*3]=-1.0e100l;
  }
  mF[0]=0.0l;
  mF[1]=-as.gapopen;
  mF[2]=-as.gapopen;

  uint64_t ca,cb;
  for (i=1; i<w; ++i){
    ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    for (j=1; j<h; ++j) {
      cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      updateF(i,j,h,mF,mc_score[ca^cb],as);
    }
  }

  if (adata.s1==-1 || adata.s1>pa) adata.s1=pa;
  if (adata.s2==-1 || adata.s2>pb) adata.s2=pb;
  if (adata.e1==-1 || adata.e1<ea) adata.e1=ea;
  if (adata.e2==-1 || adata.e2<eb) adata.e2=eb;
  adata._score+=MAX3(mF[(w-1)*h*3+(h-1)*3],mF[(w-1)*h*3+(h-1)*3+1],mF[(w-1)*h*3+(h-1)*3+2]);

  aligned=0;
  matches=0;

  long ti=ea-1;
  i=w-1; j=h-1;
  while (i>0 && j>0){
    ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    switch (maxv3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2])){
      case 0: ++aligned; if (ca==cb) { ++matches; adata.profile.add(AT_MATCH); } else adata.profile.add(AT_MISS); --j; --i; break;
      case 1: --i; ++gaps; adata.profile.add(AT_INS); break;
      case 2: --j; ++gaps; adata.profile.add(AT_DEL); break;
     default:
      ldie("unknown code");
    }
  }
  for (;i>0; --i){
    ++gaps;
    adata.profile.add(AT_INS);
  }
  for (;j>0; --j){
    ++gaps;
    adata.profile.add(AT_DEL);
  }

  adata.matches+=matches;
  adata.mismatches+=aligned-matches;
  adata.gaps+=gaps;
}



void seqcalign_global_noleftedgegap(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const ealignscore& as)
{
  uint64_t *psa=reinterpret_cast<uint64_t*>(a.seq._str),*psb=reinterpret_cast<uint64_t*>(b.seq._str);
  float mc_score[4]={as.match,-as.mismatch,-as.mismatch,-as.mismatch};
  float score=0.0;
  int aligned=0,gaps=0,matches=0,mismatches=0;
  long maxsize=MAX(ea-pa,eb-pb);
  long minsize=MIN(ea-pa,eb-pb);
  if (minsize==0) { // we set the leftmost part to the rightmost part of the sequence (in case we are aligning 0bp against Xbp long sequence)
    if (adata.s1==-1 || adata.s1>ea) adata.s1=ea;
    if (adata.s2==-1 || adata.s2>eb) adata.s2=eb;
//    if (adata.e1==-1 || adata.e1<ea) adata.e1=ea;
//    if (adata.e2==-1 || adata.e2<eb) adata.e2=eb;
    return;
  }

  int w=(ea-pa)+1;
  int h=(eb-pb)+1;

  int i,j;
  ws.reserve(w*h*3);
  double *&mF(ws.mF);

  // no gap penalty for beginning gaps
  for (i=0; i<w; ++i){
    mF[i*h*3+1]=0.0l;
    mF[i*h*3+2]=-1.0e100l;
    mF[i*h*3]=-1.0e100l;
  }
  for (j=0; j<h; ++j){
    mF[j*3+1]=-1.0e100l;
    mF[j*3+2]=0.0l;
    mF[j*3]=-1.0e100l;
  }
  mF[0]=0.0l;
  mF[1]=0.0l;
  mF[2]=0.0l;

  uint64_t ca,cb;
  for (i=1; i<w; ++i){
    ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    for (j=1; j<h; ++j) {
      cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      updateF(i,j,h,mF,mc_score[ca^cb],as);
    }
  }

  aligned=0;
  matches=0;


/*
  estr as1,as2;
  as1.reserve(h+w-2);
  as2.reserve(h+w-2);
*/

  long ti=ea-1;
  i=w-1; j=h-1;
  adata._score+=MAX3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2]);
  if (adata.aln){
    while (i>0 && j>0){
      ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
      cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      switch (maxv3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2])){
  /*
        case 0: as1 += nuc2chr(ca); as2 += nuc2chr(cb); countgaps=1; ++aligned; if (ca==cb) { adata.aln[ti]=0x01u; ++matches; } else adata.aln[ti]=0x02u; --j; --i; --ti; break;
        case 1: as1 += nuc2chr(ca); as2 += '-'; adata.aln[ti]=0x00u; --i; --ti; gaps+=countgaps; break;
        case 2: as1 += '-'; as2 += nuc2chr(cb); --j; gaps+=countgaps; break;
  */
        case 0: ++aligned; if (ca==cb) { adata.aln[ti]=0x01u; ++matches; } else adata.aln[ti]=0x02u; --j; --i; --ti; break;
        case 1: adata.aln[ti]=0x00u; --i; --ti; ++gaps; break;
        case 2: --j; ++gaps; break;
       default:
        ldie("unknown code");
      }
    }
    for (;j==0 && i>0;--i,--ti){
      adata.aln[ti]=0x00u;
  //    as2+='-'; as1+=nuc2chr(ca);
    }
  }else{
    while (i>0 && j>0){
      ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
      cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      switch (maxv3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2])){
        case 0: ++aligned; if (ca==cb) { ++matches; adata.profile.add(AT_MATCH); } else adata.profile.add(AT_MISS); --j; --i; break;
        case 1: --i; ++gaps; adata.profile.add(AT_INS); break;
        case 2: --j; ++gaps; adata.profile.add(AT_DEL); break;
       default:
        ldie("unknown code");
      }
    }
  }


/*
  invertstr(as1);
  invertstr(as2);
  cout << as1 << endl;
  cout << as2 << endl;
*/

  if (adata.s1==-1 || adata.s1>pa+i) adata.s1=pa+i;
  if (adata.s2==-1 || adata.s2>pb+j) adata.s2=pb+j;
//  if (adata.s1==-1 || adata.s1>pa) adata.s1=pa;
//  if (adata.s2==-1 || adata.s2>pb) adata.s2=pb;
//  if (adata.e1==-1 || adata.e1<ea) adata.e1=ea;
//  if (adata.e2==-1 || adata.e2<eb) adata.e2=eb;
  adata.matches+=matches;
  adata.mismatches+=aligned-matches;
  adata.gaps+=gaps;
}

void seqcalign_global_norightedgegap(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const ealignscore& as)
{
  uint64_t *psa=reinterpret_cast<uint64_t*>(a.seq._str),*psb=reinterpret_cast<uint64_t*>(b.seq._str);
  float mc_score[4]={as.match,-as.mismatch,-as.mismatch,-as.mismatch};
  float score=0;
  int aligned=0,gaps=0,matches=0,mismatches=0;
  long maxsize=MAX(ea-pa,eb-pb);
  long minsize=MIN(ea-pa,eb-pb);
  if (minsize==0) {
    if (adata.e1==-1 || adata.e1<pa) adata.e1=pa;
    if (adata.e2==-1 || adata.e2<pb) adata.e2=pb;
    return;
  }

  int w=(ea-pa)+1;
  int h=(eb-pb)+1;

  int i,j;
  ws.reserve(w*h*3);
  double *&mF(ws.mF);

  // no gap penalty for beginning gaps
  for (i=0; i<w; ++i){
    mF[i*h*3+1]=-as.gapopen-i*as.gapext;
    mF[i*h*3+2]=-1.0e100l;
    mF[i*h*3]=-1.0e100l;
  }
  for (j=0; j<h; ++j){
    mF[j*3+1]=-1.0e100l;
    mF[j*3+2]=-as.gapopen-j*as.gapext;
    mF[j*3]=-1.0e100l;
  }
  mF[0]=0.0l;
  mF[1]=-as.gapopen;
  mF[2]=-as.gapopen;

  uint64_t ca,cb;
  for (i=1; i<w; ++i){
    ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    for (j=1; j<h; ++j) {
      cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      updateFNEG(i,j,w,h,mF,mc_score[ca^cb],as);
    }
  }

/*
  for (i=1; i<w-1; ++i){
    ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    for (j=1; j<h-1; ++j) {
      cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      updateFNEG2(i,j,w,h,mF,mc_score[ca^cb],as);
    }
    j=h-1;
    cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    updateFNEGh(i,h-1,w,h,mF,mc_score[ca^cb],as);
  }
  i=w-1;
  ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
  for (j=1; j<h-1; ++j) {
    cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    updateFNEGw(i,j,w,h,mF,mc_score[ca^cb],as);
  }
  j=h-1;
  cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
  updateFNEGwh(i,j,w,h,mF,mc_score[ca^cb],as);
*/

  aligned=0;
  matches=0;

  int tmpsa=0,tmpsb=0;
  long ti=ea-1;
  i=w-1; j=h-1;
  adata._score+=MAX3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2]);
  int countgaps=0;
  if (adata.aln){
    while (i>0 && j>0){
      ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
      cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      switch (maxv3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2])){
        case 0: if (countgaps==0) { tmpsa=i; tmpsb=j; } countgaps=1; ++aligned; if (ca==cb) { adata.aln[ti]=0x01u; ++matches; } else adata.aln[ti]=0x02u; --j; --i; --ti; break;
        case 1: adata.aln[ti]=0x00u; --i; --ti; gaps+=countgaps; break;
        case 2: --j; gaps+=countgaps; break;
       default:
        ldie("unknown code");
      }
    }
    for (;j==0 && i>0;--i,--ti)
      adata.aln[ti]=0x00u;
  }else{
    while (i>0 && j>0){
      ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
      cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      switch (maxv3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2])){
        case 0: if (countgaps==0) { tmpsa=i; tmpsb=j; } countgaps=1; ++aligned; if (ca==cb) { ++matches; adata.profile.add(AT_MATCH); } else adata.profile.add(AT_MISS); --j; --i; break;
        case 1: --i; gaps+=countgaps; adata.profile.add(AT_INS); break;
        case 2: --j; gaps+=countgaps; adata.profile.add(AT_DEL); break;
       default:
        ldie("unknown code");
      }
    }
    for (;i>0; --i){
      ++gaps;
      adata.profile.add(AT_INS);
    }
    for (;j>0; --j){
      ++gaps;
      adata.profile.add(AT_DEL);
    }
  }

  if (adata.e1==-1 || adata.e1<pa+tmpsa) adata.e1=pa+tmpsa;
  if (adata.e2==-1 || adata.e2<pb+tmpsb) adata.e2=pb+tmpsb;

//  if (adata.s1==-1 || adata.s1>pa) adata.s1=pa;
//  if (adata.s2==-1 || adata.s2>pb) adata.s2=pb;
//  if (adata.e1==-1 || adata.e1<ea) adata.e1=ea;
//  if (adata.e2==-1 || adata.e2<eb) adata.e2=eb;
  adata.matches+=matches;
  adata.mismatches+=aligned-matches;
  adata.gaps+=gaps;
}




void seqcalign_local_leftext(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const ealignscore& as)
{
  uint64_t *psa=reinterpret_cast<uint64_t*>(a.seq._str),*psb=reinterpret_cast<uint64_t*>(b.seq._str);
  float mc_score[4]={as.match,-as.mismatch,-as.mismatch,-as.mismatch};
  float score=0.0;
  int aligned=0,gaps=0,matches=0,mismatches=0;
  long maxsize=MAX(ea-pa,eb-pb);
  long minsize=MIN(ea-pa,eb-pb);
  if (minsize==0) { // we set the leftmost part to the rightmost part of the sequence (in case we are aligning 0bp against Xbp long sequence)
    if (adata.s1==-1 || adata.s1>ea) adata.s1=ea;
    if (adata.s2==-1 || adata.s2>eb) adata.s2=eb;
//    if (adata.e1==-1 || adata.e1<ea) adata.e1=ea;
//    if (adata.e2==-1 || adata.e2<eb) adata.e2=eb;
    return;
  }

  int w=(ea-pa)+1;
  int h=(eb-pb)+1;
//  cout << "pa: " << pa << " ea: " << ea << " pb: " << pb << " eb: " << eb << " w: " << w << " h: " << h << endl;

  int i,j;
  ws.reserve(w*h*3);
  double *&mF(ws.mF);

/*
  for (i=0; i<w; ++i){
    for (j=0; j<h; ++j){
      mF[i*h*3+j*3]=-1.0e100l;
      mF[i*h*3+j*3+1]=-1.0e100l;
      mF[i*h*3+j*3+2]=-1.0e100l;
    }
  }
*/

  for (i=0; i<w; ++i){
    mF[i*h*3]=-1.0e100l;
    mF[i*h*3+1]=-as.gapopen-i*as.gapext;
    mF[i*h*3+2]=-1.0e100l;
  }
  for (j=0; j<h; ++j){
    mF[j*3]=-1.0e100l;
    mF[j*3+1]=-1.0e100l;
    mF[j*3+2]=-as.gapopen-j*as.gapext;
  }
  mF[0]=0.0l;
  mF[1]=-1.0e100l;
  mF[2]=-1.0e100l;
  double bscore=0.0l,tmpbs;
  long bi=0l,bj=0l;
  long j0=1l,j1=h;

  uint64_t ca,cb;
  for (i=1; i<w; ++i){
    ca=(psa[(ea-i)/32u]>>(((ea-i)%32u)*2ul))&0x03u;
    if (j0>1) { mF[i*h*3+(j0-1)*3]=-1.0e100l; mF[i*h*3+(j0-1)*3+1]=-1.0e100l; mF[i*h*3+(j0-1)*3+2]=-1.0e100l; }
    for (j=j0; j<h; ++j) {
      cb=(psb[(eb-j)/32u]>>(((eb-j)%32u)*2ul))&0x03u;
      if (j>j1) { mF[(i-1)*h*3+j*3]=-1.0e100l; mF[(i-1)*h*3+j*3+1]=-1.0e100l; mF[(i-1)*h*3+j*3+2]=-1.0e100l; }
      tmpbs=updateFNEG2(i,j,w,h,mF,mc_score[ca^cb],as);
      if (tmpbs<bscore-as.dropoff){  // dropoff algorithm
        if (j==j0) {
          ++j0; // hit left edge, skip this position in next row
//          cout << " j0: " << j0 << endl;
        } else {
          j1=j;
//          cout << " j1: " << j1 << endl;
          break; // hit right edge, skip next cells in this row
        }
      } else if (tmpbs>bscore) { // found better score
//        cout << "i: " << i << " j: " << j << " j0: " << j0 << " j1: " << j1 << " bscore: " << bscore << " tmpbs: "<<tmpbs << endl;
        bscore=tmpbs;
        bi=i; bj=j;
      }
    }
    if (j0>j1+1) break; // no more computation needed
  }


  adata._score+=bscore;
  if (adata.s1==-1 || adata.s1>ea-bi) adata.s1=ea-bi;
  if (adata.s2==-1 || adata.s2>eb-bj) adata.s2=eb-bj;

  aligned=0;
  matches=0;


  ealignprofile tmpprof;
  int ti=ea-(w-bi);
  i=bi; j=bj;
  int countgaps=0;
  while (i>0 && j>0){
    ca=(psa[(ea-i)/32u]>>(((ea-i)%32u)*2ul))&0x03u;
    cb=(psb[(eb-j)/32u]>>(((eb-j)%32u)*2ul))&0x03u;
    switch (maxv3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2])){
      case 0: countgaps=1; ++aligned; if (ca==cb) { ++matches; tmpprof.add(AT_MATCH); } else tmpprof.add(AT_MISS); --j; --i; break;
      case 1: --i; gaps+=countgaps; tmpprof.add(AT_INS); break;
      case 2: --j; gaps+=countgaps; tmpprof.add(AT_DEL); break;
     default:
      ldie("unknown code");
    }
  }
  for (;i>0 && countgaps; --i){
    ++gaps;
    adata.profile.add(AT_INS);
  }
  for (;j>0 && countgaps; --j){
    ++gaps;
    adata.profile.add(AT_DEL);
  }

  adata.matches+=matches;
  adata.mismatches+=aligned-matches;
  adata.gaps+=gaps;
  adata.profile.addinv(tmpprof); // add inverted profile
}

void seqcalign_local_rightext(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const ealignscore& as)
{
  uint64_t *psa=reinterpret_cast<uint64_t*>(a.seq._str),*psb=reinterpret_cast<uint64_t*>(b.seq._str);
  float mc_score[4]={as.match,-as.mismatch,-as.mismatch,-as.mismatch};
  float score=0;
  int aligned=0,gaps=0,matches=0,mismatches=0;
  long maxsize=MAX(ea-pa,eb-pb);
  long minsize=MIN(ea-pa,eb-pb);
  if (minsize==0) {
    if (adata.e1==-1 || adata.e1<pa) adata.e1=pa;
    if (adata.e2==-1 || adata.e2<pb) adata.e2=pb;
    return;
  }

  int w=(ea-pa)+1;
  int h=(eb-pb)+1;

  int i,j;
  ws.reserve(w*h*3);
  double *&mF(ws.mF);

/*
  for (i=0; i<w; ++i){
    for (j=0; j<h; ++j){
      mF[i*h*3+j*3]=-1.0e100l;
      mF[i*h*3+j*3+1]=-1.0e100l;
      mF[i*h*3+j*3+2]=-1.0e100l;
    }
  }
*/

  for (i=0; i<w; ++i){
    mF[i*h*3]=-1.0e100l;
    mF[i*h*3+1]=-as.gapopen-i*as.gapext;
    mF[i*h*3+2]=-1.0e100l;
  }
  for (j=0; j<h; ++j){
    mF[j*3]=-1.0e100l;
    mF[j*3+1]=-1.0e100l;
    mF[j*3+2]=-as.gapopen-j*as.gapext;
  }
  mF[0]=0.0l;
  mF[1]=-1.0e100l;
  mF[2]=-1.0e100l;
  double bscore=0.0l,tmpbs;
  long bi=0l,bj=0l;
  long j0=1l,j1=h;

  uint64_t ca,cb;
  for (i=1; i<w; ++i){
    ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    if (j0>1) { mF[i*h*3+(j0-1)*3]=-1.0e100l; mF[i*h*3+(j0-1)*3+1]=-1.0e100l; mF[i*h*3+(j0-1)*3+2]=-1.0e100l; }
    for (j=j0; j<h; ++j) {
      cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      if (j>j1) { mF[(i-1)*h*3+j*3]=-1.0e100l; mF[(i-1)*h*3+j*3+1]=-1.0e100l; mF[(i-1)*h*3+j*3+2]=-1.0e100l; j1=j; }
      tmpbs=updateFNEG2(i,j,w,h,mF,mc_score[ca^cb],as);
      if (tmpbs<bscore-as.dropoff){  // dropoff algorithm   //TODO: need to initialize border when skipping rows and columns
        if (j==j0) {
          ++j0; // hit left edge, skip this postion in next row
        } else {
          j1=j;
          break; // hit right edge, skip next cells in this row
        }
      } else 
      if (tmpbs>bscore) { // found better score
        bscore=tmpbs;
        bi=i; bj=j;
      }
    }
    if (j0>j1+1) break; // no more computation needed
  }

//  cout << "# bscore: " << bscore << " bi: " << bi << " bj: " << bj << endl;

  if (adata.e1==-1 || adata.e1<pa+bi) adata.e1=pa+bi;
  if (adata.e2==-1 || adata.e2<pb+bj) adata.e2=pb+bj;
  adata._score+=bscore;

  aligned=0;
  matches=0;

//  int tmpsa=0,tmpsb=0;
  int ti=ea-(w-bi);
  i=bi; j=bj;
  int countgaps=0;
  while (i>0 && j>0){
    ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    switch (maxv3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2])){
      case 0: countgaps=1; ++aligned; if (ca==cb) { ++matches; adata.profile.add(AT_MATCH); } else adata.profile.add(AT_MISS); --j; --i; break;
      case 1: --i; gaps+=countgaps; adata.profile.add(AT_INS); break;
      case 2: --j; gaps+=countgaps; adata.profile.add(AT_DEL); break;
     default:
      ldie("unknown code");
    }
  }
  for (;i>0 && countgaps>0; --i){
    ++gaps;
    adata.profile.add(AT_INS);
  }
  for (;j>0 && countgaps>0; --j){
    ++gaps;
    adata.profile.add(AT_DEL);
  }

//  if (adata.s1==-1 || adata.s1>pa) adata.s1=pa;
//  if (adata.s2==-1 || adata.s2>pb) adata.s2=pb;
//  if (adata.e1==-1 || adata.e1<ea) adata.e1=ea;
//  if (adata.e2==-1 || adata.e2<eb) adata.e2=eb;
  adata.matches+=matches;
  adata.mismatches+=aligned-matches;
  adata.gaps+=gaps;
}





ealignws::ealignws(): size(0),mF(0x00) {}
ealignws::~ealignws() { if (mF) delete mF; }
void ealignws::reserve(int _size)
{
  if (mF && size < _size){
    delete mF;
    mF=0x00;
  }
  if (mF==0x00){
    size=_size;
    mF=new double[size];
  }
}

ealignelem::ealignelem(): type(AT_NONE),count(0) {}
ealignelem::ealignelem(ealignelemtype _type,int _count): type(_type),count(_count) {}


bool ealignelem::operator!=(const ealignelem& e) const
{
  if (e.type!=type || e.count!=count) return(true);
  return(false);
}
bool ealignelem::operator==(const ealignelem& e) const
{
  if (e.type!=type || e.count!=count) return(false);
  return(true);
}



bool ealignprofile::operator==(const ealignprofile& p) const
{
  if (elm.size()!=p.elm.size()) return(false);
  for (int i=0; i<elm.size(); ++i)
    if (elm[i]!=p.elm[i]) return(false);
  return(true);
}

void ealignprofile::add(const ealignprofile& p)
{
  for (int j=0; j<p.elm.size(); ++j)
    add(p.elm[j].type,p.elm[j].count);
}

void ealignprofile::addinv(const ealignprofile& p)
{
  for (int j=p.elm.size()-1; j>=0; --j)
    add(p.elm[j].type,p.elm[j].count);
}

void ealignprofile::inv()
{
  ealignprofile tmp;
  tmp.addinv(*this);
  elm.clear();
  add(tmp);
}

void ealignprofile::add(ealignelemtype type,int count)
{
  if (elm.size()==0 || elm[elm.size()-1].type!=type)
    elm.add(ealignelem(type,count));
  else
    elm[elm.size()-1].count+=count;
}

estr ealignprofile::str()
{
  estr tmpstr;
  for (int i=0; i<elm.size(); ++i){
    switch (elm[i].type){
      case AT_NONE: tmpstr+='N'; break;
      case AT_INS: tmpstr+='I'; break;
      case AT_DEL: tmpstr+='D'; break;
      case AT_MATCH: tmpstr+='M'; break;
      case AT_MISS: tmpstr+='S'; break;
      case AT_LEFT: tmpstr+='<'; break;
      case AT_RIGHT: tmpstr+='>'; break;
      case AT_ID: tmpstr+='='; break;
      case AT_ALIGN: tmpstr+='|'; break;
      case AT_PAIR: tmpstr+='/'; break;
    }
    if (elm[i].count>1)
      tmpstr+=elm[i].count;
  }
  return(tmpstr);
}



ostream& operator<<(ostream& stream,const ealignprofile& p)
{
  char tmp;
  for (int i=0; i<p.elm.size(); ++i){
    switch (p.elm[i].type){
      case AT_NONE: tmp='N'; break;
      case AT_INS: tmp='I'; break;
      case AT_DEL: tmp='D'; break;
      case AT_MATCH: tmp='M'; break;
      case AT_MISS: tmp='S'; break;
      case AT_LEFT: tmp='<'; break;
      case AT_RIGHT: tmp='>'; break;
      case AT_ID: tmp='='; break;
      case AT_ALIGN: tmp='|'; break;
      case AT_PAIR: tmp='/'; break;
    }
    stream << tmp;
    if (p.elm[i].count>1)
      stream << p.elm[i].count;
  }
  return(stream);
}

int getnum(const estr& astr,int& i)
{
  int j;
  for (j=i+1; j<astr.len() && isdigit(astr[j]); ++j);
  if (j==i+1) return(1);
   
  int count=astr.substr(i+1,j-i-1).i();
  i=j-1;
  return(count);
}

estr sali_decompress(const estr& astr,const eseq& seq2)
{
  estr str;
  int count;
  uint64_t *ps2=reinterpret_cast<uint64_t*>(seq2.seq._str);
  int i2=0;
  for (int i=0; i<astr.len(); ++i){
    switch(astr[i]){
      case 'a':
      case 't':
      case 'g':
      case 'c':
      case 'A':
      case 'T':
      case 'G':
      case 'C':
        ++i2;
        str+=astr[i];
       break;
      case 'I':
        count=getnum(astr,i);
        for (int j=i+1; j<astr.len() && j<i+1+count; ++j)
          str+=astr[j];
        i+=count;
//        cout << "I: " << count << endl;
       break;
      case 'D':
        count=getnum(astr,i);
        i2+=count;
//        cout << "D: " << count << endl;
       break;
      case 'M':
        count=getnum(astr,i);
        for (int j2=i2; i2<seq2.seqlen && j2<i2+count; ++j2)
          str+=cnuc2chru((ps2[j2/32u]>>((j2%32u)*2ul))&0x03u);
        i2+=count;       
//        cout << "M: " << count << endl;
       break;
      default:
       lerror("unexpected character: "+estr(astr[i]));
       return(estr());
    }
  }
  return(str);
}

estr ealigndata::compress(const eseq& seq1)
{
  uint64_t *ps1=reinterpret_cast<uint64_t*>(seq1.seq._str);
  ldieif(s1<0 || s1>seq1.seqlen || e1<0 || e1>seq1.seqlen,"start/end mismatch with sequence: "+estr(s1)+","+e1+" seqlen: "+seq1.seqlen);
  eseq srev(seq1);
  long ts1=s1,te1=e1;
  if (revcompl){
    ts1=seq1.seqlen-e1; te1=seq1.seqlen-s1;
    srev.revcompl();
//    srev.setrevcompl(seq1,0,seq1.seqlen);
    ps1=reinterpret_cast<uint64_t*>(srev.seq._str);
  }
  int i=ts1;
  estr str1;
  if (i>0){
    str1+="I";
    if (i>1)
      str1+=i;
    for (int l=0; l<i; ++l)
      str1+=cnuc2chru((ps1[l/32u]>>((l%32u)*2ul))&0x03u);
  }
  if (s2>0){
    str1+="D";
    if (s2>1)
      str1+=s2;
  }

  for (int k=0; k<profile.elm.size(); ++k){
    ealignelem& e(profile.elm[k]);
    switch (e.type){
      case AT_MATCH:
        str1+="M";
        if (e.count>1)
          str1+=e.count;
        i+=e.count;
       break;
      case AT_MISS:
        for (int l=0; l<e.count; ++l){
          str1+=cnuc2chru((ps1[i/32u]>>((i%32u)*2ul))&0x03u);
          ++i;
        }
       break;
      case AT_DEL:
        str1+="D";
        if (e.count>1)
          str1+=e.count;
       break;
      case AT_INS:
        str1+="I";
        if (e.count>1)
          str1+=e.count;
        for (int l=0; l<e.count; ++l){
          str1+=cnuc2chru((ps1[i/32u]>>((i%32u)*2ul))&0x03u);
          ++i;
        }
       break;
    }
  }
  if (i<seq1.seqlen){
    str1+="I";
    if (i>1)
      str1+=i;
    for (int l=0; l<i; ++l)
      str1+=cnuc2chru((ps1[l/32u]>>((l%32u)*2ul))&0x03u);
  }
  return(str1);
}



estr ealigndata::align_str(const eseq& seq1,const eseq& seq2)
{
  uint64_t *ps1=reinterpret_cast<uint64_t*>(seq1.seq._str),*ps2=reinterpret_cast<uint64_t*>(seq2.seq._str);
  ldieif(s1<0 || s1>seq1.seqlen || e1<0 || e1>seq1.seqlen,"start/end mismatch with sequence: "+estr(s1)+","+e1+" seqlen: "+seq1.seqlen);
  eseq srev(seq1);
  long ts1=s1,te1=e1;
  if (revcompl){
    ts1=seq1.seqlen-e1; te1=seq1.seqlen-s1;
    srev.revcompl();
//    srev.setrevcompl(seq1,0,seq1.seqlen);
    ps1=reinterpret_cast<uint64_t*>(srev.seq._str);
  }
  int i=ts1,j=s2;
  estr str1,strs,str2;
  for (int k=0; k<profile.elm.size(); ++k){
    ealignelem& e(profile.elm[k]);
    switch (e.type){
      case AT_MATCH:
        for (int l=0; l<e.count; ++l){
          str1+=cnuc2chr((ps1[i/32u]>>((i%32u)*2ul))&0x03u);
          strs+='|';
          str2+=cnuc2chr((ps2[j/32u]>>((j%32u)*2ul))&0x03u);
          ++i; ++j;
        }
       break;
      case AT_MISS:
        for (int l=0; l<e.count; ++l){
          str1+=cnuc2chr((ps1[i/32u]>>((i%32u)*2ul))&0x03u);
          strs+=' ';
          str2+=cnuc2chr((ps2[j/32u]>>((j%32u)*2ul))&0x03u);
          ++i; ++j;
        }
       break;
      case AT_DEL:
        for (int l=0; l<e.count; ++l){
          str1+='-';
          strs+=' ';
          str2+=cnuc2chr((ps2[j/32u]>>((j%32u)*2ul))&0x03u);
          ++j;
        }
       break;
      case AT_INS:
        for (int l=0; l<e.count; ++l){
          str1+=cnuc2chr((ps1[i/32u]>>((i%32u)*2ul))&0x03u);
          strs+=' ';
          str2+='-';
          ++i;
        }
       break;
    }
  }
  return(str1+"\n"+strs+"\n"+str2);
}








void seqcalign_global_noedgegap(const eseq& a,long pa,long ea,const eseq& b,long pb,long eb,ealigndata& adata,ealignws& ws,const ealignscore& as)
{
  uint64_t *psa=reinterpret_cast<uint64_t*>(a.seq._str),*psb=reinterpret_cast<uint64_t*>(b.seq._str);
  float mc_score[4]={as.match,-as.mismatch,-as.mismatch,-as.mismatch};
  float score=0;
  int aligned=0,gaps=0,matches=0,mismatches=0;
  long maxsize=MAX(ea-pa,eb-pb);
  if (maxsize==0) return;

  int w=(ea-pa)+1;
  int h=(eb-pb)+1;

  int i,j;
//  double *mF=new double[w*h*3];
  ws.reserve(w*h*3);
  double *&mF(ws.mF);

  // no gap penalty for beginning gaps
  for (i=0; i<w; ++i){
    mF[i*h*3+1]=0.0l;
    mF[i*h*3+2]=-1.0e100l;
    mF[i*h*3]=-1.0e100l;
  }
  for (j=0; j<h; ++j){
    mF[j*3+1]=-1.0e100l;
    mF[j*3+2]=0.0l;
    mF[j*3]=-1.0e100l;
  }
  mF[0]=0.0l;
  mF[1]=-as.gapopen;
  mF[2]=-as.gapopen;

  uint64_t ca,cb;
  for (i=1; i<w; ++i){
    ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    for (j=1; j<h; ++j) {
      cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      updateFNEG(i,j,w,h,mF,mc_score[ca^cb],as);
    }
  }

  aligned=0;
  matches=0;

  int end1=-1;
  int end2=-1;

  int ti=ea-1;
  int tj=eb-1;
  i=w-1; j=h-1;
  adata._score+=MAX3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2]);
  int countgaps=0;
  if (adata.aln){
    while (i>0 && j>0){
      ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
      cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      switch (maxv3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2])){
        case 0: if (end2==-1) end2=tj; if (end1==-1) end1=ti; countgaps=1; ++aligned; if (ca==cb) { adata.aln[ti]=0x01u; ++matches; adata.profile.add(AT_MATCH); } else { adata.aln[ti]=0x02u; adata.profile.add(AT_MISS); } --j; --i; --ti; --tj; adata.profile.add(AT_MATCH); break;
        case 1: adata.aln[ti]=0x00u; --i; --ti; gaps+=countgaps; adata.profile.add(AT_INS); break;
        case 2: --j; --tj; gaps+=countgaps; adata.profile.add(AT_DEL); break;
       default:
        ldie("unknown code");
      }
    }
    for (;j==0 && i>0;--i,--ti)
      adata.aln[ti]=0x00u;
  }else{
    while (i>0 && j>0){
      ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
      cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      switch (maxv3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2])){
        case 0: if (end2==-1) end2=tj; if (end1==-1) end1=ti; countgaps=1; ++aligned; if (ca==cb) { ++matches; adata.profile.add(AT_MATCH); } else adata.profile.add(AT_MISS); --j; --i; --ti; --tj; break;
        case 1: --i; --ti; gaps+=countgaps; adata.profile.add(AT_INS); break;
        case 2: --j; --tj; gaps+=countgaps; adata.profile.add(AT_DEL); break;
       default:
        ldie("unknown code");
      }
    }
  }

  if (adata.s1==-1 || adata.s1>pa+i) adata.s1=pa+i;
  if (adata.s2==-1 || adata.s2>pb+j) adata.s2=pb+j;
  if (adata.e1==-1 || adata.e1<end1) adata.e1=end1;
  if (adata.e2==-1 || adata.e2<end2) adata.e2=end2;
  adata.matches+=matches;
  adata.mismatches+=aligned-matches;
  adata.gaps+=gaps;
}

void seqcalign_nb_global_noedgegap(const eseq& a,int pa,int ea,const eseq& b,int pb,int eb,ealigndata& adata,unsigned char *tncount,ealignws& ws,const ealignscore& as)
{
  uint64_t *psa=reinterpret_cast<uint64_t*>(a.seq._str),*psb=reinterpret_cast<uint64_t*>(b.seq._str);
  float mc_score[4]={as.match,-as.mismatch,-as.mismatch,-as.mismatch};
  float score=0;
  int aligned=0,gaps=0,matches=0,mismatches=0;
  int maxsize=MAX(ea-pa,eb-pb);
  if (maxsize==0) return;

  int w=(ea-pa)+1;
  int h=(eb-pb)+1;

  int i,j;
//  double *mF=new double[w*h*3];
  ws.reserve(w*h*3);
  double *&mF(ws.mF);

  // no gap penalty for beginning gaps
  for (i=0; i<w; ++i){
    mF[i*h*3+1]=0.0l;
    mF[i*h*3+2]=-1.0e100l;
    mF[i*h*3]=-1.0e100l;
  }
  for (j=0; j<h; ++j){
    mF[j*3+1]=-1.0e100l;
    mF[j*3+2]=0.0l;
    mF[j*3]=-1.0e100l;
  }
  mF[0]=0.0l;
  mF[1]=-as.gapopen;
  mF[2]=-as.gapopen;

  uint64_t ca,cb;
  for (i=1; i<w; ++i){
    ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    for (j=1; j<h; ++j) {
      cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
      updateFNEG(i,j,w,h,mF,mc_score[ca^cb],as);
    }
  }

  aligned=0;
  matches=0;

  int end1=-1;
  int end2=-1;

  int ti=ea-1;
  int tj=eb-1;
  i=w-1; j=h-1;
  int countgaps=0;
  adata._score=MAX3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2]);
  while (i>0 && j>0){
    ca=(psa[(i-1+pa)/32u]>>(((i-1+pa)%32u)*2ul))&0x03u;
    cb=(psb[(j-1+pb)/32u]>>(((j-1+pb)%32u)*2ul))&0x03u;
    switch (maxv3(mF[i*h*3+j*3],mF[i*h*3+j*3+1],mF[i*h*3+j*3+2])){
      case 0: if (end2==-1) end2=tj; if (end1==-1) end1=ti; if (ca==cb) ++matches; countgaps=1; ++aligned; ++tncount[4*ti+cb]; --j; --i; --ti; --tj; break;
      case 1: --i; --ti; gaps+=countgaps; break;
      case 2: --j; --tj; gaps+=countgaps; break;
     default:
      ldie("unknown code");
    }
  }
  adata.s1=pa+i;
  adata.s2=pb+j;
  adata.e1=end1;
  adata.e2=end2;
  adata.matches=matches;
  adata.mismatches=aligned-matches;
  adata.gaps=gaps;
}


