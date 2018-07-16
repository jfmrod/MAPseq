#include "eseq.h"

//#include "kmerseqtables-data.h"
extern unsigned char seq_comp_table[];
extern unsigned char prot_comp_table[];

#define PKMERSIZE 9ul
#define PKMERBITS (PKMERSIZE*2ul)
#define PMAXSIZE (1ul<<PKMERBITS)
#define PKMERMAX (1ul<<PKMERBITS)
#define PKMERMASK (PMAXSIZE-1ul)

const uint32_t revnuc[]={0x1u,0x0u,0x3u,0x2u};
const unsigned long safe_shift[32]={0x0ul,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful};

char cnuc2chr(uint32_t cc)
{
  const char larr[]={'a','t','g','c'};
  return(larr[cc&0x3]);
}

char cnuc2chru(uint32_t cc)
{
  const char larr[]={'A','T','G','C'};
  return(larr[cc&0x3]);
}

unsigned long seqpkmer(const eseq& s,long p1)
{
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s.seq._str);
  return(((pstr1[p1/32u]>>(2u*(p1%32u)))|((pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u]))&PKMERMASK);
}

void eseq::serial(estr& sstr)
{
  seriallong(seqlen,sstr);
  seq.serial(sstr);
}

long eseq::unserial(const estr& sstr,long i)
{
  i=unseriallong(seqlen,sstr,i);
  if (i==-1) return(i);
  i=seq.unserial(sstr,i);
  return(i);
}

bool eseq::operator<(const eseq& s)
{
  return(seqlen>s.seqlen);
}

eseq::eseq(): seqlen(0),seqstart(0),prot(false) {}

eseq::eseq(const estr& ucseq): seqlen(0),seqstart(0),prot(false)
{
  setseq(ucseq);
}

void eseq::setprot(estr& ucseq)
{
  prot=true;
  long i;
  uint32_t tmp;
//  estr useq;
//  useq.reserve(int((ucseq.len()+3)/4)*4);
//  seq.reserve(int((ucseq.len()+3)/4));
  seq.clear();
  seqlen=ucseq.len()*3;
  if (ucseq.len()==0) return;

  seq.reserve(((ucseq.len()*3+31u)>>4)<<2);
  memset(seq._str,0x0u,((ucseq.len()*3+31u)>>4)<<2); 

  //TODO: skip spaces and tabs
/*
  for (i=0; i<ucseq.len(); ++i){
    tmp=prot_comp_table[ucseq._str[i]] | (prot_comp_table[ucseq._str[i+1]]<<6u) |
        (prot_comp_table[ucseq._str[i+2]]<<12u) | (prot_comp_table[ucseq._str[i+3]]<<18u) | 
        (prot_comp_table[ucseq._str[i+4]]<<24u) | (prot_comp_table[ucseq._str[i+5]]<<30u);
    
  }
*/
 
  tmp=(prot_comp_table[ucseq._str[0]]<<6u);
  tmp>>=3*2u;
  for (; i<ucseq.len(); ++i){
    tmp|=(prot_comp_table[ucseq._str[i]]<<6u);
    seq._str[i*6/8]=(tmp>>((6-(i*6)%8)))&0xff;
    tmp>>=3*2u;
  }
//  seq._strlen=int((ucseq.len()+3)/4);
  seq._strlen=((ucseq.len()*3+31u)>>4)<<2;
}

void eseq::setseq(const estr& ucseq)
{
  prot=false;
  long i;
  uint32_t tmp;
//  estr useq;
//  useq.reserve(int((ucseq.len()+3)/4)*4);
//  seq.reserve(int((ucseq.len()+3)/4));
  seq.clear();
  npos.clear();
  seq.reserve(((ucseq.len()+31u)>>4)<<2);

/*
  ucseq.reserve(int((ucseq.len()+3)/4)*4);
  useq=ucseq;
  useq.lowercase();
  useq.replace("u","t");
*/

/*
  unsigned long pstr=reinterpret_cast<unsigned long*>(seq._str);
  unsigned long sstr=reinterpret_cast<unsigned long*>(ucseq._str);
  for (i=0; i<ucseq.len()-4; i+=4,++pstr){
    unsigned long tmp1=sstr[i/8];
    unsigned long tmp2=sstr[i/8+1];
    unsigned long tmp3=sstr[i/8+2];
    unsigned long tmp4=sstr[i/8+3];
    *pstr=seq_comp_table[tmp1&0xFFFFul]|(seq_comp_table[(tmp1>>16u)&0xFFFFul]<<4u)|(seq_comp_table[(tmp1>>32u)&0xFFFFul]<<8u)|(seq_comp_table[(tmp1>>48u)&0xFFFFul]<<12u);
    *pstr|=(seq_comp_table[tmp2&0xFFFFul]<<16u)|(seq_comp_table[(tmp2>>16u)&0xFFFFul]<<20u)|(seq_comp_table[(tmp2>>32u)&0xFFFFul]<<24u)|(seq_comp_table[(tmp2>>48u)&0xFFFFul]<<28u);
    *pstr|=(seq_comp_table[tmp3&0xFFFFul]<<32u)|(seq_comp_table[(tmp3>>16u)&0xFFFFul]<<36u)|(seq_comp_table[(tmp3>>32u)&0xFFFFul]<<40u)|(seq_comp_table[(tmp3>>48u)&0xFFFFul]<<44u);
    *pstr|=(seq_comp_table[tmp4&0xFFFFul]<<48u)|(seq_comp_table[(tmp4>>16u)&0xFFFFul]<<52u)|(seq_comp_table[(tmp4>>32u)&0xFFFFul]<<56u)|(seq_comp_table[(tmp4>>48u)&0xFFFFul]<<60u);
  }
*/

  //TODO: skip spaces and tabs

  for (i=0; i/4<ucseq.len()/4; i+=4){
    tmp=*(uint32_t*)(&ucseq._str[i]);
    seq._str[i/4]=seq_comp_table[tmp&0xffffu]|(seq_comp_table[(tmp>>16)&0xffffu]<<4);
//    for (int j=0; j<4; ++j,tmp>>=8u)
    for (int j=i; j<i+4; ++j)
      if (ucseq[j]=='N')
        npos.add(j);
  }
  switch (ucseq.len()%4){
    case 3: 
      tmp=*(uint32_t*)(&ucseq._str[i]);
      seq._str[i/4]=seq_comp_table[tmp&0xffffu]|(seq_comp_table[(tmp>>16)&0x00ffu]<<4);
      for (int j=i; j<i+3; ++j)
        if (ucseq[j]=='N')
          npos.add(j);
     break;
    case 2: 
      tmp=*(uint32_t*)(&ucseq._str[i]);
      seq._str[i/4]=seq_comp_table[tmp&0xffffu];
      for (int j=i; j<i+2; ++j)
        if (ucseq[j]=='N')
          npos.add(j);
     break;
    case 1: 
      tmp=*(uint32_t*)(&ucseq._str[i]);
      seq._str[i/4]=seq_comp_table[tmp&0x00ffu];
      if (ucseq[i]=='N')
        npos.add(i);
     break;
  }
//  seq._strlen=int((ucseq.len()+3)/4);
  seq._strlen=((ucseq.len()+31u)>>4)<<2;
  seqlen=ucseq.len();
}

void eseq::revcompl()
{
  estr tmpseq;
  tmpseq.reserve(seq.len());
  tmpseq._strlen=seq._strlen;
  
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(tmpseq._str);
  unsigned long v1,v2;
  int p1=seqlen-1;
  int p2=0;

  for (; p1>=0; --p1,++p2){
    v1=(pstr1[p1/32u]>>(2u*(p1%32u)))&0x3ul;
    pstr2[p2/32u]=(pstr2[p2/32u]&(~(0x3ul<<(2u*(p2%32u))))) | (((long int)(revnuc[v1]))<<(2u*(p2%32u)));
  }
  seq=tmpseq;
  if (npos.size()){
    eintarray tmpnpos;
    for (int j=npos.size()-1; j>=0; --j)
      tmpnpos.add(seqlen-npos[j]);
    npos=tmpnpos;
  }
}


void eseq::setrevcompl(const eseq& rseq,long i,long e)
{
  seq.reserve(((e-i+31u)>>4)<<2);
  seq._strlen=((e-i+31u)>>4)<<2;
  seqlen=e-i;

  unsigned long *pstr1=reinterpret_cast<unsigned long*>(rseq.seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(seq._str);
  unsigned long v1,v2;
  int p1=e-1;
  int p2=0;

  for (; p1>=i; --p1,++p2){
    v1=(pstr1[p1/32u]>>(2u*(p1%32u)))&0x3ul;
    pstr2[p2/32u]=(pstr2[p2/32u]&(~(0x3ul<<(2u*(p2%32u))))) | (((long int)(revnuc[v1]))<<(2u*(p2%32u)));
  }
}

eseq eseq::subseq(int i,int e) const
{
  ldieif(e>seqlen,"end of subsequence past sequence: "+estr(e)+" seqlen: "+seqlen);
  ldieif(e<0,"end of subsequence is negative: "+estr(e)+" seqlen: "+seqlen);
  ldieif(i<0,"start of subsequence is negative: "+estr(i)+" seqlen: "+seqlen);
  ldieif(e<i,"end of subsequence before beginning of subsequence: "+estr(i)+" "+estr(e)+" seqlen: "+seqlen);
  eseq tmps;
  tmps.seq.reserve(((e-i+31u)>>4)<<2);
  tmps.seq._strlen=((e-i+31u)>>4)<<2;
  tmps.seqlen=e-i;
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(tmps.seq._str);
  unsigned long v1;
  for (; i<e; ++pstr2,i+=32){
    v1=pstr1[i/32u]>>(2u*(i%32u));
    v1|=(pstr1[i/32u+1u]<<(64u-2u*(i%32u)))&safe_shift[i%32u];
    *pstr2=v1;
  }
  return(tmps);
}

estr eseq::print_pseq(int i,int e) const
{
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(seq._str);
  unsigned long v1;
  if (e==-1) e=seqlen;
  ldieif(i>e,"error start pos larger than end pos");
  estr tmp;
  tmp.reserve(e-i);
  for (; i<e; i+=3){
    v1=(i%32u==0?pstr1[i/32u]:pstr1[i/32u]>>(2u*(i%32u)))&((0x1u<<6u)-1u);
    switch(v1){
      case 0x09u: tmp+='*'; break;
      case 0x0Eu: tmp+='A'; break;
      case 0x10u: tmp+='N'; break;
      case 0x19u: tmp+='C'; break;
      case 0x0Fu: tmp+='P'; break;
      case 0x12u: tmp+='D'; break;
      case 0x03u: tmp+='Q'; break;
      case 0x02u: tmp+='E'; break;
      case 0x0Bu: tmp+='R'; break;
      case 0x15u: tmp+='F'; break;
      case 0x0Du: tmp+='S'; break;
      case 0x0Au: tmp+='G'; break;
      case 0x0Cu: tmp+='T'; break;
      case 0x13u: tmp+='H'; break;
      case 0x06u: tmp+='V'; break;
      case 0x04u: tmp+='I'; break;
      case 0x29u: tmp+='W'; break;
      case 0x00u: tmp+='K'; break;
      case 0x11u: tmp+='Y'; break;
      case 0x07u: tmp+='L'; break;
      case 0x24u: tmp+='M'; break;
     default:
      tmp+="!";
    }
  }
  return(tmp);
}

estr eseq::print_seq(int i,int e) const
{
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(seq._str);
  unsigned long v1;
  if (e==-1) e=seqlen;
  ldieif(i>e,"error start pos larger than end pos");
  estr tmp;
  tmp.reserve(e-i);
  for (int ti=i; ti<e; ++ti){
    v1=(ti%32u==0?pstr1[ti/32u]:pstr1[ti/32u]>>(2u*(ti%32u)))&0x3u;
    switch(v1){
      case 0: tmp+='A'; break;
      case 1: tmp+='T'; break;
      case 2: tmp+='G'; break;
      case 3: tmp+='C'; break;
    }
  }
  for (int j=0; j<npos.size(); ++j){
    if (npos[j]<i) continue;
    else if (npos[j]>=e) break;
    tmp[npos[j]-i]='N';
  }
    
  return(tmp);
}

ostream& operator<<(ostream& stream,const eseq& seq)
{
  for (int i=0; i<seq.seq.len(); ++i){
    unsigned char t=seq.seq[i];
    for (int k=0; k<4; ++k,t>>=2u){
      switch (0x03u&t){
        case 0x00u: stream << "a"; break;
        case 0x01u: stream << "u"; break;
        case 0x02u: stream << "g"; break;
        case 0x03u: stream << "c"; break;
      }
    }
  }
  return(stream);
}


