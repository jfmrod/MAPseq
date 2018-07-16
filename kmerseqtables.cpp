#include <eutils/emain.h>
#include <eutils/efile.h>

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

#define PKMERSIZE 9ul
#define PKMERBITS (PKMERSIZE*2ul)
#define PMAXSIZE (1ul<<PKMERBITS)
#define PKMERMAX (1ul<<PKMERBITS)
#define PKMERMASK (PMAXSIZE-1ul)

#define IKMERSIZE 8ul
#define IKMERBITS (IKMERSIZE*2ul)
#define IKMERMAX (1ul<<IKMERBITS)
#define IKMERMASK (IKMERMAX-1ul)

#include "data/prot_trans_table.h"
const uint32_t aacount=21u;
char aa[]={'A','N','C','P','D','Q','E','R','F','S','G','T','H','V','I','W','K','Y','L','M','*'};
uint32_t aacomp[]={0x0u,0x1u,0x2u,0x3u,0x4u,0x5u,0x6u,0x7u,0x8u,0x9u,0xAu,0xBu,0xCu,0xDu,0xEu,0xFu,0x10u,0x11u,0x12u,0x13u,0x14u};

const int prot_blosum62_size=aacount*aacount;
int8_t prot_blosum62[prot_blosum62_size];



uint32_t nuccount=12u;
uint32_t nuc[]={'a','t','g','c','A','T','u','U','G','C','n','N'};
uint32_t compnuc[]={0x0u,0x1u,0x2u,0x3u,0x0u,0x1u,0x1u,0x1u,0x2u,0x3u,0x0u,0x0u};
uint32_t revnuc[]={0x1u,0x0u,0x3u,0x2u};
uint32_t consnuc[]={0x0001u,0x0010u,0x0100u,0x1000u};

const unsigned int codon2prot_size=1u<<8u;
uint8_t codon2prot[codon2prot_size];

unsigned char seq_match_table[IKMERMAX];
unsigned char seq_comp_table[IKMERMAX];
unsigned char prot_comp_table[1u<<8u];
unsigned char seq_ident_table[IKMERMAX];
unsigned int kmer_prot_lt[PKMERMAX];
unsigned int kmer_protrev_lt[PKMERMAX];
unsigned long kmer_rev_lt[KMERMAX];
unsigned long kmer_rev_lt2[KMERMAX2];
uint64_t seq_alignment_lt[IKMERMAX];

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

void char2int(char *tmpstr,char c)
{
  sprintf(tmpstr,"%hhi",c);
}

void uchr2hex(char *tmpstr,unsigned char c)
{
  tmpstr[0]='0';
  tmpstr[1]='x';
  tmpstr[2]=oct2hex((c>>4u)&0xf);
  tmpstr[3]=oct2hex(c&0xf);
  tmpstr[4]=0x00;
}

void int2hex(char *tmpstr,unsigned int c)
{
  tmpstr[0]='0';
  tmpstr[1]='x';
  for (unsigned int k=0; k<sizeof(c)*8u/4u; ++k)
    tmpstr[2+k]=oct2hex((c>>(sizeof(c)*8u-4u-k*4u))&0xf);
//  tmpstr[2]=oct2hex((c>>12u)&0xf);
//  tmpstr[3]=oct2hex((c>>8u)&0xf);
//  tmpstr[4]=oct2hex((c>>4u)&0xf);
//  tmpstr[5]=oct2hex(c&0xf);
  tmpstr[2+sizeof(c)*8u/4u]=0x00;
}

void uint2hex(char *tmpstr,unsigned int c)
{
  tmpstr[0u]='0';
  tmpstr[1u]='x';
  for (unsigned int k=0; k<sizeof(c)*8u/4u; ++k)
    tmpstr[2u+k]=oct2hex((c>>(sizeof(c)*8u-4u-k*4u))&0xf);
  tmpstr[2u+sizeof(c)*8u/4u]='u';
  tmpstr[2u+sizeof(c)*8u/4u+1u]=0x00;
}


void long2hex(char *tmpstr,unsigned long c)
{
  tmpstr[0u]='0';
  tmpstr[1u]='x';
  for (unsigned int k=0; k<sizeof(c)*8u/4u; ++k)
    tmpstr[2u+k]=oct2hex((c>>(sizeof(c)*8u-4u-k*4u))&0xf);
  tmpstr[2u+sizeof(c)*8u/4u]='l';
  tmpstr[2u+sizeof(c)*8u/4u+1u]=0x00;
}

void ulong2hex(char *tmpstr,unsigned long c)
{
  tmpstr[0u]='0';
  tmpstr[1u]='x';
  for (unsigned int k=0; k<sizeof(c)*8u/4u; ++k)
    tmpstr[2u+k]=oct2hex((c>>(sizeof(c)*8u-4u-k*4u))&0xf);
  tmpstr[2u+sizeof(c)*8u/4u]='u';
  tmpstr[2u+sizeof(c)*8u/4u+1u]='l';
  tmpstr[2u+sizeof(c)*8u/4u+2u]=0x00;
}
void initMatchTable()
{
  unsigned char tmp[4];
  tmp[0x00u]=1u;
  tmp[0x01u]=0u;
  tmp[0x02u]=0u;
  tmp[0x03u]=0u;
  
  for (uint32_t i=0u; i<IKMERMAX; ++i){
    seq_match_table[i]=0u;
    for (unsigned int k=0u; k<16u; k+=2u)
      seq_match_table[i]+=tmp[0x03u&(i>>k)];
  }
}

void initProtCompressionTable()
{
  uint32_t i;
  for (i=0; i<(1u<<8u); ++i)
    prot_comp_table[i]=0x00u;
  for (int j=0; j<aacount; ++j){
    prot_comp_table[aa[j]]=aacomp[j];
    if (aa[j]!='*')
      prot_comp_table[tolower(aa[j])]=aacomp[j];
  }

/*
  prot_comp_table['*']=0x09u;
  prot_comp_table['A']=prot_comp_table['a']=0x0Eu;
  prot_comp_table['N']=prot_comp_table['n']=0x10u;
  prot_comp_table['C']=prot_comp_table['c']=0x19u;
  prot_comp_table['P']=prot_comp_table['p']=0x0Fu;
  prot_comp_table['D']=prot_comp_table['d']=0x12u;
  prot_comp_table['Q']=prot_comp_table['q']=0x03u;
  prot_comp_table['E']=prot_comp_table['e']=0x02u;
  prot_comp_table['R']=prot_comp_table['r']=0x0Bu;
  prot_comp_table['F']=prot_comp_table['f']=0x15u;
  prot_comp_table['S']=prot_comp_table['s']=0x0Du;
  prot_comp_table['G']=prot_comp_table['g']=0x0Au;
  prot_comp_table['T']=prot_comp_table['t']=0x0Cu;
  prot_comp_table['H']=prot_comp_table['h']=0x13u;
  prot_comp_table['V']=prot_comp_table['v']=0x06u;
  prot_comp_table['I']=prot_comp_table['i']=0x04u;
  prot_comp_table['W']=prot_comp_table['w']=0x29u;
  prot_comp_table['K']=prot_comp_table['k']=0x00u;
  prot_comp_table['Y']=prot_comp_table['y']=0x11u;
  prot_comp_table['L']=prot_comp_table['l']=0x07u;
  prot_comp_table['M']=prot_comp_table['m']=0x24u;
*/
  for (int i=0; i<prot_blosum62_size; ++i)
    prot_blosum62[i]=0x00;

  efile f;
  f.open("data/BLOSUM62.mat","r");
  estrarray arr;
  estrarray head;
  while (!f.eof() && f.readarr(arr)){
    if (arr.size()==0 || arr[0][0]=='#') continue;
    if (head.size()==0) { head=arr; continue; }
    for (int i=1; i<arr.size(); ++i){
//      cout << arr[0][0] << " " << head[i][0] << " " << arr[i].i() << endl;
      prot_blosum62[prot_comp_table[arr[0][0]]*aacount+prot_comp_table[head[i][0]]]=arr[i].i();
    }
  }
}

void initProtTable()
{
  efile f;
  f.open("data/transl_table","r");

  int il=0;
  estr line[5];
  while (!f.eof() && il<5 && f.readln(line[il])) ++il;
  ldieif(il!=5,"not enough lines in prot_trans_table: "+estr(il));

  for (int i=0; i<codon2prot_size; ++i)
    codon2prot[i]=0;
  
  for (int i=0; i<line[0].len(); ++i){
    unsigned int pkmer=((seq_comp_table[line[4][i]]&0x3u)<<4u) | ((seq_comp_table[line[3][i]]&0x3u)<<2u) | seq_comp_table[line[2][i]]&0x3u;
    codon2prot[pkmer]=prot_comp_table[line[0][i]];
//    cout << line[2][i] << line[3][i] << line[4][i] << " " << line[0][i] << endl;
  }

  const unsigned long codonmask=((0x1u<<6u) - 1u);
  for (unsigned long i=0; i<PKMERMAX; ++i){
    unsigned long tmp=i;
    unsigned long res=0x0ul;
    kmer_prot_lt[i]=codon2prot[i&codonmask] | ((unsigned long)(codon2prot[(i>>6u)&codonmask])<<6u) | ((unsigned long)(codon2prot[(i>>12u)&codonmask])<<12u);
  }
  for (unsigned long i=0; i<PKMERMAX; ++i){
    unsigned long tmp=i;
    unsigned long res=0x0ul;
    for (int j=0; j<PKMERSIZE; ++j){
      res=(res<<2u)|(unsigned long)(revnuc[0x3ul&tmp]);
      tmp>>=2u;
    }
    kmer_protrev_lt[i]=kmer_prot_lt[res];
  }
}

void initCompressionTable()
{
  uint32_t i;
  for (i=0; i<IKMERMAX; ++i)
    seq_comp_table[i]=0x00u;

  for (i=0x0000u; i<=0xff00u; i+=0x0100u){
    for (int j=0; j<nuccount; ++j)
      seq_comp_table[i | nuc[j]]|=compnuc[j];
  }
  for (i=0; i<nuccount; ++i){
    for (uint32_t j=0x0u; j<=0xffu; ++j)
      seq_comp_table[(nuc[i]<<8u) | j]|=(compnuc[i]<<2u);
  }
}


void initCompressionTableStatus()
{
  uint32_t i;
  for (i=0; i<IKMERMAX; ++i)
    seq_comp_table[i]=0x10u; // left 4bits are status code 0x1X=error by default, right 4bits are compressed 2 nucleotides

  for (i=0x0000u; i<=0xff00u; i+=0x0100u){
    for (int j=0; j<nuccount; ++j)
      seq_comp_table[i | nuc[j]]|=compnuc[j];
  }
  for (i=0; i<nuccount; ++i){
    for (uint32_t j=0x0u; j<=0xffu; ++j)
      seq_comp_table[(nuc[i]<<8u) | j]|=(compnuc[i]<<2u);
  }
  for (int j=0; j<nuccount; ++j){
    seq_comp_table[(uint32_t(' ')<<8u) | nuc[j]]=0x20u|compnuc[j];
    seq_comp_table[(uint32_t('\t')<<8u) | nuc[j]]=0x20u|compnuc[j];
  }
  for (int i=0; i<nuccount; ++i){
    seq_comp_table[(nuc[i]<<8u)|' ']=0x20u|compnuc[i];
    seq_comp_table[(nuc[i]<<8u)|'\t']=0x20u|compnuc[i];
  }
  seq_comp_table[(' '<<8u)|' ']=0x30u; // 0x30: both missing nucleotides, skip both spaces
  seq_comp_table[(' '<<8u)|'\t']=0x30u; // 0x30: both missing nucleotides, skip both spaces
  seq_comp_table[('\t'<<8u)|' ']=0x30u; // 0x30: both missing nucleotides, skip both spaces
  seq_comp_table[('\t'<<8u)|'\t']=0x30u; // 0x30: both missing nucleotides, skip both spaces
  
  for (i=0; i<nuccount; ++i){
    for (int j=0; j<nuccount; ++j)
      seq_comp_table[(nuc[i]<<8u) | nuc[j]]=0x00|(compnuc[i]<<2u)|compnuc[j]; // when both characters are nucleotides reset the status code to 0x0X = noerror
  }
}
/*
void initProtTable()
{
  const unsigned long codonmask=((0x1u<<6u) - 1u);
  for (unsigned long i=0; i<PKMERMAX; ++i){
    unsigned long tmp=i;
    unsigned long res=0x0ul;
    kmer_prot_lt[i]=nuc2prot[i&codonmask] | ((unsigned long)(nuc2prot[(i>>6u)&codonmask])<<6u) | ((unsigned long)(nuc2prot[(i>>12u)&codonmask])<<12u);
  }
  for (unsigned long i=0; i<PKMERMAX; ++i){
    unsigned long tmp=i;
    unsigned long res=0x0ul;
    for (int j=0; j<PKMERSIZE; ++j){
      res=(res<<2u)|(unsigned long)(revnuc[0x3ul&tmp]);
      tmp>>=2u;
    }
    kmer_protrev_lt[i]=kmer_prot_lt[res];
  }
}
*/

void initRevComplementTable()
{
  for (unsigned long i=0; i<KMERMAX; ++i){
    unsigned long tmp=i;
    unsigned long res=0x0ul;
    for (int j=0; j<KMERSIZE; ++j){
      res=(res<<2u)|(unsigned long)(revnuc[0x3ul&tmp]);
      tmp>>=2u;
    }
    kmer_rev_lt[i]=res;
  }
  for (unsigned long i=0; i<KMERMAX2; ++i){
    unsigned long tmp=i;
    unsigned long res=0x0ul;
    for (int j=0; j<KMERSIZE2; ++j){
      res=(res<<2u)|(unsigned long)(revnuc[0x3ul&tmp]);
      tmp>>=2u;
    }
    kmer_rev_lt2[i]=res;
  }
}


void initSeqIdent()
{
  for (uint32_t i=0; i<IKMERMAX; ++i){
    for (uint32_t j=0; j<IKMERMAX; ++j){
      seq_ident_table[i^j]=0u;
      uint32_t ti=i; uint32_t tj=j;
      for (int k=0; k<IKMERSIZE; ++k,tj>>=2u,ti>>=2u){
        if ((0x03u&ti)==(0x03u&tj))
          ++seq_ident_table[i^j];
      } 
    }
  }
}

uint32_t seqatgc_count_lt[1u<<16u];
uint64_t seqatgc_conv_lt[1u<<8u];

void initATGCCount()
{
  uint32_t ti;
  uint32_t a=0u,t=0u,g=0u,c=0u;
  for (uint32_t i=0; i<(1u<<16u); ++i){
    a=0u; t=0u; g=0u; c=0u;
    ti=i;
    for (int k=0; k<8u; ++k,ti>>=2u){
      a>>=1u; t>>=1u; g>>=1u; c>>=1u;
      switch (0x03u&ti) {
        case 0x0u: a|=0x80u; break;
        case 0x1u: t|=0x80u; break;
        case 0x2u: g|=0x80u; break;
        case 0x3u: c|=0x80u; break;
      }
    }
    seqatgc_count_lt[i]=(a<<24u)|(t<<16u)|(g<<8u)|c;
  }
  uint64_t r;
  for (uint32_t i=0; i<(1u<<8u); ++i){
    ti=i;
    r=0ul;
    for (int k=0u; k<8u; ++k,ti>>=1u){
      r>>=8u;
      r|=(ti&0x1u?0x0100000000000000ul:0ul);
    }
    seqatgc_conv_lt[i]=r;
  }
}

void initAlignment()
{
  uint64_t id;
  for (uint32_t i=0; i<IKMERMAX; ++i){
    for (uint32_t j=0; j<IKMERMAX; ++j){
      uint32_t ti=i; uint32_t tj=j;
      id=0u;
      for (int k=0; k<IKMERSIZE; ++k,tj>>=2u,ti>>=2u){
        id>>=8u;
        if ((0x03u&ti)==(0x03u&tj))
          id|=0x0100000000000000ul;
        else
          id|=0x0200000000000000ul;
      } 
      seq_alignment_lt[i^j]=id;
    }
  }
}




int emain()
{
  char tmpstr[3+sizeof(long)*2u];
  efile f,f2;
  f.open("kmerseqtables-data.cpp","w");
  f2.open("kmerseqtables-data.h","w");
  f.write("#include <stdint.h>\n\n");

//  f2.write("extern const uint32_t aacount;\n\n");
  f2.write("const uint32_t aacount=");
  uint2hex(tmpstr,aacount);
  f2.write(tmpstr);
  f2.write(";\n\n");


  initMatchTable();
  f2.write("extern unsigned char seq_match_table[];\n\n");
  f.write("unsigned char seq_match_table[]={\n");
  uchr2hex(tmpstr,seq_match_table[0]);
  f.write(tmpstr);
  for (int i=1; i<(1u<<16u); ++i){
    f.write(",");
    uchr2hex(tmpstr,seq_match_table[i]);
    f.write(tmpstr);
    if (i%20==0) f.write("\n");
  }
  f.write("\n};\n\n");

  
  initCompressionTable();
  f2.write("extern unsigned char seq_comp_table[];\n\n");
  f.write("unsigned char seq_comp_table[]={\n");
  uchr2hex(tmpstr,seq_comp_table[0]);
  f.write(tmpstr);
  for (int i=1; i<(1u<<16u); ++i){
    f.write(",");
    uchr2hex(tmpstr,seq_comp_table[i]);
    f.write(tmpstr);
    if (i%20==0) f.write("\n");
  }
  f.write("\n};\n\n");

  initProtCompressionTable();
  f2.write("extern unsigned char prot_comp_table[];\n\n");
  f.write("unsigned char prot_comp_table[]={\n");
  uchr2hex(tmpstr,prot_comp_table[0]);
  f.write(tmpstr);
  for (int i=1; i<(1u<<8u); ++i){
    f.write(",");
    uchr2hex(tmpstr,prot_comp_table[i]);
    f.write(tmpstr);
    if (i%20==0) f.write("\n");
  }
  f.write("\n};\n\n");

  initSeqIdent();
  f2.write("extern unsigned char seq_ident_table[];\n\n");
  f.write("unsigned char seq_ident_table[]={\n");
  uchr2hex(tmpstr,seq_ident_table[0]);
  f.write(tmpstr);
  for (int i=1; i<IKMERMAX; ++i){
    f.write(",");
    uchr2hex(tmpstr,seq_ident_table[i]);
    f.write(tmpstr);
    if (i%20==0) f.write("\n");
  }
  f.write("\n};\n\n");

  initProtTable();
  f2.write("extern unsigned int kmer_prot_lt[];\n\n");
  f.write("unsigned int kmer_prot_lt[]={\n");
  int2hex(tmpstr,kmer_prot_lt[0]);
  f.write(tmpstr);
  for (int i=1; i<PKMERMAX; ++i){
    f.write(",");
    int2hex(tmpstr,kmer_prot_lt[i]);
    f.write(tmpstr);
    if (i%20==0) f.write("\n");
  }
  f.write("\n};\n\n");

  f2.write("extern unsigned int kmer_protrev_lt[];\n\n");
  f.write("unsigned int kmer_protrev_lt[]={\n");
  int2hex(tmpstr,kmer_protrev_lt[0]);
  f.write(tmpstr);
  for (int i=1; i<PKMERMAX; ++i){
    f.write(",");
    int2hex(tmpstr,kmer_protrev_lt[i]);
    f.write(tmpstr);
    if (i%20==0) f.write("\n");
  }
  f.write("\n};\n\n");


  initRevComplementTable();
  f2.write("extern unsigned int kmer_rev_lt[];\n\n");
  f.write("unsigned int kmer_rev_lt[]={\n");
  int2hex(tmpstr,kmer_rev_lt[0]);
  f.write(tmpstr);
  for (int i=1; i<KMERMAX; ++i){
    f.write(",");
    int2hex(tmpstr,kmer_rev_lt[i]);
    f.write(tmpstr);
    if (i%20==0) f.write("\n");
  }
  f.write("\n};\n\n");

  f2.write("extern unsigned int kmer_rev_lt2[];\n\n");
  f.write("unsigned int kmer_rev_lt2[]={\n");
  int2hex(tmpstr,kmer_rev_lt2[0]);
  f.write(tmpstr);
  for (int i=1; i<KMERMAX2; ++i){
    f.write(",");
    int2hex(tmpstr,kmer_rev_lt2[i]);
    f.write(tmpstr);
    if (i%20==0) f.write("\n");
  }
  f.write("\n};\n\n");

  initATGCCount();
  f2.write("extern uint32_t seqatgc_count_lt[];\n\n");
  f.write("uint32_t seqatgc_count_lt[]={\n");
  uint2hex(tmpstr,seqatgc_count_lt[0]);
  f.write(tmpstr);
  for (int i=1; i<(1u<<16u); ++i){
    f.write(",");
    uint2hex(tmpstr,seqatgc_count_lt[i]);
    f.write(tmpstr);
    if (i%20==0) f.write("\n");
  }
  f.write("\n};\n\n");

  f2.write("extern uint64_t seqatgc_conv_lt[];\n\n");
  f.write("uint64_t seqatgc_conv_lt[]={\n");
  ulong2hex(tmpstr,seqatgc_conv_lt[0]);
  f.write(tmpstr);
  for (int i=1; i<(1u<<8u); ++i){
    f.write(",");
    ulong2hex(tmpstr,seqatgc_conv_lt[i]);
    f.write(tmpstr);
    if (i%20==0) f.write("\n");
  }
  f.write("\n};\n\n");

  initAlignment();
  f2.write("extern uint64_t seq_alignment_lt[];\n\n");
  f.write("uint64_t seq_alignment_lt[]={\n");
  ulong2hex(tmpstr,seq_alignment_lt[0]);
  f.write(tmpstr);
  for (int i=1; i<IKMERMAX; ++i){
    f.write(",");
    ulong2hex(tmpstr,seq_alignment_lt[i]);
    f.write(tmpstr);
    if (i%20==0) f.write("\n");
  }
  f.write("\n};\n\n");

  f2.write("extern uint8_t codon2prot[];\n\n");
  f.write("uint8_t codon2prot[]={\n");
  uchr2hex(tmpstr,codon2prot[0]);
  f.write(tmpstr);
  for (int i=1; i<codon2prot_size; ++i){
    f.write(",");
    uchr2hex(tmpstr,codon2prot[i]);
    f.write(tmpstr);
  }
  f.write("\n};\n\n");

  f2.write("extern int8_t prot_blosum62[];\n\n");
  f.write("int8_t prot_blosum62[]={\n(int8_t)");
  char2int(tmpstr,prot_blosum62[0]);
  f.write(tmpstr);
  for (int i=1; i<prot_blosum62_size; ++i){
    f.write(",(int8_t)");
    char2int(tmpstr,prot_blosum62[i]);
    f.write(tmpstr);
  }
  f.write("\n};\n\n");
  f.close();

  f2.close();

  return(0);
}



