#ifndef EMSEQ_DEFS_H
#define EMSEQ_DEFS_H

#include <eutils/ebasicarray.h>
#include <eutils/eblockarray.h>


#include <deque>

using namespace std;

#include "mapseq-config.h"
#include "eseqali.h"
#include "kmerseqtables-data.h"

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


typedef ebasicarray<unsigned int> euintarray;
//typedef eblockarray<uint32_t> ekmerarray;  // leads to high vmem usage
typedef deque<uint32_t> ekmerarray;

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




const uint32_t nuc[]={'a','t','g','c','u'};
const uint32_t compnuc[]={0x0u,0x1u,0x2u,0x3u,0x1u};
const uint32_t consnuc[]={0x0001u,0x0010u,0x0100u,0x1000u};
const unsigned long safe_shift[32]={0x0ul,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful,0xfffffffffffffffful};

//const ealignscore as(1.0,2.0,10.0,1.0);
//const ealignscore as(2.0,1.0,5.0,1.0,30.0);
const ealignscore as(1.0,1.0,5.0,1.0,20.0);
const epalignscore pas(prot_blosum62,aacount,4.0,10.0,2.0,30.0);
//const ealignscore as(2.0,1.0,5.0,1.0,30.0);

const int NCOUNT_MAXLEN=1600; // IMPORTANT: must be disivible by 8 to be 64bit aligned

const uint32_t BMASK31=(1u<<31u)-1u;
const uint32_t BMASK16=(1u<<16u)-1u;
const uint32_t BMASK11=(1u<<11u)-1u;
const uint32_t BMASK10=(1u<<10u)-1u;

const int SEQSEGSIZE=2200;
//const int SEQSEGSIZE=1500;


#endif

