#include "eseqdb.h"

#include <eutils/etimer.h>
#include <eutils/efile.h>
#include <eutils/estr.h>
#include <eutils/eheap.h>
#include <deque>
#include <map>

using namespace std;

float swmin=0.05;
float swmax=0.3;
int minid1=1;
int minid2=1;
//int topotus=10;
//int tophits=20;
//int minscore=30;
//int otulim=50;

float sweight=30.0;
float sweightabs=0.1;
float cfthres=0.5;

int minlen=75;
float minqual=0.05;


// unsigned int akmers[MAXSIZE]; // 6bases

emtdata mtdata;



estr outfmt_simple(const etax& tax,const earrayof<double,int>& ptax,const efloatarray& mcfarr)
{
  estr res;
  for (int l=0; l<ptax.size() && l<tax.names.size(); ++l){
    if (mcfarr[l]<cfthres) break;
    res+=(ptax.keys(l)==-1?estr("N/A"):tax.names[l].at(ptax.keys(l)))+";";
  }
  res.del(-1);
  return(res);
}

estr outfmt_confidences(const etax& tax,const earrayof<double,int>& ptax,const efloatarray& mcfarr)
{
  estr res;
  for (int l=0; l<ptax.size() && l<tax.names.size(); ++l)
    res+=(ptax.keys(l)==-1?estr("N/A"):tax.names[l].at(ptax.keys(l)))+"\t"+mcfarr[l]+"\t"+ptax.values(l)+"\t";
  res.del(-1);
  return(res);
}



void taskSearchPaired()
{
  estr outstr;
  esearchws searchws;
  searchws.initPaired(*mtdata.seqdb);
/*
//  searchws.seqkmers.init(MAXSIZE,-1);
//  searchws.revseqkmers.init(MAXSIZE,-1);
  searchws.kmerbitmask=new uint64_t[MAXSIZE/64+1];
  searchws.bitmask=new uint64_t[(mtdata.seqdb->otus.size()*2)/64+1];
  searchws.otukmerpos.init(mtdata.seqdb->otus.size()*2,0);
  searchws.idcount.init(mtdata.seqdb->otus.size()*2,0);
  searchws.idcount2.init(mtdata.seqdb->otus.size()*2,0);
  searchws.kmerpos2.init(MAXSIZE,0u);

  searchws.kmerpos.init(MAXSIZE,0u);
  searchws.kmerposrev.init(MAXSIZE,0u);
  searchws.kmerpos3.init(MAXSIZE,0u);
  searchws.kmerposrev3.init(MAXSIZE,0u);
  searchws.offset=1u;
  searchws.offset2=1u;
  searchws.offset3=1u;
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

    for (int i=0; i+1<pbuf->size(); i+=2){
      eseq& s(pbuf->values(i));
      eseq& s2(pbuf->values(i+1));

      earray<epredinfo> pinfoarr;
  //    ebasicarray<ealigndata> matchcounts;
      mtdata.seqdb->seqsearchpair(pbuf->keys(i),s,s2,pinfoarr,searchws);
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
  
//        outstr+=pbuf->keys(i)+"\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+"\t"+pinfo.tophit.score()+"\t"+bid+"\t"+pinfo.tophit.matches+"\t"+pinfo.tophit.mismatches+"\t"+pinfo.tophit.gaps+"\t"+(s.seqstart+pinfo.tophit.s1)+"\t"+(s.seqstart+pinfo.tophit.e1)+"\t"+pinfo.tophit.s2+"\t"+pinfo.tophit.e2+"\t"+taxcutoffmin+"\t"+pinfo.matchcounts.size()+"\t";
        outstr+=pbuf->keys(i)+"\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+"\t"+pinfo.tophit.score()+"\t"+pinfo.tophit.identity()+"\t"+pinfo.tophit.matches+"\t"+pinfo.tophit.mismatches+"\t"+pinfo.tophit.gaps+"\t"+(s.seqstart+pinfo.tophit.s1)+"\t"+(s.seqstart+pinfo.tophit.e1)+"\t"+pinfo.tophit.s2+"\t"+pinfo.tophit.e2+"\t"+(pinfo.tophit.revcompl?"-":"+")+"\t";
       
        if (mtdata.seqdb->taxa.size()==0){
          double topscore=pinfo.tophit.score();
          double tscore=0.0;
          for (int t=0; t<pinfo.matchcounts.size(); ++t)
//            tscore+=exp((1.0l-topscore/pinfo.matchcounts[t].score())*sweight);
            tscore+=exp((1.0l-topscore/pinfo.matchcounts[t].score())*sweight);
          outstr+="\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+"\t"+estr(1.0/tscore);
        }

        for (int t=0; t<mtdata.seqdb->taxa.size(); ++t){
          etax& tax(mtdata.seqdb->taxa.at(t));
      
          earrayof<double,int> ptax;
          edoublearray taxscores;
          efloatarray mcfarr;
          taxScoreSum(taxscores,pinfo,tax,searchws.taxcounts,s.seqlen);
          taxScore(ptax,mcfarr,pinfo.tophit,pinfo,taxscores,tax,s.seqlen);
      
          outstr+="\t";
          outstr+=(*mtdata.outfmt)(tax,ptax,mcfarr);
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
          estr salistr=pinfo.tophit.compress(s);
          outstr+=salistr+ "\n";
          outstr+=sali_decompress(salistr,mtdata.seqdb->seqs.values(pinfo.tophit.seqid))+"\n";
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
    }

    mtdata.m.lock();
    cout << outstr;
    mtdata.sbuffer.add(pbuf);
    mtdata.sbufferSignal.signal();
    mtdata.m.unlock();
  }
//  delete searchws.bitmask;
}

void taskSearch()
{
  estr outstr;
  esearchws searchws(*mtdata.seqdb);

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
      if (mtdata.galign)
        mtdata.seqdb->seqsearch_global(pbuf->keys(i),s,pinfoarr,searchws);
      else
        mtdata.seqdb->seqsearch(pbuf->keys(i),s,pinfoarr,searchws);

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
      outstr+=pbuf->keys(i)+(mtdata.print_hits?"\t0":"")+"\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+(mtdata.print_hits?estr("\t")+pinfo.tophit.kmercount:estr())+"\t"+pinfo.tophit.score()+"\t"+pinfo.tophit.identity()+"\t"+pinfo.tophit.matches+"\t"+pinfo.tophit.mismatches+"\t"+pinfo.tophit.gaps+"\t"+(s.seqstart+pinfo.tophit.s1)+"\t"+(s.seqstart+pinfo.tophit.e1)+"\t"+pinfo.tophit.s2+"\t"+pinfo.tophit.e2+"\t"+(pinfo.tophit.revcompl?"-":"+")+"\t";
     
      if (mtdata.seqdb->taxa.size()==0){
        double topscore=pinfo.tophit.score();
        double tscore=0.0;
        for (int t=0; t<pinfo.matchcounts.size(); ++t)
//          tscore+=exp((1.0l-topscore/pinfo.matchcounts[t].score())*sweight);
          tscore+=exp((pinfo.matchcounts[t].score()-topscore)*sweightabs);
        outstr+="\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+"\t"+estr(1.0/tscore);
      }

      earray<edoublearray> taxscores;
      taxscores.init(mtdata.seqdb->taxa.size());
      for (int t=0; t<mtdata.seqdb->taxa.size(); ++t){
        etax& tax(mtdata.seqdb->taxa.at(t));
    
        efloatarray mcfarr;
        earrayof<double,int> ptax;
        efloatarray tmpmcfarr;
        earrayof<double,int> tmptax;
        taxScoreSumE(taxscores[t],pinfo,tax,searchws.taxcounts,s.seqlen);
//        taxScoreSum(taxscores[t],pinfo,tax,searchws.taxcounts,s.seqlen);
        earrayof<double,int> mixedtax; // best mixed prediction over all taxonomies
        efloatarray mixedmcfarr;

        taxScoreE(ptax,mcfarr,pinfo.tophit,pinfo,taxscores[t],tax,s.seqlen);
//        taxScore(ptax,mcfarr,pinfo.tophit,pinfo,taxscores[t],tax,s.seqlen);
        mixedtax=ptax;
        mixedmcfarr=mcfarr;
        for (int l=pinfo.matchcounts.size()-2; l>=0; --l){
          ealigndata& adata(pinfo.matchcounts[l]);
          taxScoreE(tmptax,tmpmcfarr,adata,pinfo,taxscores[t],tax,s.seqlen);
//          taxScore(tmptax,tmpmcfarr,adata,pinfo,taxscores[t],tax,s.seqlen);
          for (int tl=0; tl<tmptax.size(); ++tl){ // choose top hit in list which is not always best aligned (when including evidence or indirect taxonomy)
            if (tmpmcfarr[tl]>mixedmcfarr[tl]){ // should check per level or just for top?
//            if (tmptax.values(tl)>ptax.values(tl)){ // previous code selecting only on alignment score, but id cutoff is important too
              ldieif(mixedtax.size()!=tmptax.size(),"size mismatch: "+estr(mixedtax.size())+" "+tmptax.size());
              mixedtax.values(tl)=tmptax.values(tl); 
              mixedtax.keys(tl)=tmptax.keys(tl); 
              mixedmcfarr[tl]=tmpmcfarr[tl];
            }
            if (tmpmcfarr[tl]>mcfarr[tl]){ // should check per level or just for top?
//            if (tmptax.values(tl)>ptax.values(tl)){ // previous code selecting only on alignment score, but id cutoff is important too
              ptax=tmptax; 
              mcfarr=tmpmcfarr;
            }
          }
        }

        outstr+="\t";
//        outstr+=(*mtdata.outfmt)(tax,ptax,mcfarr);
        outstr+=(*mtdata.outfmt)(tax,mixedtax,mixedmcfarr);
        outstr+="\t";
      }
      outstr+="\n";
      if (mtdata.print_hits){
//        etax& tax(mtdata.seqdb->taxa.at(0));
        for (int l=pinfo.matchcounts.size()-2; l>=0; --l){
          ealigndata& adata(pinfo.matchcounts[l]);
          outstr+=pbuf->keys(i)+"\t"+(pinfo.matchcounts.size()-l-1)+"\t"+mtdata.seqdb->seqs.keys(adata.seqid)+"\t"+adata.kmercount+"\t"+adata.score()+"\t"+adata.identity()+"\t"+adata.matches+"\t"+adata.mismatches+"\t"+adata.gaps+"\t"+(s.seqstart+adata.s1)+"\t"+(s.seqstart+adata.e1)+"\t"+adata.s2+"\t"+adata.e2+"\t"+(adata.revcompl?"-":"+")+"\t";
          for (int t=0; t<mtdata.seqdb->taxa.size(); ++t){
            etax& tax(mtdata.seqdb->taxa.at(t));
        
            efloatarray mcfarr;
            earrayof<double,int> ptax;
            taxScoreE(ptax,mcfarr,adata,pinfo,taxscores[t],tax,s.seqlen);
//            taxScore(ptax,mcfarr,adata,pinfo,taxscores[t],tax,s.seqlen);
            outstr+="\t";
            outstr+=(*mtdata.outfmt)(tax,ptax,mcfarr);
            outstr+="\t";
          }
          outstr+="\n";
        }
      }
      if (mtdata.print_kmerhits){
        ldieif(searchws.idcount.size()!=searchws.seqids.size(),"mismatch idcount and seqids: "+estr(searchws.idcount.size())+" "+searchws.seqids.size());
        for (int l=0; l<searchws.seqids.size(); ++l)
          outstr+="#"+pbuf->keys(i)+"\t"+searchws.idcount[l]+"\t"+mtdata.seqdb->seqs.keys(searchws.seqids[l]%mtdata.seqdb->seqs.size())+"\n";
      }
//      outstr+=pinfo.tophit.profile.str() + "\n";
      if (mtdata.print_align){
        pinfo.tophit.profile.inv();
        outstr+=pinfo.tophit.profile.str() + "\n";
        estr salistr=pinfo.tophit.compress(s);
        outstr+=salistr+ "\n";
        outstr+=sali_decompress(salistr,mtdata.seqdb->seqs.values(pinfo.tophit.seqid))+"\n";
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
//  delete searchws.bitmask;
}


void taskSearchReturn()
{
  estr outstr;
  esearchws searchws(*mtdata.seqdb);

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
      if (mtdata.galign)
        mtdata.seqdb->seqsearch_global(pbuf->keys(i),s,pinfoarr,searchws);
      else
        mtdata.seqdb->seqsearch(pbuf->keys(i),s,pinfoarr,searchws);

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
      outstr+=pbuf->keys(i)+(mtdata.print_hits?"\t0":"")+"\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+"\t"+pinfo.tophit.score()+"\t"+pinfo.tophit.identity()+"\t"+pinfo.tophit.matches+"\t"+pinfo.tophit.mismatches+"\t"+pinfo.tophit.gaps+"\t"+(s.seqstart+pinfo.tophit.s1)+"\t"+(s.seqstart+pinfo.tophit.e1)+"\t"+pinfo.tophit.s2+"\t"+pinfo.tophit.e2+"\t"+(pinfo.tophit.revcompl?"-":"+")+"\t";
     
      if (mtdata.seqdb->taxa.size()==0){
        double topscore=pinfo.tophit.score();
        double tscore=0.0;
        for (int t=0; t<pinfo.matchcounts.size(); ++t)
//          tscore+=exp((1.0l-topscore/pinfo.matchcounts[t].score())*sweight);
          tscore+=exp((pinfo.matchcounts[t].score()-topscore)*sweightabs);
        outstr+="\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+"\t"+estr(1.0/tscore);
      }

      earray<edoublearray> taxscores;
      taxscores.init(mtdata.seqdb->taxa.size());
      for (int t=0; t<mtdata.seqdb->taxa.size(); ++t){
        etax& tax(mtdata.seqdb->taxa.at(t));
    
        efloatarray mcfarr;
        earrayof<double,int> ptax;
        efloatarray tmpmcfarr;
        earrayof<double,int> tmptax;
        taxScoreSum(taxscores[t],pinfo,tax,searchws.taxcounts,s.seqlen);

        taxScore(ptax,mcfarr,pinfo.tophit,pinfo,taxscores[t],tax,s.seqlen);
        for (int l=pinfo.matchcounts.size()-2; l>=0; --l){
          ealigndata& adata(pinfo.matchcounts[l]);
          taxScore(tmptax,tmpmcfarr,adata,pinfo,taxscores[t],tax,s.seqlen);
          for (int tl=0; tl<tmptax.size(); ++tl){
            if (tmptax.values(tl)>ptax.values(tl)){
              ptax=tmptax; 
              mcfarr=tmpmcfarr;
            }
          }
        }

        outstr+="\t";
        outstr+=(*mtdata.outfmt)(tax,ptax,mcfarr);
        outstr+="\t";
      }
      outstr+="\n";
      if (mtdata.print_hits){
        etax& tax(mtdata.seqdb->taxa.at(0));
        for (int l=pinfo.matchcounts.size()-2; l>=0; --l){
          ealigndata& adata(pinfo.matchcounts[l]);
          outstr+=pbuf->keys(i)+"\t"+(pinfo.matchcounts.size()-l-1)+"\t"+mtdata.seqdb->seqs.keys(adata.seqid)+"\t"+adata.score()+"\t"+adata.identity()+"\t"+adata.matches+"\t"+adata.mismatches+"\t"+adata.gaps+"\t"+(s.seqstart+adata.s1)+"\t"+(s.seqstart+adata.e1)+"\t"+adata.s2+"\t"+adata.e2+"\t"+(adata.revcompl?"-":"+")+"\t";
          for (int t=0; t<mtdata.seqdb->taxa.size(); ++t){
            etax& tax(mtdata.seqdb->taxa.at(t));
        
            efloatarray mcfarr;
            earrayof<double,int> ptax;
            taxScore(ptax,mcfarr,adata,pinfo,taxscores[t],tax,s.seqlen);
            outstr+="\t";
            outstr+=(*mtdata.outfmt)(tax,ptax,mcfarr);
            outstr+="\t";
          }
          outstr+="\n";
        }
      }
//      outstr+=pinfo.tophit.profile.str() + "\n";
      if (mtdata.print_align){
        pinfo.tophit.profile.inv();
        outstr+=pinfo.tophit.profile.str() + "\n";
        estr salistr=pinfo.tophit.compress(s);
        outstr+=salistr+ "\n";
        outstr+=sali_decompress(salistr,mtdata.seqdb->seqs.values(pinfo.tophit.seqid))+"\n";
        outstr+=pinfo.tophit.align_str(s,mtdata.seqdb->seqs.values(pinfo.tophit.seqid));
        outstr+="\n";
      }
      mtdata.m.lock();
      mtdata.output.add(outstr); outstr.clear();
//      cout << outstr; outstr.clear();
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
    mtdata.output.add(outstr); outstr.clear();
    mtdata.sbuffer.add(pbuf);
    mtdata.sbufferSignal.signal();
    mtdata.m.unlock();
  }
//  delete searchws.bitmask;
}



void taskProtSearch()
{
  estr outstr;
  esearchws searchws;

  searchws.initProt(*mtdata.seqdb);
/*
//  searchws.seqkmers.init(MAXSIZE,-1);
//  searchws.revseqkmers.init(MAXSIZE,-1);

  searchws.kmerbitmask=new uint64_t[MAXSIZE/64+1];
  searchws.bitmask=new uint64_t[(mtdata.seqdb->otus.size()*2)/64+1];
  searchws.otukmerpos.init(mtdata.seqdb->otus.size()*2,0);
  searchws.idcount.init(mtdata.seqdb->otus.size()*2,0);
  searchws.idcount2.init(mtdata.seqdb->otus.size()*2,0);
  searchws.kmerpos.init(MAXSIZE,0u);
  searchws.kmerpos2.init(PMAXSIZE,0u);
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
  //    ebasicarray<ealigndata> matchcounts;
      mtdata.seqdb->pseqsearch(pbuf->keys(i),s,pinfoarr,searchws);
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
      outstr+=pbuf->keys(i)+(mtdata.print_hits?"\t0":"")+"\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+"\t"+pinfo.tophit.score()+"\t"+pinfo.tophit.identity()+"\t"+pinfo.tophit.matches+"\t"+pinfo.tophit.mismatches+"\t"+pinfo.tophit.gaps+"\t"+(s.seqstart+pinfo.tophit.s1)+"\t"+(s.seqstart+pinfo.tophit.e1)+"\t"+pinfo.tophit.s2+"\t"+pinfo.tophit.e2+"\t"+(pinfo.tophit.revcompl?"-":"+")+"\t";
     
      if (mtdata.seqdb->taxa.size()==0){
        double topscore=pinfo.tophit.score();
        double tscore=0.0;
        for (int t=0; t<pinfo.matchcounts.size(); ++t)
//          tscore+=exp((1.0l-topscore/pinfo.matchcounts[t].score())*sweight);
          tscore+=exp((pinfo.matchcounts[t].score()-topscore)*sweightabs);
        outstr+="\t"+mtdata.seqdb->seqs.keys(pinfo.tophit.seqid)+"\t"+estr(1.0/tscore);
      }

      earray<edoublearray> taxscores;
      taxscores.init(mtdata.seqdb->taxa.size());
      for (int t=0; t<mtdata.seqdb->taxa.size(); ++t){
        etax& tax(mtdata.seqdb->taxa.at(t));
    
        efloatarray mcfarr;
        earrayof<double,int> ptax;
        efloatarray tmpmcfarr;
        earrayof<double,int> tmptax;
        taxScoreSum(taxscores[t],pinfo,tax,searchws.taxcounts,s.seqlen);

        taxScore(ptax,mcfarr,pinfo.tophit,pinfo,taxscores[t],tax,s.seqlen);
        for (int l=pinfo.matchcounts.size()-2; l>=0; --l){
          ealigndata& adata(pinfo.matchcounts[l]);
          taxScore(tmptax,tmpmcfarr,adata,pinfo,taxscores[t],tax,s.seqlen);
          for (int tl=0; tl<tmptax.size(); ++tl){
            if (tmptax.values(tl)>ptax.values(tl)){
              ptax=tmptax; 
              mcfarr=tmpmcfarr;
            }
          }
        }

        outstr+="\t";
        outstr+=(*mtdata.outfmt)(tax,ptax,mcfarr);
        outstr+="\t";
      }
      outstr+="\n";
      if (mtdata.print_hits){
        etax& tax(mtdata.seqdb->taxa.at(0));
        for (int l=pinfo.matchcounts.size()-2; l>=0; --l){
          ealigndata& adata(pinfo.matchcounts[l]);
          outstr+=pbuf->keys(i)+"\t"+(pinfo.matchcounts.size()-l-1)+"\t"+mtdata.seqdb->seqs.keys(adata.seqid)+"\t"+adata.score()+"\t"+adata.identity()+"\t"+adata.matches+"\t"+adata.mismatches+"\t"+adata.gaps+"\t"+(s.seqstart+adata.s1)+"\t"+(s.seqstart+adata.e1)+"\t"+adata.s2+"\t"+adata.e2+"\t"+(adata.revcompl?"-":"+")+"\t";
          for (int t=0; t<mtdata.seqdb->taxa.size(); ++t){
            etax& tax(mtdata.seqdb->taxa.at(t));
        
            efloatarray mcfarr;
            earrayof<double,int> ptax;
            taxScore(ptax,mcfarr,adata,pinfo,taxscores[t],tax,s.seqlen);
            outstr+="\t";
            outstr+=(*mtdata.outfmt)(tax,ptax,mcfarr);
            outstr+="\t";
          }
          outstr+="\n";
        }
      }
//      outstr+=pinfo.tophit.profile.str() + "\n";
      if (mtdata.print_align){
        pinfo.tophit.profile.inv();
        outstr+=pinfo.tophit.profile.str() + "\n";
        outstr+=pinfo.tophit.palign_str(s,mtdata.seqdb->seqs.values(pinfo.tophit.seqid));
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
//  delete searchws.bitmask;
}

char qlt[256];

unsigned char lt[256];
unsigned char lt16[1u<<16u];

void eseqdb::printSearchHeader()
{
  cout << "#query\tdbhit\tbitscore\tidentity\tmatches\tmismatches\tgaps\tquery_start\tquery_end\tdbhit_start\tdbhit_end\tstrand\t";
  if (mtdata.outfmt==&outfmt_simple){
    for (int i=0; i<taxa.size(); ++i)
      cout << "\t" << taxa[i].name << "\t";
  }else if (mtdata.outfmt==&outfmt_confidences){
    for (int i=0; i<taxa.size(); ++i){
      etax& tax(taxa[i]);
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
  }else{
    ldie("unknown output format chosen");
  }
  cout << endl;
}


void fasta_filter(const estr& id,estr& seq){
  if (id.len()==0) return;
  int lastbad=0;
  int bad=0;
  int e=seq.len();
  for (int j=0; j<seq.len(); ++j){
//    ldieif(qlt[qual[j]]==-1,"unknown quality char: "+id+" '"+estr(qual[j])+"'");
    if (seq[j]=='N'){
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
 
  if (float(bad)/e<=minqual && e>=minlen){
//    cout << '>' << id << endl;
    if (e<seq.len())
      seq.del(e);
//    cout << (e==seq.len()?seq:seq.substr(0,e)) << endl;
  }else
    seq.clear();
}


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


int eseqdb::processQueryFASTQ(const estr& fname,void (*taskfunc)(),ethreads& t){
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
//  ethreads t;
  emutex m;
  
  ldieif(fname=="-","reading from stdin not supported with paired end data");
 
  f.open(fname,"r");

  mtdata.seqdb=this;
  mtdata.finished=false;
  mtdata.output.clear();

  t1.reset();
  t.run(taskfunc,evararray());

  mtdata.sbuffer.reserve(t.threads.size()*2);
  for (int i=0; i<t.threads.size()*2; ++i){
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
  return(0);
}


int eseqdb::processQueryFASTA(const estr& fname,void (*taskfunc)(),ethreads& t)
{
  etimer t1;
  t1.reset();

  egzfile f;
  estr line;
  estrarray args;

  t1.reset();
//  ethreads t;
  emutex m;
 
  if (fname=="-")
    f.open(stdin,"r");
  else
    f.open(fname,"r");

  mtdata.seqdb=this;
  mtdata.finished=false;
  mtdata.output.clear();

  t1.reset();
  t.run(taskfunc,evararray());

  mtdata.sbuffer.reserve(t.threads.size()*2);
  for (int i=0; i<t.threads.size()*2; ++i){
    estrarrayof<eseq> *sarr=new estrarrayof<eseq>;
    for (int j=0; j<100; ++j) sarr->add(estr(),eseq());
    mtdata.sbuffer.add(sarr);
  }

  int reterr=0;

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
    if(str2id.len()==0 || str2id[0]!='>') {
      lerror("Unexpected line: "+str2id);
      reterr=-1;
      break;
    }
    str2id.del(0,1);
    int i=str2id.findchr(" \t");
    if (i!=-1l) str2id.del(i); // only keep id up to first white space

    str2seq.clear();
    //TODO: read a limited number of nucleotides each time, so the code does not depend on line breaks. For every nucleotide "chunk", compress the nucleotides and add to sequence
    while (f.readln(line) && line.len() && line[0]!='>'){
      if (str2seq.len()+line.len()>MAXSEQLEN){
        cbuf->keys(cbufind)=str2id;
        eseq& s(cbuf->values(cbufind));
        int slen=str2seq.len();
        fasta_filter(str2id,str2seq);
        if (str2seq.len()==0) { seqstart+=slen; continue; }
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
    fasta_filter(str2id,str2seq);
    if (str2seq.len()==0) { ++seqcount; continue; }
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
  return(reterr);
}


int eseqdb::processQueryPairend(const estr& fname,const estr& fname2,void (*taskfunc)(),ethreads& t)
{
  etimer t1;
  t1.reset();

  egzfile f,f2;
  estr line,line2;
  estrarray args;

  t1.reset();
  emutex m;
  
  int reterr=0;

  ldieif(fname=="-","reading from stdin not supported with paired end data");
 
  f.open(fname,"r");
  f2.open(fname2,"r");

  mtdata.seqdb=this;
  mtdata.finished=false;
  mtdata.output.clear();

  t1.reset();
  t.run(taskfunc,evararray());

  mtdata.sbuffer.reserve(t.threads.size()*2);
  for (int i=0; i<t.threads.size()*2; ++i){
    estrarrayof<eseq> *sarr=new estrarrayof<eseq>;
    for (int j=0; j<100; ++j) sarr->add(estr(),eseq());
    mtdata.sbuffer.add(sarr);
  }

//  cerr << "# processing input... ";
  long seqcount=0;
  fprintf(stderr,"# processing input... ");
  fflush(stderr);
  lerrorif(!f.readln(line),"unable to read query file");
  lerrorif(!f2.readln(line2),"unable to read query file2");
  estrarrayof<eseq> *cbuf=0x00;

  int cbufind=0;
  const int MAXSEQLEN=100000;
  long seqstart,seqstart2;
  estr seq,id;
  estr seq2,id2;
  seq.reserve(MAXSEQLEN);
  seq2.reserve(MAXSEQLEN);
  while (!f.eof() && !f2.eof()){
    if (cbuf==0x00){
      mtdata.m.lock();
      while (mtdata.sbuffer.size()==0) mtdata.sbufferSignal.wait(mtdata.m);
      cbuf=mtdata.sbuffer[mtdata.sbuffer.size()-1];
      mtdata.sbuffer.erase(mtdata.sbuffer.size()-1);
      mtdata.m.unlock();
      cbufind=0;
    }

    seqstart=0;
    seqstart2=0;
    id=line; id2=line2;

    if (id.len()==0 || id[0]!='>'){
      lerror("Unexpected line: "+id);
      reterr=-1;
      break;
    }
    if (id2.len()==0 || id2[0]!='>'){
      lerror("Unexpected line: "+id2);
      reterr=-1;
      break;
    }

    id.del(0,1);
    id2.del(0,1);

    int i=id.findchr(" \t"),i2=id2.findchr(" \t");
    if (i!=-1l) id.del(i); // only keep id up to first white space
    if (i2!=-1l) id2.del(i2); // only keep id up to first white space

    ldieif(i!=i2,"paired end read ids do not match or reads are not in sync: "+id+" != "+id2);

    seq.clear();
    seq2.clear();
    //TODO: read a limited number of nucleotides each time, so the code does not depend on line breaks. For every nucleotide "chunk", compress the nucleotides and add to sequence
    while (f.readln(line) && line.len() && line[0]!='>'){
      ldieif(seq.len()+line.len()>MAXSEQLEN,"paired end reads larger than "+estr(MAXSEQLEN)+" not supported");
      seq+=line;
    }
    while (f2.readln(line2) && line2.len() && line2[0]!='>'){
      ldieif(seq2.len()+line2.len()>MAXSEQLEN,"paired end reads larger than "+estr(MAXSEQLEN)+" not supported");
      seq2+=line2;
    }

    cbuf->keys(cbufind)=id;
    eseq& s(cbuf->values(cbufind));
    fasta_filter(id,seq);
    s.setseq(seq);
    s.seqstart=seqstart;
    ++cbufind;

    cbuf->keys(cbufind)=id2;
    eseq& s2(cbuf->values(cbufind));
    fasta_filter(id2,seq2);
    s2.setseq(seq2);
    s2.seqstart=seqstart2;
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
      fprintf(stderr,"\r# processing input... paired end reads: %li",seqcount);
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
  return(reterr);
}





void randomize(ernd& prnd,eintarray& seqids)
{
  for (int i=seqids.size()-1; i>0; --i){
    int r=int(prnd.uniform()*(i+1));
    if (r!=i) seqids.swap(r,i);
  }
}


void calcConfidence(earrayof<double,int>& ptax,efloatarray &mcfarr,const ealigndata& adata,const epredinfo& pinfo,const etax& tax)
{
  float bid=adata.identity();
  if (tax.seqs[pinfo.tophit.seqid]!=0x00){
    eseqtax &tmptaxhit(*tax.seqs[adata.seqid]);
    // adjust id to closest gold hit
    if (tmptaxhit.bid>0.0) {
      bid=bid*tmptaxhit.bid;
      for (int l=0; l<ptax.size(); ++l)
//        ptax.values(l)=ptax.values(l)*tmptaxhit.tl[l].cf;
        ptax.values(l)=MIN(ptax.values(l),tmptaxhit.tl[l].cf);
    }
  }

  float lastmcf=0.0;
  //  adjust computed confidences
  for (int l=MIN(ptax.size(),tax.names.size())-1; l>=0; --l){
    ldieif(ptax.keys(l)!=-1 && ptax.keys(l)>=tax.names[l].size(),"key out of tax: "+estr(l)+" "+ptax.keys(l)+" "+tax.names[l].size());
    mcfarr[l]=ptax.values(l);

    if (tax.cutoff.size()>0){ // if only fixed id threshold exists
      float ncf=(bid-tax.cutoff[l]+0.02)/tax.cutoffcoef[l];
      mcfarr[l]=ncf<ptax.values(l)?ncf:ptax.values(l);
    }

    if (mcfarr[l]>1.0) mcfarr[l]=1.0; else if (mcfarr[l]<0.0) mcfarr[l]=0.0;
    if (mcfarr[l]<lastmcf) mcfarr[l]=lastmcf;  // do not let confidences get smaller
    lastmcf=mcfarr[l];
  }
}



void taxScoreE(earrayof<double,int>& ptax,efloatarray& mcfarr,ealigndata& adata,epredinfo& pinfo,edoublearray& taxscores,etax& tax,int slen)
{
  ptax.clear();
//  while (ptax.size()<tax.names.size())
//    ptax.add(-1,0.0);
  for (int i=0; i<tax.names.size(); ++i)
    ptax.add(-1,0.0);


//  float sw=sweightabs; //MAX(0.025,MIN(0.3,(0.3-MAX(0.0,(log(slen)-log(80.0)))*(0.3-0.05)/(log(500.0)-log(80.0)))));
//  float sw=MAX(0.025,MIN(0.3,(0.3-MAX(0.0,(log(slen)-log(80.0)))*(0.3-0.05)/(log(500.0)-log(80.0)))));

  float sw=MAX(swmin,MIN(swmax,(swmax-MAX(0.0,(log(slen)-log(80.0)))*(swmax-swmin)/(log(500.0)-log(80.0)))));

  mcfarr.init(ptax.size(),0.0);
  if(tax.seqs[adata.seqid]==0x00)
    return;

//  ldieif(tax.seqs[adata.seqid]==0x00,"Missing taxonomy for sequence");
  float cfw;
  eseqtax &seqtax(*tax.seqs[adata.seqid]);
  for (int l=0; l<seqtax.tl.size(); ++l){
    ldieif(seqtax.tl[l].tid>=long(tax.names[l].size()),estr("key out of tax: ")+adata.seqid+" "+estr(l)+" "+seqtax.tl[l].tid+" "+tax.names[l].size());
    cfw=1.0;
    if (seqtax.bid>0.0)
      cfw=0.5*seqtax.tl[l].cf+0.5;
    ptax.keys(l)=seqtax.tl[l].tid;
//    ptax.values(l)=exp((1.0l-pinfo.tophit.score()/adata.score())*sweight)*cfw/taxscores[l];
//    ptax.values(l)=exp((adata.score()-pinfo.tophit.score())*sweightabs)*cfw/taxscores[l];
    if (ptax.keys(l)==-1)
      ptax.values(l)=0.0;
    else
      ptax.values(l)=exp((adata.score()-pinfo.tophit.score())*sw)*cfw*seqtax.tl[l].evidence/taxscores[l];
  }
  calcConfidence(ptax,mcfarr,adata,pinfo,tax);
}

void taxScoreSumE(edoublearray& taxscores,epredinfo& pinfo,etax& tax,ebasicarray<eintarray>& taxcounts,int slen)
{
  taxscores.init(tax.names.size(),0.0);
  int tophitl=-1;

  for (int l=pinfo.matchcounts.size()-1; l>=0; --l){
    int sbest=pinfo.matchcounts[l].seqid;
    if (tax.seqs[sbest]!=0x00) { tophitl=l; break; }
  }
  if (tophitl==-1) return;

  ealigndata &tophit(pinfo.matchcounts[tophitl]);
  pinfo.tophit=tophit; 
  ldieif(tax.seqs[tophit.seqid]==0x00,"top hit does not have taxonomy: "+estr(tophit.seqid));

  double topscore=tophit.score();

  int taxid=0;
  ebasicarray<efloatarray> taxmaxscores;
  taxmaxscores.init(tax.names.size());
  for (int i=0; i<taxmaxscores.size(); ++i)
    taxmaxscores[i].init(tax.names[i].size(),0.0);
//  taxcounts.init(tax.names.size());
//  for (int i=0; i<taxcounts.size(); ++i)
//    taxcounts[i].init(tax.names[i].size(),-1);

//  float sw=MAX(0.025,MIN(0.3,(0.3-MAX(0.0,(log(slen)-log(80.0)))*(0.3-0.05)/(log(500.0)-log(80.0)))));
//  float sw=MAX(0.05,MIN(0.3,(0.3-MAX(0.0,(log(slen)-log(80.0)))*(0.3-0.05)/(log(500.0)-log(80.0)))));
  float sw=MAX(swmin,MIN(swmax,(swmax-MAX(0.0,(log(slen)-log(80.0)))*(swmax-swmin)/(log(500.0)-log(80.0)))));
//  float sw=sweightabs; // MAX(0.025,MIN(0.3,(0.3-MAX(0.0,(log(slen)-log(80.0)))*(0.3-0.05)/(log(500.0)-log(80.0)))));

  for (int l=pinfo.matchcounts.size()-1; l>=0; --l){
    int sbest=pinfo.matchcounts[l].seqid;
    if (tax.seqs[sbest]==0x00) continue;
//    if (pinfo.matchcounts[l].score()<=0.0) break; // do not use alignments with less than or zero score
    eseqtax &taxhit(*tax.seqs[sbest]);

    float cfw,w;
    for (int k=0; k<taxhit.tl.size(); ++k){
//      if (taxcounts[k][taxhit.tl[k].tid]==taxid) continue;
      if (taxhit.tl[k].tid==-1) continue;
      cfw=1.0;
      if (taxhit.bid>0.0)
        cfw=0.5*taxhit.tl[k].cf+0.5;
//      taxcounts[k][taxhit.tl[k].tid]=taxid;
//      w=exp((1.0l-topscore/pinfo.matchcounts[l].score())*sweight)*cfw;
//      w=exp((pinfo.matchcounts[l].score()-topscore)*sweightabs)*cfw;
      w=exp((pinfo.matchcounts[l].score()-topscore)*sw)*taxhit.tl[k].evidence*cfw;
      if (w>taxmaxscores[k][taxhit.tl[k].tid]){
//        cout << "l: " << l << " k: " << k << " w: " << w << " score: " << pinfo.matchcounts[l].score() << " topscore: " << topscore << endl;
        taxscores[k]+=w-taxmaxscores[k][taxhit.tl[k].tid];
        taxmaxscores[k][taxhit.tl[k].tid]=w;
      }
//      taxscores[k]+=exp(log(pinfo.matchcounts[l].identity()*100.0)*30.0);
    }
  }
//  for (int i=0; i<taxscores.size(); ++i)
//    cout << " i: " << i << " maxscore: " << taxscores[i] << endl;
}




void taxScore(earrayof<double,int>& ptax,efloatarray& mcfarr,ealigndata& adata,epredinfo& pinfo,edoublearray& taxscores,etax& tax,int slen)
{
  ptax.clear();
//  while (ptax.size()<tax.names.size())
//    ptax.add(-1,0.0);
  for (int i=0; i<tax.names.size(); ++i)
    ptax.add(-1,0.0);


//  float sw=sweightabs; //MAX(0.025,MIN(0.3,(0.3-MAX(0.0,(log(slen)-log(80.0)))*(0.3-0.05)/(log(500.0)-log(80.0)))));
//  float sw=MAX(0.025,MIN(0.3,(0.3-MAX(0.0,(log(slen)-log(80.0)))*(0.3-0.05)/(log(500.0)-log(80.0)))));

  float sw=MAX(swmin,MIN(swmax,(swmax-MAX(0.0,(log(slen)-log(80.0)))*(swmax-swmin)/(log(500.0)-log(80.0)))));

  mcfarr.init(ptax.size(),0.0);
  if(tax.seqs[adata.seqid]==0x00)
    return;

//  ldieif(tax.seqs[adata.seqid]==0x00,"Missing taxonomy for sequence");
  float cfw;
  eseqtax &seqtax(*tax.seqs[adata.seqid]);
  for (int l=0; l<seqtax.tl.size(); ++l){
    ldieif(seqtax.tl[l].tid>=long(tax.names[l].size()),estr("key out of tax: ")+adata.seqid+" "+estr(l)+" "+seqtax.tl[l].tid+" "+tax.names[l].size());
    cfw=1.0;
    if (seqtax.bid>0.0)
      cfw=0.5*seqtax.tl[l].cf+0.5;
    ptax.keys(l)=seqtax.tl[l].tid;
//    ptax.values(l)=exp((1.0l-pinfo.tophit.score()/adata.score())*sweight)*cfw/taxscores[l];
//    ptax.values(l)=exp((adata.score()-pinfo.tophit.score())*sweightabs)*cfw/taxscores[l];
    ptax.values(l)=exp((adata.score()-pinfo.tophit.score())*sw)*cfw/taxscores[l];
  }
  calcConfidence(ptax,mcfarr,adata,pinfo,tax);
}

void taxScoreSum(edoublearray& taxscores,epredinfo& pinfo,etax& tax,ebasicarray<eintarray>& taxcounts,int slen)
{
  taxscores.init(tax.names.size(),0.0);
  int tophitl=-1;

  for (int l=pinfo.matchcounts.size()-1; l>=0; --l){
    int sbest=pinfo.matchcounts[l].seqid;
    if (tax.seqs[sbest]!=0x00) { tophitl=l; break; }
  }
  if (tophitl==-1) return;

  ealigndata &tophit(pinfo.matchcounts[tophitl]);
  pinfo.tophit=tophit; 
  ldieif(tax.seqs[tophit.seqid]==0x00,"top hit does not have taxonomy: "+estr(tophit.seqid));

  double topscore=tophit.score();

  int taxid=0;
  ebasicarray<efloatarray> taxmaxscores;
  taxmaxscores.init(tax.names.size());
  for (int i=0; i<taxmaxscores.size(); ++i)
    taxmaxscores[i].init(tax.names[i].size(),0.0);
//  taxcounts.init(tax.names.size());
//  for (int i=0; i<taxcounts.size(); ++i)
//    taxcounts[i].init(tax.names[i].size(),-1);

//  float sw=MAX(0.025,MIN(0.3,(0.3-MAX(0.0,(log(slen)-log(80.0)))*(0.3-0.05)/(log(500.0)-log(80.0)))));
//  float sw=MAX(0.05,MIN(0.3,(0.3-MAX(0.0,(log(slen)-log(80.0)))*(0.3-0.05)/(log(500.0)-log(80.0)))));
  float sw=MAX(swmin,MIN(swmax,(swmax-MAX(0.0,(log(slen)-log(80.0)))*(swmax-swmin)/(log(500.0)-log(80.0)))));
//  float sw=sweightabs; // MAX(0.025,MIN(0.3,(0.3-MAX(0.0,(log(slen)-log(80.0)))*(0.3-0.05)/(log(500.0)-log(80.0)))));

  for (int l=pinfo.matchcounts.size()-1; l>=0; --l){
    int sbest=pinfo.matchcounts[l].seqid;
    if (tax.seqs[sbest]==0x00) continue;
//    if (pinfo.matchcounts[l].score()<=0.0) break; // do not use alignments with less than or zero score
    eseqtax &taxhit(*tax.seqs[sbest]);

    float cfw,w;
    for (int k=0; k<taxhit.tl.size(); ++k){
//      if (taxcounts[k][taxhit.tl[k].tid]==taxid) continue;
      if (taxhit.tl[k].tid==-1) continue;
      cfw=1.0;
      if (taxhit.bid>0.0)
        cfw=0.5*taxhit.tl[k].cf+0.5;
//      taxcounts[k][taxhit.tl[k].tid]=taxid;
//      w=exp((1.0l-topscore/pinfo.matchcounts[l].score())*sweight)*cfw;
//      w=exp((pinfo.matchcounts[l].score()-topscore)*sweightabs)*cfw;
      w=exp((pinfo.matchcounts[l].score()-topscore)*sw)*cfw;
      if (w>taxmaxscores[k][taxhit.tl[k].tid]){
//        cout << "l: " << l << " k: " << k << " w: " << w << " score: " << pinfo.matchcounts[l].score() << " topscore: " << topscore << endl;
        taxscores[k]+=w-taxmaxscores[k][taxhit.tl[k].tid];
        taxmaxscores[k][taxhit.tl[k].tid]=w;
      }
//      taxscores[k]+=exp(log(pinfo.matchcounts[l].identity()*100.0)*30.0);
    }
  }
//  for (int i=0; i<taxscores.size(); ++i)
//    cout << " i: " << i << " maxscore: " << taxscores[i] << endl;
}






inline void sumcounts(unsigned char* tncount,int p1,uint64_t tmpnuc)
{
  *reinterpret_cast<uint64_t*>(&tncount[p1])|=tmpnuc;
/*
  for (int i=0; i<8; ++i,tmpnuc>>=8u)
    tncount[p1+i]|=(tmpnuc&0xff);
*/
//  cout << uint2hex(tmpstr,*p) << " " << uint2hex(tmpstr2,tmpnuc) << " " << uint2hex(tmpstr3,tncount[p1]) << endl;
}

unsigned long seqkmer(const eseq& s,long p1)
{
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s.seq._str);
  return(((pstr1[p1/32u]>>(2u*(p1%32u)))|((pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u]))&KMERMASK);
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

void pseqncount(const eseq& s1,long p1,long e1,const eseq& s2,long p2,long e2,ealigndata& adata,const epalignscore& pas)
{
//  cout << "s2seqlen: " << s2.seqlen << endl;
  int k;
//  for (int i=0; i+PKMERSIZE<=e1-p1; i+=PKMERSIZE){
//    cout << p1+i << ": " << ckmer2str(seqpkmer(s1,p1+i)) << " = p2: " << p2+i << ": " << pkmer2str(seqpkmer(s2,p2+i)) << endl;
//    cout << " = p1: " << p1+i-1 << ": " << ckmer2str(seqpkmer(s1,p1+i-1)) << endl;
//    cout << " = p1: " << p1+i+1 << ": " << ckmer2str(seqpkmer(s1,p1+i+1)) << endl;
//  }
  adata.matches+=(e1-p1)/3;
  adata._score+=(e1-p1)*pas.pmatch/3;
  adata.profile.add(AT_MATCH,(e1-p1)/3);
}




void find_diags(const eseq& s1,const eseq& s2,long s1_start,long s1_end,euintarray& kmerpos1,int offset1,esearchws& sws,ebasicarray<ediag> &diags)
{
  long p1,p2;
  unsigned long *pstr1=reinterpret_cast<unsigned long*>(s1.seq._str);
  unsigned long *pstr2=reinterpret_cast<unsigned long*>(s2.seq._str);
  unsigned long v1,v2;
  long lastkmerlen,lastkmerpos,lastndelta;
  long s1_len=s1_end-s1_start;
  int k;

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
              if (lastkmerpos-lastkmerlen+KMERSIZE>s1.seqlen){
                cerr << "p1: " << lastkmerpos-lastkmerlen+KMERSIZE<< "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
                exit(-1);
              }
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
              if (lastkmerpos-lastkmerlen+KMERSIZE>s1.seqlen){
                cerr << "p1: " << lastkmerpos-lastkmerlen+KMERSIZE<< "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
                exit(-1);
              }
            diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE,lastndelta,lastkmerlen));
          }
          lastkmerlen=KMERSIZE;
          lastndelta=kpos2-p1-k;
        }
        lastkmerpos=p1+k;
      }
    }
    if (lastkmerlen){
              if (lastkmerpos-lastkmerlen+KMERSIZE>s1.seqlen){
                cerr << "p1: " << lastkmerpos-lastkmerlen+KMERSIZE<< "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
                exit(-1);
              }

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
        if (kpos1>offset1){
          kpos1+=s1_start-offset1;
          if (kpos1>s1.seqlen){
            cerr << endl << "#kmerpos1 out of bounds: " << endl;
            exit(0);
          }
          long d=p2+k-lastkmerpos;
          if (d<=KMERSIZE && ((lastndelta==p2+k-kpos1) || (p2+k-lastndelta+KMERSIZE<=s1_end && (v2&KMERMASK)==seqkmer(s1,p2+k-lastndelta)))){
            lastkmerlen+=d;
          }else{
            if (lastkmerlen){
//              cout << "p2: " << lastkmerpos-lastkmerlen+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
              if (lastkmerpos-lastkmerlen+KMERSIZE-lastndelta>s1.seqlen){
                cerr << "p2: " << lastkmerpos-lastkmerlen+KMERSIZE<< "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
                exit(-1);
              }
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
      if (kpos1>offset1){
        kpos1+=s1_start-offset1;
          if (kpos1>s1.seqlen){
            cerr << endl << "#kmerpos1 out of bounds: " << endl;
            exit(0);
          }
        long d=p2+k-lastkmerpos;
        if (d<=KMERSIZE && ((lastndelta==p2+k-kpos1) || (p2+k-lastndelta+KMERSIZE<=s1_end && (v2&KMERMASK)==seqkmer(s1,p2+k-lastndelta)))){
          lastkmerlen+=d;
        }else{
          if (lastkmerlen){
              if (lastkmerpos-lastkmerlen+KMERSIZE-lastndelta>s1.seqlen){
                cerr << "p2: " << lastkmerpos-lastkmerlen+KMERSIZE<< "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
                exit(-1);
              }
            diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE-lastndelta,lastndelta,lastkmerlen));
          }
          lastkmerlen=KMERSIZE;
          lastndelta=p2+k-kpos1;
        }
        lastkmerpos=p2+k;
      }
    }
    if (lastkmerlen){
              if (lastkmerpos-lastkmerlen+KMERSIZE-lastndelta>s1.seqlen){
                cerr << "p2: " << lastkmerpos-lastkmerlen+KMERSIZE<< "," << lastkmerpos+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
                exit(-1);
              }
//      cout << "end p2: " << lastkmerpos-lastkmerlen+KMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
      diags.add(ediag(lastkmerpos-lastkmerlen+KMERSIZE-lastndelta,lastndelta,lastkmerlen));
    }
  }
}

ediag *sparse_dynamic_prog(ebasicarray<ediag> &diags,const ealignscore& as)
{
  if (diags.size()==0) return(0x00);

  // sparse dynamic programming to find longest chain of segments
  multimap<long,ediag*> segments;

  ediag *bestseg=&diags[0];
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
  return(bestseg);
}

#define align(a,b) (a-(a)%b)

void pseqident_local(const eseq& s1,euintarray& kmerpos1,unsigned int offset1,const eseq& s2,ealigndata& adata,esearchws& sws,const epalignscore& pas,long s1_start=0,long s1_end=0)
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

//  t2.reset();
    lastkmerpos=s1_start-2*long(PKMERSIZE);
    lastkmerlen=0;
    lastndelta=-s1_end-s2.seqlen;
    for (p1=s1_start; p1<s1_end-32; p1+=k){
      v1=pstr1[p1/32u]>>(2u*(p1%32u));
      v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
      for (k=0; k<32u-PKMERSIZE; k+=3,v1>>=6u){
        long d=p1+k-lastkmerpos;
        if (d<=PKMERSIZE && p1+k+lastndelta>=0 && p1+k+lastndelta+PKMERSIZE<=s2.seqlen && kmer_prot_lt[(v1&PKMERMASK)]==seqpkmer(s2,p1+k+lastndelta))
        {
          lastkmerlen+=d;
          lastkmerpos=p1+k;
        }else{
          long kpos2=sws.kmerpos2[kmer_prot_lt[v1&PKMERMASK]];
          if (kpos2>sws.offset2){
            kpos2-=sws.offset2;
          
            if (d<=PKMERSIZE && (lastndelta==kpos2-p1-k)){
              lastkmerlen+=d;
            }else{
              if (lastkmerlen){
                cout << "p1: " << lastkmerpos-lastkmerlen+PKMERSIZE<< "," << lastkmerpos+PKMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << " kpos2: " << kpos2 << " p1+k: " << p1+k << endl;
                diags.add(ediag(lastkmerpos-lastkmerlen+PKMERSIZE,lastndelta,lastkmerlen));
              }
              lastkmerlen=PKMERSIZE;
              lastndelta=kpos2-p1-k;
            }
            lastkmerpos=p1+k;
          }
        }
      }
    }
    v1=pstr1[p1/32u]>>(2u*(p1%32u));
    v1|=(pstr1[p1/32u+1u]<<(64u-2u*(p1%32u)))&safe_shift[p1%32u];
    for (k=0; p1+k+PKMERSIZE<s1_end; k+=3,v1>>=6u){
      long kpos2=sws.kmerpos2[kmer_prot_lt[v1&PKMERMASK]];
      if (kpos2>sws.offset2){
        kpos2-=sws.offset2;
        long d=p1+k-lastkmerpos;
        if (d<=PKMERSIZE && ((lastndelta==kpos2-p1-k) || (p1+k+lastndelta+PKMERSIZE<=s2.seqlen && kmer_prot_lt[(v1&PKMERMASK)]==seqpkmer(s2,p1+k+lastndelta)))){
          lastkmerlen+=d;
        }else{
          if (lastkmerlen){
            diags.add(ediag(lastkmerpos-lastkmerlen+PKMERSIZE,lastndelta,lastkmerlen));
          }
          lastkmerlen=PKMERSIZE;
          lastndelta=kpos2-p1-k;
        }
        lastkmerpos=p1+k;
      }
    }
    if (lastkmerlen){
      diags.add(ediag(lastkmerpos-lastkmerlen+PKMERSIZE,lastndelta,lastkmerlen));
      cout << "end p1: " << lastkmerpos-lastkmerlen+PKMERSIZE << "," << lastkmerpos+PKMERSIZE << " len: " << lastkmerlen << " ndelta: " << lastndelta << endl;
    }
//  tdp1=tdp1*0.99+t2.lap()*0.01;

  if (diags.size()==0){
//    lderror("no shared segments found");
    return;
  }

  // sparse dynamic programming to find longest chain of segments
  ediag *bestseg=&diags[0];
  multimap<long,ediag*> segments;
  for (long j=0; j<diags.size(); ++j){
    ediag& diag(diags.at(j));
    diag.V=(diag.j-diag.i)*pas.pmatch/3; // should compute protein match score
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
      double tmpbestscore=(sit->second->j-sit->second->i)*pas.pmatch/3;
      for (slb=Lsegments.begin(); slb!=Lsegments.end(); ++slb){
        if (slb->first >= sit->second->i2) break;
//      sit->second->V=sit->second->j-sit->second->i+slb->second->V;
        long gapdiff=abs((sit->second->i-sit->second->i2)-(slb->second->i-slb->second->i2));
        double tmpscore=(sit->second->j-sit->second->i)*pas.pmatch/3 - gapdiff*pas.gapext/3 -(gapdiff>0?pas.gapopen:0) + slb->second->V;
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
//  tdp2=tdp2*0.99+t2.lap()*0.01;

  int rf=s1_start%3;
  float tmpflanks=0.0;
  ldieif(s2.seqlen<=bestseg->j2,"segment start larger than sequence length: "+estr(s2.seqlen)+" < "+estr(bestseg->j2));
  long minright=MIN(align(s1.seqlen-bestseg->j,3)+rf,s2.seqlen-bestseg->j2);
//  seqcalign_global_norightedgegap(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,alignws);
//  cout << "rf: " << rf << " j: " << bestseg->j << " j2: "<<bestseg->j2 << endl;
//  ldieif(bestseg->j2%3!=0,"not protein aligned");
//  ldieif(bestseg->j%3!=s1_start%3,"not equal codon alignment: "+estr(bestseg->j)+" "+estr(s1_start));
  pseqcalign_local_rightext(s1,bestseg->j,align(MIN(s1.seqlen,bestseg->j+2*minright),3)+rf,s2,bestseg->j2,align(MIN(s2.seqlen,bestseg->j2+2*minright),3),adata,sws.alignws,pas);
//  tmpflanks+=t2.lap();
  adata.profile.add(AT_RIGHT);
  LDEBUG(D_SEQIDENT,cout << "rightid: " << bestseg->j << "," << s1_end << " l: " << s1_end-bestseg->j << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
  LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,s1_start,s1_end));
  while (bestseg!=0x00 && bestseg->bestseg != 0x00){
    pseqncount(s1,bestseg->i,bestseg->j,s2,bestseg->i2,bestseg->j2,adata,pas);
    adata.profile.add(AT_ID);
    LDEBUG(D_SEQIDENT,cout << "seg: "<< bestseg->i << "," << bestseg->j << " l: " << (bestseg->j-bestseg->i) << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
    LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,s1_start,s1_end));
    long gapdiff=abs((bestseg->bestseg->i-bestseg->bestseg->i2)-(bestseg->i-bestseg->i2));
    pseqcalign_global(s1,bestseg->bestseg->j,bestseg->i,s2,bestseg->bestseg->j2,bestseg->i2,adata,sws.alignws,pas);
    adata.profile.add(AT_ALIGN);
    LDEBUG(D_SEQIDENT,cout << "segid: " << bestseg->bestseg->j << "," << bestseg->i << " l: " << bestseg->i-bestseg->bestseg->j << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << " gapdiff: " << gapdiff << endl);
    LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,s1_start,s1_end));
//    cout << "segid: " << bestseg->i << "," << bestseg->j << " l: " << bestseg->j-bestseg->i << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl;
    bestseg=bestseg->bestseg;
  }
  pseqncount(s1,bestseg->i,bestseg->j,s2,bestseg->i2,bestseg->j2,adata,pas);
  adata.profile.add(AT_ID);
  LDEBUG(D_SEQIDENT,cout << "seg: "<< bestseg->i << "," << bestseg->j << " l: " << (bestseg->j-bestseg->i) << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
  LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,s1_start,s1_end));


//  tdpmd=tdpmd*0.99+t2.lap()*0.01;
//  seqident_seg_left_local(s1,bestseg->i,s2,bestseg->i2,adata,tncounts);
  long minleft=MIN(align(bestseg->i,3)+rf,bestseg->i2);
  pseqcalign_local_leftext(s1,rf+MAX(0,align(bestseg->i-2*minleft,3)),bestseg->i,s2,MAX(0,align(bestseg->i2-2*minleft,3)),bestseg->i2,adata,sws.alignws,pas);
  adata.profile.add(AT_LEFT);
//  tmpflanks+=t2.lap();
//  tdpfl=tdpfl*0.99+tmpflanks*0.01;

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


void seqident_local_leftext(const estr& s1id,const estr& s2id,const eseq& s1,euintarray& kmerpos1,unsigned int offset1,const eseq& s2,ealigndata& adata,esearchws& sws,const ealignscore& as,long s1_start=0,long s1_end=0)
{
  if (s1_end==0) s1_end=s1.seqlen;
  long s1_len=s1_end-s1_start;

  adata.matches=0; adata.mismatches=0; adata.gaps=0; adata._score=0.0;
  adata.s1=-1; adata.e1=-1;
  adata.s2=-1; adata.e2=-1;

  char tmp[5];
  char itmp;
  ebasicarray<ediag> diags;

  LDEBUG(D_PROFILE,t2.reset());
  find_diags(s1,s2,s1_start,s1_end,kmerpos1,offset1,sws,diags);
  LDEBUG(D_PROFILE,tdp1=tdp1*0.99+t2.lap()*0.01);

  if (diags.size()==0){
    return;
  }

  ediag *bestseg=sparse_dynamic_prog(diags,as);
  LDEBUG(D_PROFILE,tdp2=tdp2*0.99+t2.lap()*0.01);


  float tmpflanks=0.0;
  ldieif(s1.seqlen<bestseg->j,"segment start larger than sequence s1 length: "+estr(s1.seqlen)+" < "+estr(bestseg->j)+" s1 id: "+s1id+" s2 id: "+s2id+" i: "+bestseg->i);
  ldieif(s2.seqlen<bestseg->j2,"segment start larger than sequence s2 length: "+estr(s2.seqlen)+" < "+estr(bestseg->j2)+" s1 id: "+s1id+" s2 id: "+s2id);
  long minright=MIN(s1.seqlen-bestseg->j,s2.seqlen-bestseg->j2);
//  seqcalign_global_norightedgegap(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,alignws);
  seqcalign_global_norightedgegap(s1,bestseg->j,MIN(s1.seqlen,bestseg->j+2*minright),s2,bestseg->j2,MIN(s2.seqlen,bestseg->j2+2*minright),adata,sws.alignws,as);
//  seqcalign_local_rightext(s1,bestseg->j,MIN(s1.seqlen,bestseg->j+2*minright),s2,bestseg->j2,MIN(s2.seqlen,bestseg->j2+2*minright),adata,sws.alignws,as);
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
//    if (gapdiff==0) {
//      seqident_seg_nogaps(s1,bestseg->bestseg->j,bestseg->i,s2,bestseg->bestseg->j2,bestseg->i2,adata,as);
//      adata.profile.add(AT_ALIGNF);
//    } else { 
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
//  seqcalign_global_noleftedgegap(s1,MAX(0,bestseg->i-2*minleft),bestseg->i,s2,MAX(0,bestseg->i2-2*minleft),bestseg->i2,as,as);
//  seqcalign_global_noleftedgegap(s1,0,bestseg->i,s2,0,bestseg->i2,as1,as2);
  adata.profile.add(AT_LEFT);
  LDEBUG(D_PROFILE,tmpflanks+=t2.lap());
//  tdpfl=tdpfl*0.99+tmpflanks*0.01;

//  estr as1,as2;
//  cout << endl;
//  cout << as1 << endl << as2 << endl;
  LDEBUG(D_SEQIDENT,cout << "leftid: "<< bestseg->i << " l: " << (bestseg->i<bestseg->i2?bestseg->i:bestseg->i2) << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
  LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,s1_start,s1_end));
//  cout << "partial alignment" << endl;
//  print_seqali(pas1,pas2);

  return;
}

void seqident_local_rightext(const estr& s1id,const estr& s2id,const eseq& s1,euintarray& kmerpos1,unsigned int offset1,const eseq& s2,ealigndata& adata,esearchws& sws,const ealignscore& as,long s1_start=0,long s1_end=0)
{
  if (s1_end==0) s1_end=s1.seqlen;
  long s1_len=s1_end-s1_start;

  adata.matches=0; adata.mismatches=0; adata.gaps=0; adata._score=0.0;
  adata.s1=-1; adata.e1=-1;
  adata.s2=-1; adata.e2=-1;

  char tmp[5];
  char itmp;
  ebasicarray<ediag> diags;

  LDEBUG(D_PROFILE,t2.reset());
  find_diags(s1,s2,s1_start,s1_end,kmerpos1,offset1,sws,diags);
  LDEBUG(D_PROFILE,tdp1=tdp1*0.99+t2.lap()*0.01);

  if (diags.size()==0){
    return;
  }

  ediag *bestseg=sparse_dynamic_prog(diags,as);
  LDEBUG(D_PROFILE,tdp2=tdp2*0.99+t2.lap()*0.01);


  float tmpflanks=0.0;
  ldieif(s2.seqlen<bestseg->j2,"segment start larger than sequence length: "+estr(s2.seqlen)+" < "+estr(bestseg->j2)+" s1 id: "+s1id+" s2 id: "+s2id);
  long minright=MIN(s1.seqlen-bestseg->j,s2.seqlen-bestseg->j2);
//  seqcalign_global_norightedgegap(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,alignws);
//  seqcalign_global_norightedgegap(s1,bestseg->j,MIN(s1.seqlen,bestseg->j+2*minright),s2,bestseg->j2,MIN(s2.seqlen,bestseg->j2+2*minright),adata,sws.alignws);
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
//    if (gapdiff==0) {
//      seqident_seg_nogaps(s1,bestseg->bestseg->j,bestseg->i,s2,bestseg->bestseg->j2,bestseg->i2,adata,as);
//      adata.profile.add(AT_ALIGNF);
//    } else { 
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
//  seqcalign_local_leftext(s1,MAX(0,bestseg->i-2*minleft),bestseg->i,s2,MAX(0,bestseg->i2-2*minleft),bestseg->i2,adata,sws.alignws,as);
  seqcalign_global_noleftedgegap(s1,MAX(0,bestseg->i-2*minleft),bestseg->i,s2,MAX(0,bestseg->i2-2*minleft),bestseg->i2,adata,sws.alignws,as);
//  seqcalign_global_noleftedgegap(s1,0,bestseg->i,s2,0,bestseg->i2,as1,as2);
  adata.profile.add(AT_LEFT);
  LDEBUG(D_PROFILE,tmpflanks+=t2.lap());
//  tdpfl=tdpfl*0.99+tmpflanks*0.01;

//  estr as1,as2;
//  cout << endl;
//  cout << as1 << endl << as2 << endl;
  LDEBUG(D_SEQIDENT,cout << "leftid: "<< bestseg->i << " l: " << (bestseg->i<bestseg->i2?bestseg->i:bestseg->i2) << " d: " << bestseg->i2-bestseg->i << " m: " << adata.matches << " miss: " << adata.mismatches << " gaps: " << adata.gaps << endl);
  LDEBUG(D_SEQALIGNMENT,print_tncount(adata.aln,s1_start,s1_end));
//  cout << "partial alignment" << endl;
//  print_seqali(pas1,pas2);

  return;
}


void seqident_local(const estr& s1id,const estr& s2id,const eseq& s1,euintarray& kmerpos1,unsigned int offset1,const eseq& s2,ealigndata& adata,esearchws& sws,const ealignscore& as,long s1_start=0,long s1_end=0)
{
  if (s1_end==0) s1_end=s1.seqlen;
  long s1_len=s1_end-s1_start;

//  for (int i=0; i<11; ++i)
//    match[i]=0;
  adata.matches=0; adata.mismatches=0; adata.gaps=0; adata._score=0.0;
  adata.s1=-1; adata.e1=-1;
  adata.s2=-1; adata.e2=-1;

  char tmp[5];
  char itmp;
  ebasicarray<ediag> diags;

  LDEBUG(D_PROFILE,t2.reset());
  find_diags(s1,s2,s1_start,s1_end,kmerpos1,offset1,sws,diags);
  LDEBUG(D_PROFILE,tdp1=tdp1*0.99+t2.lap()*0.01);

  if (diags.size()==0){
//    lderror("no shared segments found");
    return;
  }

  ediag *bestseg=sparse_dynamic_prog(diags,as);
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
//  tdpfl=tdpfl*0.99+tmpflanks*0.01;

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


void seqident_global(const estr& s1id,const estr& s2id,const eseq& s1,euintarray& kmerpos1,unsigned int offset1,const eseq& s2,ealigndata& adata,esearchws& sws,const ealignscore& as,long s1_start=0,long s1_end=0)
{
  if (s1_end==0) s1_end=s1.seqlen;
  long s1_len=s1_end-s1_start;

//  for (int i=0; i<11; ++i)
//    match[i]=0;
  adata.matches=0; adata.mismatches=0; adata.gaps=0; adata._score=0.0;
  adata.s1=-1; adata.e1=-1;
  adata.s2=-1; adata.e2=-1;

  char tmp[5];
  char itmp;
  ebasicarray<ediag> diags;

  LDEBUG(D_PROFILE,t2.reset());
  find_diags(s1,s2,s1_start,s1_end,kmerpos1,offset1,sws,diags);
  LDEBUG(D_PROFILE,tdp1=tdp1*0.99+t2.lap()*0.01);

  if (diags.size()==0){
//    lderror("no shared segments found");
    return;
  }

  ediag *bestseg=sparse_dynamic_prog(diags,as);
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
  seqcalign_global(s1,bestseg->j,s1.seqlen,s2,bestseg->j2,s2.seqlen,adata,sws.alignws,as);
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
  seqcalign_global(s1,0,bestseg->i,s2,0,bestseg->i2,adata,sws.alignws,as);
  adata.profile.add(AT_LEFT);
  LDEBUG(D_PROFILE,tmpflanks+=t2.lap());
//  tdpfl=tdpfl*0.99+tmpflanks*0.01;

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




void setpkmerpos(euintarray& kmerpos,eseq& s,unsigned int offset)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  int k;
  long p;
  for (p=0; p+32<s.seqlen; p+=k){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k+PKMERSIZE<32; k+=3u,v>>=6u)
//      kmerpos[kmer_prot_lt[v&PKMERMASK]]=offset+p+k;
      kmerpos[v&PKMERMASK]=offset+p+k;
  }
  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k+PKMERSIZE<s.seqlen; k+=3u,v>>=6u)
//    kmerpos[kmer_prot_lt[v&PKMERMASK]]=offset+p+k;
    kmerpos[v&PKMERMASK]=offset+p+k;
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



void kmercount_single(ebasicarray<ekmerarray>& otukmers,eseq& s,ebasicarray<uint32_t>& idcount,eintarray& kmerpos,const euintarray& otukmersize)
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
        ekmerarray &otuind(otukmers[v&KMERMASK]);
        unsigned int otuindsize=otukmersize[v&KMERMASK];
        for (unsigned int l=0; l<otuindsize; ++l){
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
      ekmerarray &otuind(otukmers[v&KMERMASK]);
      unsigned int otuindsize=otukmersize[v&KMERMASK];
      for (unsigned int l=0; l<otuindsize; ++l){
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


void kmercount_single(ebasicarray<ekmerarray>& otukmers,eseq& s,ebasicarray<uint32_t>& idcount,eintarray& kmerpos)
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
        ekmerarray &otuind(otukmers[v&KMERMASK]);
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
      ekmerarray &otuind(otukmers[v&KMERMASK]);
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

void kmercount_both_nopos2_skip(int scount,ebasicarray<ekmerarray>& otukmers,eseq& s,unsigned long start,unsigned long end,eintarray& idcount,eintarray& kmerpos,unsigned int *akmers,unsigned long akmask,int skipn)
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
        ekmerarray &otuind(otukmers[v&KMERMASK]);
        for (unsigned int l=0; l<otuind.size(); ++l)
          ++idcount[(otuind[l]&BMASK31)];
//          idcount[(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
      }
      if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
//      if (otukmers[kmer_rev_lt[v&KMERMASK]]){
        ekmerarray &otuind(otukmers[kmer_rev_lt[v&KMERMASK]]);
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
      ekmerarray &otuind(otukmers[v&KMERMASK]);
      for (unsigned int l=0; l<otuind.size(); ++l)
        ++idcount[(otuind[l]&BMASK31)];
//        idcount[(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
    }
//    if (otukmers[kmer_rev_lt[v&KMERMASK]]){
    if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
      ekmerarray &otuind(otukmers[kmer_rev_lt[v&KMERMASK]]);
      for (unsigned int l=0; l<otuind.size(); ++l)
        ++idcount[scount+(otuind[l]&BMASK31)];
//        idcount[scount+(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
    }
  }
}


void kmercount_both_prot_skip(int scount,ebasicarray<ekmerarray>& otukmers,eseq& s,unsigned long start,unsigned long end,eintarray& idcount,eintarray& kmerpos,unsigned int *akmers,unsigned long akmask,int skipn)
{
  unsigned long *pstr=reinterpret_cast<unsigned long*>(s.seq._str);
  unsigned long v;
  int k;
  long p=start;
  for (; p+32<end; p+=32-PKMERSIZE){
    v=pstr[p/32u]>>(2u*(p%32u));
    v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
    for (k=0; k<32-PKMERSIZE; ++k,v>>=2u){
//      cout << p+k << " : " << kmer2str(v&KMERMASK) << " -- " << pkmer2str(kmer_prot_lt[v&PKMERMASK]) << endl;
//      if ((p+k)%skipn) continue;
//      cout << kmer2str(v&KMERMASK) << " " << kmer2str(kmer_rev_lt[v&KMERMASK]) << endl;
//      if (otukmers[v&KMERMASK]){
//      if (akmers[v&akmask]){
      {
        ekmerarray &otuind(otukmers[kmer_prot_lt[v&PKMERMASK]]);
        for (unsigned int l=0; l<otuind.size(); ++l)
          ++idcount[scount*((p+k)%3)*2+(otuind[l]&BMASK31)];
//          idcount[scount*((p+k)%3)*2+(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
      }
//      }
//      if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
//      if (otukmers[kmer_rev_lt[v&KMERMASK]]){
      {
        ekmerarray &otuind(otukmers[kmer_protrev_lt[v&PKMERMASK]]);
        for (unsigned int l=0; l<otuind.size(); ++l)
          ++idcount[scount*((p+k)%3)*2+scount+(otuind[l]&BMASK31)];
//          idcount[scount*((p+k)%3)*2+scount+(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
      }
//      }
    }
  }

//  cout << "p: " << p << " seqlen: " << s.seqlen << endl;

  v=pstr[p/32u]>>(2u*(p%32u));
  v|=(pstr[p/32u+1u]<<(64u-2u*(p%32u)))&safe_shift[p%32u];
  for (k=0; p+k+PKMERSIZE<end; ++k,v>>=2u){
//    if ((p+k)%skipn) continue;
//    if (otukmers[v&KMERMASK]){
//    if (akmers[v&akmask]){
    {
      ekmerarray &otuind(otukmers[kmer_prot_lt[v&PKMERMASK]]);
      for (unsigned int l=0; l<otuind.size(); ++l)
        ++idcount[scount*((p+k)%3)*2+(otuind[l]&BMASK31)];
//        idcount[scount*((p+k)%3)*2+(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
    }
//    }
//    if (otukmers[kmer_rev_lt[v&KMERMASK]]){
//    if (akmers[kmer_rev_lt[v&KMERMASK]&akmask]){
    {
      ekmerarray &otuind(otukmers[kmer_protrev_lt[v&PKMERMASK]]);
      for (unsigned int l=0; l<otuind.size(); ++l)
        ++idcount[scount*((p+k)%3)*2+scount+(otuind[l]&BMASK31)];
//        idcount[scount*((p+k)%3)*2+scount+(otuind[l]&BMASK31)]+=idc_lt[otuind[l]>>31u];
    }
//    }
  }
}




void otukmeradd(ebasicarray<ekmerarray>& otukmers,int i,eseq& s,eintarray& tmpkmers,int ti,unsigned int *akmers,unsigned long akmask,euintarray &kmermask)
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
      if (kmermask[v&KMERMASK]==1 || tmpkmers[v&KMERMASK]==ti) continue;
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
    if (kmermask[v&KMERMASK]==1 || tmpkmers[v&KMERMASK]==ti) continue;
    tmpkmers[v&KMERMASK]=ti;
//    if (otukmers[v&KMERMASK]==0x00){
//      otukmers[v&KMERMASK]=new deque<int>();
//      otukmers[v&KMERMASK]->reserve(1000);
//    }
    otukmers[v&KMERMASK].push_back(i);
  }
}

void otukmeradd(ebasicarray<ekmerarray>& otukmers,int i,eseq& s,eintarray& tmpkmers,int ti,unsigned int *akmers,unsigned long akmask)
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


esearchws::esearchws() {}

esearchws::esearchws(const eseqdb& seqdb)
{
  init(seqdb);
}

void esearchws::init(const eseqdb& seqdb)
{
//  searchws.kmerbitmask=new uint64_t[MAXSIZE/64+1];
//  searchws.bitmask=new uint64_t[(mtdata.seqdb->otus.size()*2)/64+1];
  otukmerpos.init(seqdb.otus.size()*2,0);
  idcount.init(seqdb.otus.size()*2,0);
  idcount2.init(seqdb.otus.size()*2,0);
  kmerpos.init(MAXSIZE,0u);
  kmerpos2.init(MAXSIZE,0u);
  kmerposrev.init(MAXSIZE,0u);
  offset=1u;
  offset2=1u;
  kmermask.init(MAXSIZE,0u);
  maskid=1u;
/*
  // needed for paired end search
  searchws.kmerpos3.init(MAXSIZE,0u);
  searchws.kmerposrev3.init(MAXSIZE,0u);
  searchws.offset3=1u;
*/
}

void esearchws::initPaired(const eseqdb& seqdb)
{
//  searchws.kmerbitmask=new uint64_t[MAXSIZE/64+1];
//  searchws.bitmask=new uint64_t[(mtdata.seqdb->otus.size()*2)/64+1];
  otukmerpos.init(seqdb.otus.size()*2,0);
  idcount.init(seqdb.otus.size()*2,0);
  idcount2.init(seqdb.otus.size()*2,0);
  kmerpos.init(MAXSIZE,0u);
  kmerpos2.init(MAXSIZE,0u);
  kmerposrev.init(MAXSIZE,0u);
  offset=1u;
  offset2=1u;
  kmermask.init(MAXSIZE,0u);
  maskid=1u;

  // needed for paired end search
  kmerpos3.init(MAXSIZE,0u);
  kmerposrev3.init(MAXSIZE,0u);
  offset3=1u;
}

void esearchws::initProt(const eseqdb& seqdb)
{
//  searchws.kmerbitmask=new uint64_t[MAXSIZE/64+1];
//  searchws.bitmask=new uint64_t[(mtdata.seqdb->otus.size()*2)/64+1];
  otukmerpos.init(seqdb.otus.size()*2,0);
  idcount.init(seqdb.otus.size()*2,0);
  idcount2.init(seqdb.otus.size()*2,0);
  kmerpos.init(MAXSIZE,0u);
  kmerpos2.init(PMAXSIZE,0u);
//  searchws.kmerpos2.init(MAXSIZE,0u);
  kmerposrev.init(MAXSIZE,0u);
  offset=1u;
  offset2=1u;
  kmermask.init(MAXSIZE,0u);
  maskid=1u;

/*
  // needed for paired end search
  searchws.kmerpos3.init(MAXSIZE,0u);
  searchws.kmerposrev3.init(MAXSIZE,0u);
  searchws.offset3=1u;
*/
}

//eseqdb::eseqdb(): minscore(30),tophits(20),topotus(10),otulim(50)
// changed minscore to 60 to reduce false classifications in WGS sequence data
eseqdb::eseqdb(): minscore(60),tophits(20),topotus(10),otulim(50)
{
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
}


void eseqdb::seqsearch(const estr& str2id,eseq& s,earray<epredinfo>& pinfoarr,esearchws& sws)
{
//  t1.reset();
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
  
    sws.otukmerpos.init(otus.size()*2,0);
    sws.idcount.init(otus.size()*2,0);

//    cout << "# counting" << endl;
    kmercount_both_nopos2_skip(otus.size(),otukmers,s,s_start,s_end,sws.idcount,sws.otukmerpos,akmers,0x0Fu,step);
//    ti=ti*0.99+t1.lap()*0.01;
  
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
//    ts=ts*0.99+t1.lap()*0.01;

    eintarray &seqids(sws.seqids);
    seqids.clear();
    for (int l=0; l<best.size() && l<topotus; ++l){
      int ibest=best[l];
      if (ibest<otus.size()){
  //      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest][0]) << "\t" << idcount[ibest] << endl;
        for (int l2=0; l2<otus[ibest].size() && (otulim==0 || l2<otulim); ++l2)
          seqids.add(otus[ibest][l2]);
      }else{
  //      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest-otus.size()][0]) << "\t" << idcount[ibest] << " reversed" << endl;
        ibest-=otus.size();
        for (int l2=0; l2<otus[ibest].size() && (otulim==0 || l2<otulim); ++l2)
          seqids.add(seqs.size()+otus[ibest][l2]);
      }
    }

    eintarray otureps;
  
    // add representatives from all non-chosen OTUS (may improve confidence and novel otu estimation)
    for (int l=topotus; l<topotus+3 && l<best.size(); ++l){
      int ibest=best[l];
      if (ibest<otus.size()){
        if (otus[ibest].size()==0) continue;
        otureps.add(otus[ibest][0]);
      } else {
  //      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest-otus.size()][0]) << "\t" << idcount[ibest] << " reversed" << endl;
        ibest-=otus.size();
        if (otus[ibest].size()==0) continue;
        otureps.add(seqs.size()+otus[ibest][0]);
      }
    }
  
//    cout << "# 2nd counting" << endl;

    randomize(sws.rng,seqids);
    sws.idcount.init(seqids.size(),0);
  
    if (sws.maskid+1u<sws.maskid){
      sws.maskid=1u;
      sws.kmermask.init(KMERSIZE,0u);
    }
    ++sws.maskid;
  
//    ti2=ti2*0.99+t1.lap()*0.01;
  //  memset(sws.kmerbitmask,0,(MAXSIZE2/64+1)*sizeof(uint64_t));
  //  setkmermask(sws.kmerbitmask,s,akmers,0xFul);
  //  kmercount_mask(db.seqs,seqids,sws.kmerbitmask,sws.maskid,sws.idcount,akmers,0xFul);
//    cout << "# 2nd counting -- setkmermask" << endl;
    setkmermask(sws.kmermask,s,sws.maskid,akmers,0xFul,s_start,s_end);
//    cout << "# 2nd counting -- kmercount_mask" << endl;
    kmercount_mask(seqs,seqids,sws.kmermask,sws.maskid,sws.idcount);
//    ts2=ts2*0.99+t1.lap()*0.01;
  
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

    for (int i=0; i<otureps.size(); ++i){
      best.add(seqids.size());
      seqids.add(otureps[i]);
      sws.idcount.add(-1); // make sure that idcount and seqids match
    }
//    best+=otureps;
   
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
      int sbest=seqids[best[l]]%seqs.size();
      eseq &sdb(seqs.values(sbest));
      if (sws.offset2+(unsigned int)(seqs.values(sbest).seqlen)<sws.offset2){  // need an unsigned int here otherwise the comparison is made in long and the offset is not correctly reset, signed int overflows are undefined so this cannot be done with signed ints either
        sws.offset2=1u;
        sws.kmerpos2.init(MAXSIZE,0);
      }
      setkmerpos(sws.kmerpos2,seqs.values(sbest),sws.offset2);
  
      ealigndata adata;
      adata.seqid=sbest;
      adata.revcompl=(seqids[best[l]]>=seqs.size());
      adata.kmercount=sws.idcount[best[l]];
    
      if (seqids[best[l]]<seqs.size()){
//        cout << "# forward align" << endl;
        seqident_local(str2id,seqs.keys(sbest),s,sws.kmerpos,sws.offset,sdb,adata,sws,as,s_start,s_end);
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
        seqident_local(str2id,seqs.keys(sbest),srev,sws.kmerposrev,sws.offset,sdb,adata,sws,as);
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
//        adata._eval=sdb.seqlen*exp(-lambda*adata.score()); // for K-A stats we need (*s.seqlen) but this is constant
        pinfo.matchcounts.add(adata);
      }
  
      sws.offset2+=seqs.values(sbest).seqlen;
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
//  ta=ta*0.99+t1.lap()*0.01;
}

void eseqdb::seqsearch_global(const estr& str2id,eseq& s,earray<epredinfo>& pinfoarr,esearchws& sws)
{
  //t1.reset();
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
  
    sws.otukmerpos.init(otus.size()*2,0);
    sws.idcount.init(otus.size()*2,0);

//    cout << "# counting" << endl;
    kmercount_both_nopos2_skip(otus.size(),otukmers,s,s_start,s_end,sws.idcount,sws.otukmerpos,akmers,0x0Fu,step);
//    ti=ti*0.99+t1.lap()*0.01;
  
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
//    ts=ts*0.99+t1.lap()*0.01;

    eintarray seqids;
    for (int l=0; l<best.size() && l<topotus; ++l){
      int ibest=best[l];
      if (ibest<otus.size()){
  //      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest][0]) << "\t" << idcount[ibest] << endl;
        for (int l2=0; l2<otus[ibest].size() && (otulim==0 || l2<otulim); ++l2)
          seqids.add(otus[ibest][l2]);
      }else{
  //      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest-otus.size()][0]) << "\t" << idcount[ibest] << " reversed" << endl;
        ibest-=otus.size();
        for (int l2=0; l2<otus[ibest].size() && (otulim==0 || l2<otulim); ++l2)
          seqids.add(seqs.size()+otus[ibest][l2]);
      }
    }

    eintarray otureps;
  
    // add representatives from all non-chosen OTUS (may improve confidence and novel otu estimation)
    for (int l=topotus; l<topotus+3 && l<best.size(); ++l){
      int ibest=best[l];
      if (ibest<otus.size()){
        if (otus[ibest].size()==0) continue;
        otureps.add(otus[ibest][0]);
      } else {
  //      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest-otus.size()][0]) << "\t" << idcount[ibest] << " reversed" << endl;
        ibest-=otus.size();
        if (otus[ibest].size()==0) continue;
        otureps.add(seqs.size()+otus[ibest][0]);
      }
    }
  
//    cout << "# 2nd counting" << endl;

    randomize(sws.rng,seqids);
    sws.idcount.init(seqids.size(),0);
  
    if (sws.maskid+1u<sws.maskid){
      sws.maskid=1u;
      sws.kmermask.init(KMERSIZE,0u);
    }
    ++sws.maskid;
  
//    ti2=ti2*0.99+t1.lap()*0.01;
  //  memset(sws.kmerbitmask,0,(MAXSIZE2/64+1)*sizeof(uint64_t));
  //  setkmermask(sws.kmerbitmask,s,akmers,0xFul);
  //  kmercount_mask(db.seqs,seqids,sws.kmerbitmask,sws.maskid,sws.idcount,akmers,0xFul);
//    cout << "# 2nd counting -- setkmermask" << endl;
    setkmermask(sws.kmermask,s,sws.maskid,akmers,0xFul,s_start,s_end);
//    cout << "# 2nd counting -- kmercount_mask" << endl;
    kmercount_mask(seqs,seqids,sws.kmermask,sws.maskid,sws.idcount);
//    ts2=ts2*0.99+t1.lap()*0.01;
  
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

    for (int i=0; i<otureps.size(); ++i){
      best.add(seqids.size());
      seqids.add(otureps[i]);
    }
//    best+=otureps;
   
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
      int sbest=seqids[best[l]]%seqs.size();
      eseq &sdb(seqs.values(sbest));
      if (sws.offset2+(unsigned int)(seqs.values(sbest).seqlen)<sws.offset2){  // need an unsigned int here otherwise the comparison is made in long and the offset is not correctly reset, signed int overflows are undefined so this cannot be done with signed ints either
        sws.offset2=1u;
        sws.kmerpos2.init(MAXSIZE,0);
      }
      setkmerpos(sws.kmerpos2,seqs.values(sbest),sws.offset2);
  
      ealigndata adata;
      adata.seqid=sbest;
      adata.revcompl=(seqids[best[l]]>=seqs.size());
      adata.kmercount=sws.idcount[best[l]];
    
      if (seqids[best[l]]<seqs.size()){
//        cout << "# forward align" << endl;
        seqident_global(str2id,seqs.keys(sbest),s,sws.kmerpos,sws.offset,sdb,adata,sws,as,s_start,s_end);
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
        seqident_global(str2id,seqs.keys(sbest),srev,sws.kmerposrev,sws.offset,sdb,adata,sws,as);
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
      adata.globalalign(as);
  //    LDEBUG(D_SEQALIGNMENT,print_tncount(&tncounts[l*NCOUNT_MAXLEN],0,s.seqlen));
      if (adata.matches+adata.mismatches>0 && adata.score()>=minscore){
//        adata._eval=sdb.seqlen*exp(-lambda*adata.score()); // for K-A stats we need (*s.seqlen) but this is constant
        pinfo.matchcounts.add(adata);
      }
  
      sws.offset2+=seqs.values(sbest).seqlen;
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
//  ta=ta*0.99+t1.lap()*0.01;
}


void eseqdb::seqalign_global(const estr& str2id,eseq& s,earray<epredinfo>& previnfoarr,earray<epredinfo>& pinfoarr,esearchws& sws)
{
  pinfoarr.clear();

  epredinfo pinfo;
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

  for (int l=0; l<previnfoarr[0].matchcounts.size(); ++l){
//    epinfo &pinfo(previnfoarr[l]);
    ealigndata &prevadata(previnfoarr[0].matchcounts[l]);

    int sbest=prevadata.seqid;
    eseq &sdb(seqs.values(sbest));
    if (sws.offset2+(unsigned int)(seqs.values(sbest).seqlen)<sws.offset2){  // need an unsigned int here otherwise the comparison is made in long and the offset is not correctly reset, signed int overflows are undefined so this cannot be done with signed ints either
      sws.offset2=1u;
      sws.kmerpos2.init(MAXSIZE,0);
    }
    setkmerpos(sws.kmerpos2,seqs.values(sbest),sws.offset2);

    ealigndata adata;
    adata.seqid=sbest;
    adata.revcompl=prevadata.revcompl;
//    adata.kmercount=sws.idcount[best[l]];
  
    if (!prevadata.revcompl){
      seqident_global(str2id,seqs.keys(sbest),s,sws.kmerpos,sws.offset,sdb,adata,sws,as,s_start,s_end);
    }else{
      seqident_global(str2id,seqs.keys(sbest),srev,sws.kmerposrev,sws.offset,sdb,adata,sws,as);
      int tmp=srev.seqlen-adata.s1+s_start; adata.s1=srev.seqlen-adata.e1+s_start; adata.e1=tmp; 
      adata.revcompl=true;
    }
    adata.globalalign(as);
//    if (adata.matches+adata.mismatches>0 && adata.score()>=minscore){
      pinfo.matchcounts.add(adata);
//    }

    sws.offset2+=seqs.values(sbest).seqlen;
  }
  sws.offset+=s_len;

//  if (pinfo.matchcounts.size()>0){
//    heapsort(pinfo.matchcounts);
    pinfo.tophit=pinfo.matchcounts[pinfo.matchcounts.size()-1];
    pinfo.seqid=pinfo.tophit.seqid;
    pinfoarr.add(pinfo);
//  }
}


void eseqdb::seqsearchpair(const estr& id,eseq& s,eseq& srev2,earray<epredinfo>& pinfoarr,esearchws& sws)
{
  //t1.reset();
  pinfoarr.clear();
//  memset(sws.bitmask,0,((db.otus.size()*2)/64+1)*sizeof(uint64_t));
//  kmercount_both(db.otus.size(),db.otukmers,s,sws.idcount2,sws.bitmask,akmers,0x0Fu);
//  kmercount_both2(db.otus.size(),db.otukmers,s,sws.idcount,sws.otukmerpos,akmers,0x0Fu);
//  kmercount_both_nopos2(db.otus.size(),db.otukmers,s,sws.idcount,sws.otukmerpos,akmers,0x0Fu);
  ldieif(s.seqlen>SEQSEGSIZE || srev2.seqlen>SEQSEGSIZE,"paired end reads longer than "+estr(SEQSEGSIZE)+" not supported: "+id);

  long segi=0;
    epredinfo pinfo;
    long s_start=0 , s_start2=0;
    long s_end=s.seqlen,  s_end2=srev2.seqlen;
    long s_len=s.seqlen,  s_len2=srev2.seqlen;
    
//    int step=1;
    int step=(s_len+s_len2)/150+1;
//    if (s_len<150) step=
    eseq srev,s2;
    srev.setrevcompl(s,s_start,s_end);
    s2.setrevcompl(srev2,s_start2,s_end2);
//    cout << "# " << str2id << " start: " << s_start << " end: " << s_end << " len: " << s_len << " seqlen: " << s.seqlen << " srev: " << srev.seqlen << endl;
  
    sws.otukmerpos.init(otus.size()*2,0);
    sws.idcount.init(otus.size()*2,0);

//    cout << "# counting" << endl;
    kmercount_both_nopos2_skip(otus.size(),otukmers,s,s_start,s_end,sws.idcount,sws.otukmerpos,akmers,0x0Fu,step);
    kmercount_both_nopos2_skip(otus.size(),otukmers,s2,s_start2,s_end2,sws.idcount,sws.otukmerpos,akmers,0x0Fu,step);
//    ti=ti*0.99+t1.lap()*0.01;
  
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
    if (best.size()==0)
      return;
//    ts=ts*0.99+t1.lap()*0.01;

    eintarray seqids;
    for (int l=0; l<best.size() && l<topotus; ++l){
      int ibest=best[l];
      if (ibest<otus.size()){
  //      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest][0]) << "\t" << idcount[ibest] << endl;
        for (int l2=0; l2<otus[ibest].size() && (otulim==0 || l2<otulim); ++l2)
          seqids.add(otus[ibest][l2]);
      }else{
  //      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest-otus.size()][0]) << "\t" << idcount[ibest] << " reversed" << endl;
        ibest-=otus.size();
        for (int l2=0; l2<otus[ibest].size() && (otulim==0 || l2<otulim); ++l2)
          seqids.add(seqs.size()+otus[ibest][l2]);
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
  
//    ti2=ti2*0.99+t1.lap()*0.01;
  //  memset(sws.kmerbitmask,0,(MAXSIZE2/64+1)*sizeof(uint64_t));
  //  setkmermask(sws.kmerbitmask,s,akmers,0xFul);
  //  kmercount_mask(db.seqs,seqids,sws.kmerbitmask,sws.maskid,sws.idcount,akmers,0xFul);
//    cout << "# 2nd counting -- setkmermask" << endl;
    setkmermask(sws.kmermask,s,sws.maskid,akmers,0xFul,s_start,s_end);
//    cout << "# 2nd counting -- kmercount_mask" << endl;
    kmercount_mask(seqs,seqids,sws.kmermask,sws.maskid,sws.idcount);

    if (sws.maskid+1u<sws.maskid){
      sws.maskid=1u;
      sws.kmermask.init(KMERSIZE,0u);
    }
    ++sws.maskid;

    setkmermask(sws.kmermask,s2,sws.maskid,akmers,0xFul,s_start2,s_end2);
//    cout << "# 2nd counting -- kmercount_mask" << endl;
    kmercount_mask(seqs,seqids,sws.kmermask,sws.maskid,sws.idcount);

 
//    ts2=ts2*0.99+t1.lap()*0.01;
  
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
    if (best.size()==0)
      return;

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
   
    if (sws.offset3+(unsigned int)(s_len2)<sws.offset3){ // need an int here otherwise the comparison is made in long and the offset is not correctly reset
      sws.offset3=1u;
      sws.kmerpos3.init(MAXSIZE,0u);
      sws.kmerposrev3.init(MAXSIZE,0u);
    }
//    cout << "# 2nd counting -- setkmerpos" << endl;
    setkmerpos(sws.kmerpos3,s2,sws.offset3,s_start2,s_end2);
//    cout << "# 2nd counting -- setkmerposrev" << endl;
    setkmerpos(sws.kmerposrev3,srev2,sws.offset3);
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
      int sbest=seqids[best[l]]%seqs.size();
      eseq &sdb(seqs.values(sbest));
      if (sws.offset2+(unsigned int)(seqs.values(sbest).seqlen)<sws.offset2){  // need an unsigned int here otherwise the comparison is made in long and the offset is not correctly reset, signed int overflows are undefined so this cannot be done with signed ints either
        sws.offset2=1u;
        sws.kmerpos2.init(MAXSIZE,0);
      }
      setkmerpos(sws.kmerpos2,seqs.values(sbest),sws.offset2);
  
      ealigndata adata,adata1,adata2;
      adata.seqid=sbest;
      adata.revcompl=(seqids[best[l]]>=seqs.size());
      adata.kmercount=sws.idcount[best[l]];
    
      if (seqids[best[l]]<seqs.size()){
        seqident_local_leftext(id,seqs.keys(sbest),s,sws.kmerpos,sws.offset,sdb,adata1,sws,as,s_start,s_end);
        seqident_local_rightext(id,seqs.keys(sbest),s2,sws.kmerpos3,sws.offset3,sdb,adata2,sws,as,s_start2,s_end2);

//        seqident_local(id,db.seqs.keys(sbest),s,sws.kmerpos,sdb,adata1,sws,as,s_start,s_end);
//        seqident_local(id,db.seqs.keys(sbest),s2,sws.kmerpos3,sdb,adata2,sws,as,s_start2,s_end2);
      }else{
        seqident_local_rightext(id,seqs.keys(sbest),srev,sws.kmerposrev,sws.offset,sdb,adata1,sws,as);
        seqident_local_leftext(id,seqs.keys(sbest),srev2,sws.kmerposrev3,sws.offset3,sdb,adata2,sws,as);
//        seqident_local(id,db.seqs.keys(sbest),srev,sws.kmerposrev,sdb,adata1,sws,as);
//        seqident_local(id,db.seqs.keys(sbest),srev2,sws.kmerposrev3,sdb,adata2,sws,as);

        // flip 
        int tmp;
        tmp=srev.seqlen-adata1.s1+s_start; adata1.s1=srev.seqlen-adata1.e1+s_start; adata1.e1=tmp; 
        tmp=srev2.seqlen-adata2.s1+s_start2; adata2.s1=srev2.seqlen-adata2.e1+s_start2; adata2.e1=tmp; 
        adata.revcompl=true;
      }
      adata.s1=adata1.s1<adata2.s1?adata1.s1:adata2.s1;
      adata.s2=adata1.s2<adata2.s2?adata1.s2:adata2.s2;
      adata.e1=adata1.e1>adata2.e1?adata1.e1:adata2.e1;
      adata.e2=adata1.e2>adata2.e2?adata1.e2:adata2.e2;
      adata.matches=adata1.matches+adata2.matches;
      adata.mismatches=adata1.mismatches+adata2.mismatches;
      adata.gaps=adata1.gaps+adata2.gaps;
      adata._score=adata1._score+adata2._score;
//      cout << "pair end: " << adata1.s2 << " " << adata1.e2 << " - " << adata2.s2 << " " << adata2.e2 << " " << (adata.revcompl?"-":"+") << " " << seqids[best[l]] << " " << seqs.size() << endl;
      if (adata1.e2 > adata2.s2 && adata2.e2 > adata1.s2){ // some overlap, need to subtract overlap scores and stats
//        lwarn("overlap: "+estr(adata1.s2)+" "+adata1.e2+" - "+adata2.s2+" "+adata2.e2);
        int tmpscore,tmpmatches,tmpmismatches,tmpgaps;
        if (!adata.revcompl && adata2.s2 < adata1.s2 || adata.revcompl && adata1.s2 < adata2.s2){ // incorrect order of paired ends, take highest scoring of both alignments
//          lerror("2nd pair end before 1st?");
          adata._score=0; // this is not a correct match, skip this alignment
          if (adata1._score>adata2._score){
            adata.s1=adata1.s1;
            adata.s2=adata1.s2;
            adata.e1=adata1.e1;
            adata.e2=adata1.e2;
            adata.matches=adata1.matches;
            adata.mismatches=adata1.mismatches;
            adata.gaps=adata1.gaps;
            adata._score=adata1._score;
          }else{
            adata.s1=adata2.s1;
            adata.s2=adata2.s2;
            adata.e1=adata2.e1;
            adata.e2=adata2.e2;
            adata.matches=adata2.matches;
            adata.mismatches=adata2.mismatches;
            adata.gaps=adata2.gaps;
            adata._score=adata2._score;
          }
        }else{
          if (!adata.revcompl)
            adata2.partscore(adata2.s2-adata1.s2+1,adata2.e1,tmpscore,tmpmatches,tmpmismatches,tmpgaps,as);
          else
            adata1.partscore(adata1.s2-adata2.s2+1,adata1.e1,tmpscore,tmpmatches,tmpmismatches,tmpgaps,as);

          adata._score-=tmpscore;
          adata.matches-=tmpmatches;
          adata.mismatches-=tmpmismatches;
          adata.gaps-=tmpgaps;
        }
      }

  //    LDEBUG(D_SEQALIGNMENT,print_tncount(&tncounts[l*NCOUNT_MAXLEN],0,s.seqlen));
      if (adata.matches+adata.mismatches>0 && adata.score()>=minscore){
//        adata._eval=sdb.seqlen*exp(-lambda*adata.score()); // for K-A stats we need (*s.seqlen) but this is constant
        pinfo.matchcounts.add(adata);
      }
  
      sws.offset2+=seqs.values(sbest).seqlen;
    }
    sws.offset+=s_len;
    sws.offset3+=s_len2;
 
    if (pinfo.matchcounts.size()>0){
      heapsort(pinfo.matchcounts);
    //  for (int i=0; i<pinfo.matchcounts.size(); ++i)
    //    cout << "#best: " << i << " " << seqs.keys(pinfo.matchcounts[i].seqid) << " " << pinfo.matchcounts[i].score() << endl;
      pinfo.tophit=pinfo.matchcounts[pinfo.matchcounts.size()-1];
      pinfo.seqid=pinfo.tophit.seqid;
      pinfoarr.add(pinfo);
    }
//  ta=ta*0.99+t1.lap()*0.01;
}



void eseqdb::pseqsearch(const estr& str2id,eseq& s,earray<epredinfo>& pinfoarr,esearchws& sws)
{
//  t1.reset();
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
    eseq srev;
    srev.setrevcompl(s,s_start,s_end);
//    cout << "# " << str2id << " start: " << s_start << " end: " << s_end << " len: " << s_len << " seqlen: " << s.seqlen << " srev: " << srev.seqlen << endl;
  
//    sws.otukmerpos.init(db.otus.size()*6,0);
    sws.idcount.init(otus.size()*6,0);

//    cout << "# counting" << endl;
//    kmercount_both_nopos2_skip(db.otus.size(),db.otukmers,s,s_start,s_end,sws.idcount,sws.otukmerpos,akmers,0x0Fu,s_len/100+1);
    kmercount_both_prot_skip(otus.size(),otukmers,s,s_start,s_end,sws.idcount,sws.otukmerpos,akmers,0x0Fu,1);
//    ti=ti*0.99+t1.lap()*0.01;
  
    eintarray best;
  //  ebasicarray<ealigndata> matchcounts;
  
    // choosing sequences for kmercounting step
  //  eintarray& best(sws.best);
  //  eintarray bestcount;
  //  best.clear();
    int l;
    for (l=0; l<sws.idcount.size(); ++l){
  //    uint16_t bt=(sws.bitmask[l/64u]>>(l%64u))&0x1ul;
      if (sws.idcount[l]>0) { best.add(l); break; }
    }
    for (; l<sws.idcount.size(); ++l){
  //    uint16_t bt=(sws.bitmask[l/64u]>>(l%64u))&0x1ul;
      if (sws.idcount[l]==0 || (sws.idcount[l]<sws.idcount[best[best.size()-1]] && best.size()==3*topotus)) continue; // worst id than the bottom of the list
  
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
//    ts=ts*0.99+t1.lap()*0.01;
    for (int ti=0; ti<3 && ti<best.size(); ++ti)
      cout << ti << " " << str2id << " score: " << sws.idcount[best[ti]] << " best: "<< seqs.keys(best[ti]%otus.size()) << " " << best[ti]/otus.size() << endl;
//    continue;

    eintarray seqids; 
    for (int l=0; l<best.size() && l<topotus; ++l){
      int ibest=best[l];
//      if ((ibest/db.otus.size())%2==0){
  //      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest][0]) << "\t" << idcount[ibest] << endl;
        for (int l2=0; l2<otus[ibest%otus.size()].size(); ++l2)
          seqids.add(otus[ibest%otus.size()][l2]+seqs.size()*(ibest/otus.size()));
//      }else{
//  //      cout << "#\t" << str2id << "\t" << seqs.keys(otus[ibest-otus.size()][0]) << "\t" << idcount[ibest] << " reversed" << endl;
//        ibest=ibest%db.otus.size();
//        for (int l2=0; l2<db.otus[ibest].size(); ++l2)
//          seqids.add(db.seqs.size()+db.otus[ibest][l2]);
//      }
    }
  
 
//    cout << "# 2nd counting" << endl;
    randomize(sws.rng,seqids);

/*
    sws.idcount.init(seqids.size(),0);
  
    if (sws.maskid+1<sws.maskid){
      sws.maskid=-1;
      sws.kmermask.init(KMERSIZE,-1);
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
  
 
    best.clear();
    for (l=0; l<sws.idcount.size(); ++l){
      if (sws.idcount[l]>0){ best.add(l); break; }
    }
    for (; l<sws.idcount.size(); ++l){
  //    if (idcount[l]<idcount[best[best.size()-1]] && (idcount[l]<0.8*idcount[ibest] && best.size()>=20 || best.size()==100)) continue; // worse id than the bottom of the list
      if (sws.idcount[l]==0 || (sws.idcount[l]<sws.idcount[best[best.size()-1]] && best.size()>=tophits)) continue; // zero counts or worse id than the bottom of the list
  
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
  
    if (sws.offset+s_len<sws.offset){
      sws.offset=0;
      sws.kmerpos.init(PMAXSIZE,-1);
      sws.kmerposrev.init(PMAXSIZE,-1);
    }
//    cout << "# 2nd counting -- setkmerpos" << endl;
//    cout << "# 2nd counting -- setkmerposrev" << endl;

//    for (int l=0; l<6; ++l)
//    setpkmerpos(sws.kmerpos,s,sws.offset,s_start,s_end);


  
  //  eintarray seqboth;
  //  seqboth.init(db.seqs.size(),-1);
//    cout << "# aligning" << endl;
*/
  
    for (int l=0; l<seqids.size(); ++l){
      int sbest=seqids[l]%seqs.size();
      eseq &sdb(seqs.values(sbest));
      if (sws.offset2+(unsigned int)(seqs.values(sbest).seqlen)<sws.offset2){
        sws.offset2=1u;
        sws.kmerpos2.init(PMAXSIZE,0);
      }
      setpkmerpos(sws.kmerpos2,sdb,sws.offset2);
  
      ealigndata adata;
      adata.seqid=sbest;
      adata.revcompl=(seqids[l]/seqs.size()%2==1);
//      adata.kmercount=sws.idcount[best[l]];
      
//      cout << "s_start: " << s_start << " " << s_start+seqids[l]/db.seqs.size()/2 << endl;
      if ((seqids[l]/seqs.size())%2==0){
        cout << "# forward align: " << seqs.keys(sbest) << " slen: " << sdb.seqlen << endl;
        pseqident_local(s,sws.kmerpos,sws.offset,sdb,adata,sws,pas,s_start+seqids[l]/seqs.size()/2,s_end);
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
        cout << "# reverse align: " << seqs.keys(sbest) << " slen: " << sdb.seqlen << endl;
//        pseqident_local(srev,sws.kmerposrev,sdb,adata,sws,as,s_start+seqids[l]/db.seqs.size()/2,s_end);
        int fs=(s_end - seqids[l]/seqs.size()/2)%3;
        pseqident_local(srev,sws.kmerposrev,sws.offset,sdb,adata,sws,pas,fs,srev.seqlen);
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
      if (adata.matches+adata.mismatches>0 && adata.score()>=20){
        pinfo.matchcounts.add(adata);
//        adata.profile.inv();
//        cout << "l: " << l << endl << adata.palign_str(s,db.seqs.values(seqids[l]%db.seqs.size())) << endl;
      }
  
      sws.offset2+=sdb.seqlen;
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
//  ta=ta*0.99+t1.lap()*0.01;
}


void eseqdb::loadTaxFormat1(efile& f,etax& tax)
{
  estr line;
  estrarray parts,parts2;
  while (!f.eof()){
    if (line.len()==0) continue;
    if (line[0]=='#'){
      f.readln(line);
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

/*  simple taxonomy */
    eseqtax *newtax=new eseqtax();
    newtax->tl.reserve(parts.size()-1);
    for (int i=1; i<parts.size(); ++i)
      newtax->tl.add(eseqtaxlevel(parts[i].i()));

/* mseq data with predicted confidences        
    eseqtax *newtax=new eseqtax(parts[1].f());
    newtax->tl.reserve((parts.size()-2)/2);
    for (int i=2; i<parts.size(); i+=2)
      newtax->tl.add(eseqtaxlevel(parts[i].i(),parts[i+1].f()));
*/
    int si=parts[0].i();
    ldieif(si>=seqs.size(),"sequence number in taxonomy file does not exist in database, please make sure to recreate the taxonomy file everytime the fasta database changes");
    tax.seqs[si]=newtax;
  }
}


void eseqdb::loadTaxonomy(const estr& fname)
{
  efile f;
  estr line;
  estrarray parts,parts2;
  int notfound=0;
  f.open(fname,"r");
  etax& tax(taxa.add(etax()));
  tax.seqs.init(seqs.size(),0x00);
  int taxind=0;
  int taxlevels=-1;
  while (!f.eof()){
    f.readln(line);
    if (line.len()==0) continue;
    if (line[0]=='#'){
      if (line=="#taxformat1"){
        loadTaxFormat1(f,tax);
        break;
      }
      parts=line.explode(" ");
      if (parts[0]=="#cutoff:"){
        ldieif(tax.cutoff.size()>0,"duplicate #cutoff lines!!");
        taxlevels=parts.size()-1;
        for (int i=1; i<parts.size(); ++i){
          parts2=parts[i].explode(":");
          ldieif(parts2.size()<2,"not enough parts on #cutoff line, i.e.: 0.97:0.02");
          tax.cutoff.add(parts2[0].f());
          tax.cutoffcoef.add(parts2[1].f());
        }
      }else if (parts[0]=="#sweight:" && parts.size()>1){
        if (parts.size()>1)
           swmin=parts[1].f();
        if (parts.size()>2)
           swmax=parts[2].f();
      }else if (parts[0]=="#name:" && parts.size()>1){
        tax.name=parts[1];
      }else if (parts[0]=="#levels:" && parts.size()>1){
        ldieif(taxlevels>0 && taxlevels!=parts.size()-1,"Taxlevels mismatch: "+estr(taxlevels)+" != "+estr(parts.size()-1)+" in #levels cutoffs");
        for (int i=1; i<parts.size(); ++i)
          tax.levels.add(parts[i]);
      }
      continue; 
    }
    parts=line.explode("\t");
    ldieif(parts.size()<2,"loading taxonomy, not enough fields in line: "+line);
    if (seqind.exists(parts[0])){
      if (parts.size()==2 || parts.size()==3){ // simple taxonomy file
        parts2=parts[1].explode(";");
        ldieif(taxlevels>0 && taxlevels!=parts2.size(),"Taxlevels mismatch: "+estr(taxlevels)+" != "+estr(parts2.size())+" in taxonomy for sequence: "+parts[0]+" :: "+parts[1]);
        ldieif(parts2.size()==0,"empty tax: "+line);
        eseqtax *newtax=new eseqtax();
        newtax->tl.reserve(parts2.size());
        for (int i=0; i<parts2.size(); ++i){
          if (i>=tax.ind.size() || i>=tax.names.size()){
            tax.names.add(earray<estr>());
            tax.ind.add(estrhashof<int>());
          }
          if (parts2[i].len()==0)
            newtax->tl.add(eseqtaxlevel(-1));
          else {
            if (!tax.ind[i].exists(parts2[i])){
              tax.ind[i].add(parts2[i],tax.names[i].size());
              tax.names[i].add(parts2[i]);
            }
            if (parts.size()==3) // with evidence weights: SEQID \t SEQTAX \t EVIDENCEWEIGHT
              newtax->tl.add(eseqtaxlevel(tax.ind[i][parts2[i]],1.0,parts[2].f()));
            else
              newtax->tl.add(eseqtaxlevel(tax.ind[i][parts2[i]],1.0));
          }
        }
        tax.seqs[seqind[parts[0]]]=newtax;
      }else if (parts.size()>11) { // indirect taxonomy file with confidences
        
        if (taxind==0){
          for (int i=1; i<parts.size(); ++i)
            if (parts[i].len()==0)
              { taxind=i+1; break; }
        }
        int dbhit=(seqind.exists(parts[1])?seqind[parts[1]]:-1);
//        lerrorif(dbhit==-1,"gold hit not found: "+parts[1]);
        eseqtax *newtax=new eseqtax(parts[3].f());
        int taxfields=(parts.size()-taxind)/3;
        newtax->tl.reserve(taxfields);
//        cout << parts[0] << " " << parts[3] << " " << taxfields;
        for (int i=0; i<taxfields; ++i){
          int field=i*3+taxind;
          if (i>=tax.ind.size() || i>=tax.names.size()){
            tax.names.add(earray<estr>());
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
        tax.seqs[seqind[parts[0]]]=newtax;
      }else
        ldie("unrecognized tax format: "+line);
    }else{
      ++notfound;
    }
  }
  if (tax.name.len()==0)
    tax.name=basename(fname);
  lwarnif(notfound>0,"loading taxonomy, "+estr(notfound)+" sequences not found in sequence database");
  f.close();
}

void taskClusterDB()
{
//  ebasicarray<ekmerarray> otukmers;
//  otukmers.init(MAXSIZE);

  eseqdb &db(*mtdata.seqdb);

  euintarray otukmersize;
  eintarray seqotu;
  euintarray idcount;
  eintarray kmerpos;
  eintarray tmpmapped;  

  idcount.reserve(db.seqs.size());
  kmerpos.reserve(db.seqs.size());

  eintarray &len_si(mtdata.len_si);

  mtdata.m.lock();
  seqotu=db.seqotu;
  mtdata.m.unlock();


  int dbotus=0;
  int ipos=0,ilen=0;
  int mapped=0;
  while (1) {
    mtdata.m.lock();
    while (mtdata.ilen==0) mtdata.seqsSignal.wait(mtdata.m);
    if (mtdata.ilen==-1) { mtdata.m.unlock(); break; }
    ipos=mtdata.ipos+mtdata.ilen*mtdata.threadi;
    ilen=mtdata.ilen;
    dbotus=mtdata.dbotus;

    int totalt=(db.seqs.size()-mtdata.ipos)/mtdata.ilen;
    if (totalt==0) break;
    if (mtdata.threadi+1==totalt){
      cout << endl << mtdata.threadi << " " << ipos << " " << ipos+ilen << " " << totalt << endl;
    }
    ldieif(ipos<0,"ipos is negative");
  
    mtdata.threadi=(mtdata.threadi+1)%totalt;
    otukmersize=mtdata.otukmersize;
/*
    for (int i=0; i<db.otukmers.size(); ++i){
      for (int j=db.otukmers[i].size(); j<db.otukmers[i].size(); ++j)
        otukmers[i].add(db.otukmers[i][j]);
    }
*/
    for (int i=mapped; i<mtdata.mapped.size(); ++i)
      seqotu[mtdata.mapped[i]]=db.seqotu[mtdata.mapped[i]];
    mapped=mtdata.mapped.size();
//    cout << endl << mtdata.threadi << " " << ipos << " " << ipos+ilen << endl;
    mtdata.m.unlock();
    for (int i=ipos; i<ipos+ilen && i<len_si.size(); ++i){
      if (seqotu[len_si[i]]!=-1) continue;

//      cout << endl << i << endl;
      eseq& s(db.seqs.values(len_si[i]));
      idcount.init(dbotus,0);
      kmerpos.init(dbotus,0u);

      kmercount_single(db.otukmers,s,idcount,kmerpos,otukmersize);

      int ibest=0;
      for (int l=1; l<idcount.size(); ++l){
        if (idcount[l]>idcount[ibest]) ibest=l;
      }
      if (float(idcount[ibest])/s.seqlen >= 0.80){
        tmpmapped.add(len_si[i]);
        seqotu[len_si[i]]=ibest;
      }
    }
    mtdata.m.lock();
//    cout << endl << "mapped: " << tmpmapped.size() << endl;
    for (int i=mapped; i<mtdata.mapped.size(); ++i)
      seqotu[mtdata.mapped[i]]=db.seqotu[mtdata.mapped[i]];
    for (int i=0; i<tmpmapped.size(); ++i){
      db.seqotu[tmpmapped[i]]=seqotu[tmpmapped[i]];
      mtdata.mapped.add(tmpmapped[i]);
    }
    mapped=mtdata.mapped.size();
    mtdata.m.unlock();
    tmpmapped.clear();
   }
}


void eseqdb::makeClusterMT(ethreads& t) //,const estr& cfile)
{
  euintarray kmermask;

  kmermask.init(MAXSIZE,0);
//  eseqdb &db(*mtdata.seqdb);
  mtdata.ipos=0;
  mtdata.ilen=0;
  mtdata.dbotus=0;
  mtdata.threadi=0;
  mtdata.otukmersize.init(MAXSIZE,0);
//  db.otukmers.init(MAXSIZE);
  for (uint32_t i=0; i<MAXSIZE; ++i){
    if (akmers[i&0xFu]==0u) continue;
//    otukmers[i].reserve(seqs.size());
  }

  eintarray &len_si(mtdata.len_si);

  len_si=iheapsort(seqs);
  otus.reserve(seqs.size()); // worst case scenario
  otus.add(eintarray(len_si[0]));
  seqotu.init(seqs.size(),-1);
  seqotu[len_si[0]]=0;

  eintarray tmpkmers;
  euintarray idcount;
  eintarray kmerpos;

  tmpkmers.init(MAXSIZE,-1);
  otukmeradd(otukmers,0,seqs.values(otus[0][0]),tmpkmers,0,akmers,0xFul);
  idcount.reserve(seqs.size());
  kmerpos.reserve(seqs.size());
  idcount.init(otus.size(),0u);
  kmerpos.init(otus.size(),0u);
 
  const int seqcount=10000;
  const int updcount=100;
  eintarray tmpseqotu;
  tmpseqotu=seqotu;

//  ethreads t;
//  t.run(taskClusterDB,evararray(),nthreads);
  t.run(taskClusterDB,evararray());

  int mapped=0;

  for (int i=0; i<seqs.size(); ++i){
    if (i%updcount==0 && i>=1000) {
      mtdata.m.lock();
      if (i==1000){
        int count=0;
        for (int l=0; l<otukmers.size(); ++l){
          if (otukmers[l].size()>0.5*otus.size()){
            kmermask[l]=1;
            otukmers[l].clear();
            ++count;
          }
        }
        cout << "fcount: " << count << " db.otus: " << otus.size() << endl;
      }
      mtdata.dbotus=otus.size();
      for (int l=0; l<otukmers.size(); ++l)
        mtdata.otukmersize[l]=otukmers[l].size();

      for (int l=mapped; l<mtdata.mapped.size(); ++l)
        tmpseqotu[mtdata.mapped[l]]=seqotu[mtdata.mapped[l]];
      mapped=mtdata.mapped.size();
      mtdata.ipos=i+seqcount*2;
      mtdata.ilen=seqcount;
      mtdata.seqsSignal.broadcast();
      fprintf(stderr,"\rseq: %i otus: %li mapped: %li",i,otus.size(),mtdata.mapped.size());
      mtdata.m.unlock();
    }
    if (tmpseqotu[len_si[i]]!=-1) continue;

    eseq& s(seqs.values(len_si[i]));
    idcount.init(otus.size(),0);
    kmerpos.init(otus.size(),0u);

    kmercount_single(otukmers,s,idcount,kmerpos);

    int ibest=0;
    for (int l=1; l<idcount.size(); ++l){
      if (idcount[l]>idcount[ibest]) ibest=l;
    }
    if (float(idcount[ibest])/s.seqlen >= 0.80) {
    } else {
      otus.add(eintarray(len_si[i]));
      ibest=otus.size()-1;
      tmpseqotu[len_si[i]]=ibest;
//      if (i<1000)
      otukmeradd(otukmers,ibest,s,tmpkmers,ibest,akmers,0xFul,kmermask);
    }
  }

  t.wait();

  fprintf(stderr,"\rseq: %li clusters: %li\n",seqs.size(),(long)otus.size());
/*
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
*/
 
}


void eseqdb::loadSequencesBinary(const estr& filename)
{
  efile f;
  estr tmpstr;
  f.open(filename,"r");
  long si=0,tsi;
  unsigned int ver;
  const int bufsize=10000;
  f.read(tmpstr,bufsize);
  si=unserialuint(ver,tmpstr,si);
  ldieif(si==-1,"loading binary database");
  ldieif(ver!=0,"unknown database version: "+estr(ver));
  eseq s;
  estr sid;
  while (!f.eof() || si < tmpstr.len()){
    tsi=sid.unserial(tmpstr,si);
    ldieif(tsi==-1 && f.eof(),"Unexpected end of file");
    if (tsi==-1){
      tmpstr=tmpstr.substr(si);
      f.read(tmpstr,bufsize);
      si=0;
      tsi=sid.unserial(tmpstr,si);
      ldieif(tsi==-1,"Buffer not long enough: "+estr(seqs.size())+" "+estr(tmpstr.len()));
    }
    si=tsi;

    tsi=s.unserial(tmpstr,si);
    ldieif(tsi==-1 && f.eof(),"Unexpected end of file");
    if (tsi==-1){
      tmpstr=tmpstr.substr(si);
      f.read(tmpstr,bufsize);
      si=0;
      tsi=s.unserial(tmpstr,si);
      ldieif(tsi==-1,"Buffer not long enough: "+estr(seqs.size())+" "+estr(tmpstr.len()));
    }
    si=tsi;
    seqs.add(sid,s);
  }
}


void eseqdb::saveSequences(const estr& filename)
{
  efile f;
  estr tmpstr;
  f.open(filename,"w");
  serialuint(0,tmpstr);
  f.write(tmpstr);
  tmpstr.clear();
  for (int i=0; i<seqs.size(); ++i){
    seqs.keys(i).serial(tmpstr);
    seqs.values(i).serial(tmpstr);
    f.write(tmpstr);
    tmpstr.clear();
  }
  f.close();
}


void eseqdb::loadSequences(const estr& dbfile)
{
  estr str2id,str2seq,line;
  efile f;
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

    seqind.add(str2id,seqs.size());
    seqs.add(str2id,eseq(str2seq));
  }
  f.close();
}

void eseqdb::init(const estr& dbfile,bool nocluster,outfmt_fd ofmt)
{
  mtdata.outfmt=ofmt;

  seqotu.init(seqs.size(),-1);
  otukmers.init(MAXSIZE);
  estr cfile=dbfile+".mscluster";
  if (nocluster){
    eintarray tmpkmers;
    otus.reserve(seqs.size());
    otus.add(eintarray(0));
    tmpkmers.init(MAXSIZE,-1);
    otukmeradd(otukmers,0,seqs.values(otus[0][0]),tmpkmers,0,akmers,0xFul);
    for (int i=1; i<seqs.size(); ++i){
      eseq& s(seqs.values(i));
      otus.add(eintarray(i));
      int ibest=otus.size()-1;
      seqotu[i]=ibest;
      otukmeradd(otukmers,ibest,s,tmpkmers,ibest,akmers,0xFul);
    }
//  } else if (dbfilter.len()==0 && efile(cfile).exists()){ // load clustering
  } else if (efile(cfile).exists()){ // load clustering
    loadCluster(cfile);
    ldieif(otus.size()==0,"no clusters in database");
  }else{ // perform clustering
    fprintf(stderr,"# no clustering file found, performing clustering\n");
    makeCluster(cfile);
  }

  int fcount=0;
  int okmercount=0;
  for (int i=0; i<otukmers.size(); ++i){
    if (otukmers[i].size()>0)
      ++okmercount;
    if (otukmers[i].size()>1000 && otukmers[i].size()>otus.size()*0.50){
      otukmers[i].clear();
      ++fcount;
    }
  }
}

void eseqdb::makeCluster(const estr& cfile)
{
    // phase 1: find all cluster seeds, but do not add non-seeds in first phase
  eintarray tmpkmers;
  euintarray idcount;
  eintarray kmerpos;
  int i;

  eintarray len_si(iheapsort(seqs));
  otus.reserve(seqs.size()); // worst case scenario
  otus.add(eintarray(len_si[0]));
  seqotu[len_si[0]]=0;

  tmpkmers.init(MAXSIZE,-1);
  otukmeradd(otukmers,0,seqs.values(otus[0][0]),tmpkmers,0,akmers,0xFul);
  idcount.reserve(seqs.size());
  kmerpos.reserve(seqs.size());
  idcount.init(otus.size(),0u);
  kmerpos.init(otus.size(),0u);
  
  for (i=1; i<seqs.size(); ++i){
    eseq& s(seqs.values(len_si[i]));
    if (i%100==0) 
      fprintf(stderr,"\rseq: %i otus: %li",i,otus.size());
    idcount.init(otus.size(),0);
    kmerpos.init(otus.size(),0u);

    kmercount_single(otukmers,s,idcount,kmerpos);

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
      otus.add(eintarray(len_si[i]));
//        otukmercount.add(kmercount);
      ibest=otus.size()-1;
      seqotu[len_si[i]]=ibest;
      otukmeradd(otukmers,ibest,s,tmpkmers,ibest,akmers,0xFul);
    }
//      ++sp;
  }
  fprintf(stderr,"\rseq: %i clusters: %li\n",i,otus.size());

  fprintf(stderr,"phase 2:\n");
  // phase 2: after finding all seeds, add all non-seeds to the most similar clusters
  for (i=0; i<seqotu.size(); ++i){
    eseq& s(seqs.values(i));
    if (i%100==0)
      fprintf(stderr,"\rseq: %i",i);
    if (seqotu[i]!=-1) continue;

//      if (sp+1<sp){
      idcount.init(otus.size(),0);
//        sp=0;
//      }
//      if (kp+s.seqlen<kp){
      kmerpos.init(otus.size(),0u);
//        kp=8u;
//      }

//      kmercount_single3(db.otukmers,s,sp,kp,idcount,kmerpos);
//      kmercount_single2(db.otukmers,s,kp,idcount,kmerpos);
    kmercount_single(otukmers,s,idcount,kmerpos);
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
    otus[ibest].add(i);
    seqotu[i]=ibest;
//      ++sp;
  }
//    return(0);
//    db.otus.init(dbotus.size());
//    for (int i=0; i<db.seqotu.size(); ++i)
//      db.otus[db.seqotu[i]].add(i);


//  if (dbfilter.len()==0){ // do not save cluster file for a filtered db
    fprintf(stderr,"\rseq: %li clusters: %li\n",seqs.size(),(long)otus.size());
    if (cfile.len()){
      efile f(cfile,"w");
      for (int i=0; i<otus.size(); ++i){
        f.write(estr(i));
        for (int j=0; j<otus[i].size(); ++j)
          f.write(" "+estr(otus[i][j]));
        f.write("\n");
      }
      f.close();
    }
//  }
}


void eseqdb::loadCluster(const estr& cfile)
{
  efile f;
  f.open(cfile,"r");
  estr line;
  eintarray tmpkmers;
  estrarray parts;
  int i;

  tmpkmers.init(MAXSIZE,0);
  while (!f.eof() && f.readarr(line,parts)){
    if (line.len()==0) continue;
    ldieif(parts.size()<2,"less than 2 parts: "+line);

    int tmpi,repid=-1;
    for (i=1; i<parts.size(); ++i){
      tmpi=parts[i].i();
// //        if (ignseqs.size() && !ignseqs.exists(seqs.keys(tmpi)) || seqs.values(tmpi).tax!=0x00) {repid=tmpi; break;}
//        if (ignseqs.size() && ignseqs.exists(db.seqs.keys(tmpi))) continue;
      repid=tmpi;
      break;
    }
    if (repid==-1) continue; // no sequence found with taxonomic annotation

    otus.add(eintarray());
    eintarray &tmpo(otus[otus.size()-1]);
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
    eseq &s(seqs.values(tmpi));
//      kmercount_single_nopos(db.otukmers,s,bitmask,idcount,tmpkmers,db.otus.size(),kmerlen,ibest,bcount,akmers,0xFul);
//      if (idcount.size()==0){
    otukmeradd(otukmers,otus.size()-1,s,tmpkmers,otus.size(),akmers,0xFul);
//        tmpkmers.init(MAXSIZE,0);
//        kmerhash(db.otukmers,db.kmerlast,s,tmpkmers,ti++,db.seqs.size(),akmers,0xFul);
//      }else
//        otuaddkmerdiff(db.otukmers,db.seqs.values(db.otus[ibest][0]),s,tmpkmers,i,db.otus.size()-1,akmers,0xFul);
//      idcount.add(0);

      

    for (i=1; i<parts.size(); ++i){
      tmpi=parts[i].i();
//   //        if (ignseqs.size() && !ignseqs.exists(seqs.keys(tmpi)) || seqs.values(tmpi).tax==0x00) continue; // do not add sequence if there is no taxonomic annotation
//        if (ignseqs.size() && ignseqs.exists(db.seqs.keys(tmpi))) continue;
      ldieif(tmpi>=seqs.size(),"cluster file has more sequence ids than original file, please remove cluster file: "+cfile);
      tmpo.add(tmpi);
      seqotu[tmpi]=otus.size()-1;
    }
  }
}


