#include "eseqdb.h"

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
        for (int i=1; i<parts.size(); ++i)
          tax.levels.add(parts[i]);
      }
      continue; 
    }
    parts=line.explode("\t");
    ldieif(parts.size()<2,"loading taxonomy, not enough fields in line: "+line);
    if (seqind.exists(parts[0])){
      if (parts.size()==2){ // simple taxonomy file
        parts2=parts[1].explode(";");
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
            newtax->tl.add(eseqtaxlevel(tax.ind[i][parts2[i]]));
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
      }
    }else{
      ++notfound;
    }
  }
  if (tax.name.len()==0)
    tax.name=basename(fname);
  lwarnif(notfound>0,"loading taxonomy, "+estr(notfound)+" sequences not found in sequence database");
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

void eseqdb::init(const estr& dbfile)
{
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


