#!/bin/bash

gawk '
{ arr[NR]=$0; }

END {
  order["A"]=0;
  order["T"]=1;
  order["G"]=2;
  order["C"]=3;
  b2h["0000"]="0";
  b2h["0001"]="1";
  b2h["0010"]="2";
  b2h["0011"]="3";
  b2h["0100"]="4";
  b2h["0101"]="5";
  b2h["0110"]="6";
  b2h["0111"]="7";
  b2h["1000"]="8";
  b2h["1001"]="9";
  b2h["1010"]="A";
  b2h["1011"]="B";
  b2h["1100"]="C";
  b2h["1101"]="D";
  b2h["1110"]="E";
  b2h["1111"]="F";
  b2i["0000"]=0;
  b2i["0001"]=1;
  b2i["0010"]=2;
  b2i["0011"]=3;
  b2i["0100"]=4;
  b2i["0101"]=5;
  b2i["0110"]=6;
  b2i["0111"]=7;
  b2i["1000"]=8;
  b2i["1001"]=9;
  b2i["1010"]=10;
  b2i["1011"]=11;
  b2i["1100"]=12;
  b2i["1101"]=13;
  b2i["1110"]=14;
  b2i["1111"]=15;
  nuc["A"]="00";
  nuc["T"]="01";
  nuc["G"]="10";
  nuc["C"]="11";
  for (i=1; i<=length(arr["1"]); ++i){
    if (!(substr(arr[1],i,1) in aa))
      aa[substr(arr[1],i,1)]=length(aa);
    codon=substr(arr[3],i,1) substr(arr[4],i,1) substr(arr[5],i,1);
    ntrans[codon]=substr(arr[1],i,1);
    htrans[codon]="0x" b2h["00" nuc[substr(arr[5],i,1)]] b2h[nuc[substr(arr[4],i,1)] nuc[substr(arr[3],i,1)]];
    itrans[b2i["00" nuc[substr(arr[5],i,1)]]*16 + b2i[nuc[substr(arr[4],i,1)] nuc[substr(arr[3],i,1)]]]=codon;
    trans["0x" b2h["00" nuc[substr(arr[5],i,1)]] b2h[nuc[substr(arr[4],i,1)] nuc[substr(arr[3],i,1)]]]=substr(arr[1],i,1);
  }
  for (i in ntrans){
    if (!(ntrans[i] in listaa) || order[substr(i,3,1)] < order[substr(listaa[ntrans[i]],3,1)])
      listaa[ntrans[i]]=i;
  }
  for (i in listaa)
    print i,htrans[listaa[i]];
  str="";
  for (i=0; i in itrans; ++i)
    str=str "," htrans[listaa[ntrans[itrans[i]]]];
  print "unsigned char nuc2prot[]={"substr(str,2) "};";
}' transl_table

