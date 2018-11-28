#!/bin/bash

FA=example/700016012.V35.fasta.10000
DB=data/mapseqref.fasta
TAX=ncbitax
TAXFILE=$DB.$TAX
OUT=tests/$(basename $FA).$TAX.mseq

#echo ./mapseq -nthreads 1 $FA data/mapref-2.2b.fna data/mapref-2.2b.fna.ncbitax \> $OUT 2\> /dev/null
#./mapseq -nthreads 1 $FA data/mapref-2.2b.fna data/mapref-2.2b.fna.ncbitax > $OUT 2> /dev/null

if [ ! -f "$OUT" ]; then echo "FAIL: output file missing $OUT"; exit 99; fi # output file not found

gawk '
function ABS(a){ return(a>=0.0?a:-a);}
BEGIN {
  FS="\t";
  while (getline<ARGV[1] > 0){
    if (/^#/) continue;
    arr[$1]=$0;
    arrscore[$1]=$3;
    arrid[$1]=$4;
  }
  while (getline<ARGV[2] > 0){
    if (/^#/) continue;
    arr2[$1]=$0;
    arrscore2[$1]=$3;
    arrid2[$1]=$4;
  }
#  if (length(arr)!=length(arr2)) { print "FAIL: number of entries different " length(arr) " != " length(arr2); exit(99); }
  for (i in arr){
#    if (!(i in arr2)) { print "FAIL: entry not present in both files " i; exit(99); }
    if (!(i in arr2)) { print "FAIL: entry not present in both files " i; }
#    if (ABS(arrscore[i]-arrscore2[i])/arrscore[i]>=0.01){ print "FAIL: difference between scores larger than 1% " i " :: " arrscore[i] " " arrscore2[i]; exit(99); }
    if (ABS(arrscore[i]-arrscore2[i])/arrscore[i]>=0.01){ print "FAIL: difference between scores larger than 1% " i " :: " arrscore[i] " " arrscore2[i]; }
    if (ABS(arrid[i]-arrid2[i])>=0.01){ print "FAIL: difference between identities larger than 1% " i " :: " arrid[i] " " arrid2[i]; print arr[i]; print arr2[i]; }
#    if (ABS(arrid[i]-arrid2[i])>=0.01){ print "FAIL: difference between identities larger than 1% " i " :: " arrid[i] " " arrid2[i]; print arr[i]; print arr2[i]; exit(99); }
  }
}' $OUT $OUT.ref

#if [ "$OTUCOUNT" != "$OTUCOUNT2" -o "$LASTDIST" != "$LASTDIST2" ]; then exit 99; fi

exit 0
