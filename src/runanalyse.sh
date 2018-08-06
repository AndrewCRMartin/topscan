#!/bin/bash
#
# Assumes the PDB file is in the current directory named pdbXXXX.ent
#
# Run the program as:
#    runanalyse.sh XXXX C.A.T
# where XXXX is the PDB code and C.A.T is the CAT assignment

topdir=/home/amartin/topscan
export dbname=cath

topscan -s -p -e 4 -h 4 pdb$1.ent $topdir/cathsn_44.top > $1.e4h4.out
sort -n +1 $1.e4h4.out | tail -50 | $topdir/analyse.perl > $1.e4h4.anal50
sort -n +1 $1.e4h4.out | tail -100 | $topdir/analyse.perl > $1.e4h4.anal100
echo "$1 $2 e4h4 50"
grep $2 $1.e4h4.anal50 | wc
echo "$1 $2 e4h4 100"
grep $2 $1.e4h4.anal100 | wc


topscan -s -p -e 3 -h 3 pdb$1.ent $topdir/cathsn_33.top > $1.e3h3.out
sort -n +1 $1.e3h3.out | tail -50 | $topdir/analyse.perl > $1.e3h3.anal50
sort -n +1 $1.e3h3.out | tail -100 | $topdir/analyse.perl > $1.e3h3.anal100
echo "$1 $2 e3h3 50"
grep $2 $1.e3h3.anal50 | wc
echo "$1 $2 e3h3 100"
grep $2 $1.e3h3.anal100 | wc


topscan -s -p -e 3 -h 4 pdb$1.ent $topdir/cathsn_e3h4.top > $1.e3h4.out
sort -n +1 $1.e3h4.out | tail -50 | $topdir/analyse.perl > $1.e3h4.anal50
sort -n +1 $1.e3h4.out | tail -100 | $topdir/analyse.perl > $1.e3h4.anal100
echo "$1 $2 e3h4 50"
grep $2 $1.e3h4.anal50 | wc
echo "$1 $2 e3h4 100"
grep $2 $1.e3h4.anal100 | wc


topscan -s -p -e 4 -h 3 pdb$1.ent $topdir/cathsn_e4h3.top > $1.e4h3.out
sort -n +1 $1.e4h3.out | tail -50 | $topdir/analyse.perl > $1.e4h3.anal50
sort -n +1 $1.e4h3.out | tail -100 | $topdir/analyse.perl > $1.e4h3.anal100
echo "$1 $2 e4h3 50"
grep $2 $1.e4h3.anal50 | wc
echo "$1 $2 e4h3 100"
grep $2 $1.e4h3.anal100 | wc


