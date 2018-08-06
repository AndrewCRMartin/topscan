#!/bin/bash

topdir=/home/amartin/topscan
export dbname=cath

echo -n "$1 $2 "
sort -r -n +1 $1 | $topdir/src/analyse.perl | $topdir/src/findline.perl $2
