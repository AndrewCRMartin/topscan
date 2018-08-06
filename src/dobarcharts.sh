#!/bin/bash

for file in *.out
do
    awk '{print $2}' $file >bars/tmp.bar
    bars -n 20 bars/tmp.bar >bars/$file.bar
done
rm bars/tmp.bar

for file in bars/*.bar
do
    echo "barchart" >bars/tmp.bar
    echo "shrfirst" >>bars/tmp.bar
    echo "boxed" >>bars/tmp.bar
    echo "centxlab" >>bars/tmp.bar
    echo "ticks 5 100" >>bars/tmp.bar
    echo "bounds 0 100 0 1000" >>bars/tmp.bar
    cat $file >>bars/tmp.bar
    mv -f bars/tmp.bar $file
    amplot $file >$file.ps
    rm -f $file
done

