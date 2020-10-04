#!/bin/sh

for p in $(ls ../problems/*/example.bc | sed 's/\/example\.bc//g')
do
        echo '===================================='
        cmd="$p/biqcrunch $p/example.bc $p/biq_crunch.param"
        echo $cmd
        $cmd | grep 'value' > $p/example.out
        if [ -z "`diff -q $p/example.out $p/example.ans`" ]
        then echo `cat $p/example.out` ... OK
        else echo `cat $p/example.out` ... FAIL
        fi
        rm -f $p/example.bc.output* $p/example.out
done
echo '===================================='

