#!/bin/sh

make

echo '===================================='
./kc2bc kc_example.txt 2 -w 0
../problems/k-cluster/biqcrunch kc_example_k2.bc ../problems/k-cluster/biq_crunch.param 
echo '===================================='
./kc2bc kcw_example.dat -w 1
../problems/k-cluster/biqcrunch kcw_example_k2.bc ../problems/k-cluster/biq_crunchW.param 
echo '===================================='
./lp2bc.py lp_example.txt > lp_example.bc
../problems/generic/biqcrunch lp_example.bc ../problems/generic/biq_crunch.param 
echo '===================================='
./mc2bc.py mc_example.txt > mc_example.bc
../problems/max-cut/biqcrunch mc_example.bc ../problems/max-cut/biq_crunch.param 
echo '===================================='
./qp2bc qp_example.txt > qp_example.bc
../problems/max-cut/biqcrunch qp_example.bc ../problems/max-cut/biq_crunch.param 
echo '===================================='

