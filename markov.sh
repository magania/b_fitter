#! /bin/bash

. /work/elsanto-clued0/Z/root_SL4_i386_gcc43/bin/thisroot.sh 
cd /work/elsanto-clued0/Z/Bs/fit/BDT20_Dmfix
root -b -q -l client_markov.C+
