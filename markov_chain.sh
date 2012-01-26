#!/bin/bash

ID=0
DIR=/work/elsanto-clued0/Z/Bs/fit/CUT
EVENTS=10000
BETA=-0.2
DG=-0.1
ROOT=ws_CUT_fit.root
NUMCPU=8


#CUTHERE


die()
{
echo "$@"
exit 1
}

if [ ! -n "$PBS_JOBID" ]; then
  echo "EE: Workdir not found."
  exit 1
fi


cd /scratch/$PBS_JOBID || die "EE: CD Workdir failed."

kbatch || die "EE: kbatch."
scp elsanto-clued0:/work/elsanto-clued0/Z/root_tgz/root_v5.26.00.Linux-slc4-gcc3.4.tar.gz .  || die "EE: Copy root failed."
tar xzf root_v5.26.00.Linux-slc4-gcc3.4.tar.gz || die "EE: Untar failed."
. root/bin/thisroot.sh || die "EE: Setup root failed."

cp /lib/libssl.so.6 libssl.so.4 && cp /lib/libcrypto.so.6 libcrypto.so.4 && export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.  || die "EE: Library Fix failed."

scp elsanto-clued0:/work/elsanto-clued0/Z/Bs/fit/b_fitter/markov_chain.C . || die "EE: Copy script failed."
scp elsanto-clued0:/work/elsanto-clued0/Z/Bs/fit/b_fitter/markov_chain_C.so . || die "EE: Copy script failed."
mkdir lib || die "EE: Creating lib dir."
scp elsanto-clued0:${DIR}/lib/libBFitter.so lib/ || die "EE: Copy lib failed."
scp elsanto-clued0:${DIR}/${ROOT} ws_data_fit.root || die "EE: Copy root file failed."

time root -b -q -l markov_chain.C+\(\"ws_data_fit.root\",${EVENTS},$NUMCPU,$BETA,$DG\) || die "EE: markov chain failed."

kbatch || die "EE: kbatch2"
scp mcmc.root elsanto-clued0:${DIR}/mcmc_${ID}_${EVENTS}_${BETA}_${DG}.root || die "Copy mcmc.root back."
