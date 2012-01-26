#!/bin/bash

#DIR=/work/elsanto-clued0/Z/Bs/fit/BDT20
DIR=/work/elsanto-clued0/Z/Bs/fit/CUT_fs0
LOCK=CUT_fs0
#ROOT=ws_BDT20_fit.root
ROOT=ws_CUT_fit.root
CHILDREN=100

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

cp /lib/libssl.so.6 libssl.so.4 && cp /lib/libcrypto.so.6 libcrypto.so.4 && export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.

scp elsanto-clued0:${DIR}/server_markov.C . || die "EE: Copy script failed 1."
scp elsanto-clued0:${DIR}/server_markov_C.so . || die "EE: Copy script failed 2."
mkdir lib || die "EE: Creating lib dir."
scp elsanto-clued0:${DIR}/lib/libBFitter.so lib/ || die "EE: Copy lib failed."
scp elsanto-clued0:${DIR}/${ROOT} workspace.root || die "EE: Copy script failed."

touch ~/${LOCK} || die "EE: Can't create lock."
while [ -f ~/${LOCK} ]
do
  sleep 1;
done

time root -b -q -l server_markov.C+\(${CHILDREN}\) || die "EE: markov chain failed."
kbatch || die "EE: kbatch2."
scp mcmc.root elsanto-clued0:${DIR}/  || die "EE: Copy mcmc.root back."
