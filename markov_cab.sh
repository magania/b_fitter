#!/bin/bash

#DIR=/work/elsanto-clued0/Z/Bs/fit/BDT20
DIR=/work/elsanto-clued0/Z/Bs/fit/CUT_fs0
SERVER_FILE=server.txt

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

scp elsanto-clued0:${DIR}/client_markov.C . || die "EE: Copy script failed."
scp elsanto-clued0:${DIR}/client_markov_C.so . || die "EE: Copy script failed."
mkdir lib || die "EE: Creating lib dir."
scp elsanto-clued0:${DIR}/lib/libBFitter.so lib/ || die "EE: Copy lib failed."

while [ ! -f ~/${SERVER_FILE} ]
do
  sleep 30;
done

SERVER_NAME=`cat ~/${SERVER_FILE}`

time root -b -q -l client_markov.C+\(\"${SERVER_NAME}\"\) || die "EE: markov chain failed."
