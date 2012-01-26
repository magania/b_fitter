#!/bin/bash

DG=0.1
BETA=0.1
SAMPLE=BDT20









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

#scp elsanto-clued0:/work/elsanto-clued0/Z/Bs/b_fitter/ws_${SAMPLE}_fit.root . || die "EE: Copy fit failed."
scp elsanto-clued0:/work/elsanto-clued0/Z/Bs/b_fit_${SAMPLE}/ws_${SAMPLE}_fit.root . || die "EE: Copy fit failed."
scp elsanto-clued0:/work/elsanto-clued0/Z/Bs/b_fit_${SAMPLE}/ScanBs . || die "EE: Copy ScanBs failed."

./ScanBs -d $DG -b $BETA -i ws_${SAMPLE}_fit.root -o ws_${SAMPLE}_fit_${DG}_${BETA}.root || die "EE: ScanBs failed."

scp ws_${SAMPLE}_fit_${DG}_${BETA}.root elsanto-clued0:/work/elsanto-clued0/Z/Bs/scan/${SAMPLE}/ || die "Copy root file back failed."
