#!/bin/bash

TOTALEVENTS=1000000
NJOBS=100
ROOT=ws_CUT_fit.root

sendJob()
{
echo "#!/bin/bash" > /tmp/job.sh
echo "ID=${ID}" >> /tmp/job.sh
echo "DIR=`pwd`" >> /tmp/job.sh 
echo "EVENTS=$[TOTALEVENTS/NJOBS]" >> /tmp/job.sh
echo "BETA=${BETA}" >> /tmp/job.sh 
echo "DG=${DG}" >> /tmp/job.sh
echo "ROOT=${ROOT}" >> /tmp/job.sh
echo "NUMCPU=2" >> /tmp/job.sh

awk '$0=="#CUTHERE" {doPrint=1}; {if (doPrint==1) print $0 }' markov_chain.sh >> /tmp/job.sh

#cat /tmp/job.sh
#echo ----------------------------------------------------------------------------------------
qsub -q sam_hi@d0cabsrv2 /tmp/job.sh
}


BETA=0.2
DG=0.1
for i in `seq 1 $[NJOBS/4]`; do
 ID=a${i}
 echo $ID, $BETA, $DG
 sendJob 
done

BETA=-0.2
DG=0.1
for i in `seq 1 $[NJOBS/4]`; do
 ID=o${i}
 echo $ID, $BETA, $DG
 sendJob 
done

BETA=0.2
DG=-0.1
for i in `seq 1 $[NJOBS/4]`; do
 ID=e${i}
 echo $ID, $BETA, $DG
 sendJob 
done

BETA=-0.2
DG=-0.1
for i in `seq 1 $[NJOBS/4]`; do
 ID=u${i}
 echo $ID, $BETA, $DG
 sendJob 
done
