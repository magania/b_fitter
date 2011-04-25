#!/bin/bash

DIR=`pwd`

for beta in `seq -1.4 0.20 1.4`; do
	for dg in `seq -0.30 0.06 0.31`; do	 
              echo "#!/bin/bash" > /tmp/job.sh
              echo "DG=$dg" >> /tmp/job.sh
              echo "BETA=$beta" >> /tmp/job.sh 
              echo "SAMPLE=CUT" >> /tmp/job.sh
              tail -n 25 fit.sh >> /tmp/job.sh
              echo send $beta $dg
              /usr/bin/qsub -q medium@d0cabsrv1 -l nodes=1 -k oe -m ae /tmp/job.sh 
	done
done
