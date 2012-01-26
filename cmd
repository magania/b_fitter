cat angleCUT.out | awk 'function abs(value){return (value<0?-value:value);}; abs($3)<$4 {print "ws->var(\""$2"\")->setConstant(kTRUE); //"$1}'

----------------------------------------------------------------------------
rm /RunII/home/magania/BDT20
sleep 30
echo d0cs3060.fnal.gov > /RunII/home/magania/server.txt

for t in `seq 1 10`; do qsub markov_cab.sh -q sam_hi@d0cabsrv1; done
qsub markov_scab.sh -q sam_hi@d0cabsrv1 -k oe

rm /RunII/home/magania/server.txt
----------------------------------------------------------------------------
