cat angleCUT.out | awk 'function abs(value){return (value<0?-value:value);}; abs($3)<$4 {print "ws->var(\""$2"\")->setConstant(kTRUE); //"$1}'
