#!/bin/sh
END=99
for ((i=1;i<=END;i++)); do
    export A=`head -$i runs_trial1.dat| tail -1`
    echo $A
    echo "/data/gem/xrayscan/COMPASS_4GEM_trial1/${A}"
    ./rawdata ${A} ${i}
done
