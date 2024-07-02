#!/bin/bash

for datedir in /xnilu_wrk/users/dgho/gfas/archive/mirrorMARS/0001/2015*/
do
    echo $datedir
    for i in {0030..2330..100}
    do
        if ! (test -d $datedir$i); then
            mkdir $datedir$i
            echo $datedir$i
        fi
        for filename in $datedir/1200/*97.grb.gz
        do
            ln -sr $filename /$datedir$i/
        done
        echo $i
    done
done
