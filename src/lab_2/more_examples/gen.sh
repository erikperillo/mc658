#!/bin/bash

for n in 10 50 100 200 300 350 400 450 500 550 600 700 800 900 1000; do
    echo "in n = $n"
    for i in $(seq 3); do
        d=$(($RANDOM%(n/2) + 1))
        B=$((2*($RANDOM%(n/3)) + n/3))
        sh ./gen_input.sh $n $d $B > $n"_"$d"_"$B".in"
    done
done
