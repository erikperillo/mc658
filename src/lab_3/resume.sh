#!/bin/bash

dir="./examples2"
t=15

echo "arquivo,tMax(s),custo sol,tempo(s),timeout"
#for _f in $(ls $dir/*.in | sort -n); do
for _f in ./2001_1940.in; do
    f=out
    ./transportation.e -i $_f -t $t > $f
    name=$(basename $_f)
    [[ -z $(grep "Time limit reached" $f) ]] && timeout=0 || timeout=1
    sol_value=$(cat $f | grep "cost =" | cut -f3 -d" ")
    time_elapsed=$(cat $f | grep " iterations and " | head -n1 |\
        rev | cut -f2 -d" " | rev)
    echo "$name,15,$sol_value,$time_elapsed,$timeout"
done
