#!/bin/bash

dir=$1

echo "name,time_limit,sol_value,time_elapsed,timeout"
for _f in $(ls $dir/*.in | sort -n); do
    f=$(echo $_f | cut -f1 -d.)_bb.stats
    name=$(basename $_f)
    timeout=$(cat $f | grep timeout | cut -f2 -d=)
    sol_value=$(cat $f | grep sol_value | cut -f2 -d=)
    time_elapsed=$(cat $f | grep "time=" | cut -f2 -d=)
    echo "$name,15,$sol_value,$time_elapsed,$timeout"
done
