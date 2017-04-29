#!/bin/sh

#read n
#read d
#read B
n=$1
d=$2
B=$3

classes=10
echo $n $d $B
for i in `seq 1 $n`
do
	p=$(($RANDOM%40))
	w=$(($RANDOM%80))
	c=$(($RANDOM%10))
	echo $(($p+1)) $(($w+1)) $(($c+1))
done
