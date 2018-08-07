#!/bin/bash

if [ "$#" -ne 5 ]; then
    echo "usage : $0 <reference.fa> <numparts> <output_index> <minimap_executable> <minimap_profile>"
        exit
fi
REFERENCE=$1
NUMPARTS=$2
OUTIDX=$3
MINIMAP=$4
PROFILE=$5

BASEDIR=$(dirname "$0")

test -e stat.csv && rm *.csv

set -e

echo "Compiling"
gcc -Wall -O3 $BASEDIR/divide.c -o divide 

echo "Running divider"
/usr/bin/time -v ./divide $REFERENCE $NUMPARTS  2> >(tee divider.log) 

echo "Minimap2 Indexing"

LIST="";

NUMPARTS2=$(echo "$NUMPARTS-1" | bc )
for i in $(seq 0 $NUMPARTS2)
do
	echo "Index for part $i"
	/usr/bin/time -v $MINIMAP -a -x $PROFILE -d part$i.idx part$i.fa 2> >(tee part$i.log)
	LIST="$LIST"" part$i.idx"
	rm part$i.fa
done


max_mem=$(grep "Maximum resident set size" *.log | awk '{print $NF/1024}' | sort -n | tail -1);
echo "Max mem usage : $max_mem MB"

	
echo "concatenating $LIST"	
cat $LIST > $OUTIDX
rm $LIST

for i in $(seq 0 $NUMPARTS2)
do
	rm part$i.log
done