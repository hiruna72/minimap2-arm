#!/bin/bash


../divide reads.fasta 16

for file in exp/*.fa
do
	diff  $file $(basename $file) -q || echo "Diff failed!"
done
rm *.csv
rm part*.fa