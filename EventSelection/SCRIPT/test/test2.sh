#!/bin/bash


root -l -b -q 'test.C' >&$1.log&

#ii=0
#for aaa in 1 2 3; do
#	for bbb in 100 200 300; do
#		ii=$(echo "$ii + 1" | bc)
#		xxx=$(echo "$aaa + $bbb" | bc)
#		echo $xxx
#		echo $ii
#	done
#done

