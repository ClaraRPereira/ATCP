#!/bin/sh

ncols=$(sed -e "1q" "$1" | grep -o " " | wc -l)
echo $ncols $((ncols-1))
if [ "$2" = "nrg" ]
then
	start=$ncols
	end=$ncols
else
	start=2
	end=$((ncols-1))
fi
gnuplot -e "filename=\"$1\";start=$start;end=$end" plot.gp
