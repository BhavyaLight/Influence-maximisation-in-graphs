#!/bin/bash

# input
# $1 - Graph file part1 
# $2 - Graph file part2
# $3 - targetfile - NO NEW LINES IN THIS
# $4 - number of seed nodes required
# $5 - max number of nodes in graph
# $6 - number of colors required
# $7 - value of R

# Outputs
# part1_output.txt
# part2_output.txt
# samplecplor.txt
# part2_query.txt
# statistics.txt -> contains summary of part 1 and part 2 outputs

echo "This is your input: $1 $2 $3 $4 $5 $6 $7"

echo "Running code to find seed nodes..."
./inf_ic_bfs_celf $1 part1_output.txt $3 samplecolor.txt $4 $5

echo "Converting output to correct format for part 2 of the code..."
# Result modification for part 2
awk '{ if($1 != "Time") {printf "%s ",$2}}' part1_output.txt > part2_query.txt
echo -n ";" >> part2_query.txt
cat $3 >> part2_query.txt

echo "Storing current results..."
awk '{ if($1 != "Time") {printf "%s ",$2} else {print $1,$2 }}' part1_output.txt >> statistics.txt

echo "Running code to find top-k colors..."
./colrel.out h rel $2 part2_query.txt part2_output.txt $6 $7

echo "Converting output to correct format for re-iteration of part1"
awk 'BEGIN { FS=";"}; {print $3}' part2_output.txt >samplecolor.txt
echo "Storing current results..."
awk 'BEGIN { FS=";"}; {print $3,$6}' part2_output.txt >>statistics.txt

echo "" >> statistics.txt