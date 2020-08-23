#!/bin/bash
beds_dir=$1
ref_beds_dir=$2
intersections=$3
output=$4

while IFS=$'\t' read -r -a myArray
do
    ref_id=${myArray[0]};
    ref_bed="$ref_beds_dir/$ref_id.bed";
    tool_id=${myArray[1]};
    tool_bed="$beds_dir/$tool_id.bed";
    seqpare $ref_bed $tool_bed | awk 'BEGIN{OFS="\t"} NR==1 {s=$3/($1+$2-$3); print "'$ref_id'","'$tool_id'",s}'
done < $intersections > $output;
#
# awk 'NF==4 {print $4}' $workdir/pairings.tsv | Rscript -e 'summary (as.numeric (readLines ("stdin")))' >> "$workdir/summary.txt"
