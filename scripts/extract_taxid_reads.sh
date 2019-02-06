#!/usr/bin/env bash

if [ $# -ne 6 ]; then 
    echo "$#"
    echo "USAGE: `basename $0` input.krak taxid reads_1.fq reads_2.fq out_1.fq out_2.fq"
    echo "Currently only gets reads directly assigned to taxid, not further classified"
    echo "To get reads for a genus, run for the genus and all species underneath"
    echo "Stores read identifiers in kraken_taxid_reads.txt in the current working directory"
    echo "Requires filterbyname.sh from the bbmap package to work"
    exit 85
fi

in_krak=$1
taxid=$2
in_1=$3
in_2=$4
out_1=$5
out_2=$6

awk -v taxid="$taxid" '$3==taxid { print }' "$in_krak" | cut -f 2 > kraken_"$taxid"_reads.txt
# only do this if there are some lines in the file, elese exit with error
if [[ $(wc -l < kraken_"$taxid"_reads.txt) -ge 0 ]]
then
    filterbyname.sh in="$in_1" in2="$in_2" names=kraken_"$taxid"_reads.txt include=true out="$out_1" out2="$out_2"
else
    echo "No reads found with $taxid"
    exit 1
fi 

