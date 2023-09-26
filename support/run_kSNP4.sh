#!/bin/bash -l

gdir=$1 #"Genomes_all"
filein=$2 #"V_vulnificus_genomes_all"
outdir=$3 #"kSNP4_all_output"
MakeKSNP4infile -indir $gdir -outfile $filein
Kchooser4 -in $filein
value=$(grep "The optimum value of k is" Kchooser4_$filein.report | sed s/"The optimum value of k is "// )
kSNP4 -in $filein -k $value -outdir $outdir -CPU $4 $5 $6
