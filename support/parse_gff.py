#!/usr/bin/env python3
import os
import argparse
import re

usage = 'python parse_gff.py [options]'
description = 'This program adds Orthologues number to genomes gff files'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i',dest='i',help='directory containing the gff files', required=True)
parser.add_argument('-l',dest='l',help='Orthofinder file, Orthologue groups.tsv', required=True)
parser.add_argument('-o', dest='o', help='output directory', required=True)
args = parser.parse_args()


files = [f for f in os.listdir(args.i) if f.endswith(".gff")]

def Ortho_groups(fi):

    GeC={}
    first_line=True
    with open(fi, "r") as fin:
        for line in fin:
          line=line.rstrip()
          line=line.split("\t")

          if first_line:
              genomes=[w.strip() for w in line[1:] ]
              first_line =False

          else:

            Orth=line[0]

            lp=[ y.split(",") for y in line[1:] ]  #list of proteins

            for i in range(0,len(line[1:])) :
                genome=genomes[i]

                if len(lp[i]) == 1 and lp[i] != "":
                    p= lp[i][0]

                    if genome in GeC:
                        GeC[genome][p]= Orth

                    else:
                        GeC[genome]={p: Orth}

                if len(lp[i])>1: #if there is several genes in the same orthologue group

                    proteins=[q.strip() for q in lp[i] ]

                    for p in proteins:
                        if p != "":
                            if genome in GeC:
                                GeC[genome][p]= Orth

                            else:
                                GeC[genome]={p: Orth }


    return GeC

OG=Ortho_groups(args.l)


if not os.path.exists(args.o):
    os.makedirs(args.o)

#print(OG.keys())
with open(os.path.join(args.o, "orthologues_gff.tsv" ), "a") as fout:
    print("Genome\tOrthologue\tContig\tType\tStart\tStop\tDirection\tDescription", file=fout)
    for f in files:
        Geno=f.split(".gff")[0]
        if Geno in OG.keys():
            with open(os.path.join(args.i,f), "r") as fin:
                for line in fin:
                    line=line.rstrip()
                    if not line.startswith("#"):
                        line=line.split("\t")
                         #PDNR01000754.1	Prodigal:002006	CDS	993	1613	.	+	0	ID=MPFGEDNF_00001;inference=ab initio prediction:Prodigal:002006;locus_tag=MPFGEDNF_00001;product=hypothetical protein
                        contig=line[0]
                        type=line[2]
                        start=line[3]
                        stop=line[4]
                        direction=line[6]
                        ID=line[8]
                        prot=ID.split(";")[0].split("=")[1]

                        if prot in OG[Geno].keys():
                            O=OG[Geno][prot]
                        elif type != "CDS":
                            O=type
                        else:
                            O="NA"

                        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(Geno,O,contig,type,start,stop,direction,ID), file=fout)
        else:
            print(Geno)  # it seems can be removed 
