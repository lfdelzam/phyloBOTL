#!/usr/bin/env python3
import os
import argparse
import re
from collections import OrderedDict


usage = 'python parse_Uniprot.py [options]'
description = 'This program parses the Orthologues Uniprot annotation table and prints out the annotation per cluster of othologue group'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i',dest='i',help='Uniprot annotations ', required=True)
parser.add_argument('-c',dest='c',help='Cluster file ', required=True)
parser.add_argument('-o', dest='o', help='output directory', required=True)


args = parser.parse_args()


def get_cluster_groups(filein):
    cl=dict()
    with open(filein, "r") as fin:
        for line in fin:
            line=line.rstrip()
            if line.startswith(">"):
                cluster=line[1:]
            else:
                if cluster in cl:
                    cl[cluster].append(line)
                else:
                    cl[cluster]=[line]
    return(cl)

def get_uniprot_annotations(filein):
    #Ortholog	Reference_in_db	Description	Frequency (%)	Max_Score
    Ort_annot=dict()
    with open(filein, "r") as fin:
        for line in fin:
            line=line.rstrip()
            if not line.startswith("#"):
                line=line.split("\t")
                if line[1] != "No hits":
                    Ort_annot[line[0]]=line[2]
                else:
                    Ort_annot[line[0]]=line[1]
    return(Ort_annot)




if not os.path.exists(args.o):
    os.makedirs(args.o)

clusters=get_cluster_groups(args.c)
annotations=get_uniprot_annotations(args.i)

with open(os.path.join(args.o,"cluster_most_frequent_annotations.tsv"), "w") as fout:
    fout.write("#Cluster\tOrthologue\tAnnotation\n")
    for c in clusters:
        for o in clusters[c]:
            if o in annotations:
                anota=annotations[o]
            else:
                anota="No hit"
            fout.write("{}\t{}\t{}\n".format(c,o,anota))

