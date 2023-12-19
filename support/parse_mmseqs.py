#!/usr/bin/env python3
import os
import argparse
import re
from collections import OrderedDict


usage = 'python parse_mmseqs.py [options]'
description = 'This program parses the mmseqs easy-search output tables and prints out the hit with the higest score and more frequent per Ortholog'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i',dest='i',help='directory containing the mmseqs tables', required=True)
parser.add_argument('-r',dest='r',help='path to tax_specific_protein_annotation', required=True)
parser.add_argument('-p',dest='p',help='prefix', required=True)
parser.add_argument('-o', dest='o', help='output directory', required=True)

args = parser.parse_args()


files = [f for f in os.listdir(args.i) if f.endswith("hits.txt")]


def get_ref_annotations(file):
    ref_annot={}
    with open(file, "r") as fin:
        for line in fin:
          line=line.rstrip()
          if line.startswith(">"): #it should follow this estructure
#UniRef90_A0A023NA98 MARTX n=9 Tax=Vibrio TaxID=662 RepID=A0A023NA98_VIBVL
            line=line.split()
            id=line[0][1:]
            ref_annot[id]=" ".join(line[1:])
    return(ref_annot)


def get_annotations(file, ref_annot):
#The file is formatted as a tab-separated list with 12 columns: (1,2) identifiers for query and target sequences/profiles, (3) sequence identity, (4) alignment length, (5) number of mismatches, (6) num- ber of gap openings, (7-8, 9-10) domain start and end-position in query and in target, (11) E-value, and (12) bit score
#get the hit with highest score and more frequent

  scoremax=0
  refer_freq={}
  refer_scores={}
  refer_max_score_per_prot={}

  with open(file, "r") as fin:
    counter=0
    for line in fin:
      line=line.rstrip()
      if not line.startswith("#"):
        line=line.split("\t")
        if line[1] in refer_freq:
            refer_freq[line[1]] += 1
        else:
            refer_freq[line[1]] = 1
        if int(line[11]) in refer_scores:
            refer_scores[line[11]].add(set([line[1]]))
        else:
            refer_scores[line[11]]=set([line[1]])
        if int(line[11]) > scoremax:
            scoremax=int(line[11])
        if line[1] in refer_max_score_per_prot:
            if int(line[11]) > refer_max_score_per_prot[line[1]] :
                refer_max_score_per_prot[line[1]] = int(line[11])
        else:
            refer_max_score_per_prot[line[1]] = int(line[11])
      counter += 1

  list_best_hits=list(refer_scores[str(scoremax)])
  sorted_refer_freq = sorted(refer_freq.items(), key=lambda x:x[1], reverse=True)
  sorted_refer_freq_dict = dict(sorted_refer_freq)
  best_hit_to_find=True
  most_freq_hit_to_find=True
  ortho=str(os.path.basename(file).split(".")[0])
  with open(os.path.join(args.o,ortho+"_annotations.txt"), "w") as fout:
      fout.write("#Ortholog\tReference_in_db\tDescription\tFrequency (%)\tMax_Score\n")
      for o in sorted_refer_freq_dict:
          if most_freq_hit_to_find:
              thehit_most_freq=[o,ref_annot[o],str(round(float(sorted_refer_freq_dict[o]/counter)*100,2)),str(refer_max_score_per_prot[o])]
              most_freq_hit_to_find=False

          if o in list_best_hits and best_hit_to_find:
              thehit=[o,ref_annot[o],str(round(float(sorted_refer_freq_dict[o]/counter)*100,2)),str(scoremax)]
              best_hit_to_find=False


          fout.write("{}\t{}\t{}\t{}\t{}\n".format(ortho,o,ref_annot[o],str(round(float(sorted_refer_freq_dict[o]/counter)*100,2)),str(refer_max_score_per_prot[o])))
      return(thehit, thehit_most_freq)





uniref_annot=get_ref_annotations(args.r)


if not os.path.exists(args.o):
    os.makedirs(args.o)

with open(os.path.join(args.o, args.p+"_orthologs_best_annotation.tsv"), "w") as fout, open(os.path.join(args.o,args.p+"_orthologs_most_frequent_annotation.tsv"), "w") as fout2:
    fout.write("#Ortholog\tReference_in_db\tDescription\tFrequency (%)\tMax_Score\n")
    fout2.write("#Ortholog\tReference_in_db\tDescription\tFrequency (%)\tMax_Score\n")
    for f in files:
        full_path_f=os.path.join(args.i,f)
        ort=f.split(".")[0]
        if os.path.exists(full_path_f) and os.path.getsize(full_path_f) > 0:
            thehit, thehit_most_freq=get_annotations(full_path_f, uniref_annot)
            fout.write("{}\t{}\n".format(ort,"\t".join(thehit)))
            fout2.write("{}\t{}\n".format(ort,"\t".join(thehit_most_freq)))
        else:
            fout.write("{}\tNo hits\n".format(ort))
            fout2.write("{}\tNo hits\n".format(ort))

