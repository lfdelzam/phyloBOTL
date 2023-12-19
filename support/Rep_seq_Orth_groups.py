#!/usr/bin/env python3
import os
import argparse
import re
import random


usage = 'python Rep_seq_Orth_groups.py [options]'
description = 'This program select one sequence from a Ortholog group'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i',dest='i',help='directory containing the Orthologue group sequences', required=True)
parser.add_argument('-c',dest='c',help='selection criterium. Options longest, random', default="longest")
parser.add_argument('-o', dest='o', help='output directory', required=True)
args = parser.parse_args()


files = sorted([f for f in os.listdir(args.i) if f.endswith(".fa")])

def select_Ortho_seqs(fi, criterium):
    seq=""
    sel=0
    counter=0
    num = len([1 for line in open(fi) if line.startswith(">")])
    if criterium == "random":
        sel=random.randint(1,num)


    with open(fi, "r") as fin:
            for line in fin:
              line=line.rstrip()

              if line.startswith(">"):
                  if seq != "":
                        counter +=1
                        if criterium == "longest" and len(seq) > sel:
                            selected=(id, seq, num)
                            sel=len(seq)
                        if criterium == "random" and counter == sel:
                            selected=(id, seq, num)
                            return selected
                            #break

                  id=line
                  seq=""
              else:
                   seq+=line

            if seq != "":
                counter +=1
                if criterium == "longest" and len(seq) > sel:
                    selected=(id, seq, num)
                if criterium == "random" and counter == sel:
                    selected=(id, seq, num)

    return selected


pathout=os.path.join(args.o,"selection")
if not os.path.exists(pathout):
    os.makedirs(pathout, exist_ok=True)

with open(os.path.join(args.o, "report.txt"), "w") as foutall:
    print("#criterium "+str(args.c), file=foutall)
    print("#Orthologue_group\tN.seqs_in_the_group\tID_selected_sequence\tlength", file=foutall)
    for fi in files:

        pathfile=os.path.join(args.i,fi)
        sequence=select_Ortho_seqs(pathfile, args.c)
        sel_orth=sequence[0]
        len_seq=len(sequence[1])
        with open(os.path.join(args.o, "selection",fi), "w") as fout:
            print("{}\n{}".format(sel_orth, sequence[1]), file=fout)

        name=fi.split(".fa")[0]
        print("{}\t{}\t{}\t{}".format(name,sequence[2],sel_orth[1:],len_seq), file=foutall)

