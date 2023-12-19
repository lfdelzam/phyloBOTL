#!/usr/bin/env python3
import os
import argparse
import re
from collections import OrderedDict


usage = 'python parse_eggNOGs.py [options]'
description = 'This program parses the Orthologues EggNOG/KEGG annotations and prints out the unique eggNOG_OGs, KEGG_Pathway and KEGG_Module annotation per othologue'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i',dest='i',help='directory containing the EggNOG annotations', required=True)
parser.add_argument('-o', dest='o', help='output directory', required=True)
parser.add_argument('-p', dest='p', help='List of Kegg pathways, .tsv')
parser.add_argument('-m', dest='m', help='List of Kegg Modules, .tsv')


args = parser.parse_args()


files = [f for f in os.listdir(args.i) if f.endswith(".annotations")]

def KEGG_Pathway_list(file):
    KEGGP={}
    with open (file, "r") as fin:
        for line in fin:
            line=line.rstrip()
            if line.startswith("map"):
                line=line.split("\t")
                KEGGP[line[0]]=line[1]
    return KEGGP


def KEGG_Module_list(file):
    KEGGM={}
    with open (file, "r") as fin:
        for line in fin:
            line=line.rstrip()
            if line.startswith("M0"):
                line=line.split("\t")
                KEGGM[line[0]]=line[1]
    return KEGGM

def update_annot2(dict, co, param):
    if param != "":
        tempo=set(param.split(","))
        if co in dict:
            dict[co].update(tempo) 
        else:
            dict[co]=tempo
    return dict

def update_paths(dict, co, param):
    if param != "":
        tempo=set([ m for m in param.split(",") if m.startswith("map")])
        if co in dict:
            dict[co].update(tempo) 
        else:
            dict[co]=tempo
    return dict

def get_annotations2(file):
#query  seed_ortholog   evalue  score   eggNOG_OGs      max_annot_lvl   COG_category    Description     Preferred_name  GOs     EC      
#KEGG_ko KEGG_Pathway    KEGG_Module     KEGG_Reaction   KEGG_rclass     BRITE       KEGG_TC CAZy    BiGG_Reaction   PFAMs

###query_name seed_eggNOG_ortholog seed_ortholog_evalue seed_ortholog_score best_tax_level Preferred_name GOs EC KEGG_ko KEGG_Pathway KEGG_Module 
###KEGG_Reaction KEGG_rclass BRITE KEGG_TC CAZy BiGG_Reaction taxonomic scope eggNOG OGs best eggNOG OG COG Functional cat. eggNOG free text desc.


  cog_d=list()
  Pathway={}
  Module={}
  description={}
  ec={}
  gene={}
  reaction={}
  rc={}
  go={}
  bt={}
  with open(file, "r") as fin:
    for line in fin:
      line=line.rstrip()
      
      if not line.startswith("#"):
        line=line.split("\t")
        cs=line[18].split(",")  #4
        cogs=set()
        noncogs=set()
        for ca in cs:
            c=ca.split("@")[0]
            if c != "":
                if not c.startswith("COG"):
                    noncogs.update([c])
                elif c.startswith("COG"):
                    cogs.update([c])
        if len(cogs) == 0:
            cogs=noncogs

        cogs=sorted(cogs)
        tcogs=",".join(cogs)
        cog_d.extend([tcogs])

        description=update_annot2(description,tcogs, line[len(line)-1]) #7
        gene=update_annot2(gene,tcogs, line[5]) #8
        go=update_annot2(go,tcogs, line[6]) #9
        ec=update_annot2(ec,tcogs, line[7]) #10

        Pathway=update_paths(Pathway,tcogs, line[9]) #12
        Module=update_annot2(Module,tcogs, line[10]) #13
        reaction=update_annot2(reaction,tcogs,line[11]) #14
        rc=update_annot2(rc,tcogs,line[12]) #15
        bt=update_annot2(bt,tcogs,line[13]) #16

  return Pathway,Module,go,ec,gene,reaction,rc,description, cog_d, bt

def tex_to(DK):
        texto=",".join(DK)
        return texto

def update_dict(DK1, D2, nam):
    for g in DK1:
        if g in D2:
            D2[g].append(nam)
        else:
            D2[g]=[nam]
    texto= tex_to(DK1)
    D2= OrderedDict(sorted(D2.items(), key=lambda x: len(x[1]), reverse=True))


    return D2, texto

o_modules=dict()
o_paths=dict()
o_react=dict()
o_rc=dict()
o_bt=dict()
o_go=dict()

KEGG_P=KEGG_Pathway_list(args.p)
KEGG_M=KEGG_Module_list(args.m)

if not os.path.exists(args.o):
    os.makedirs(args.o)

with open(os.path.join(args.o,"Kegg_annotations.tsv"), "w") as fout:
    fout.write("#Ortholog\teggNOG_OGs\tCOG most Frequent (%)\tKEGG_Pathway\tKEGG_Module\tGOs\tEC\tGene\tKEGG_Reaction\tKEGG_rclass\tBRITE\tDescription\n") #, file=fout)
    for f in files:

        KEGG_Pathway,KEGG_Module, GO,EC,Gene,Reaction,RC,Description, COG_L, BT=get_annotations2(os.path.join(args.i,f))

        name=f.split(".")[0]
        if len(COG_L)>0:
            max_COG=max(set(COG_L), key=COG_L.count)
            perct_COGmax=round(COG_L.count(max_COG)*100/len(COG_L),2)
            tp=""
            tm=""
            tr=""
            trc=""
            tgo=""
            tec=""
            tg=""
            td=""
            tbt=""
            if max_COG in KEGG_Pathway:
                o_paths, tp=update_dict(KEGG_Pathway[max_COG], o_paths, name)
            if max_COG in KEGG_Module:
                o_modules, tm=update_dict(KEGG_Module[max_COG], o_modules, name)
            if max_COG in Reaction:
                o_react, tr=update_dict(Reaction[max_COG], o_react, name)
            if max_COG in RC:
                o_rc, trc=update_dict(RC[max_COG], o_rc, name)
            if max_COG in BT:
                o_bt, tbt=update_dict(BT[max_COG], o_bt, name)
            if max_COG in GO:
                o_go, tgo=update_dict(GO[max_COG], o_go, name)
            if max_COG in EC:
                tec=tex_to(EC[max_COG])
            if max_COG in Gene:
                tg=tex_to(Gene[max_COG])
            if max_COG in Description:
                td=tex_to(Description[max_COG])

            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(name, " | ".join(set(COG_L)),str(max_COG)+" ("+str(perct_COGmax)+")",tp, tm,tgo,tec,tg,tr,trc,tbt,td)) #, file=fout)

        else:
            print("INFO: Orthologue "+name+" has not EGGNOG hit")



with open(os.path.join(args.o,"orthologues_sharing_Kegg_modules.tsv"), "w") as fout2:
    fout2.write("#Module\tOrthologues\tDescription\n") #, file=fout2)
    for m,o in o_modules.items():
        if len(o) >0:
            if m in KEGG_M:
                fout2.write("{}\t{}\t{}\n".format(m,",".join(o), KEGG_M[m])) #, file=fout2)
            else:
                fout2.write("{}\t{}\t{}\n".format(m,",".join(o), "")) #, file=fout2)


with open(os.path.join(args.o,"orthologues_sharing_Kegg_pathway.tsv"), "w") as fout3:
    fout3.write("#Pathway\tOrthologues\tDescription\n") #, file=fout3)
    for map,o in o_paths.items():
        if len(o) >0:
            if map in KEGG_P:
                fout3.write("{}\t{}\t{}\n".format(map,",".join(set(o)),KEGG_P[map]) ) #,file=fout3)
            else:
                fout3.write("{}\t{}\t{}\n".format(map,",".join(set(o)),"") ) #,file=fout3)

with open(os.path.join(args.o,"orthologues_sharing_Kegg_reaction.tsv"), "w") as fout4:
    fout4.write("#Reaction\tOrthologues\n") #, file=fout4)
    for rc,o in o_react.items():
        if len(o) >1:
            fout4.write("{}\t{}\n".format(rc,",".join(set(o)))) #, file=fout4)

with open(os.path.join(args.o,"orthologues_sharing_Kegg_rclass.tsv"), "w") as fout5:
    fout5.write("#rClass\tOrthologues\n") #, file=fout5)
    for rn,o in o_rc.items():
        if len(o) >1:
            fout5.write("{}\t{}\n".format(rn,",".join(set(o))) ) #, file=fout5)

with open(os.path.join(args.o,"orthologues_sharing_GOs.tsv"), "w") as fout6:
    fout6.write("#GOs\tOrthologues\n") #, file=fout6)
    for g,o in o_go.items():
        if len(o) >1:
            fout6.write("{}\t{}\n".format(g,",".join(set(o)))) #, file=fout6)
