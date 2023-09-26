#!/usr/bin/env python3
import os
import argparse
import re

usage = 'python parse_gbk.py [options]'
description = 'This program filters the Enriched/Depleted Ortholgues from the gdk files'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i',dest='i',help='directory containing the gdk files', required=True)
parser.add_argument('-e',dest='e',help='orthologue gff table', required=True)
parser.add_argument('-s',dest='s',help='Enriched/Depleted orthologue list', required=True)
parser.add_argument('-o', dest='o', help='output directory', required=True)
args = parser.parse_args()


files = [f for f in os.listdir(args.i) if f.endswith(".gbk")]

def list_enriched(fe):
    enrd=[]
    with open(fe) as fin:
        for line in fin:
          line=line.rstrip()
          if not line.startswith("#"):
              enrd.append(line)
    return enrd


def selected_gff_lines(fi,fen):
    enriched=list_enriched(fen)
    #print(enriched)
    GeC={}
    C_id={}
    O_enriched={}
    with open(fi, "r") as fin:
        for line in fin:
          line=line.rstrip()
          if not "genome" in line:
            line=line.split()
        #    print(line)
            genome=line[0]
            Orth=line[1]
            if Orth in enriched:
                contig=line[2]
                ID='"'+str(line[7].split(";")[0][3:])+'"' #Protein
                if genome in GeC:
                    if contig in GeC[genome]:
                        GeC[genome][contig][ID]= Orth
                    else:
                        GeC[genome][contig]={ID: Orth}
                else:
                    GeC[genome]={contig:{ID: Orth}}

    return GeC

GC=selected_gff_lines(args.e, args.s)

if not os.path.exists(args.o):
    os.makedirs(args.o)

non_in_enrcihed=set()
for f in files:
    Geno=f.split(".gbk")[0]
    if Geno in GC.keys():
#        print("INFO: Filtering "+Geno)
        with open(os.path.join(args.i,f), "r") as fin, open(os.path.join(args.o,f), "w") as fout :
            copy=False
            save=False
            save_origine=False
            save_CDS=False
            for line in fin:
                line=line.rstrip()
                if line.startswith("LOCUS"):
                    if save:
                        for s in saved_lines:
                            print(s, file= fout)
                            save=False
                            save_origine=False
                    LOC=line.split()[1]
                    TE=[i for i in GC[Geno].keys() if re.match(i, LOC) ]
                    if len(TE) >0:
                        LOC=TE[0]
#                    if LOC in GC[Geno]:
                        saved_lines=[]
                        saved_lines.append(line)
                        copy=True
                        save_origine=True
                        copy_extra=True


                    else:
                        copy=False

                if not line.startswith("LOCUS") and copy:
                    if line.startswith("DEFINITION"):
                        saved_lines.append("Vibrio vulnificus "+Geno+" "+Geno+".")
                    elif line.startswith("SOURCE"):
                        saved_lines.append("Vibrio vulnificus")
                    elif line.startswith("KEYWORDS"):
                        nl=",".join( [ v for v in GC[Geno][LOC].values() ] )
                        l=re.sub("\.",nl, line)
                        saved_lines.append(l)
                    elif "ORGANISM" in line:
                        l=re.sub("Genus species", "Vibrio vulnificus", line)
                        saved_lines.append(l)
                    elif "/organism" in line:
                        l=re.sub("Genus species", "Vibrio vulnificus", line)
                        saved_lines.append(l)
                    elif "/strain" in line:
                        l=re.sub('"strain"', '"'+Geno+'"', line)
                        saved_lines.append(l)

                    elif "CDS" in line:
                        if save:
                            saved_lines.extend(saved_lines_CDS)
                            save_origine=True

                        save_CDS=True
                        copy_extra=False
                        saved_lines_CDS=[line]

                    elif "/locus_tag" in line:
                        prot=line.split("=")[1]

                        if prot in GC[Geno][LOC]:
                            save=True
                            #saved_lines_CDS.append(line)
                            newline=str(line.split("=")[0])+'="'+str(GC[Geno][LOC][prot])+'"'
                            saved_lines_CDS.append(newline)
                        else:
                            save_CDS=False
                            saved_lines_CDS=[]

                    elif save_CDS:
                        saved_lines_CDS.append(line)
                    elif copy_extra:
                        saved_lines.append(line)

                if line.startswith("ORIGIN") and save_origine:
                    if save_CDS:
                        saved_lines.extend(saved_lines_CDS)
                        save_CDS=False
                    else:
                        saved_lines.append(line)
                    copy_extra=True
                    save=True

#Printing out the last line of the locus with enriched orthologue(s)
            if save:
                for s in saved_lines:
                    print(s, file= fout)
    else:
        non_in_enrcihed.add(Geno)

with open(os.path.join(args.o,"Genomes_without_selected_orthologs.txt"), "w") as fout2 :
        print("#Genomes without enriched orthologues:", file=fout2)
        for e in non_in_enrcihed:
            print(e, file=fout2)
