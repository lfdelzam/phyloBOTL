#!/usr/bin/env python3
import os
import argparse
import re
import pathlib


usage = 'python parse_genomad.py [options]'
description = 'This program print out a table indicating plamid contigs and when present, Enriched Orthologs from Isolates'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i',dest='i',help='directory Genomad_output', required=True)
parser.add_argument('-f',dest='f',help='Annotations/Orthologues/orthologues_gff.tsv', required=True)
parser.add_argument('-e',dest='e',help='Enriched orthologs list, Vv_Enriched.txt', required=True)
parser.add_argument('-o', dest='o', help='output directory', required=True)
parser.add_argument('-c', dest='c', help='score cutoff', default=0.7)
args = parser.parse_args()


def Enriched_list(filein):
    enriched=set()
    with open(filein, "r") as fin:
        for line in fin:
            line=line.rstrip()
            if not line.startswith("#"):
                enriched.add(line)
    return(enriched)


def get_contigs(filein, region):
    contigs=set()
    size=0
    with open(filein, "r") as fin:
    #seq_name	length	topology	n_genes	genetic_code	plasmid_score	fdr	n_hallmarks	marker_enrichment	conjugation_genes	amr_genes
    #JZES01000051.1	5654	No terminal repeats	7	11	0.9633	NA	0.7537	NA	NF033693
        for line in fin:
            line=line.rstrip()
            line=line.split("\t")
            if region == "plasmid" and not line[0].startswith("seq_name"):
                if float(line[5]) > args.c:
                    contigs.add(line[0])
                    size += int(line[1])
            if region == "virus" and not line[0].startswith("seq_name"):
                if float(line[6]) > args.c:
                    contigs.add(line[0].split("|")[0])
                    size += int(line[1])
    return(contigs, size)


def Orthologues(filein):
    Orthos={}
    with open(filein, "r") as fin:
#Genome	Orthologue	Contig	Type	Start	Stop	Direction	Description
#NV1	OG0004874	JZET01000001.1	CDS	152	397	+	ID=LMCCGNMF_00001;inference=ab initio prediction:Prodigal:002006;locus_tag=LMCCGNMF_00001;product=hypothetical protein
        for line in fin:
            line=line.rstrip()
            line=line.split("\t")
            if not line[0].startswith("Genome"):
                Genome=line[0]
                Ortho=line[1]
                contig=line[2]
                if Genome in Orthos:
                    if contig in Orthos[Genome]: #and Ortho.startswith("OG"):
                        Orthos[Genome][contig].add(Ortho)
                    #elif Ortho.startswith("OG"):
                    else:
                        Orthos[Genome][contig] = set([Ortho])
                else:
                #elif Ortho.startswith("OG"):
                    Orthos[Genome]={ contig: set([Ortho]) }

    return(Orthos)

def print_out(contig_fil, Orthologs, Enriched_Orthologs_list, region):
    Genomes_with_region={}
    Region_size={}
    for fl in contig_fil:
        first_line=True
        first_line_e= True
        fl_bs= os.path.basename(fl)
        Genome=fl_bs.split("_"+region+"_summary.tsv")[0]
        plasmid_contigs, total_size=get_contigs(fl, region)
        Genomes_with_region[Genome]=[c for c in plasmid_contigs ]
        Region_size[Genome]=total_size
#        for c in plasmid_contigs:
#                if first_line:
#                    with open(os.path.join(args.o, Genome+"_"+region+"_contigs_orthologues.tsv"), "w") as fout:
#                        print("#Genome\tContig\torthologs_in_contig", file=fout)
#                    first_line=False

#                with open(os.path.join(args.o, Genome+"_"+region+"_contigs_orthologues.tsv" ), "a") as fout:
#                        print("{}\t{}\t{}".format(Genome, c,",".join(Orthologs[Genome][c])), file=fout)
#
#                enr_o=[ e for e in Orthologs[Genome][c] if e in Enriched_Orthologs_list ]
#
#                if len(enr_o) > 0:
#                    if first_line_e:
#                        with open(os.path.join(args.o, Genome+"_"+region+"_contigs_Enriched_orthologues.tsv" ), "w") as fout:
#                            print("#Genome\tContig\tEnriched_orthologs_in_contig", file=fout)
#                        first_line_e=False

#                    with open(os.path.join(args.o, Genome+"_"+region+"_contigs_Enriched_orthologues.tsv" ), "a") as fout:
#                        print("{}\t{}\t{}".format(Genome, c,",".join(enr_o)), file=fout)


    with open(os.path.join(args.o,"Genomes_with_"+region+".tsv" ), "w") as fout2, open(os.path.join(args.o,"Genomes_contigs_without_orthologues_"+region+".tsv" ), "w") as fout3 :
        print("#Genome\tN_"+region+"_Contigs\tTotal_size_bs\tEnriched_orthologs_in_"+region+"\t"+region+"_Contigs", file=fout2)
        for g,cs in Genomes_with_region.items():
            enr_ot=[]
            for c in cs:
                if c in Orthologs[g]:
                    check=[ e for e in Orthologs[g][c] if e in Enriched_Orthologs_list ]
                    if len(check)>0:
                        enr_ot.extend(check)
                else:
                    print("Genome:{}\tContig:{}\tregion:{}".format(g,c,region), file=fout3)


            print("{}\t{}\t{}\t{}\t{}".format(g, len(cs),Region_size[g], ",".join(enr_ot),",".join(cs)), file=fout2)




Orthologs=Orthologues(args.f)
Enriched_Orthologs_list=Enriched_list(args.e)

contig_files = [f for f in pathlib.Path(args.i).glob("Genomad_output_*/*_summary/*_plasmid_summary.tsv") ]
contig_files_virus = [f for f in pathlib.Path(args.i).glob("Genomad_output_*/*_summary/*_virus_summary.tsv") ]

if not os.path.exists(args.o):
    os.makedirs(args.o)

print_out(contig_files, Orthologs, Enriched_Orthologs_list, "plasmid")
print_out(contig_files_virus, Orthologs, Enriched_Orthologs_list,"virus")
