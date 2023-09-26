import os
import argparse

usage = 'python DNA_selected_orthologues.py -i -d -g -o'
description = 'This program extracts the DNA sequences of selected orthologues per annotation'


parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i',dest='i',help='input file of selected orthologues, .tsv', default='Selected_orthologues_annotation.tsv')
parser.add_argument('-d',dest='d',help='genes directory', default='Genes')
parser.add_argument('-g',dest='g',help='Orthogroups file .tsv', default='phyloglm_input/Orthogroups.tsv')
parser.add_argument('-o', dest='o', help='Output directory', default='Selected_Orthologues_DNA_sequences')

args = parser.parse_args()

def extract_orthologues(fi):
    set_ort=set()
    with open(fi, "r") as fin:
      for line in fin:
        line=line.rstrip()
        if not line.startswith("#"):
          set_ort.add(line)
    return set_ort

def genes_per_ortholog(fils, mydict):
#    og={}
    go={}
    with open(fils, "r") as fin:
        for line in fin:
            line=line.rstrip()
            if not line.startswith("Orthogroup"):
                #OG0000000       PJFELKFM_02705, PJFELKFM_02712  IMDEMLBJ_01706, IMDEMLBJ_01707, IMDEMLBJ_01709, IMDEMLBJ_01712  PONCNELG_00563  JIJPHLGC_03615  NHPIGIFA_00284, NHPIGIFA_00717, NHPIGIFA_00721, NHPIGIFA_00722, NHP
                line=line.split("\t")
                if line[0] in mydict:
                    genes_list=[o for o in line[1:] if o != ""]
                #    og[line[0]]=genes_list
                    for g in genes_list:
                        go[g]=line[0]
    #return og, go
    return go

def get_sequences(fil, gene_dict,ort_N_seqs):

    with open(fil, "r") as fin:
        tocopy=False
        for line in fin:
            line=line.rstrip()
            if line.startswith(">"):
                if tocopy:
                    with open(fileout, "a") as fo:
                        print("{}\n{}".format(header, seqs), file=fo)
                        tocopy= False
                        if orthol in ort_N_seqs:
                            ort_N_seqs[orthol]+=1
                        else:
                            ort_N_seqs[orthol]=1


                #>KADHOGOG_00001 hypothetical protein
                copyline=line
                line=line.split(" ")
                gene=line[0][1:]
                if gene in gene_dict:
                    #annt=" ".join(line[1:])
                    orthol=gene_dict[gene]
                    #annotations=ann_dict[orthol]
                    #for a in annotations:
                    fileout= os.path.join(args.o, str(orthol+".fna"))
                    header=copyline
                    tocopy=True
                    seqs=""
            else:
                if tocopy:
                    seqs+=line



        if tocopy:
            with open(fileout, "a") as fo:
                print("{}\n{}".format(header, seqs), file=fo)
                #tocopy= False
                if orthol in ort_N_seqs:
                    ort_N_seqs[orthol]+=1
                else:
                    ort_N_seqs[orthol]=1
    return ort_N_seqs

##main
file_list=[os.path.join(args.d, f) for f in os.listdir(args.d) if f.endswith(".ffn")]
Ort_Annot=extract_orthologues(args.i)
#Ort_genes, gene_ort=genes_per_ortholog(args.g, Ort_Annot)
gene_ort=genes_per_ortholog(args.g, Ort_Annot)

if not os.path.exists(args.o):
    os.mkdir(args.o)

O_Seqs={}
for g in file_list:
    O_Seqs=get_sequences(g, gene_ort,O_Seqs)

with open(os.path.join(args.o, "report.txt"), "w") as fou:
    print("#Orthologue\tNumber_of_sequences", file=fou)
    for k,v in O_Seqs.items():
        print("{}\t{}".format(k,v), file=fou)
