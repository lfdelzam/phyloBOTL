rm(list = ls())
options(warn=-1)

list_of_packages <- c("ape","tidyverse","gggenomes", "igraph","tidytree","gridExtra","future","readr","ggtree", "ggnewscale", "randomcoloR")

suppressPackageStartupMessages(library(argparser))

# arguments
p <- arg_parser("synteny_visual.R")
p <- add_argument(p, "-a", help="Path to orthologues_gff.tsv", default="Annotations/Orthologues/orthologues_gff.tsv")
p <- add_argument(p, "-m", help="Minimum Weigth in graph to be included", default=4)#4 esta bien
p <- add_argument(p, "-t", help="Path to cleaned tree, .txt", default="Results_415_with_IQTREE_rooted/tree/Vv_Cleaned_tree.txt") #must bee the cleaned one
p <- add_argument(p, "-s", help="loci size to look for co-location", default=40000) #40K es optimo
p <- add_argument(p, "-f", help="max number of genomes per figure", default=22)
p <- add_argument(p, "-o", help="output directory", default="Results_415_with_IQTREE_rooted/BORRELO")
p <- add_argument(p, "-l", help="Label of Selected Orthologues in figures", default=" Enriched")
p <- add_argument(p, "-e", help="Path to enriched/depleted orthologues list", default="Results_415_with_IQTREE_rooted/Vv_Enriched.txt")
p <- add_argument(p, "-i", help="genomes metadata", default="GENOME_LIST")
p <- add_argument(p, "-l1", help="Group 1 label", default="C")
p <- add_argument(p, "-g1", help="Group 1 name", default="Clinical")
p <- add_argument(p, "-g2", help="Group 2 name", default="Environmental")

argv <- parse_args(p)

argv <- parse_args(p)

#libraries
for (s in list_of_packages) { suppressPackageStartupMessages(library(s, character.only = TRUE))}
plan(multisession) #multithread


#functions

visualizing_ortholgues_in_cluster<-function(dfP, clster, prefixT) {

  SIZE=3
  palette <- distinctColorPalette(length(unique(dfP$Cluster))-1)

  p1<-ggtree(tree, layout="fan", open.angle = 0, ladderize=TRUE, branch.length="branch.length") +
    geom_tiplab(size=1.4, align=T, linesize=.05, offset=0.003)+
    geom_text2(aes(subset = !isTip, label=label), size=SIZE,angle = 0,vjust = "inward", hjust = "inward",check_overlap = TRUE)+
    geom_nodepoint(color="black", alpha=1/3, size=0.5 )+ geom_rootedge(rootedge = 0.005)

  NA_vector=rep(NA,(length(p1[["data"]][["label"]]))-length(dfP$Cluster))
  new_labels=c(gsub("Absent", "", dfP$Cluster), NA_vector)
  Pa <- p1 +geom_tiplab2(size=1.5, align=T,geom="text",linesize=0.0, offset=0.0021,aes(label=new_labels))

  P2<-gheatmap(Pa, dfP[,1, drop=FALSE] , offset=0, width=0.02, colnames=F) +scale_fill_manual(
    values = c("black","blue"), name="Isolate")

  P3 <- P2 + new_scale_fill()

  P4<-gheatmap(P3, dfP[,2, drop=FALSE], offset=0.001, width=0.02, colnames=F) +
    scale_fill_manual(values = c(palette,"white"),
                      name=paste0("Cluster ",clster, " subgroups")) #+theme(legend.position="bottom")

  ggsave(paste(prefixT,paste("cluster",clster, "tree_fan_ladderized_branchlength.pdf", sep="_"), sep="/"), P4, width = 42, height = 60, units = "cm", device = "pdf")

}


region_in_genome<-function(geno,tDF) {
  #geno="SE_VV_18_05"  #maybe if a convert tDF to data.frame works with allocated memory
  #tDF=O2 names(O2)
  contig_where_loci=unique(tDF$seq_id[tDF$bin_id == geno])
  outdf=as.data.frame(matrix(nrow=0, ncol=ncol(O1)))

  for (c in contig_where_loci) {

    #c="NZ_JALQSE010000018.1"
    sudtDF=subset(tDF,tDF$bin_id == geno & tDF$seq_id == c )
    start_p=min(tDF$start[tDF$bin_id == geno & tDF$seq_id == c ]) #-5000
    end_p=max(tDF$end[tDF$bin_id == geno & tDF$seq_id == c ]) #+5000
    Orts_IN_genome=subset(O1, O1$bin_id == geno)
    Orts_in_contig=subset(Orts_IN_genome, Orts_IN_genome$seq_id == c )
    Orts_in_contigs_where_loci_is=subset(Orts_in_contig, Orts_in_contig$start >= start_p & Orts_in_contig$end <= end_p)
    outdf=rbind(outdf,Orts_in_contigs_where_loci_is)
  }
  return(outdf)
}

plot_synteny<-function(c1,loci1,O1d, texto, pout) {
#  c1=1
#  loci1=cls_leiden[[c1]]
#  O1d=O1
  #pout=FALSE
  if (length(loci1) >0) {

    cat("Orthologues in cluster", c1, ":", loci1, "\n", file = paste(argv$o,paste("Leiden_clusters_using_max_loci_size",max_loci,"orthologues.txt" ,sep="_"), sep="/"),
        append = T)
    if (pout) {
      cat(">",c1,"\n", file = paste(argv$o,paste("Leiden_clusters_using_max_loci_size",max_loci,"orthologues.ft",sep="_"), sep="/"),
          append = T, sep="")
      for (k in loci1) {
      cat(k,"\n", file = paste(argv$o,paste("Leiden_clusters_using_max_loci_size",max_loci,"orthologues.ft",sep="_"), sep="/"),
          append = T)
      }
    }

    O2d=subset(O1d, O1d$feat_id %in% loci1)

    averd=O2d[names(O2d) %in% c("bin_id","feat_id", "strand", "seq_id")]

    list_of_genomesd=as.vector(unique(O2d$bin_id))

    references=rep("0", length(list_of_genomesd))
    similars=vector("list", length(list_of_genomesd)) #list()
    names(similars)<-list_of_genomesd
    k=1
    while ( length(list_of_genomesd) >0 ) {

      aver1d=region_in_genome(list_of_genomesd[1],O2d)

      order1=as.vector(aver1d$feat_id)
      strand1=as.vector(aver1d$strand)
      strand1=gsub("\\+",1,strand1)
      strand1=gsub("\\-", -1, strand1)
      strand1=as.numeric(strand1)
      cluster1=list_of_genomesd[1]
      if (length(list_of_genomesd) >1) {
      for (i in 2:length(list_of_genomesd)) {
        #i=2
        if (!list_of_genomesd[i] %in% cluster1 && !is.na(list_of_genomesd[i]) ) {

          gr2=region_in_genome(list_of_genomesd[i],O2d)
          order2=as.vector(gr2$feat_id)
          strand2=as.vector(gr2$strand)
          strand2=gsub("\\+",1,strand2)
          strand2=gsub("\\-", -1, strand2)
          strand2=as.numeric(strand2)
          if ( (identical(order1, order2) && identical(strand1, strand2)) || ( (identical(rev(order1), order2)) && identical(rev(strand1*(-1)), strand2)) ) {
            cluster1=c(cluster1,list_of_genomesd[i])
          }
        }
      }
      }
      r=cluster1[1]
      references[k]=r
      k=k+1
      similars[[r]] = cluster1

      list_of_genomesd=list_of_genomesd[!list_of_genomesd %in% cluster1]

    } #end while
    references=references[references != "0"]
    similars=similars[lapply(similars,length)>0]

    sink(paste(argv$o,paste("group_of_genomes",texto,"in_figure.txt", sep="_"),sep="/"))
    print(similars)
    sink()

    phyloorder=Y$label[1:ngen][!is.na(match(Y$label[1:ngen], references))]

    list_of_bins=vector("character",length = length(phyloorder))

    fdfd=data.frame(matrix(nrow=0,ncol=ncol(O1d)))
    emale_nagF=data.frame(matrix(nrow=0,ncol=ncol(emale_nag)))
    for (j in 1:length(phyloorder)){
      #j=1
      g=phyloorder[j]
      tempdf=region_in_genome(g,O2d)
      ngeno=length(similars[[g]])
      tempdf$bin_id=paste(tempdf$bin_id," (",ngeno,")", sep=" ")
      fdfd=rbind(fdfd,tempdf)
      temnag=subset(emale_nag, emale_nag$bin_id == g)
      temnag$bin_id=paste(temnag$bin_id," (",ngeno,")", sep=" ")
      emale_nagF=rbind(emale_nagF,temnag)
      list_of_bins[j]=paste(g," (",ngeno,")", sep=" ")
    }

    fdfd$bin_id=factor(fdfd$bin_id, levels = unique(fdfd$bin_id))
    fdfd$feat_id <- factor(fdfd$feat_id, levels = unique(fdfd$feat_id))
    fdfd$strand <- factor(fdfd$strand, levels = unique(fdfd$strand))
    emale_nagF$bin_id<- factor(emale_nagF$bin_id, levels = unique(emale_nagF$bin_id))
    ### if N genomes > min_genome_per_fig
    if (length(references) >ste) {
      f_numb=1
      for (counter in seq(1,length(references),ste)) {
        end_v=(counter+ste-1)
        if (end_v > ngen) {end_v = length(references) }
        # counter=1
        #  end_v=10
        flocsub <- gggenomes(genes=fdfd[fdfd$bin_id %in% list_of_bins[counter:end_v] ,],spacing = 0.05, theme ="clean") %>%
          add_feats(ngaros=emale_nagF)+
          geom_seq() +  geom_bin_label(size=3,stat="identity") +
          geom_gene(aes(fill=feat_id),position="strand",stat="identity",  na.rm=TRUE, show.legend = F) +
          geom_feat_tag(aes(label=feat_id, color=type), angle = 20, nudge_y=0.2, check_overlap = TRUE, size=2, na.rm=TRUE,
                        key_glyph = "rect")+
          scale_color_manual(values=c("darkblue","grey", "red", "black", "purple"),
                             aesthetics = "colour",
                             name="Orthologue label type")+
          theme(legend.position="bottom",legend.key.size = unit(0.2, 'cm'))

        ggsave(paste(argv$o,paste(texto,f_numb,"by_ref_genomes","figure.pdf", sep="_"),sep="/"), flocsub, width = 32, height = 60, units = "cm", device = "pdf")
        f_numb=f_numb+1
      }
    } else {

        floc=gggenomes(genes=fdfd,spacing = 0.05, theme ="clean") %>%
        add_feats(ngaros=emale_nagF)+
        geom_seq() +  geom_bin_label(size=3,stat="identity") +
        geom_gene(aes(fill=feat_id),position="strand",stat="identity",  na.rm=TRUE, show.legend = F) +
        geom_feat_tag(aes(label=feat_id, color=type), key_glyph = "rect", angle = 20, nudge_y=0.2, check_overlap = TRUE, size=2, na.rm=TRUE)+
        scale_color_manual(values=c("darkblue","grey", "red", "black", "purple"),aesthetics = "colour",name="Orthologue label type")+
        theme(legend.position="bottom",legend.key.size = unit(0.2, 'cm'))

        ggsave(paste(argv$o,paste(texto,"by_ref_genomes","figure.pdf", sep="_"),sep="/"), floc, width = 32, height = 60, units = "cm", device = "pdf")
       }
  }
  return(similars)
}

find_co_locations <- function(genomesDF) {
  #genomesDF=genomeDFall
  co_locations=vector("list", nrow(genomesDF)) #list()
  bloq=1
  for (g in unique(genomesDF$genome)) {  #same genome
    #g=unique(genomesDF$genome)[1]
    sdf_g=subset(genomesDF, genomesDF$genome == g)
    cat("INFO: finding synteny in genome:",g, "\n      Contings:\n")
    for (c in unique(sdf_g$Contig)) { #same contig
      #c=unique(sdf_g$Contig)[1]
      cat(c,",", sep="")
      sdf_c=subset(sdf_g, sdf_g$Contig == c)
      if (length(sdf_c$Orthologue) >1) {
        i=1
        while (i < length(sdf_c$Orthologue)) {
          o=sdf_c$Orthologue[i]
          end1=max(sdf_c$end[sdf_c$Orthologue == o])
          sdf_bloq=subset(sdf_c, sdf_c$end >= end1 & sdf_c$end <= max_loci + end1) #una sola direction -->
          #  ultimo=sdf_bloq$Orthologue[nrow(sdf_bloq)]

          if (length(sdf_bloq$Orthologue) >1) {
            co_locations[[bloq]] = unique(sdf_bloq$Orthologue)
            bloq=bloq +1
          } #if more than one orthologue in 40K
          i=i+1 #max(which(sdf_c$Orthologue == ultimo))+1  #maybe I need to modify this, and go directly to the next orthologue
        } #close while
      } #close if two orthologues in contig

    } #close contig lopp
    cat("\n")
  } # close genome loop
  co_locations=co_locations[lapply(co_locations,length)>0]
}

create_network_table <-function(co_locations){
  #co_locations=co_location
  co_tables=data.frame(matrix(ncol=2, nrow=length(co_locations)*15*2))
  names(co_tables)=c("Ortholog1","Ortholog2")
  k=1
  for (n in 1:length(co_locations) ) {
    #n=1
    for (i in 1:(length(co_locations[[n]])-1)) {
      c1=co_locations[[n]][i]
      for (i2 in (i+1):length(co_locations[[n]])) {
        c2=co_locations[[n]][i2]
        #co_tables[nrow(co_tables)+1,]=c(c1,c2)
        co_tables[k,]=c(c1,c2)
        k=k+1
      }
    }

  }
  co_tables=na.omit(co_tables)
  return(co_tables)
}

create_graph_input <- function(co_tabl) {
  #co_tabl=co_table
#  co_table_counts=data.frame(matrix(ncol=3, nrow=0))
  co_table_counts=data.frame(matrix(ncol=3, nrow=50*50))
  names(co_table_counts)=c("Ortholog1","Ortholog2", "weigth")
  k=1
  cat("INFO: Generating graph input\n")
  for (u1 in unique(co_tabl$Ortholog1)) {
    for (u2 in unique(co_tabl$Ortholog2)) {
      count_coloc=nrow(subset(co_tabl, co_tabl$Ortholog1 == u1 & co_tabl$Ortholog2 == u2))
      if (count_coloc >argv$m) {
        #co_table_counts[nrow(co_table_counts)+1,]=c(u1,u2,count_coloc)
        co_table_counts[k,]=c(u1,u2,count_coloc)
        k=k+1
      }
    }
  }
  co_table_counts=na.omit(co_table_counts)
  return(co_table_counts)
}

#Reading info
list_enriched=read.csv(argv$e)
names(list_enriched)="Enriched"
list_enriched$Enriched=gsub(" ","",list_enriched$Enriched)

genomeDFWall <- read_tsv(argv$a, col_names = TRUE, show_col_types = FALSE) #read.table(argv$a, header=T, sep="\t") doesn't work
genomeDFall <- genomeDFWall[genomeDFWall$Orthologue %in% list_enriched$Enriched,]
genomeDFall <- genomeDFall %>% rename(genome=Genome, strand=Direction,start=Start, end=Stop, attributes=Description, type=Type)
genomeDFall <- genomeDFall[,!names(genomeDFall) %in% c("type")]

dir.create(argv$o, showWarnings = FALSE)
#### Find co-location of enriched orthologues
max_loci=argv$s
if (!file.exists(paste(argv$o,paste(max_loci,"graph.rds", sep="_"), sep="/"))) {

  co_location <- find_co_locations(genomeDFall)

  cat("INFO: Creating network input table\n")

  co_table %<-% {create_network_table(co_location)}
  co_table <- as.data.frame(co_table)

  co_table_count %<-% {create_graph_input(co_table)}
  co_table_count <- as.data.frame(co_table_count)


  cat("INFO: Generating graph and Leiden clustering\n")
  gra=graph_from_data_frame(co_table_count, directed = FALSE)
  saveRDS(gra, file = paste(argv$o,paste(max_loci,"graph.rds", sep="_"), sep="/"))
} else {
gra=readRDS(paste(argv$o,paste(max_loci,"graph.rds", sep="_"), sep="/"))
}
set.seed(1)
cls_leiden=cluster_leiden(
  gra,objective_function="modularity", # objective_function = c("CPM", "modularity"), #  weights = NULL,
  resolution_parameter = 1, beta = 0.01, n_iterations = 100, # vertex_weights = NULL, #  initial_membership = NULL
)


pdf(paste(argv$o,paste("Leiden_clusters_using_max_loci_size",max_loci, sep="_"), sep="/"))
plot(cls_leiden, gra,layout=layout_nicely,
     vertex.label.cex = 0.6,
     vertex.size=8,vertex.label.dist=0.5, edge.arrow.size=0.5, rescale=TRUE,
     main = paste(argv$l,"Orthologues"), sub= "Co-localization clusters")

dev.off()
####
for (c in 1:cls_leiden$nb_clusters) {
  if (length(cls_leiden[[c]]) > 9) {
          sg = induced_subgraph(gra, cls_leiden[[c]])
          sub=cluster_leiden(sg,objective_function="modularity")
          pdf(paste(argv$o,paste("Leiden_clusters_using_max_loci_size",max_loci,"subgraph_cluster",c, sep="_"), sep="/"))
          plot(sub,sg, cex=0.6)
          dev.off()

  }
}
###
tree=read.tree(argv$t)
Y=as_tibble(tree)


emale_genes_all <- genomeDFWall %>% rename(seq_id=Contig, bin_id=Genome, feat_id=Orthologue, strand=Direction,
                                           start=Start, end=Stop, attributes=Description, type=Type)

emale_nag=emale_genes_all
emale_nag$type[emale_genes_all$feat_id %in% unique(genomeDFall$Orthologue) ]=argv$l #" Enriched"

emale_genes_all=emale_genes_all[,!names(emale_genes_all) %in% c("Type","attributes")]
ngen=length(unique(emale_genes_all$bin_id)) #415
ste=argv$f

###-- Form here is new
cat("INFO: Reading genome/sample information \n")
genome_table <- read.csv(argv$i, sep = ",", na.strings = "NA", header = FALSE, col.names = c("Strain", "File","Type","Special"))
genome_table %<-% {mutate_all(genome_table, function(x) gsub("\\s+"," ",x))}
#Clinical Isolates names
PATHOGENIC_STRAINS=genome_table[,1][genome_table[,3] == argv$l1]
y=c()
for (n in tree$tip.label) {
  if (n %in% PATHOGENIC_STRAINS) y=c(y,1) else y=c(y,0)
}

if(length(y[y>0]) != length(PATHOGENIC_STRAINS)) stop("Strain names and tree labels don't match, please check spelling in GENOME_LIST file")

tipo<-y
tipo[tipo==0]<-argv$g2
tipo[tipo==1]<-argv$g1
DFP <- data.frame(Isolates=as.character(tipo))
rownames(DFP) <- tree$tip.label
###-- up to here is new


cat("INFO: Co-location - Output per similar genomes\n")

O1=subset(emale_genes_all, emale_genes_all$bin_id %in% Y$label[1:ngen])
O1$feat_id <- factor(O1$feat_id, levels = unique(O1$feat_id))
O1$strand <- factor(O1$strand, levels = unique(O1$strand))
O1$bin_id <- factor(O1$bin_id, levels = unique(O1$bin_id))
for (c in 1:cls_leiden$nb_clusters) {
  #c=1
  cat("INFO: Cluster -",c,"\n")

  loci=cls_leiden[[c]]
  #plot_synteny(c,loci,O1,paste("grouped_cluster_leiden", c, sep="_"), pout=TRUE)
  ###-- Form here is new
  cluster_group=plot_synteny(c,loci,O1,paste("grouped_cluster_leiden", c, sep="_"), pout=TRUE)
  DFP_c=DFP
  DFP_c$Cluster="Absent"
  for (z in 1:length(cluster_group)  ) {
    DFP_c$Cluster[ row.names(DFP_c) %in% cluster_group[[z]] ]=z
  }
  visualizing_ortholgues_in_cluster(DFP_c, c, vis_outdir)

  ###-- up to here is new
}

cat("INFO: co-location - all genomes\n")

dirall=paste(argv$o, "All_genomes", sep="/")
dir.create(dirall, showWarnings = FALSE)

#ploting all the genomes
for (c in 1:cls_leiden$nb_clusters) {
  #c=1
  loci=cls_leiden[[c]]
  if (length(loci) >1){
    texto=paste("cluster_leiden", c, sep="_")
    fnumb=1
    for (counter in seq(1,ngen,ste)) {
      end_v=(counter+ste-1)
      if (end_v > ngen) {end_v = ngen}
      #  counter=1
      # end_v=10
      O1=subset(emale_genes_all, emale_genes_all$bin_id %in% Y$label[counter:end_v])
      O1$feat_id <- factor(O1$feat_id, levels = unique(O1$feat_id))
      O2=subset(O1, O1$feat_id %in% loci)

      if (nrow(O2)>0) {
        phyloorder=Y$label[counter:end_v]

        #----

        fdfa=data.frame(matrix(nrow=0,ncol=length(names(O1))))
        for (g in phyloorder){
          fdfa=rbind(fdfa,region_in_genome(g,O2))
        }

        fdfa$bin_id=factor(fdfa$bin_id, levels = unique(fdfa$bin_id))
        fdfa$feat_id <- factor(fdfa$feat_id, levels = unique(fdfa$feat_id))
        fdfa$strand <- factor(fdfa$strand, levels = unique(fdfa$strand))
        #---

        floca=gggenomes(genes=fdfa,spacing = 0.05, theme ="clean") %>%
          add_feats(ngaros=emale_nag)+
          geom_seq() +  geom_bin_label(size=3,stat="identity") +
          geom_gene(aes(fill=feat_id),position="strand",stat="identity",  na.rm=TRUE, show.legend = F) +
          geom_feat_tag(aes(label=feat_id, color=type, key_glyph = "rect"), angle = 20, nudge_y=0.2, check_overlap = TRUE, size=2, na.rm=TRUE)+
          scale_color_manual(values=c("darkblue","grey", "red", "black", "purple"),
                             name="Orthologue label type")+
          theme(legend.position="bottom", legend.key.size = unit(1, 'cm'))

        ggsave(paste(dirall,paste("genomes",texto,fnumb,"figure.pdf", sep="_"),sep="/"), floca, width = 32, height = 60, units = "cm", device = "pdf")
        fnumb=fnumb+1
      }

    } # close  for (counter in seq(1,ngen,ste))

  }
}
