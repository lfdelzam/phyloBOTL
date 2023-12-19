rm(list = ls())
options(warn=-1)

suppressPackageStartupMessages(library(argparser))

# arguments
p <- arg_parser("phylocomgen")
p <- add_argument(p, "-t", help="phylogenetic tree file", default="Tree_rooted_boot.treefile")
p <- add_argument(p, "-g", help="gene/orthologue count file", default="Orthogroups.GeneCount.tsv")
p <- add_argument(p, "-i", help="genomes metadata", default="GENOME_LIST")
p <- add_argument(p, "-e", help="Special group", default="Baltic Sea")
p <- add_argument(p, "-s", help="Orthogroups list file", default="phyloglm_input/Orthogroups.tsv")
p <- add_argument(p, "-m", help="list of protein annotation from prokka file", default="phyloglm_input/Annotations.txt")
p <- add_argument(p, "-l", help="outfiles prefix", default="Vv")
p <- add_argument(p, "-o", help="output directory", default="OUTPUT")
p <- add_argument(p, "-r", help="FILTER: ratio of orthologue presence[absence] in the genome dataset,
                  if higher than this value, the orthologues is not considered in the analysis", default=0.95)
p <- add_argument(p, "-b", help="phyloglm Bootnumber", default=0)
p <- add_argument(p, "-a", help="phyloglm btol number", default=10)
p <- add_argument(p, "-q", help="p.adj.value cutoff", default=0.05)
p <- add_argument(p, "-c", help="number of clusters", default=2)
p <- add_argument(p, "-N", help="FastANI output file", default="")
p <- add_argument(p, "-k", help="Path to gff files", default="Annotations/GFF_files")
p <- add_argument(p, "-y", help="Group 1 name", default="Clinical")
p <- add_argument(p, "-w", help="Group 1 label", default="C")
p <- add_argument(p, "-z", help="Group 2 name", default="Environmental")
p <- add_argument(p, "-v", help="Group 2 label", default="E")


argv <- parse_args(p)

##extra parameter - on development
core_fraction=0.97
Pathogenicity_probability=F
special_nodes = "none" #c(418, 419, 630, 638, 773)
#labels_nodes=c("root","L5", "L2","L3","L1","L4")
#type_name_special_nodes="Lineage"
labels_nodes_colors=c("gray","purple","black","red","blue","darkgreen", "orange", "yellow")
##

#libraries
list_of_packages <- c("dplyr","ggtree","ggplot2", "phylolm", "factoextra", "cluster",
                      "FactoMineR", "future","umap","tidyr","ape","phangorn","pheatmap",
                      "vegan","ggrepel","tidyverse", "data.table", "RColorBrewer")
for (s in list_of_packages) { suppressPackageStartupMessages(library(s, character.only = TRUE))}

plan(multisession)

#functions

clustersout<-function(clust, vec, list_pathog_str) {
  cat("\n ---------------- \n")
  cat("******** Cluster", clust, "\n")
  PAM_patho=names(vec$clustering[vec$clustering == clust])
  unknown_in_cluster=rep("Na", length(PAM_patho))
  patho_in_cluster=rep("Na", length(PAM_patho))
  cat("Total strains in cluster",clust, ":" ,length(PAM_patho), "\n")
  for (w in PAM_patho) {
    if (w %in% list_pathog_str) {
      patho_in_cluster=c(patho_in_cluster,w)
    } else {
      unknown_in_cluster=c(unknown_in_cluster,w)
    }
  }
  patho_in_cluster=patho_in_cluster[patho_in_cluster != "Na"]
  unknown_in_cluster=unknown_in_cluster[unknown_in_cluster != "Na"]

  cp=length(patho_in_cluster)
  Np=length(unknown_in_cluster)
  cat(paste0("Known ",argv$y, ": "), cp, "->", cp*100/length(list_pathog_str), paste0("% of the known ",argv$y, " strains"), "\n List:", patho_in_cluster,"\n")
  cat("Unknown:", Np, "\n List ",unknown_in_cluster,"\n")
}

print_out_annot<-function(vect_features, Orthogroups, annots) {
  #vect_features=present_markets
  #Orthogroups=Orthogroups_g
  #annots=annots_g
  anot_table=data.frame()
  #anot_table=as.data.table(matrix(nrow=length(vect_features), ncol=ncol(Orthogroups)+1))

  for (s in vect_features) {
    lista=unlist(Orthogroups[Orthogroups[[1]] == s,])
    lista=lista[! lista %in% c("",s)]
    lista=lista[!is.na(lista)]
    Annotations=vector(length = length(lista))
    for (i in lista) {
      Annotations[i]=annots[[2]][annots[[1]] == i]
    }
    Annotations=Annotations[ Annotations != FALSE ]
    if (length(Annotations) != 0 ) {

    anot_table=rbind(anot_table,data.frame(Orthologue=s, as.data.frame(table(Annotations))))
    }
  }
  return(anot_table)
}

unsupervised_learning<-function(texto, features, mDF ) {
 # texto="Unsupervise"
#  features=selected_features
 # mDF=df
    dirU=paste(argv$o, texto, sep="/")
    dir.create(dirU, showWarnings = FALSE)
    prefixU=paste(dirU,argv$l, sep="/")
    test=mDF[,features, drop=FALSE]
    test[test>1]=1

    dist_mat <- dist(test, method = "binary")
    try(hclust_avg <- hclust(dist_mat, method = 'average'), silent = TRUE)
    if (exists("hclust_avg")) {

      pdf(paste(prefixU,"hclust.pdf", sep = "_"))
      plot(hclust_avg, hang = -1, cex = 0.3)
      dev.off()

     treehc=as.phylo(hclust_avg)
     SIZE=3
     p<-ggtree(treehc, layout="fan", open.angle = 0, ladderize=TRUE, branch.length="none") +   #"fan" "none"
       geom_tiplab(size=1.8, align=T, linesize=.3, offset=2)
     P1<-gheatmap(p, dfP, offset=0, width=.01,legend_title = "Isolates", colnames = FALSE)
     ggsave(paste(prefixU,"hclust_as_tree.pdf", sep = "_"), P1, width = 84, height = 120, units = "cm", device = "pdf")


     if (exists("subSELDFPATHO")) {
     p<-ggtree(treehc, layout="fan", open.angle = 5, ladderize=TRUE, branch.length="none") +   #"fan" "none"
       geom_tiplab(size=2.5, align=T, linesize=.3, offset=12)
     P2<-gheatmap(p, subSELDFPATHO, offset=0, width=0.04,legend_title = "", colnames=T, colnames_position="top",
                   colnames_angle=30,hjust=0, font.size=3) +scale_fill_manual(
                     values = c("white","blue","black","white", brewer.pal(n = 9, name = "Reds")), name="")

     ggsave(paste(prefixU,"hclust_as_tree_with_probabilities.pdf", sep = "_"), P2, width = 84, height = 120, units = "cm", device = "pdf")
     }

    }
  # Cutting tree by no. of clusters
      try(Members_in_clusters <- cutree(hclust_avg, k = 2 ), silent = TRUE)
      if (exists("Members_in_clusters")) {
      sink(paste(prefixU,"hclust_table.txt", sep = "_"))
      print(table(Members_in_clusters))
      cat("\nMembership: \n")
      print(Members_in_clusters)
      sink()
      }

  try(pheatmap(test,
           clustering_distance_rows = "binary",
           clustering_distance_cols = "binary",
           cutree_rows = 2,
           cutree_cols = 2,fontsize = 1.8, fontsize_row=2,filename= paste(prefixU,"heatmap_orthologs_binary.png", sep = "_"), main="Orthologues clustering using binary distances"
            ), silent = TRUE)


  ### Beta-diversity (Bray-Curtis):

  jaccard_dist = as.matrix(vegdist(t(test), method = "jaccard")) #Orthologues "bray"
  jaccard_distst = as.matrix(vegdist(test, method = "jaccard")) #Strains

  try(pheatmap(
    jaccard_dist,
    clustering_distance_rows = as.dist(jaccard_dist),
    clustering_distance_cols = as.dist(jaccard_dist),
  cutree_rows = 2,
  cutree_cols = 2,fontsize = 1.8, filename= paste(prefixU,"heatmap_orthologs_jaccard.png", sep = "_"), main="Orthologs clustering using jaccard dissimilarity"
  ), silent = TRUE)
  try(pheatmap(
    jaccard_distst,
    clustering_distance_rows = as.dist(jaccard_distst),
    clustering_distance_cols = as.dist(jaccard_distst),
    cutree_rows = 2,
    cutree_cols = 2,fontsize = 1.8,filename= paste(prefixU,"heatmap_strains_jaccard.png", sep = "_"), main="Strains clustering using jaccard dissimilarity"
  ), silent = TRUE)


  try(pheatmap(test,
           clustering_distance_rows = as.dist(jaccard_distst),
           clustering_distance_cols = as.dist(jaccard_dist),
           cutree_rows = 2,
           cutree_cols = 2,fontsize = 1.8, fontsize_row=1.8,filename= paste(prefixU,"heatmap_orthologs_jaccard_strain_otholg.png", sep = "_"), main="Orthologues and Strain clustering using Jaccard distances"
  ), silent = TRUE)


 try(pheatmap(test,
           clustering_distance_rows = as.dist(1-cor(t(test), method = "spearman")),
           clustering_distance_cols = as.dist(1-cor(test, method = "spearman")),
           cutree_rows = 2,
           cutree_cols = 2,fontsize = 1.8, fontsize_row=1.8,filename= paste(prefixU,"heatmap_orthologs_spearman_strain_otholg.png", sep = "_"), main="Orthologues and Strain clustering using Spearman distances"
  ), silent = TRUE)


  # run the mds algorithm
  try( mds <- metaMDS(jaccard_distst), silent = TRUE)
  if (exists("mds")) {
    color_type = rep("white", ncol(jaccard_distst))
    color_type[colnames(jaccard_distst) %in% PATHOGENIC_STRAINS]="red"

    pch_type = rep(21,  ncol(jaccard_distst)) # Circle symbols
    pch_type[colnames(jaccard_distst) %in% PATHOGENIC_STRAINS] = 23 # Diamond symbols

    df_mds=data.frame( X1=mds$points[,1], X2=mds$points[,2])
    df_mds$Isolate="Env"
    df_mds$Isolate[row.names(mds$points) %in% PATHOGENIC_STRAINS]="Clin"
    #plot the results
    set.seed(42)
    gmds=ggplot(data=df_mds, aes(x = X1, y = X2, color = Isolate))+
      geom_point()+ labs(x = "MDS1", y = "MDS2", subtitle = "MDS using Jaccard distance")+
      geom_text_repel(label=row.names(df_mds), #mds$points
        size = 2.5,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      )+
      theme_bw()
      ggsave(paste(prefixU,"MDS_sel_features.pdf", sep="_"), gmds, width = 42, height = 60, units = "cm", device = "pdf")
   }

    res.pca3 <- PCA(test,  graph = FALSE)
    A<-fviz_pca_ind(res.pca3, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE, ggtheme = theme_minimal())

    ggsave(paste(prefixU,"PCA.pdf", sep="_"), A, width = 32, height = 32, units = "cm", device = "pdf")
    ####
    #UMAP

    df.umapSEL = umap(test, n_components = 2, random_state = 15)  #df[,selected_features]
    layoutSEL <- df.umapSEL[["layout"]]
    finalSEL <- data.frame(layoutSEL)
    finalSEL$Isolate="Env"
    finalSEL$Isolate[row.names(layoutSEL) %in% PATHOGENIC_STRAINS]="Clin"

    UMFSEL=ggplot(data=finalSEL, aes(x = X1, y = X2, color = Isolate))+
      geom_point()+
      labs(x = "UMAP1", y = "UMAP2", subtitle = "UMAP plot")+
      geom_text_repel(label=row.names(finalSEL),max.overlaps = 30,
                    size = 2,
                    box.padding = 0.25,
                    point.padding = 0.25 )+theme_bw()

    ggsave(paste(prefixU,"UMAP_sel_features.pdf", sep="_"), UMFSEL, width = 32, height = 32, units = "cm", device = "pdf")

    ####
    # Control variable colors using their contributions
    B<-fviz_pca_var(res.pca3, col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
    ggsave(paste(prefixU,"PCA_variables.pdf", sep="_"), B, width = 32, height = 32, units = "cm", device = "pdf")

    cl.res <- clara(test,argv$c)
    C<-fviz_cluster(cl.res, stand=FALSE, repel = TRUE, show.clust.cent=FALSE, ggtheme = theme_bw(), main="")
    ggsave(paste(prefixU,"Clara_clusters.pdf", sep="_"), C, width = 32, height = 32, units = "cm", device = "pdf")


    sink(paste(prefixU,"Clusters.txt", sep="_"))

    for (u in 1:argv$c) { clustersout(u, cl.res, PATHOGENIC_STRAINS ) }
    sink()
}

visualizing_ortholgues<-function(feature, mdf, parte) {
#Can be imporved... keep only the best figure option
    cat("INFO: Visualization steps - ", parte, "Orthologues\n")
    SELDF=mdf[,feature, drop=FALSE]
    SELDF$Isolate=0
    SELDF$Isolate[row.names(SELDF) %in% PATHOGENIC_STRAINS]=1
    SELDF=SELDF[order(SELDF$Isolate),]
    subSELDF=SELDF[row.names(SELDF) %in% tree$tip.label,]

    for (t in names(subSELDF)) {
        if (t == "Isolate") {
            subSELDF[[t]][subSELDF[[t]] > 0] = paste(argv$y,"isolate", sep=" ")
            subSELDF[[t]][subSELDF[[t]] == 0] = paste(argv$z,"isolate", sep=" ")
        } else {
          subSELDF[[t]][subSELDF[[t]] > 0] = "Orthologue present"
          subSELDF[[t]][subSELDF[[t]] == 0] =  "Orthologue absent"
        }
    }
    if (argv$e != "None" ) {
          Special_group=genome_table[,1][genome_table[,4] == argv$e]
          subSELDF[[argv$e]]="."
          subSELDF[[argv$e]][row.names(subSELDF) %in% Special_group]=argv$e
    }



    N2<-ggtree(tree, layout="fan", open.angle = 12, ladderize=TRUE, branch.length="branch.length", size=0.5) +
      geom_tiplab(size=1.3, align=T, linesize=.05, offset =0.02*length(feature)/20 )+ #0.125
      geom_text2(aes(subset = !isTip, label=label), size=1.0,angle = 0,vjust = "inward", hjust = "inward",check_overlap = TRUE)+
      geom_nodepoint(color="black", alpha=1/3, size=0.5 )+ geom_rootedge(rootedge = 0.005)


    TOR2<-gheatmap(N2, subSELDF, offset=0, width=0.02*length(feature),legend_title = "", colnames=T, colnames_position="top",
                   colnames_angle=80,hjust=0, font.size=0.8) +scale_fill_manual(
                     values = c("white", "blue", "red", "purple", "gray", "black", "darkgreen","brown","orange","yellow"), name="")

    ggsave(paste(prefixV,"Tree_with",parte,"orthologues_fan_branch.length.pdf", sep="_"), TOR2, width = 32, height = 32, units = "cm", device = "pdf")


        #Importance

    try(fitSEL <- phyloglm(Isolate~.,phy=tree,data=SELDF, boot=argv$b, method = meth, btol=argv$a), silent = TRUE)

        if (exists("fitSEL")) {
          sink(paste(prefix,parte,"feature_importance.txt", sep="_"))
          print(summary(fitSEL))
          sink()
    }
}

clinical_prediction<-function(feature, mdf, parte) {

    cat("INFO: Pathogenicity classification - ", parte, "Orthologues\n")
    SELDF=mdf[,feature, drop=FALSE]
    SELDF$Isolate=0
    SELDF$Isolate[row.names(SELDF) %in% PATHOGENIC_STRAINS]=1
    SELDF=SELDF[order(SELDF$Isolate),]
    subSELDF=SELDF[row.names(SELDF) %in% tree$tip.label,]

    for (t in names(subSELDF)) {
      if (t == "Isolate") {
        subSELDF[[t]][subSELDF[[t]] > 0] = paste(argv$y,"isolate", sep=" ")
        subSELDF[[t]][subSELDF[[t]] == 0] =  paste(argv$z,"isolate", sep=" ")
      } else {
        subSELDF[[t]][subSELDF[[t]] > 0] = "Orthologue present"
        subSELDF[[t]][subSELDF[[t]] == 0] =  "Orthologue absent"
      }
    }
    if (argv$e != "None" ) {
      Special_group=genome_table[,1][genome_table[,4] == argv$e]
      subSELDF[[argv$e]]="."
      subSELDF[[argv$e]][row.names(subSELDF) %in% Special_group]=argv$e
    }

    mRF=tuneRF(SELDF[,!names(SELDF) %in% c("Isolate")], SELDF$Isolate, 5, ntree=2000)
    mymtry=mRF[mRF[,2]==min(mRF[,2])][1]
    rfc %<-% {randomForest(Isolate~.,data=SELDF,importance=TRUE, mtry=mymtry, keep.forest=TRUE, ntree=2000)} #, cutoff=c(0.7,0.3))
   # rfc <-randomForest(Isolate~.,data=SELDF,importance=TRUE, mtry=mymtry, keep.forest=TRUE) #, cutoff=c(0.7,0.3))


    PDF=data.frame(Strain=names(rfc$predicted), Pathogenticity_probability=rfc$predicted, Clinical=SELDF$Isolate)
    PDF=PDF[order(PDF$Pathogenticity_probability, decreasing = T ),]

    pdf(paste(prefixP, parte,"Random_Forest_strain_patthogenicity_probability.pdf", sep="_"))
    roc_object <- roc(SELDF$Isolate, rfc$predicted,plot = TRUE, print.auc = TRUE)
    dev.off()

    sink(paste(prefixP, parte, "Random_Forest_importance.txt", sep="_"))
    print(rfc$importance[order(rfc$importance[,1], decreasing = T ),])
    sink()
    fwrite(PDF, file=paste(prefixP,parte, "strain_pathogenicity_probability.tsv", sep="_"), quote=FALSE, sep='\t', row.names = FALSE)

        return(PDF)
  }

############# set input/output files
dir.create(argv$o, showWarnings = FALSE)
dirtree=paste(argv$o, "tree", sep="/")
dir.create(dirtree, showWarnings = FALSE)

prefix=paste(argv$o,argv$l, sep="/")
prefixT=paste(dirtree,argv$l, sep="/")

cat("INFO: Reading data\n")
cat("INFO: Reading gene count table\n")

gene_count <- read_tsv(argv$g, col_names =TRUE, show_col_types = FALSE)

df_A<-as.data.frame(t(gene_count))
names(df_A)<-df_A[1,]
df_A <- df_A[-1,]
df_A <- df_A[-nrow(df_A),] #Removing total count row

df_A %<-% {mutate_all(df_A, function(x) as.numeric(as.character(x)))}

###### read tree
cat("INFO: Reading phylogenetic tree\n")
if (argv$N != "") {
  cat("INFO: FASTANI output will be used to create the phylogenetic tree\n")
  table <- read_tsv(argv$N, col_names = F)
  table<-table[!(table$V1==table$V2),]
  table$V3 <- 1 - (table$V3/100)
  table$V4 <- NULL
  table$V5 <- NULL

  matrix <- spread(data = table, key = V1, value = V3 )
  row.names(matrix) <- matrix$V2

  my_tree <- as.phylo(upgma(dist(matrix)))

  my_tree$tip.label=gsub("Genomes/","",my_tree$tip.label)
  my_tree$tip.label=gsub(".fa","",my_tree$tip.label)
  my_tree$tip.label=gsub("'","",my_tree$tip.label)

  write.tree(my_tree, paste(prefixT,"Cleaned_tree.txt", sep="_"))
  tree=my_tree
} else {

  tree = read.tree(argv$t)
  ix1 = c(grep("^RS_GCF_", tree$tip.label), grep("^GB_GCA_", tree$tip.label), grep("^GCA_", tree$tip.label))
  tree= drop.tip(tree, ix1)

  tree$tip.label=gsub("'","",tree$tip.label)

   write.tree(tree, paste(prefixT,"Cleaned_tree.txt", sep="_"))

}
# read strain/genome information
cat("INFO: Reading genome/sample information \n")
genome_table <- read.csv(argv$i, sep = ",", na.strings = "NA", header = FALSE, col.names = c("Strain", "File","Type","Special"))

genome_table %<-% {mutate_all(genome_table, function(x) gsub("\\s+"," ",x))}

#Clinical Isolates names
PATHOGENIC_STRAINS=genome_table[,1][genome_table[,3] == argv$w ]

#remove orthologues present/absent in more than user defined percentage (default 95%) of the total strains in the dataset
if (argv$r == "") {
  NC=length(PATHOGENIC_STRAINS)
  NE=length(genome_table[,1][genome_table[,3] == argv$v])
  NT=NC+NE
  CUT_OFF_ORT=0.5*(1+ NC/NT)
} else { CUT_OFF_ORT=argv$r }


toremove=rep(0, ncol(df_A))
for (c in 1:ncol(df_A)) {
  if(length(df_A[[c]][df_A[[c]]!=0])/nrow(df_A) > CUT_OFF_ORT) {toremove[c]=c}
  else if(length(df_A[[c]][df_A[[c]]==0])/nrow(df_A) > CUT_OFF_ORT) {toremove[c]=c}
}
toremove=toremove[toremove>0]
df_A=df_A[,-toremove]

#Visualizing tree

y=c()
for (n in tree$tip.label) {
      if (n %in% PATHOGENIC_STRAINS) y=c(y,1) else y=c(y,0)
      }

if(length(y[y>0]) != length(PATHOGENIC_STRAINS)) stop("Strain names and tree labels don't match, please check spelling in GENOME_LIST file")

tipo<-y
tipo[tipo==0]<-argv$z
tipo[tipo==1]<-argv$y
dfP <- data.frame(Isolates=as.character(tipo))
rownames(dfP) <- tree$tip.label

if (argv$e != "None" ) {
  Special_group=genome_table[,1][genome_table[,4] == argv$e]
  dfP[[argv$e]]=" ISOLATES"
  dfP[[argv$e]][row.names(dfP) %in% Special_group]=argv$e
}

SIZE=3

OP<-ggtree(tree, layout="rectangular", open.angle = 0, ladderize=TRUE, branch.length="branch.length") +   #"fan" "none"
  geom_tiplab(size=1.8, align=TRUE, linesize=.05, offset=0.02)+
  geom_text2(aes(subset = !isTip, label=label), size=SIZE,angle = 0,vjust = "inward", hjust = "inward",check_overlap = TRUE)+
  geom_nodepoint(color="black", alpha=1/3, size=1)+ geom_rootedge(rootedge = 0.005)


ggsave(paste(prefixT,"tree_ladderized_branchlength.pdf", sep="_"), OP, width = 84, height = 120, units = "cm", device = "pdf")


if (argv$e == "None" ) {
  p2<-ggtree(tree, layout="fan", open.angle = 0, ladderize=TRUE, branch.length="branch.length") +   #"fan" "none"
    geom_tiplab(size=1.8, align=T, linesize=.05, offset=0.0015)+
    geom_text2(aes(subset = !isTip, label=label), size=SIZE,angle = 0,vjust = "inward", hjust = "inward",check_overlap = TRUE)+
    geom_nodepoint(color="black", alpha=1/3, size=1)+ geom_rootedge(rootedge = 0.005)
  P3<-gheatmap(p2, dfP, offset=0, width=.03,legend_title = "Isolates", colnames = FALSE)

}else{
  ofs=0.025*(ncol(dfP))

if (special_nodes != "none") {
  tree2 <- groupClade(tree, .node=special_nodes)
  p2<-ggtree(tree2, aes(color=group),layout="fan", open.angle = 10, ladderize=TRUE, branch.length="branch.length") +
    geom_tiplab(size=1.6, align=T, linesize=.05, offset=0.0028)+
    geom_text2(aes(subset = !isTip, label=label), size=SIZE,angle = 0,vjust = "inward", hjust = "inward",check_overlap = TRUE)+
    geom_nodepoint(color="black", alpha=1/3, size=1 )+ geom_rootedge(rootedge = 0.005)+
    scale_color_manual(values=labels_nodes_colors, labels=labels_nodes, name=type_name_special_nodes)
} else {
  p2<-ggtree(tree, layout="fan", open.angle = 0, ladderize=TRUE, branch.length="branch.length") +
    geom_tiplab(size=1.6, align=T, linesize=.05, offset=0.0028)+
    geom_text2(aes(subset = !isTip, label=label), size=SIZE,angle = 0,vjust = "inward", hjust = "inward",check_overlap = TRUE)+
    geom_nodepoint(color="black", alpha=1/3, size=1 )+ geom_rootedge(rootedge = 0.005)
}
  P3<-gheatmap(p2, dfP, offset=0, width=ofs, colnames=F) +scale_fill_manual(
                 values = c("white", "blue", "red", "purple","darkgreen", "grey", "black","brown","orange","yellow"), name=""
               )
}
ggsave(paste(prefixT,"tree_fan_ladderized_branchlength.pdf", sep="_"), P3, width = 42, height = 60, units = "cm", device = "pdf")

sink(paste(prefixT,"labels_tree.tsv", sep="_"))
print(data.frame(Pathogenicity=y,Strain=tree$tip.label), quote = FALSE, row.names = FALSE)
sink()

############# Phylogenetic Generalized Linear Model

df=df_A[match(tree$tip.label,row.names(df_A)),]

if (all(row.names(df)==tree$tip.label)) cat("INFO: Comparative genomics starts\n") else stop("Strain names in Gene count table and tree labels don't match, please check spelling in GENOME_LIST file")


dirR=paste(argv$o, "R_objects", sep="/")
if (!dir.exists(dirR)) dir.create(dirR, showWarnings = FALSE)

if (!file.exists(paste(dirR,"RES_R", sep="/")) & !file.exists(paste(dirR,"COEFs_R", sep="/")) & !file.exists(paste(dirR,"ALPHAs_R", sep="/")) ) {

    meth="logistic_IG10" #options ("poisson_GEE","logistic_MPLE","logistic_IG10")
    set.seed(123)
    RES=vector("list", ncol(df))
    names(RES)=names(df)
    COEFs=vector("list", ncol(df))
    names(COEFs)=names(df)
    ALPHAs=vector("list", ncol(df))
    names(ALPHAs)=names(df)
    pb<-txtProgressBar(min=0, max=ncol(df), style = 3, width=50, char="-")
    for (i in 1:ncol(df)) {
      #i=1
      dat = data.frame(row.names=tree$tip.label,trait = y, predictor = df[[i]])
      dat$predictor[dat$predictor>1]=1

      try(fit %<-%  {phyloglm(trait~predictor,phy=tree,data=dat, boot=argv$b, method = meth, btol=argv$a)}, silent = TRUE)

      if (exists("fit")) {
                    sta=summary(fit)
                    if (argv$b == 0) { pvalue=sta$coefficients[2,4]} else { pvalue=sta$coefficients[2,6] }

                    if (pvalue <= 1) {
                      RES[[names(df)[i]]]=pvalue
                      COEFs[[names(df)[i]]]=sta$coefficients[2,1]
                      ALPHAs[[names(df)[i]]]=sta$alpha
                    }
      }
      setTxtProgressBar(pb,i)
    }
    close(pb)

    RES=RES[lapply(RES,length)>0]
    COEFs=COEFs[lapply(COEFs,length)>0]
    ALPHAs=ALPHAs[lapply(ALPHAs,length)>0]

    saveRDS(RES, file=paste(dirR,"RES_R", sep="/"))
    saveRDS(COEFs, file=paste(dirR,"COEFs_R", sep="/"))
    saveRDS(ALPHAs, file=paste(dirR,"ALPHAs_R", sep="/"))
} else {
    RES=readRDS(file=paste(dirR,"RES_R", sep="/"))
    COEFs=readRDS(file=paste(dirR,"COEFs_R", sep="/"))
    ALPHAs=readRDS(file=paste(dirR,"ALPHAs_R", sep="/"))
    meth="logistic_IG10"
}

res2=p.adjust(RES, method = "fdr")
res3=res2[res2<argv$q]
#### Order by pvalue
res3=res3[order(res3)]
####

gene_count_clinical = subset(df, row.names(df) %in% PATHOGENIC_STRAINS)

gene_count_others = subset(df, !row.names(df) %in% PATHOGENIC_STRAINS)

Freqs_orth=list()
for (o in names(res3)) {

  w1=gene_count_clinical[[o]]
  w1[w1>0]<-1
  w2=gene_count_others[[o]]
  w2[w2>0]<-1
  
  name_g1_p=paste0(argv$y,"_Present" )
  name_g1_a=paste0(argv$y,"_Abs" )
  name_g2_p=paste0(argv$z,"_Present" )
  name_g2_a=paste0(argv$z,"_Abs" )
  
  Freqs_orth[[o]][[name_g1_p]]=length(w1[w1>0])
  Freqs_orth[[o]][[name_g1_a]]=length(w1[w1==0])
  Freqs_orth[[o]][[name_g2_p]]=length(w2[w2>0])
  Freqs_orth[[o]][[name_g2_a]]=length(w2[w2==0])

}

cat("\n ---------------- \n")


dirA=paste(argv$o, "Annotations", sep="/")
dir.create(dirA, showWarnings = FALSE)
prefixA=paste(dirA,argv$l, sep="/")

cat("Ortolog",paste0("Present in ",argv$y),paste0("Absent in ", argv$y),paste0("Present in ", argv$z),paste0("Absent in ", argv$z),"PadjValue","Coeffcient", "alpha\n", sep="\t",
    file=paste(prefixA,"Candidates_enriched_orthologues.tsv", sep="_"))

cat("Ortolog",paste0("Present in ",argv$y),paste0("Absent in ", argv$y),paste0("Present in ", argv$z),paste0("Absent in ", argv$z),"PadjValue","Coeffcient", "alpha\n", sep="\t",
        file=paste(prefixA,"Candidates_depleted_orthologues.tsv", sep="_"))

abs_markets=rep("0", length(names(res3)))
present_markets=rep("0", length(names(res3)))

k=1
for (o in names(res3)) {

  if (COEFs[[o]] < 0 ) {
    abs_markets[k]=o
    cat(o,Freqs_orth[[o]][[name_g1_p]],
        Freqs_orth[[o]][[name_g1_a]],Freqs_orth[[o]][[name_g2_p]],
        Freqs_orth[[o]][[name_g2_a]],res3[[o]], COEFs[[o]], ALPHAs[[o]], "\n", sep="\t",
        file=paste(prefixA,"Candidates_depleted_orthologues.tsv", sep="_"), append = TRUE)
  } else {
    present_markets[k]=o
    cat(o,Freqs_orth[[o]][[name_g1_p]],
        Freqs_orth[[o]][[name_g1_a]],Freqs_orth[[o]][[name_g2_p]],
        Freqs_orth[[o]][[name_g2_a]],res3[[o]], COEFs[[o]], ALPHAs[[o]], "\n", sep="\t",
        file=paste(prefixA,"Candidates_enriched_orthologues.tsv", sep="_"), append = TRUE)
  }
  k=k+1
}
abs_markets=abs_markets[abs_markets != "0" ]
present_markets=present_markets[present_markets != "0" ]

cat("Ortolog",paste0("Present in ",argv$y),paste0("Absent in ", argv$y),paste0("Present in ", argv$z),paste0("Absent in ", argv$z),"PadjValue","Coeffcient", "alpha\n", sep="\t",
    file=paste(prefixA,"core_enriched_orthologues.tsv", sep="_"))
cat("Ortolog",paste0("Present in ",argv$y),paste0("Absent in ", argv$y),paste0("Present in ", argv$z),paste0("Absent in ", argv$z),"PadjValue","Coeffcient", "alpha\n", sep="\t",
    file=paste(prefixA,"core_depleted_orthologues.tsv", sep="_"))

enriched_core=rep("0", length(present_markets))
depleted_core=rep("0", length(abs_markets))

j=1
l=1
for (o in names(res3)) {

  if (o %in% present_markets & Freqs_orth[[o]][[name_g1_p]]/(Freqs_orth[[o]][[name_g1_a]]+Freqs_orth[[o]][[name_g1_p]]) > core_fraction) {
    cat(o,Freqs_orth[[o]][[name_g1_p]],
        Freqs_orth[[o]][[name_g1_a]],Freqs_orth[[o]][[name_g2_p]],
        Freqs_orth[[o]][[name_g2_a]],res3[[o]], COEFs[[o]], ALPHAs[[o]], "\n", sep="\t",
        file=paste(prefixA,"core_enriched_orthologues.tsv", sep="_"), append = TRUE)
    enriched_core[j]=o
    j=j+1
  }
  if (o %in% abs_markets & Freqs_orth[[o]][[name_g1_a]]/(Freqs_orth[[o]][[name_g1_a]]+Freqs_orth[[o]][[name_g1_p]]) > core_fraction) {
    cat(o,Freqs_orth[[o]][[name_g1_p]],
        Freqs_orth[[o]][[name_g1_a]],Freqs_orth[[o]][[name_g2_p]],
        Freqs_orth[[o]][[name_g2_a]],res3[[o]], COEFs[[o]], ALPHAs[[o]], "\n", sep="\t",
        file=paste(prefixA,"core_depleted_orthologues.tsv", sep="_"), append = TRUE)
    depleted_core[l]=o
    l=l+1
  }


}
enriched_core=enriched_core[enriched_core != "0" ]
depleted_core=depleted_core[depleted_core != "0" ]

selected_features=c()
selected_core_features=c()
if (length(present_markets) > 0 ) {
    cat("\n List of Enriched orthologs in", argv$y, "Isolates",present_markets, "\n")
    selected_features=c(selected_features,present_markets)

    cat("#List of Enriched Orthologues\n",file=paste(prefix,"Enriched.txt", sep="_"))
    for (E in present_markets) {cat(E,"\n", file=paste(prefix,"Enriched.txt", sep="_"), append=TRUE)}

    } else { cat("No enriched Ortholog found") }

if (length(enriched_core) > 0 ) {
    cat("\n List of Enriched core orthologs in", argv$y, "Isolates",enriched_core, "\n")
    selected_core_features=c(selected_core_features,enriched_core)

    cat("#List of Enriched core Orthologues\n",file=paste(prefix,"Enriched_core.txt", sep="_"))
    for (Ec in enriched_core) {cat(Ec,"\n", file=paste(prefix,"Enriched_core.txt", sep="_"), append=TRUE)}

    } else { cat("No enriched core Ortholog found") }

if (length(abs_markets) > 0 ) {
  cat("\n List of depleted orthologs in", argv$y,"Isolates",abs_markets,"\n")
  selected_features=c(selected_features,abs_markets)

  cat("#List of depleted Orthologues\n",file=paste(prefix,"Depleted.txt", sep="_"))
  for (Dp in abs_markets) {cat(Dp,"\n", file=paste(prefix,"Depleted.txt", sep="_"), append=TRUE)}


  } else {
  cat("No depleted Ortholog found")
  }

if (length(depleted_core) > 0 ) {
    cat("\n List of depleted core orthologs in", argv$y,"Isolates",depleted_core,"\n")
    selected_core_features=c(selected_core_features,depleted_core)

    cat("#List of Depleted core Orthologues\n",file=paste(prefix,"Depleted_core.txt", sep="_"))
    for (Dpc in depleted_core) {cat(Dpc,"\n", file=paste(prefix,"Depleted_core.txt", sep="_"), append=TRUE)}

    } else {
    cat("No depleted core Ortholog found")
    }


if (length(selected_features) == 0) stop("No enriched/depleted Ortholog found")

if (Pathogenicity_probability) {

  rest_of_packages <- c("randomForest","pROC")
  for (s in rest_of_packages) { suppressPackageStartupMessages(library(s, character.only = TRUE))}

  dirpred=paste(argv$o, "Prediction", sep="/")
  dir.create(dirpred, showWarnings = FALSE)

  prefixP=paste(dirpred,argv$l, sep="/")

  Prob_patho=clinical_prediction(selected_features, df, "Selected")
  #Prob_patho=clinical_prediction(c(enriched_core, depleted_core ), df, "enriched_depleted_core")

  Prob_patho2=Prob_patho[,names(Prob_patho) %in% c("Pathogenticity_probability"), drop=FALSE]
  Prob_patho2$Pathogenticity_probability=round(Prob_patho2$Pathogenticity_probability*100,2)

  ######## Should I keep this? starts

  cat("INFO: Visualization steps - Pathogenicity probability Orthologues\n")
  SELDF=Prob_patho2
  SELDF$Isolate=0
  SELDF$Isolate[row.names(SELDF) %in% PATHOGENIC_STRAINS]=1
  SELDF=SELDF[order(SELDF$Isolate),]

  subSELDFPATHO=SELDF[row.names(SELDF) %in% tree$tip.label,]
  subSELDFPATHO[["Isolate"]][subSELDFPATHO[["Isolate"]] > 0] = paste(argv$y,"isolate", sep=" ")
  subSELDFPATHO[["Isolate"]][subSELDFPATHO[["Isolate"]] == 0] =  paste(argv$z,"isolate", sep=" ")
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability <= 1] =0
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability > 1 & subSELDFPATHO$Pathogenticity_probability <= 5] =1
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability > 5 & subSELDFPATHO$Pathogenticity_probability <= 10] =2
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability > 10 & subSELDFPATHO$Pathogenticity_probability <= 30] =3
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability > 30 & subSELDFPATHO$Pathogenticity_probability <= 50] =4
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability > 50 & subSELDFPATHO$Pathogenticity_probability <= 70] =5
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability > 70 & subSELDFPATHO$Pathogenticity_probability <= 80] =6
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability > 80 & subSELDFPATHO$Pathogenticity_probability <= 90] =7
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability > 90 ] =8

  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability == 8] = "Level 8 | Probability(%) > 90"
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability == 7] = "Level 7 | 80 < Probability(%) >= 90"
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability == 6] = "Level 6 | 70 < Probability(%) >= 80"
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability == 5] = "Level 5 | 50 < Probability(%) >= 70"
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability == 4] = "Level 4 | 30 < Probability(%) >= 50"
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability == 3] = "Level 3 | 10 < Probability(%) >= 30"
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability == 2] = "Level 2 | 5 < Probability(%) >= 10"
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability == 1] = "Level 1 | 1 < Probability(%) >= 5"
  subSELDFPATHO$Pathogenticity_probability[subSELDFPATHO$Pathogenticity_probability == 0] = "Level 0 | Probability(%) <= 1"


  if (argv$e != "None" ) {
    Special_group=genome_table[,1][genome_table[,4] == argv$e]
    subSELDFPATHO[[argv$e]]=" "
    subSELDFPATHO[[argv$e]][row.names(subSELDFPATHO) %in% Special_group]=paste(" ",argv$e,sep=" ")
  }


  N<-ggtree(tree, layout="fan", open.angle = 7, ladderize=TRUE, branch.length="none") +
    geom_tiplab(size=1.5, align=T, linesize=.05, offset =4 )
  TOR<-gheatmap(N, subSELDFPATHO, offset=0, width=0.08,legend_title = "", colnames=T, colnames_position="top",
                colnames_angle=25,hjust=0, font.size=1.5) +scale_fill_manual(
                  values = c("white","blue","black","white", brewer.pal(n = 9, name = "Reds")), name="")

  ggsave(paste(prefixP,"Tree_with_probabilities_based_on_selected_features.pdf", sep="_"), TOR, width = 32, height = 32, units = "cm", device = "pdf")
####### Should I keep this? ends
}


dirV=paste(argv$o, "Visualization", sep="/")
dir.create(dirV, showWarnings = FALSE)
prefixV=paste(dirV,argv$l, sep="/")

if (length(selected_features) >1) {
  cat("INFO: Unsupervise learning steps \n")
  unsupervised_learning("Unsupervise", selected_features, df)
  visualizing_ortholgues(selected_features, df, "Selected")
}

if (length(selected_core_features) >1) {
  cat("INFO: Unsupervise learning steps using core enriched/depleted orthologs \n")
  unsupervised_learning("Unsupervise_core", selected_core_features, df)
  visualizing_ortholgues(selected_core_features, df, "core")
  visualizing_ortholgues(enriched_core, df, "enriched_core")
  visualizing_ortholgues(depleted_core, df, "depleted_core")
}

#annotations of selected orthologs
cat("INFO: Annotation steps\n")

Orthogroups_g <- read_tsv(argv$s, col_names =FALSE, show_col_types = FALSE)

annots_g <- read_delim(argv$m, delim = ";", col_names=FALSE, show_col_types = FALSE)


if (length(present_markets) >0) {
pres_anot_table %<-% {print_out_annot(present_markets, Orthogroups_g, annots_g)}
cat("INFO: Annotation Enriched orthologues\n")
fwrite(pres_anot_table, file=paste(prefixA, "Enriched_orthologues_annotation.tsv", sep="_"), quote=FALSE, sep='\t', row.names = FALSE)
}
if (length(abs_markets) >0) {
abs_anot_table %<-% {print_out_annot(abs_markets, Orthogroups_g, annots_g)}
cat("INFO: Annotation Depleted orthologues\n")
fwrite(abs_anot_table, file=paste(prefixA, "Depleted_orthologues_annotation.tsv", sep="_"), quote=FALSE, sep='\t', row.names = FALSE)

}

if (length(enriched_core) >0) {
  enriched_core_anot_table %<-% {print_out_annot(enriched_core, Orthogroups_g, annots_g)}
  cat("INFO: Annotation Enriched core orthologues\n")
  fwrite(enriched_core_anot_table, file=paste(prefixA, "Enriched_core_orthologues_annotation.tsv", sep="_"), quote=FALSE, sep='\t', row.names = FALSE)
}

if (length(depleted_core) >0) {
  depleted_core_anot_table %<-% {print_out_annot(depleted_core, Orthogroups_g, annots_g)}
  cat("INFO: Annotation Depleted core orthologues\n")
  fwrite(depleted_core_anot_table, file=paste(prefixA, "depleted_core_orthologues_annotation.tsv", sep="_"), quote=FALSE, sep='\t', row.names = FALSE)

  }


