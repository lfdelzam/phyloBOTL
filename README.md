# Requirements and installation

1. miniconda - Follow instruction on https://docs.conda.io/projects/miniconda/en/latest/

   make sure you have the conda version 23.7.3

         conda activate base
         conda install conda=23.7.3 -y
   
3. clone the repository:

         git clone https://github.com/lfdelzam/phyloBOTL

## Installation option 1 - step by step, in case you already have some programs installed or databases downloaded
Make sure you are using the same version as indicated hereafter: 

Conda R_env:

         cd phyloBOTL
         conda env create -f conda_env/R_env.yaml -y
         conda activate R_env
         Rscript --vanilla support/Install_Rpackages.R <numberofcpus>

Snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html):

         conda create -n phylobotl_env -c bioconda snakemake=7.25.0 Python=3.11.4 -y
         conda activate phylobotl_env
         conda config --set channel_priority strict

Genomad and genomad_db (https://portal.nersc.gov/genomad/index.html)

          conda create -n genomad_env -c conda-forge -c bioconda genomad=1.6.1 -y
          conda activate genomad_env
          mkdir -p <directory_DB_path>/genomad_env/data
          genomad download-database <directory_DB_path>/genomad_env/data
          conda deactivate
    
EGGNOG: emapper.py and EGGNOG database (https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#user-content-v2112)

          conda create -n eggnog_mapper_env -c bioconda -c conda-forge eggnog-mapper=2.0.1 -y
          conda activate eggnog_mapper_env
          mkdir -p <directory_DB_path>/eggnog_mapper_env/data
          download_eggnog_data.py --data_dir <directory_DB_path>/eggnog_mapper_env/data -y
          conda deactivate

## Installation option 2 - all in one go, in case you want to install/download all the required programs and databases from scratch. 

      cd phyloBOTL
      bash install_key_envs.sh <directory_DB_path> <cpus>

<directory_DB_path> a directory where you want to download the databases. 
<cpus> number of cpus to download R packages when creating the R_env

## Optional

if you want to use gtdb or kSNP4 to build the phylogenetic tree

A. gtdbtk - source: https://ecogenomics.github.io/GTDBTk/installing/index.html#installing 
          - conda: https://ecogenomics.github.io/GTDBTk/installing/bioconda.html

    conda create -n gtdb_env -c conda-forge -c bioconda gtdbtk=2.1.1 -y
    conda activate gtdb_env
    download-db.sh
    # Set the environment variable to the directory containing the GTDB-Tk reference data
    conda env config vars set GTDBTK_DATA_PATH="/path/to/unarchived/gtdbtk/data"


B. kSNP4  - source: https://sourceforge.net/projects/ksnp/files/latest/download
          - manual: https://sourceforge.net/projects/ksnp/files/kSNP4.1%20User%20Guide.pdf/download


# Usage

Modify config file:

      nano support/config_phylobotl.json
      
      "workdir":"/abs/path/to/phyloBOTL",
      "threads": 20,   - CPUs
      "input": {
                "File_list":"GENOME_LIST", - csv file <sample name>,<abs/path/to/contig.fasta[.gz]>,<group label>,<special_group_name>
                "group_1_label": "C", 
                "group_1_name": "Clinical",
                "group_2_label": "E",
                "group_2_name": "Environmental",
                "Special_group_name": "Baltic Sea"
              },
      "tree_file": "Tree_rooted_boot.treefile", - tree file name. If you already have a tree, set the path to it here and the pipeline won't build the tree 
      "tree_using": "iqtree",
      "iqtree": {
                "iqtree_params":"-m GTR+I+G -B 1000 -bnni",
                "iqtree_rooted": "TRUE"
                },
      "FastTree_params":"-gtr -nt -gamma",
      "GTDB": {
              "gtdb_params": "--bacteria",
              "taxa_filter":"p__Proteobacteria",
              "outgroup_taxon":"p__Firmicutes"
              },
      "KSNP": {
              "kSNP4_param_phylo_method":"ML",
              "kSNP4_param_genome_fraction": "core",
              "Path_to_kSNP4pkg": "/usr/local/kSNP4pkg"
              },
      "Prokka_params": "--rawproduct --quiet",
      "Proteome_dir": "proteome", - directory containg the protein files <sample name>.faa. 
      "orthology":
                {
                "orthofinder_parameters": "-f proteome -a 4 -S diamond -og", - double check that name in -f <Proteom_dir> corresponds to "Proteome_dir" above
                "ortholog_count_table": "Orthogroups.GeneCount.tsv", 
                "ortholog_table": "Orthogroups.tsv",
                "path_to_orthologs_sequences": "Orthologues/Results_dir/Orthogroup_Sequences"
                },if you provide a full path to existing "ortholog_count_table","ortholog_table" and "path_to_orthologs_sequences" files, the pipeline won't run orthofinder  
      "PHYLOGLM": {
                  "phyloglm_Bootnumber":0,- phyloglm parameter, number of independent bootstrap replicates, 0 means no bootstrap.
                  "p_adj_value_cutoff":0.05,
                  "orthologue_ratio_in_genome_dataset":0.95,- In this case if the orthologue is present/absent in 95% of the genomes, it won't be considered in the analysis. 
                  "phyloglm_btol_number":10, - phyloglm parameter, (logistic regression only) bound on the linear predictor to bound the searching space.
                  "phyloglm_outfiles_prefix":"Vv"
                  },
      "synteny_visualization": {
                            "present_in_at_least_n_genomes": 10, - parameter used when creating the Graph, representing the weigth cutoff to be included in graph 
                            "locus_size_window": 40000, - bps, window used when looking for enriched orthologs locilised in the region
                            "max_n_genomes_per_fig": 22 - 22 representative genomes will be included in the .pdf figure
                          },
      "path_to_eggnog_db":"directory_DB_path/eggnog_mapper_env/data", - abs path to eggnog database directory. 
      "genomad":{
                  "path_to_genomad_db": "directory_DB_path/genomad_env/data/genomad_db", - abs path to genomad database directory
                  "Genomad_score_cut_off": 0.8, - Minimum score to classify a contig as plasmid (or Virus)
                  "params_genomad": "--cleanup --conservative --enable-score-calibration", - genomad parameters. The flag --enable-score-calibration will work if there are more than 1000 sequences in the fasta file (in a genome file) 
                },
      "output_dir": "Results_with_IQTREE_rooted"


      

and save CTRL+x, y

Run:

      conda activate phylobotl_env
      snakemake -s phylobotl.smk --use-conda --conda-frontend conda --cores <numberofcpus>
