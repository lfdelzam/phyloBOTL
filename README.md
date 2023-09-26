Requirements

1. miniconda - Follow instruction on https://docs.conda.io/projects/miniconda/en/latest/


2. EGGNOG - Follow instructions on https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#user-content-v2112

   Example of command-lines to install emapper.py in the file/script install_key_envs.sh

4. Conda R_env:

      conda env create -f conda_env/R_env.yaml -y
      conda activate R_env
      Rscript --vanilla Install_Rpackages.R <numberofcpus>


5. Snakemake

      conda create -n phylobotl_env -c bioconda snakemake=7.25.0 Python=3.11.4
      conda activate phylobotl_env
      conda config --set channel_priority strict


Optional
A. gtdbtk - source: https://ecogenomics.github.io/GTDBTk/installing/index.html#installing 
          - conda: https://ecogenomics.github.io/GTDBTk/installing/bioconda.html

    conda create -n gtdb_env -c conda-forge -c bioconda gtdbtk=2.1.1 -y
    conda activate gtdb_env
    download-db.sh
    # Set the environment variable to the directory containing the GTDB-Tk reference data
    conda env config vars set GTDBTK_DATA_PATH="/path/to/unarchived/gtdbtk/data"

B. kSNP4  - source: https://sourceforge.net/projects/ksnp/files/latest/download
          - manual: https://sourceforge.net/projects/ksnp/files/kSNP4.1%20User%20Guide.pdf/download


Usage

Modify config file:

      nano support/config_phylobotl.json

and save CTRL+x, y

Run:

      conda activate phylobotl_env
      snakemake -s phylobotl.smk --use-conda --conda-frontend conda --cores <numberofcpus>
