#!/bin/bash -l

conda activate phylobotl_env
snakemake --profile .config/snakemake/phylobotl_slurm -s phylobotl.smk --use-conda --conda-frontend conda $1
