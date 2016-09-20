#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -o log/
#PBS -e log/
#PBS -N bou-rda
#PBS -m ea
#PBS -l mem=125G
#PBS -l vmem=125G

cd $PBS_O_WORKDIR

SRC=$PBS_O_WORKDIR


Rscript ${SRC}/clust_kmeans_bouchard_qsub.R
