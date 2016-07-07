#!/bin/bash
#PBS -l nodes=1:ppn=5
#PBS -l walltime=20:00:00
#PBS -o log/
#PBS -e log/
#PBS -N nih-rda
#PBS -m ea
#PBS -l mem=40G
#PBS -l vmem=40G

cd $PBS_O_WORKDIR

SRC=$PBS_O_WORKDIR


Rscript ${SRC}/clust_kmeans_qsub.R
