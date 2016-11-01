#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=5:00:00
#PBS -o log/
#PBS -e log/
#PBS -N sim5
#PBS -m ea
#PBS -t 1-2
#PBS -l pmem=11700m
#PBS -A nzt-671-ab

module load foss/2015b
module load iomkl/2015b
module load R/3.3.1

cd $PBS_O_WORKDIR

SRC=$PBS_O_WORKDIR

Rscript ${SRC}/simulation5-guillimin.R ${index} $PBS_ARRAYID

