#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o log/
#PBS -e log/
#PBS -N sim6
#PBS -m ea
#PBS -t 1-200
#PBS -l pmem=7500m
#PBS -A nzt-671-ab

module load foss/2015b
module load iomkl/2015b
module load R/3.3.1

cd $PBS_O_WORKDIR

SRC=$PBS_O_WORKDIR

Rscript ${SRC}/simulation6-guillimin.R ${index} $PBS_ARRAYID

