#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=15:00:00
#PBS -o log/
#PBS -e log/
#PBS -N sim3
#PBS -m ea
#PBS -l mem=12G
#PBS -l vmem=12G
#PBS -t 1-50%4

cd $PBS_O_WORKDIR

SRC=$PBS_O_WORKDIR


Rscript ${SRC}/simulation.R ${index} $PBS_ARRAYID
