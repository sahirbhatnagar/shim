#~/bin/bash
#PBS -m ea
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -q qwork
#PBS -r n
#PBS -N sim6
#PBS -t 1-200
#PBS -o log/
#PBS -e log/
#PBS -V

module load bioinformatics/R/3.2.5

SRC=$HOME/coexpression/august2016simulation/nonlinear

Rscript $SRC/simulation6-hydra.R ${index} $PBS_ARRAYID
