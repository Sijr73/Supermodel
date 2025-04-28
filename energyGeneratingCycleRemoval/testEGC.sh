#!/bin/bash

#PBS -l walltime=23:59:00
#PBS -l select=1:ncpus=4:mem=16gb
#PBS -l place=free:group=board
#PBS -r y
#PBS -N AGC
#PBS -A GBAnonlinearSJ
#PBS -m abe
#PBS -M ghaffas@hhu.de

cd $PBS_O_WORKDIR
echo “JobID: “$PBS_JOBID

#module load intel/xe2013
#module load R/3.0.2

module load R/3.6.1

cat $PBS_NODEFILE


Rscript testEGC.R





