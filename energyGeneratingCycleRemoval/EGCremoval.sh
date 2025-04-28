#!/bin/bash

#PBS -l walltime=71:59:00
#PBS -l select=1:ncpus=4:mem=16gb
#PBS -l place=free:group=board
#PBS -J 1-101
#PBS -r y
#PBS -N AGC
#PBS -A GBAnonlinearSJ
#PBS -m abe
#PBS -M your email.gmail.com

cd $PBS_O_WORKDIR
echo “JobID: “$PBS_JOBID

module load intel/xe2013
#module load R/3.0.2

module load R/3.6.1

cat $PBS_NODEFILE


Rscript EGCremoval.R $PBS_ARRAY_INDEX





