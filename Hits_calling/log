#!/bin/tcsh

#SBATCH -J matlab
#SBATCH -p super
#SBATCH -N 1
#SBATCH -t 0-24:00:00
#SBATCH --mail-user=Gary.Hon@utsouthwestern.edu
#SBATCH --mail-type=end
#SBATCH -o preprocess.out.%j
#SBATCH -e preprocess.err.%j

echo "hello world"

matlab -nodesktop -nosplash < run1.m
