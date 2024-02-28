#!/bin/bash -l        
#SBATCH --time=24:00:00
#SBATCH --ntasks=10
#SBATCH --mem=60g
#SBATCH --tmp=60g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=mullinsj@umn.edu 
#SBATCH --output=output/bootstrap.out

cd ~/BFM
~/julia-1.9.3/bin/julia -t 10 scripts/run_bootstrap.jl

