#!/bin/bash -l        
#SBATCH --time=0:10:00
#SBATCH --ntasks=50
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=mullinsj@umn.edu 
#SBATCH --output=output/test_speeds.out

cd ~/BFM
~/julia-1.9.3/bin/julia -t 5 tests/test_speed.jl
~/julia-1.9.3/bin/julia -t 10 tests/test_speed.jl
~/julia-1.9.3/bin/julia -t 15 tests/test_speed.jl
~/julia-1.9.3/bin/julia -t 20 tests/test_speed.jl
~/julia-1.9.3/bin/julia -t 25 tests/test_speed.jl

