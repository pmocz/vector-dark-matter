#!/bin/bash
#SBATCH -p hernquist
#SBATCH -J vdm42
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -o OUTPUT.lsf
#SBATCH -e ERROR.lsf
#SBATCH -t 8640 # 4 day in min
#SBATCH --mail-user=philip.mocz@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem=12000  # 4000 for 256

module purge
module load matlab

srun -n 1 -c 16 matlab -nojvm -nodisplay -nosplash -nodesktop -r "N=200; myseed=42; vectorDM.m"



#EOF
