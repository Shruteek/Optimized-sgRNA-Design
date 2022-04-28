#!/bin/bash
#SBATCH --job-name=align_ctx_with_bowtie
#SBATCH -J download
#SBATCH --time=8:00:00
#SBATCH --account=fc_pickeringlab
#SBATCH -p savio2_htc
#SBATCH --ntasks=4
#SBATCH --mem=84g
#SBATCH --mail-type=END
#SBATCH --mail-user=shruteek@berkeley.edu
#SBATCH --output=datafiles/Outputs/%j_output.txt
#SBATCH --error=datafiles/Outputs/%j_error.txt
## Command(s) to run:
module load python
## export PATH=$PATH:/global/scratch/projects/fc_pickeringlab/tools/miniconda3/bin/
dir='/global/scratch/projects/fc_pickeringlab/Optimized-sgRNA-Design/'
module load python
source activate /global/scratch/projects/fc_pickeringlab/miniconda3/envs/bowtie
cd /global/scratch/projects/fc_pickeringlab/Optimized-sgRNA-Design/datafiles
bowtie-build 'pf3_np_subset.fasta' pf3_illumina_subset