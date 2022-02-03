#!/bin/bash
#SBATCH --job-name=analyze_tet_A
#SBATCH -J download
#SBATCH --time=72:00:00
#SBATCH --account=fc_pickeringlab
#SBATCH -p savio2_htc
#SBATCH --ntasks=4
#SBATCH --mem=84g
#SBATCH --mail-type=END
#SBATCH --mail-user=shruteek@berkeley.edu
#SBATCH --output=%j_output.txt
#SBATCH --error=%j_error.txt
## Command(s) to run:
export PATH=/global/home/users/shruteek/.local/bin:$PATH
dir='/global/scratch/projects/fc_pickeringlab/target_analysis_project/Pickering Lab/Optimized-sgRNA-Design/'
module load python
python SpacerSequenceHandler.py MTI 'datafiles/tet_A_guides.tsv' '/global/scratch/projects/fc_pickeringlab/target_analysis_project/Illumina Sequences Folder/pf3_R1_illumina_full.fasta' datafiles/analysis2.csv >& job.out