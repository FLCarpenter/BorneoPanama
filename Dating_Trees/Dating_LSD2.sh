#!/bin/bash
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=phylodate
#SBATCH --no-requeue

source activate phylosoftware

iqtree -s 5_nt_supermatrix_6ktree_4184_Binary.fasta --date lsd2_dates_NTBIN_PARS9_fig8.txt -te pars9_root.tree -m GTR2+FO+R4 --date-tip 0 -o SRAA00097,SRAA00099,GBDL01738,SRAA00098,GBDL01734,SRAA00102,SRAA00101,SRAA00103,GBDL01736 --date-ci 100 -nt 5

iqtree -s 5_nt_supermatrix_6ktree_4184_Binary_BORNEO.fasta --date lsd2_dates_NTBIN_Borneo_fig8.txt -te Reconstructed_Borneo_NTBIN_rooted.tree -m GTR2+FO+R4 --date-tip 0 --date-ci 100 -nt 5

iqtree -s 5_nt_supermatrix_6ktree_4184_Binary_PANAMA.fasta --date lsd2_dates_NTBIN_Panama_fig8.txt -te Reconstructed_Panama_NTBIN_rooted.tree -m GTR2+FO+R4 --date-tip 0 --date-ci 100 -nt 5

