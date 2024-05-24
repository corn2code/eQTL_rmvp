#!/bin/sh
#SBATCH --array=1
#SBATCH --job-name=eGWAS
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=50GB
#SBATCH --output=out.log/%x_%a.out
#SBATCH --partition=jclarke,batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vladimir.torres@unl.edu

### Initialization
# Get Array ID
i=${SLURM_ARRAY_TASK_ID}

### Parameter file to read 
pheno="mads69ID.txt"
echo "Parameter file: ${pheno}"
echo

# Read line #i from the parameter file
PHENO=$(sed "${i}q;d" ${pheno})
echo "Parameters are: ${PHENO}"
echo

ml anaconda/4.12 R/4.1 cmake/3.20
conda activate rMVP

Rscript e.GWAS.rMVP2020mads69.R ${PHENO} ${outFile} --verbose 
