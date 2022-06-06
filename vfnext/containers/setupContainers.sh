echo "| --- setupContainers.sh -----------------------------------------------|"
echo "| USAGE: /bin/bash/ setupContainers.sh"
echo "| This scripts setups containers to be used at Viralflow under nextflow "
echo "| Be sure all required singularity recipies are present at working dir."
echo "| ----------------------------------------------------------------------|"

# --- SETUP ERROR HANDLING ----------------------------------------------------
# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT
# -----------------------------------------------------------------------------

echo "@ Building singularity_bam_readcount.sif..."
singularity build singularity_bam_readcount.sif docker://mgibio/bam-readcount
echo "  > Done <"

echo "@ Building singularity_mafft.sif..."
singularity build singularity_mafft.sif docker://staphb/mafft
echo "  > Done <"

echo "@ Building singularity_bedtools.sif..."
singularity build singularity_bedtools.sif docker://staphb/bedtools
echo "  > Done <"

echo "@ Building nextclade_latest.sif..."
singularity pull docker://nextstrain/nextclade
echo "  > loading sars-cov2 nextclade dataset..."
./nextclade_latest.sif nextclade dataset get --name 'sars-cov-2' --output-dir 'nextclade_dataset/sars-cov-2'
echo "  > Done <"

echo "@ Building pangolin_latest.sif..."
singularity pull docker://staphb/pangolin
echo "  > Done <"

echo "@ Building singularity_picard.sif..."
singularity build --fakeroot singularity_picard.sif Singularity_picard
echo "  > Done <"

echo "@ Building singularity_bwa_v0717.sif..."
singularity build --fakeroot singularity_bwa_v0717.sif Singularity_bwa_v0717
echo "  > Done <"

echo "@ Building singularity_fastp_v0201.sif..."
singularity build --fakeroot singularity_fastp_v0201.sif Singularity_fastp_v0201
echo "  > Done <"

echo "@ Building singularity_bwa_v0717.sif..."
singularity build --fakeroot ivar_latest.sif Singularity_ivar
echo "  > Done <"

echo "@ Building singularity_pythonScripts.sif..."
singularity build --fakeroot singularity_pythonScripts.sif Singularity_pythonScripts
echo "  > Done <"
