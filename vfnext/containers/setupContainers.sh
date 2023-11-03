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
singularity build --fakeroot --sandbox singularity_bam_readcount.sif docker://mgibio/bam-readcount
echo "  > Done <"

echo "@ Building singularity_mafft.sif..."
singularity build --fakeroot --sandbox singularity_mafft.sif docker://staphb/mafft
echo "  > Done <"

echo "@ Building singularity_bedtools.sif..."
singularity build --fakeroot --sandbox singularity_bedtools.sif docker://staphb/bedtools
echo "  > Done <"

echo "  > loading sars-cov2 nextclade dataset..."
singularity build --fakeroot --sandbox singularity_nextclade.sif Singularity_nextclade
singularity exec -B nextclade_dataset/sars-cov-2:/tmp singularity_nextclade.sif nextclade dataset get --name 'sars-cov-2' --output-dir '/tmp'
echo "  > Done <"

echo "@ Building pangolin_latest.sif..."
#singularity pull docker://staphb/pangolin
#singularity build --fakeroot --sandbox my_sbox my_image.sif
singularity build --fakeroot --sandbox pangolin_latest.sif Singularity_pangolin
echo "  > Done <"

echo "@ Building singularity_picard.sif..."
singularity build --fakeroot --sandbox singularity_picard.sif Singularity_picard
echo "  > Done <"

echo "@ Building singularity_bwa_v0717.sif..."
singularity build --fakeroot --sandbox singularity_bwa_v0717.sif Singularity_bwa_v0717
echo "  > Done <"

echo "@ Building singularity_fastp_v0201.sif..."
singularity build --fakeroot --sandbox singularity_fastp_v0201.sif Singularity_fastp_v0201
echo "  > Done <"

echo "@ Building ivar_latest.sif..."
singularity build --fakeroot --sandbox ivar_latest.sif Singularity_ivar
echo "  > Done <"

echo "@ Building singularity_edirect.sif..."
singularity build --fakeroot --sandbox singularity_edirect.sif Singularity_edirect

echo "@ Building singularity_snpeff.sif..."
singularity build --fakeroot --sandbox singularity_snpeff.sif Singularity_snpEff

echo "@ Building singularity_samtools.sif..."
singularity build --fakeroot --sandbox singularity_samtools.sif docker://staphb/samtools:1.15

echo "@ Building vfreport.sif..."
singularity build --fakeroot --sandbox vfreport.sif Singularity_vfreport

echo "@ Building singularity_bamdash.sif..."
singularity build --fakeroot --sandbox singularity_bamdash.sif Singularity_bamdash

echo "@ Downloading snpeff database catalog..."
singularity exec ./singularity_snpeff.sif snpEff databases > snpEff_DB.catalog

echo "  > Done <"
