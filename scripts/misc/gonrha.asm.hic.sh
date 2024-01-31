#!/bin/bash
#SBATCH --account=EcoGenetics
#SBATCH --cpus-per-task=30
#SBATCH --mem=160g
#SBATCH --time=10:00:00

# Sources necessary environment
if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    #shellcheck disable=1091
    source activate genome_assembly
fi

output_dir=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/genome_assembly_and_annotation/steps/museomics/Gonepteryx_rhamni/draft_assembly/hifiasm

[ -d "$output_dir" ] || mkdir -p "$output_dir"

hifiasm \
    -t 30 \
    -o "$output_dir"/GonRha.asm \
    -l 3 \
    -s 0.1 \
    --h1 /faststorage/project/EcoGenetics/BACKUP/HiC/GonRha/Gon_Rha_L1_1.fq.gz \
    --h2 /faststorage/project/EcoGenetics/BACKUP/HiC/GonRha/Gon_Rha_L1_2.fq.gz \
    /faststorage/project/EcoGenetics/BACKUP/PacBio_HiFi/Gonepteryx_rhamni/HiFiAdapterFilt/GonRha_HMW1_m64101e_231223_171722.bc1011--bc1011.hifi_reads.filt.fastq.gz