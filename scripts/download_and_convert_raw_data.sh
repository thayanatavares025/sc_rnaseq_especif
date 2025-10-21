#!/bin/sh

sra_id=$1
echo "Recebi ID: $sra_id"

# Use pastas dentro do diretório atual (workdir)
sra_dir="./sras"
fastq_dir="./fastqs"

mkdir -p "$sra_dir" "$fastq_dir"

# Download SRA file
prefetch "$sra_id" -O "$sra_dir" -X 30G

# Convert to FASTQ
fasterq-dump "$sra_dir/$sra_id" -O "$fastq_dir" --split-files --force

mv "$fastq_dir/${sra_id}_1.fastq" "$fastq_dir/${sra_id}_S0_L001_R1_001.fastq"
mv "$fastq_dir/${sra_id}_2.fastq" "$fastq_dir/${sra_id}_S0_L001_R2_001.fastq"

# Optional cleanup
# rm -f "$sra_dir/$sra_id"
