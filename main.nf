#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


if ( !params.sra_ids )    error "No SRA ids path identified! Use --sra_ids <file>"
if ( !params.metadata_file )   error "No metadata file identified! Use --metadata_file <file>"

// ---- Input ----
Channel
  .fromPath(params.sra_ids)
  .splitText()
  .map { it.trim() }
  .filter { it }
  .set { sra_ids }

Channel
  .of(params.human_ref_link)
  .set { human_ref_link_ch }

Channel
  .fromPath(params.metadata_file)
  .set { metadata_file }


process download_and_convert_raw_data {
  label 'sra_tools'
  tag "${sra_id}"
  
  input:
  val sra_id

  output:
  tuple val(sra_id), path("fastqs")

  script:
  """
  sh ${projectDir}/scripts/download_and_convert_raw_data.sh ${sra_id}
  """
}

process download_human_genome {
  label 'align_tools'
  tag "download_genome"

  input:
  val human_ref_link

  output:
  path "refs/human_genome_ref_2024"

  script:
  """
  set -euo pipefail
  FILE=\$(basename "${human_ref_link}")
  TARGET_DIR=refs
  mkdir -p "\$TARGET_DIR"

  wget -O "\$TARGET_DIR/\$FILE" "${human_ref_link}"
  tar -xzf "\$TARGET_DIR/\$FILE" -C "\$TARGET_DIR"

  if [ -d "\$TARGET_DIR/refdata-gex-GRCh38-2024-A" ]; then
    mv "\$TARGET_DIR/refdata-gex-GRCh38-2024-A" "\$TARGET_DIR/human_genome_ref_2024"
  fi
  """
}

process align_reads_against_human_genome {
  label 'align_tools'
  tag "${sra_id}"

  input:
  tuple val(sra_id), path(fastqs), path(genome_ref_dir)

  output:
  tuple val(sra_id), path("${sra_id}_cellranger")

  script:
  def fastq_arg = (fastqs instanceof List) ? fastqs.collect{ it.toString() }.join(',') : fastqs.toString()

  """
  set -euo pipefail
  cellranger count \\
    --id=${sra_id}_cellranger \\
    --transcriptome=${genome_ref_dir} \\
    --fastqs=${fastq_arg} \\
    --create-bam=true \\
    --sample=${sra_id}
  """
}

process run_quality_control {
  label 'R_tools'
  tag "${sra_id}"

  input:
  tuple val(sra_id), path(cellranger_dir), path(metadata_file)

  output:
  tuple val(sra_id), path("${sra_id}_seurat_object.RDS", optional: true)

  script:
  """
  set -euo pipefail
  matrix_dir="${cellranger_dir}/outs/filtered_feature_bc_matrix"
  csv_metrics="${cellranger_dir}/outs/metrics_summary.csv"
  Rscript ${projectDir}/scripts/quality_control.R "\$matrix_dir" "\$csv_metrics" "$sra_id"

  """
}

process run_normalization{
  label 'R_tools'
  tag "${sra_id}"
  
  input:
  tuple val(sra_id), path(seurat_object)

  output:
  path("${sra_id}_seurat_object.RDS", optional: true)
  
  script:
  """
  set -euo pipefail
  Rscript ${projectDir}/scripts/normalize_samples.R "$sra_id" "$seurat_object"
  """
}


workflow {
  def genome_ref_ch = human_ref_link_ch | download_human_genome | first()
  def fastq_ch = sra_ids | download_and_convert_raw_data
  def aligned_ch = fastq_ch
      .combine(genome_ref_ch)              // -> (sra_id, fastqs, genome_ref_path)
      .map { sra_id, fastqs, genome_ref -> tuple(sra_id, fastqs, genome_ref) }
      | align_reads_against_human_genome   // -> (sra_id, sra_id_cellranger_dir)
  def qc_ch = aligned_ch
      .combine(metadata_file)
      .map { sra_id, cellranger_dir, metadata ->
          tuple(sra_id, cellranger_dir, metadata)
      }
      | run_quality_control                // -> (sra_id, seurat_object.RDS)
  def norm_ch = qc_ch \
      | run_normalization
}
