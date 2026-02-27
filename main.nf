#!/usr/bin/env nextflow

// ============================================================
// Ordination pipeline — Aitchison PCA + Robust PCA
// ============================================================
// Runs two ordination analyses for every counts CSV found in
// params.input_dir, one per taxonomic level:
//
//   AITCHISON_PCA  – CLR-based Aitchison PCA (01_aitchison_pca.R)
//   RPCA_DEICODE   – Robust PCA via rCLR + matrix completion (02_rpca_deicode.R)
//
// Output layout:
//   <output_dir>/<level>/aitchison_pca/
//   <output_dir>/<level>/rpca_deicode/
//
// Usage (via run_pipeline.sh):
//   ./run_pipeline.sh --input_dir data/filtered/ --metadata metadata.csv
// ============================================================

nextflow.enable.dsl = 2

log.info """
=============================================================
  O R D I N A T I O N   P I P E L I N E
=============================================================
  input_dir  : ${params.input_dir}
  metadata   : ${params.metadata}
  output_dir : ${params.output_dir}
  seed       : ${params.seed}
=============================================================
""".stripIndent()

// ============================================================
// PROCESS 1 – Aitchison PCA (CLR-based ordination)
// ============================================================
process AITCHISON_PCA {
    tag "$level"
    publishDir "${params.output_dir}/${level}/aitchison_pca", mode: 'copy'

    input:
    tuple val(level), path(counts_csv)
    path metadata_csv

    output:
    path "aitchison_scores.pdf",   emit: scores_plot
    path "aitchison_loadings.pdf", emit: loadings_plot

    script:
    """
    Rscript /workspace/01_aitchison_pca.R \
        '${counts_csv}' \
        '${metadata_csv}'
    """
}

// ============================================================
// PROCESS 2 – Robust PCA / DEICODE (rCLR + matrix completion)
// ============================================================
process RPCA_DEICODE {
    tag "$level"
    publishDir "${params.output_dir}/${level}/rpca_deicode", mode: 'copy'

    input:
    tuple val(level), path(counts_csv)
    path metadata_csv

    output:
    path "rpca_scores.pdf",   emit: scores_plot
    path "rpca_loadings.pdf", emit: loadings_plot

    script:
    """
    Rscript /workspace/02_rpca_deicode.R \
        '${counts_csv}' \
        '${metadata_csv}'
    """
}

// ============================================================
// WORKFLOW
// ============================================================
workflow {

    if (!params.input_dir) {
        error "params.input_dir must be set (directory containing counts CSVs)"
    }
    if (!params.metadata) {
        error "params.metadata must be set (path to sample metadata CSV)"
    }

    metadata_file = file(params.metadata)
    if (!metadata_file.exists()) {
        error "Metadata file not found: ${params.metadata}"
    }

    // Discover all CSV files in input_dir and extract the taxonomic level
    // from the filename.  Handles both naming conventions:
    //   filtered_ancom_count_input_<level>.csv
    //   ancom_count_input_<level>.csv
    counts_ch = Channel
        .fromPath("${params.input_dir}/*.csv")
        .map { f ->
            def level = f.baseName.replaceAll(/^.*_input_/, '')
            tuple(level, f)
        }

    AITCHISON_PCA(counts_ch, metadata_file)
    RPCA_DEICODE(counts_ch, metadata_file)
}
