#!/usr/bin/env nextflow

// ============================================================
// Community structure analysis pipeline
// ============================================================
// Runs four community structure / co-occurrence analyses in
// parallel on a counts matrix + sample metadata:
//
//   AITCHISON_PCA  – CLR-based Aitchison PCA (01_aitchison_pca.R)
//   RPCA_DEICODE   – Robust PCA via rCLR + matrix completion (02_rpca_deicode.R)
//   FASTSPAR       – FastSpar co-occurrence network (03_fastspar.R)
//   NETWORK_NULL   – Null-model validation of FastSpar network (04_network_null_model.R)
//
// Usage (via run_pipeline.sh):
//   ./run_pipeline.sh --input counts.csv --metadata meta.csv
//
// All parameters and their defaults are in nextflow.config.
// ============================================================

nextflow.enable.dsl = 2

// ── pretty-print the run parameters ─────────────────────────
log.info """
=============================================================
  C O M M U N I T Y   S T R U C T U R E   P I P E L I N E
=============================================================
  input      : ${params.input}
  metadata   : ${params.metadata}
  taxonomy   : ${params.taxonomy ?: '(none)'}
  output_dir : ${params.output_dir}
  seed       : ${params.seed}
=============================================================
""".stripIndent()

// ── helper: resolve an absolute path for a file param ───────
def resolveFile(p) {
    def f = file(p)
    if (!f.exists()) { error "Input file not found: ${p}" }
    return f
}

// ============================================================
// PROCESS 1 – Aitchison PCA (CLR-based ordination)
// ============================================================
process AITCHISON_PCA {
    publishDir "${params.output_dir}/aitchison_pca", mode: 'copy'

    input:
    path counts_csv
    path metadata_csv

    output:
    path "aitchison_scores.pdf",   emit: scores_plot
    path "aitchison_loadings.pdf", emit: loadings_plot

    script:
    """
    Rscript /workspace/01_aitchison_pca.R \\
        '${counts_csv}' \\
        '${metadata_csv}'
    """
}

// ============================================================
// PROCESS 2 – Robust PCA / DEICODE (rCLR + matrix completion)
// ============================================================
process RPCA_DEICODE {
    publishDir "${params.output_dir}/rpca_deicode", mode: 'copy'

    input:
    path counts_csv
    path metadata_csv

    output:
    path "rpca_scores.pdf",   emit: scores_plot
    path "rpca_loadings.pdf", emit: loadings_plot

    script:
    """
    Rscript /workspace/02_rpca_deicode.R \\
        '${counts_csv}' \\
        '${metadata_csv}'
    """
}

// ============================================================
// PROCESS 3 – FastSpar co-occurrence network
// ============================================================
process FASTSPAR {
    publishDir "${params.output_dir}/fastspar", mode: 'copy'

    input:
    path counts_csv
    path metadata_csv
    path taxonomy_csv  // may be a dummy file when params.taxonomy is null

    output:
    path "fastspar_edges.csv",   emit: edges
    path "fastspar_heatmap.pdf", emit: heatmap, optional: true

    script:
    // Only pass the taxonomy argument when a real file is supplied
    def tax_arg = (params.taxonomy && taxonomy_csv.name != 'NO_TAXONOMY') \
        ? "'${taxonomy_csv}'" : ""
    """
    Rscript /workspace/03_fastspar.R \\
        '${counts_csv}' \\
        '${metadata_csv}' \\
        ${tax_arg} \\
        ${task.cpus}
    """
}

// ============================================================
// PROCESS 4 – Null-model validation of FastSpar network
// ============================================================
process NETWORK_NULL {
    publishDir "${params.output_dir}/network_null", mode: 'copy'

    input:
    path edges_csv
    path counts_csv

    output:
    path "null_model_comparison.pdf", emit: null_plot

    script:
    """
    Rscript /workspace/04_network_null_model.R \\
        '${edges_csv}' \\
        '${counts_csv}' \\
        ${task.cpus}
    """
}

// ============================================================
// WORKFLOW
// ============================================================
workflow {

    if (!params.input) {
        error "params.input must be set (path to counts CSV)"
    }
    if (!params.metadata) {
        error "params.metadata must be set (path to sample metadata CSV)"
    }

    input_file    = resolveFile(params.input)
    metadata_file = resolveFile(params.metadata)

    // Optional taxonomy file: use a dummy path when not provided so
    // FASTSPAR still receives a valid path value in its input block.
    taxonomy_file = params.taxonomy
        ? resolveFile(params.taxonomy)
        : file('NO_TAXONOMY')

    AITCHISON_PCA(input_file, metadata_file)
    RPCA_DEICODE(input_file, metadata_file)
    FASTSPAR(input_file, metadata_file, taxonomy_file)
    NETWORK_NULL(FASTSPAR.out.edges, input_file)
}
