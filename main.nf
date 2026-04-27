#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================================================================
// H&E Spatial Heterogeneity Analysis Pipeline
// Wraps: he_hypothesis_testing.R
// Container: he-analysis:latest
// =============================================================================

// Help message
if (params.help) {
    log.info """
    =========================================
     H&E Spatial Heterogeneity Pipeline
    =========================================
    Usage:
        nextflow run main.nf -profile vm_docker \\
            --input_dir /path/to/h5files \\
            --output_dir /path/to/results

    Parameters:
        --input_dir   Directory containing *_Probabilities.h5 files [required]
        --output_dir  Directory for all output files  [default: ./results]
        --help        Show this message and exit
    =========================================
    """.stripIndent()
    exit 0
}

// Validate required parameters
if (!params.input_dir) {
    error "ERROR: --input_dir is required.\nRun with --help for usage."
}

// =============================================================================
// Process
// =============================================================================

process HE_ANALYSIS {

    tag "HE analysis – ${h5_files.size()} samples"

    container params.container__he_analysis

    // Copy all outputs to params.output_dir on the host after the process finishes
    publishDir params.output_dir, mode: 'copy', overwrite: true

    cpus   4
    memory '48 GB'
    time   '2h'

    input:
    // All *_Probabilities.h5 files collected into the work directory together
    path h5_files

    output:
    // Capture everything written under the output/ subdirectory
    path "results/**", emit: results

    script:
    """
    echo "=== Staged files in work dir ==="
    ls -la .
    echo "================================"

    mkdir -p results

    Rscript /app/he_hypothesis_testing.R \\
        --input_dir . \\
        --output_dir results
    """
}

// =============================================================================
// Workflow
// =============================================================================

workflow {

    log.info """
    =========================================
     H&E Spatial Heterogeneity Pipeline
    =========================================
    input_dir  : ${params.input_dir}
    output_dir : ${params.output_dir}
    =========================================
    """.stripIndent()

    // Collect all matching H5 files from input_dir into a single list
    // .collect() ensures they are all staged together into one work directory,
    // which is required because the R script processes all samples jointly.
    h5_ch = channel
        .fromPath("${params.input_dir}/*_Probabilities.h5")
        .ifEmpty { error "No *_Probabilities.h5 files found in: ${params.input_dir}" }
        .collect()

    HE_ANALYSIS(h5_ch)
}
