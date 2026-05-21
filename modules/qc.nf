// ============================================================
// modules/qc.nf
// QC, read filtering, BUSCO, MultiQC
// ============================================================

process NANOPLOT {
    tag        { sample_id }
    label      'qc'
    publishDir { "${params.outdir}/qc/nanoplot/${sample_id}" }, mode: 'copy'
    container  'quay.io/biocontainers/nanoplot:1.43.0--pyhdfd78af_1'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("nanoplot_${sample_id}/*"), emit: report

    script:
    """
    NanoPlot \\
        --threads ${task.cpus} \\
        --fastq ${reads} \\
        --outdir nanoplot_${sample_id} \\
        --tsv_stats
    """
}

process FILTER_READS {
    tag        { sample_id }
    label      'qc'
    publishDir "${params.outdir}/reads_filtered", mode: 'copy'
    container  'quay.io/biocontainers/filtlong:0.2.1--hdcf5f25_3'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.fastq.gz"), emit: reads

    script:
    """
    cat ${reads} > combined.fastq.gz
    filtlong \\
        --min_length 1000 \\
        --min_mean_q 80 \\
        combined.fastq.gz \\
        | gzip > ${sample_id}.filtered.fastq.gz
    rm combined.fastq.gz
    """
}

process BUSCO_NUCLEAR {
    tag        { sample_id }
    label      'qc_heavy'
    publishDir { "${params.outdir}/qc/busco/${sample_id}" }, mode: 'copy'
    container  'quay.io/biocontainers/busco:5.8.0--pyhdfd78af_0'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("busco_${sample_id}/*"),                  emit: full
    tuple val(sample_id), path("busco_${sample_id}/short_summary*.txt"), emit: summary

    script:
    """
    busco \\
        --in ${assembly} \\
        --out busco_${sample_id} \\
        --mode genome \\
        --lineage_dataset ${params.busco_lineage} \\
        --cpu ${task.cpus}
    """
}

process MULTIQC {
    label     'qc'
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    container 'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'

    input:
    path '*'

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}
