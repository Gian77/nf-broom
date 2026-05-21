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

process FILTER_ORGANELLE_CONTIGS {
    tag       { "${sample_id}_${compartment}" }
    label     'map'
    publishDir { "${params.outdir}/assembly/${compartment}/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/minimap2:2.30--h577a1d6_0'

    input:
    tuple val(sample_id), val(compartment), path(raw_assembly), path(reference)

    output:
    tuple val(sample_id), val(compartment), path("${sample_id}_${compartment}_filtered.fasta"), emit: filtered
    tuple val(sample_id), val(compartment), path("${sample_id}_${compartment}_alignment.paf"),  emit: paf
    tuple val(sample_id), val(compartment), path("${sample_id}_${compartment}_kept_ids.txt"),    emit: ids

    script:
    """
    # Align raw assembly contigs to organelle reference.
    # PAF columns: 1=qname 2=qlen 3=qstart 4=qend 5=strand 6=tname 7=tlen 8=tstart 9=tend 10=matches 11=alnlen
    if [ -s ${raw_assembly} ]; then
        minimap2 -cx asm10 -t ${task.cpus} ${reference} ${raw_assembly} \\
            > ${sample_id}_${compartment}_alignment.paf
    else
        : > ${sample_id}_${compartment}_alignment.paf
    fi

    # Keep contigs where:
    #   - matched bases / contig length > 0.5 (covers majority of contig)
    #   - matched bases / aln length > 0.7 (high identity)
    # Summed across multiple alignments per contig (a contig may align in pieces).
    awk 'BEGIN{OFS="\\t"}
         {
             cov[\$1] += \$10
             aln[\$1] += \$11
             qlen[\$1] = \$2
         }
         END {
             for (c in cov) {
                 frac_q   = cov[c] / qlen[c]
                 identity = cov[c] / aln[c]
                 if (frac_q > 0.5 && identity > 0.7) print c
             }
         }' ${sample_id}_${compartment}_alignment.paf | sort -u \\
         > ${sample_id}_${compartment}_kept_ids.txt

    # Extract kept contigs (using awk; no seqkit dep needed)
    if [ -s ${sample_id}_${compartment}_kept_ids.txt ]; then
        awk 'BEGIN { while ((getline line < "'${sample_id}_${compartment}_kept_ids.txt'") > 0) keep[line]=1 }
             /^>/ { name = substr(\$1, 2); flag = (name in keep); if (flag) print; next }
             { if (flag) print }' \\
             ${raw_assembly} > ${sample_id}_${compartment}_filtered.fasta
    else
        echo "No contigs passed filter for ${compartment}" >&2
        : > ${sample_id}_${compartment}_filtered.fasta
    fi
    """
}
