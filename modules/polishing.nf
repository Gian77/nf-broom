// ============================================================
// modules/polishing.nf
// Medaka polishing, purge_dups, HapDup phasing
// ============================================================

process POLISH_MEDAKA {
    tag       { sample_id }
    label     'polish'
    publishDir { "${params.outdir}/polishing/medaka/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/medaka:2.2.1--py312hc7af5e1_0'

    input:
    tuple val(sample_id), path(assembly), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_polished.fasta"), emit: assembly
    path "medaka_${sample_id}/*",                               emit: log

    script:
    """
    medaka_consensus \\
        -i ${reads} \\
        -d ${assembly} \\
        -o medaka_${sample_id} \\
        -t ${task.cpus} \\
        -m ${params.medaka_model}

    cp medaka_${sample_id}/consensus.fasta ${sample_id}_polished.fasta
    """
}

process PURGE_DUPS {
    tag       { sample_id }
    label     'polish'
    publishDir { "${params.outdir}/polishing/purge_dups/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/purge_dups:1.2.6--h577a1d6_3'

    input:
    tuple val(sample_id), path(assembly), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_purged.fasta"), emit: assembly
    path  "purge_dups_${sample_id}/*",                        emit: log

    script:
    """
    mkdir -p purge_dups_${sample_id}
    cd purge_dups_${sample_id}

    # Step 1: map reads to assembly to compute depth
    minimap2 -t ${task.cpus} -xmap-ont ../${assembly} ../${reads} \\
        | gzip -c > reads.paf.gz
    pbcstat reads.paf.gz
    calcuts PB.stat > cutoffs 2> calcuts.log

    # Step 2: self-align contigs to find duplicates
    split_fa ../${assembly} > assembly.split
    minimap2 -t ${task.cpus} -xasm5 -DP assembly.split assembly.split \\
        | gzip -c > assembly.self.paf.gz

    # Step 3: purge
    purge_dups -2 -T cutoffs -c PB.base.cov assembly.self.paf.gz > dups.bed 2> purge_dups.log
    get_seqs -e dups.bed ../${assembly}

    cd ..
    cp purge_dups_${sample_id}/purged.fa ${sample_id}_purged.fasta
    """
}


process ALIGN_FOR_HAPDUP {
    tag       { sample_id }
    label     'map'
    container 'quay.io/biocontainers/minimap2:2.31--h118bc1c_0'

    input:
    tuple val(sample_id), path(assembly), path(reads)

    output:
    tuple val(sample_id), path(assembly), path("aligned.sam"), emit: sam

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont ${assembly} ${reads} > aligned.sam
    """
}

process SORT_FOR_HAPDUP {
    tag       { sample_id }
    label     'map'
    container 'quay.io/biocontainers/samtools:1.6--h5fe306e_13'

    input:
    tuple val(sample_id), path(assembly), path(sam)

    output:
    tuple val(sample_id), path(assembly), path("aligned.bam"), path("aligned.bam.bai"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o aligned.bam ${sam}
    samtools index aligned.bam
    """
}

process HAPDUP {
    tag       { sample_id }
    label     'assemble_heavy'
    publishDir { "${params.outdir}/polishing/hapdup/${sample_id}" }, mode: 'copy'
    container 'docker://mkolmogo/hapdup:0.12'

    input:
    tuple val(sample_id), path(assembly), path(bam), path(bai)

    output:
    tuple val(sample_id), path("hapdup_${sample_id}/hapdup_dual_*.fasta"), emit: haplotypes
    path "hapdup_${sample_id}/*",                                           emit: log

    script:
    """
    hapdup \\
        --assembly ${assembly} \\
        --bam ${bam} \\
        --out-dir hapdup_${sample_id} \\
        --threads ${task.cpus} \\
        --rtype ont
    """
}
