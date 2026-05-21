// ============================================================
// modules/assembly.nf
// Flye assemblies for chloroplast, mitochondria, nuclear genome
// ============================================================

process ASSEMBLE_CP {
    tag       { sample_id }
    label     'assemble_small'
    publishDir { "${params.outdir}/assembly/chloroplast/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/flye:2.9.4--py310h2b6aa90_0'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_cp.fasta"),  emit: assembly
    path  "flye_cp_${sample_id}/*",                        emit: log

    script:
    """
    flye \\
        --nano-hq ${reads} \\
        --genome-size 150k \\
        --threads ${task.cpus} \\
        --out-dir flye_cp_${sample_id} \\
        --meta || true

    if [ -f flye_cp_${sample_id}/assembly.fasta ]; then
        cp flye_cp_${sample_id}/assembly.fasta ${sample_id}_cp.fasta
    else
        echo "Flye failed on chloroplast — emitting empty assembly"
        touch ${sample_id}_cp.fasta
    fi
    """
}

process ASSEMBLE_MT {
    tag       { sample_id }
    label     'assemble_small'
    publishDir { "${params.outdir}/assembly/mitochondria/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/flye:2.9.4--py310h2b6aa90_0'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_mt.fasta"),  emit: assembly
    path  "flye_mt_${sample_id}/*",                        emit: log

    script:
    """
    flye \\
        --nano-hq ${reads} \\
        --genome-size 500k \\
        --threads ${task.cpus} \\
        --out-dir flye_mt_${sample_id} \\
        --meta || true

    if [ -f flye_mt_${sample_id}/assembly.fasta ]; then
        cp flye_mt_${sample_id}/assembly.fasta ${sample_id}_mt.fasta
    else
        echo "Flye failed on mitochondria — emitting empty assembly"
        touch ${sample_id}_mt.fasta
    fi
    """
}

process ASSEMBLE_NUCLEAR {
    tag       { sample_id }
    label     'assemble_heavy'
    publishDir { "${params.outdir}/assembly/nuclear/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/flye:2.9.4--py310h2b6aa90_0'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_nuclear.fasta"), emit: assembly
    path  "flye_nuclear_${sample_id}/*",                       emit: log

    script:
    """
    flye \\
        --nano-hq ${reads} \\
        --genome-size ${params.genome_size} \\
        --threads ${task.cpus} \\
        --out-dir flye_nuclear_${sample_id}

    cp flye_nuclear_${sample_id}/assembly.fasta ${sample_id}_nuclear.fasta
    """
}

