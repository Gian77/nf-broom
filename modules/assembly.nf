// ============================================================
// modules/assembly.nf
// Flye assemblies (cp/mt/nuclear) and OATK alternative for organelles
// ============================================================

process ASSEMBLE_CP_FLYE {
    tag       { sample_id }
    label     'assemble_small'
    publishDir { "${params.outdir}/assembly/chloroplast/${sample_id}/raw" }, mode: 'copy'
    container 'quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_cp_raw.fasta"), emit: assembly
    path  "flye_cp_${sample_id}/*",                          emit: log

    script:
    """
    flye \\
        --nano-hq ${reads} \\
        --genome-size 150k \\
        --threads ${task.cpus} \\
        --out-dir flye_cp_${sample_id} \\
        --min-overlap 5000 \\
        --asm-coverage 50 || true

    if [ -f flye_cp_${sample_id}/assembly.fasta ]; then
        cp flye_cp_${sample_id}/assembly.fasta ${sample_id}_cp_raw.fasta
    else
        echo "Flye failed on chloroplast — emitting empty assembly" >&2
        touch ${sample_id}_cp_raw.fasta
    fi
    """
}

process ASSEMBLE_MT_FLYE {
    tag       { sample_id }
    label     'assemble_small'
    publishDir { "${params.outdir}/assembly/mitochondria/${sample_id}/raw" }, mode: 'copy'
    container 'quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_mt_raw.fasta"), emit: assembly
    path  "flye_mt_${sample_id}/*",                          emit: log

    script:
    """
    flye \\
        --nano-hq ${reads} \\
        --genome-size 500k \\
        --threads ${task.cpus} \\
        --out-dir flye_mt_${sample_id} \\
        --min-overlap 5000 || true

    if [ -f flye_mt_${sample_id}/assembly.fasta ]; then
        cp flye_mt_${sample_id}/assembly.fasta ${sample_id}_mt_raw.fasta
    else
        echo "Flye failed on mitochondria — emitting empty assembly" >&2
        touch ${sample_id}_mt_raw.fasta
    fi
    """
}

process ASSEMBLE_ORGANELLES_OATK {
    tag       { sample_id }
    label     'assemble_small'
    publishDir { "${params.outdir}/assembly/oatk/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/oatk:1.0--h4d35ad6_1'
    // Adjust tag after verifying with: quay_tags oatk

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_cp_raw.fasta"), emit: cp_assembly
    tuple val(sample_id), path("${sample_id}_mt_raw.fasta"), emit: mt_assembly
    path  "oatk_${sample_id}/*",                              emit: log

    script:
    """
    mkdir -p oatk_${sample_id}
    cd oatk_${sample_id}

    # OATK takes ALL reads (post-filtering) and detects organelle reads via HMM
    # Defaults: -k 1001, -c 30 (min coverage) — adjust if needed
    oatk \\
        -k 1001 \\
        -c 30 \\
        -t ${task.cpus} \\
        -m \$(which embryophyta_mito.fam || echo embryophyta_mito.fam) \\
        -p \$(which embryophyta_pltd.fam || echo embryophyta_pltd.fam) \\
        -o ${sample_id} \\
        ../${reads} || true

    cd ..

    # Copy outputs to canonical names (handle empty/missing gracefully)
    if [ -f oatk_${sample_id}/${sample_id}.pltd.ctg.fasta ]; then
        cp oatk_${sample_id}/${sample_id}.pltd.ctg.fasta ${sample_id}_cp_raw.fasta
    else
        echo "OATK produced no chloroplast assembly" >&2
        touch ${sample_id}_cp_raw.fasta
    fi

    if [ -f oatk_${sample_id}/${sample_id}.mito.ctg.fasta ]; then
        cp oatk_${sample_id}/${sample_id}.mito.ctg.fasta ${sample_id}_mt_raw.fasta
    else
        echo "OATK produced no mitochondrial assembly" >&2
        touch ${sample_id}_mt_raw.fasta
    fi
    """
}

// Nuclear assembly unchanged
process ASSEMBLE_NUCLEAR {
    tag       { sample_id }
    label     'assemble_heavy'
    publishDir { "${params.outdir}/assembly/nuclear/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1'

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

