// Chromosome scaffolding with RagTag.
// Runs after purge_dups: corrects chimeric contigs against the nuclear reference,
// then orders/orients them into reference-based chromosome groups, bridging
// neighbouring contigs with a 100 bp N-gap (RagTag's default scaffold gap).

process RAGTAG_SCAFFOLD {
    tag       { sample_id }
    label     'assemble_small'
    publishDir { "${params.outdir}/assembly/scaffolds/${sample_id}" }, mode: 'copy'
    container 'quay.io/biocontainers/ragtag:2.1.0--pyhdfd78af_2'

    input:
    tuple val(sample_id), path(assembly)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}_scaffold.fasta"),     emit: scaffold
    tuple val(sample_id), path("${sample_id}_scaffold.agp"),       emit: agp
    path "ragtag_${sample_id}",                                    emit: log

    script:
    """
    # 1. Correction — split contigs that disagree with the reference (chimeras).
    ragtag.py correct ${reference} ${assembly} \\
        -t ${task.cpus} \\
        -o ragtag_${sample_id}/correct

    # 2. Scaffold — order & orient corrected contigs into reference chromosome
    #    groups. RagTag inserts a 100 bp N-string between joined contigs (-g 100,
    #    the default) so scaffolds form continuous chromosome-scale sequences.
    ragtag.py scaffold ${reference} ragtag_${sample_id}/correct/ragtag.correct.fasta \\
        -t ${task.cpus} \\
        -o ragtag_${sample_id}/scaffold

    cp ragtag_${sample_id}/scaffold/ragtag.scaffold.fasta ${sample_id}_scaffold.fasta
    cp ragtag_${sample_id}/scaffold/ragtag.scaffold.agp   ${sample_id}_scaffold.agp
    """
}
