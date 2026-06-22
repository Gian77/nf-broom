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
    tuple val(sample_id), path("${sample_id}_scaffold.fasta"),                         emit: scaffold
    tuple val(sample_id), path("${sample_id}_scaffold.agp"),                           emit: agp
    tuple val(sample_id), path("ragtag_${sample_id}/scaffold/ragtag.scaffold.stats"),  emit: stats
    path "ragtag_${sample_id}",                                                         emit: log

    script:
    def do_correct   = params.ragtag_correct
    def scaffold_input = do_correct ? "ragtag_${sample_id}/correct/ragtag.correct.fasta"
                                    : "${assembly}"
    """
    mkdir -p ragtag_${sample_id}/correct

    ${do_correct ? """
    # Correction: split chimeric contigs against the reference
    ragtag.py correct ${reference} ${assembly} \\
        -t ${task.cpus} \\
        -o ragtag_${sample_id}/correct
    """ : "# correction skipped (params.ragtag_correct = false)"}

    # Scaffold: order & orient contigs into chromosome-scale groups
    ragtag.py scaffold ${reference} ${scaffold_input} \\
        -t ${task.cpus} \\
        -o ragtag_${sample_id}/scaffold

    cp ragtag_${sample_id}/scaffold/ragtag.scaffold.fasta ${sample_id}_scaffold.fasta
    cp ragtag_${sample_id}/scaffold/ragtag.scaffold.agp   ${sample_id}_scaffold.agp
    """
}
