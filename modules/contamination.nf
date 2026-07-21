// Contaminant screening (gated behind --run_kraken2).
// Kraken2 classifies each assembled contig against the shared NCBI reference DB, then its
// per-contig taxids are turned into a BlobTools "hits" file so the coverage-vs-GC blob plot
// is coloured by taxonomy — non-Viridiplantae contigs stand out as candidate contamination.
//
// The PlusPF DB (~104 GB hash) fits in RAM, so Kraken2 loads it fully (no --memory-mapping)
// — one sequential read, then fast classification. PlusPF ships its own nodes.dmp/names.dmp,
// which BlobTools uses to resolve taxids to lineage.

process KRAKEN2_CLASSIFY {
    tag           { "${sample_id}_${stage}" }
    label         'qc_heavy'   // 160 GB override set via withName in nextflow.config
    errorStrategy 'ignore'     // advisory screening — never fail the pipeline
    publishDir    { "${params.outdir}/qc/kraken2/${sample_id}" }, mode: 'copy'
    // NB: biocontainers kraken2 images use a busybox base whose hardlinked applets
    // (linuxrc, usr/bin/[) fail to unpack under this cluster's rootless apptainer
    // ("unpriv.link: operation not permitted"). staphb/kraken2 is Ubuntu-based with
    // separate coreutils binaries (no problematic hardlinks) and unpacks cleanly.
    container     'quay.io/staphb/kraken2:2.1.3'

    input:
    tuple val(sample_id), val(stage), path(assembly), path(bam), path(bai)
    path kraken2_db

    output:
    tuple val(sample_id), val(stage), path(assembly), path(bam), path(bai),
          path("${sample_id}_${stage}.kraken2.hits"),                               emit: hits
    tuple val(sample_id), val(stage), path("${sample_id}_${stage}.kraken2.report"), emit: report

    script:
    def pfx = "${sample_id}_${stage}"
    """
    # PlusPF (~104 GB) loads fully into the 160 GB RAM allocation — no --memory-mapping.
    kraken2 --db ${kraken2_db} \\
        --threads ${task.cpus} \\
        --output ${pfx}.kraken2.out \\
        --report ${pfx}.kraken2.report \\
        ${assembly}

    # BlobTools hits format: seqID <tab> taxID <tab> score. Kraken2 .out columns are
    # 1=C/U, 2=contigID, 3=taxID; keep classified contigs with a real taxid. Score is
    # unused (one hit per contig).
    awk -F'\\t' '\$1=="C" && \$3!="0" {print \$2"\\t"\$3"\\t1"}' ${pfx}.kraken2.out > ${pfx}.kraken2.hits
    """
}

process BLOBTOOLS_TAXONOMY {
    tag           { "${sample_id}_${stage}" }
    label         'qc'
    errorStrategy 'ignore'
    publishDir    { "${params.outdir}/qc/blobtools/${sample_id}" }, mode: 'copy'
    container     'quay.io/biocontainers/blobtools:1.1.1--py_1'

    input:
    tuple val(sample_id), val(stage), path(assembly), path(bam), path(bai), path(hits)
    path taxdump

    output:
    path "blobtax_${sample_id}_${stage}*", emit: results

    script:
    def pfx = "blobtax_${sample_id}_${stage}"
    """
    # Build the BlobDB with taxonomy. --nodes/--names point at the taxdump so the Kraken2
    # taxids resolve to a lineage; default taxrule (bestsum) picks the per-contig hit.
    # --db gives a WRITABLE cwd path for the built nodesDB: without it blobtools tries to
    # cache nodesDB.txt back into its read-only package dir and dies (Errno 30). The file
    # does not exist yet, so blobtools builds it here from --nodes/--names.
    blobtools create \\
        -i ${assembly} \\
        -b ${bam} \\
        -t ${hits} \\
        --nodes ${taxdump}/nodes.dmp \\
        --names ${taxdump}/names.dmp \\
        --db nodesDB.txt \\
        -o ${pfx}

    # Blob plot + table coloured by phylum; prefix outputs so they match the publish glob.
    blobtools plot -i ${pfx}.blobDB.json --rank phylum --out ${pfx}
    blobtools view -i ${pfx}.blobDB.json --rank phylum -o ${pfx}
    """
}
