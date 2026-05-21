#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ============================================================
// nf-broom plant genome assembly pipeline — ONT
// Phase 1+2: QC → organelle split → assemble → polish → purge → BUSCO
// ============================================================

// ---- Parameters ----
params.reads        = "${projectDir}/reads"
params.cp_ref       = null
params.mt_ref       = null
params.nuclear_ref  = null
params.outdir       = "results"
params.genome_size  = "720m"
params.busco_lineage = "poales_odb10"
params.medaka_model = "r1041_e82_400bps_sup_v5.0.0"
params.run_hapdup   = false
params.help = false

// ---- Help message ---- (top-level function declaration — this is allowed)
def helpMessage() {
    log.info """
    ===================================================================
     nf-broom — PLANT GENOME ASSEMBLY PIPELINE (ONT)
    ===================================================================
    Usage:
      nextflow run main.nf -profile condor --reads <dir> --cp_ref <fa> --mt_ref <fa> [options]

    Mandatory arguments:
      --reads [path]        Directory of input FASTQ (one subdir per sample)
      --cp_ref [path]       Chloroplast reference FASTA
      --mt_ref [path]       Mitochondrion reference FASTA

    Optional arguments:
      --nuclear_ref [path]  Nuclear reference FASTA               [default: none]
      --outdir [path]       Output directory                      [default: results]
      --genome_size [str]   Estimated nuclear genome size         [default: 720m]
      --busco_lineage [str] BUSCO lineage dataset                 [default: poales_odb10]
      --medaka_model [str]  Medaka model                          [default: r1041_e82_400bps_sup_v5.0.0]
      --run_hapdup [bool]   Run HapDup phasing step               [default: false]
      --help                Show this message and exit
    ===================================================================
    """.stripIndent()
}

// ---- Module imports ----
include { NANOPLOT;           FILTER_READS;       MULTIQC }       from './modules/qc.nf'
include { BUSCO_NUCLEAR }                                          from './modules/qc.nf'
include { ALIGN_TO_ORGANELLES; SORT_INDEX_BAM; EXTRACT_RAW_READSETS; DEDUP_ORGANELLE_READS; READSET_STATS } from './modules/mapping.nf'
include { ASSEMBLE_CP;       ASSEMBLE_MT;        ASSEMBLE_NUCLEAR } from './modules/assembly.nf'
include { POLISH_MEDAKA; PURGE_DUPS; ALIGN_FOR_HAPDUP; SORT_FOR_HAPDUP; HAPDUP } from './modules/polishing.nf'

// ============================================================
//  Workflow
// ============================================================
workflow {

    if (params.help) {
        helpMessage()
        return            // exits the workflow cleanly — no `exit 0` needed
    }

    // ---- Banner ----
    log.info """
        ╔══════════════════════════════════════════════════════╗
        ║                        nt-broom                      ║
        ║             PLANT GENOME ASSEMBLY PIPELINE           ║
        ║  Chloroplast · Mitochondria · Nuclear separation     ║
        ╚══════════════════════════════════════════════════════╝
        reads dir     : ${params.reads}
        cp reference  : ${params.cp_ref}
        mt reference  : ${params.mt_ref}
        nuclear ref   : ${params.nuclear_ref ?: '(none — Phase 3)'}
        output dir    : ${params.outdir}
        genome size   : ${params.genome_size}
        busco lineage : ${params.busco_lineage}
        medaka model  : ${params.medaka_model}
        run HapDup    : ${params.run_hapdup}
        """.stripIndent()

    // ---- Reference channels ----
    cp_ref_ch = Channel.value(file(params.cp_ref, checkIfExists: true))
    mt_ref_ch = Channel.value(file(params.mt_ref, checkIfExists: true))

    // ---- Sample channel ----
    raw_reads_ch = Channel
        .fromPath("${params.reads}/*/*.fastq.gz")
        .map { f -> tuple(f.parent.name, f) }
        .groupTuple()

   // 1. QC + filtering
    NANOPLOT(raw_reads_ch)
    FILTER_READS(raw_reads_ch)

    // 2a. Align filtered reads to combined cp+mt reference (minimap2)
    ALIGN_TO_ORGANELLES(FILTER_READS.out.reads, cp_ref_ch, mt_ref_ch)

    // 2b. Sort and index → BAM (samtools)
    SORT_INDEX_BAM(ALIGN_TO_ORGANELLES.out.sam)

    // 3a. Extract per-compartment FASTQ (raw cp/mt may have dup IDs from supplementary alignments)
    EXTRACT_RAW_READSETS(SORT_INDEX_BAM.out.bam)

    // 3b. Dedup organelle FASTQs by read name
    dedup_in = EXTRACT_RAW_READSETS.out.cp_raw
        .join(EXTRACT_RAW_READSETS.out.mt_raw)
    DEDUP_ORGANELLE_READS(dedup_in)

    // 3c. Per-compartment stats and estimated coverage
    readset_stats_in = DEDUP_ORGANELLE_READS.out.cp_reads
        .join(DEDUP_ORGANELLE_READS.out.mt_reads)
        .join(EXTRACT_RAW_READSETS.out.nuclear_reads)
        .combine(cp_ref_ch)
        .combine(mt_ref_ch)
    READSET_STATS(readset_stats_in)

   // 4. Organelle assemblies (use deduplicated reads)
    ASSEMBLE_CP(DEDUP_ORGANELLE_READS.out.cp_reads)
    ASSEMBLE_MT(DEDUP_ORGANELLE_READS.out.mt_reads)

    // 5. Nuclear assembly (no dedup needed — unmapped reads only appear once)
    ASSEMBLE_NUCLEAR(EXTRACT_RAW_READSETS.out.nuclear_reads)

    // 6. Polish nuclear assembly with Medaka
    polish_in = ASSEMBLE_NUCLEAR.out.assembly
        .join(EXTRACT_RAW_READSETS.out.nuclear_reads)
    POLISH_MEDAKA(polish_in)

    // 7. Purge haplotigs
    purge_in = POLISH_MEDAKA.out.assembly
        .join(EXTRACT_RAW_READSETS.out.nuclear_reads)
    PURGE_DUPS(purge_in)

    // 8. Optional: HapDup phasing
    if (params.run_hapdup) {
        // Align reads to purged assembly
        align_in = PURGE_DUPS.out.assembly
            .join(EXTRACT_RAW_READSETS.out.nuclear_reads)
        ALIGN_FOR_HAPDUP(align_in)

        // Sort and index BAM
        SORT_FOR_HAPDUP(ALIGN_FOR_HAPDUP.out.sam)

        // Run HapDup with the prepared BAM
        HAPDUP(SORT_FOR_HAPDUP.out.bam)
    }

    // 9. BUSCO — nuclear only
    BUSCO_NUCLEAR(PURGE_DUPS.out.assembly)

    // 10. MultiQC aggregation
    multiqc_in = Channel.empty()
        .mix( NANOPLOT.out.report.map      { sample, f -> f } )
        .mix( READSET_STATS.out.mqc.map    { sample, f -> f } )
        .mix( BUSCO_NUCLEAR.out.summary.map { sample, f -> f } )
        .collect()
    MULTIQC(multiqc_in)

    // ---- Completion handler ----
    def outdir = params.outdir ?: 'NA'
    
    workflow.onComplete = { meta ->
        log.info "Pipeline completed | Status: ${meta?.successful ? 'SUCCESS' : 'FAILED'} | Duration: ${meta?.duration ?: 'NA'} | Output: ${outdir}"
    }
}
