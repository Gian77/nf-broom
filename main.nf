#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ============================================================
// nf-broom plant genome assembly pipeline — ONT
// Phase 1+2: QC → organelle split → assemble → polish → purge → BUSCO
// ============================================================

// ---- Parameters ----
params.reads          = "${projectDir}/reads"
params.cp_ref         = null
params.mt_ref         = null
params.oatkdb_version = "v20230921"
params.oatkdb_commit  = "75e8db0ac4a7d508a9a518d900876003ceb70737"
params.oatk_mito_db   = "https://raw.githubusercontent.com/c-zhou/OatkDB/${params.oatkdb_commit}/${params.oatkdb_version}/embryophyta_mito.fam"
params.oatk_pltd_db   = "https://raw.githubusercontent.com/c-zhou/OatkDB/${params.oatkdb_commit}/${params.oatkdb_version}/embryophyta_pltd.fam"
params.oatkdb_recipe  = "v1-pressed"
params.outdir         = "results"
params.genome_size    = "720m"
params.busco_lineage  = "poales_odb10"
params.medaka_model   = "r1041_e82_400bps_sup_v5.0.0"
params.run_hapdup     = false
params.help           = false

params.organelle_assembler = "flye"    // "flye" (default) or "oatk"

// ---- Help message ---- (top-level function declaration — this is allowed)
def helpMessage() {
    log.info """
    ===================================================================
     nf-broom — PLANT GENOME ASSEMBLY PIPELINE (ONT)
    ===================================================================
    Usage:
      nextflow run main.nf -profile condor --reads <dir> --cp_ref <fa> --mt_ref <fa> [options]

    Mandatory arguments:
      --reads [path]                Directory of input FASTQ (one subdir per sample)
      --cp_ref [path]               Chloroplast reference FASTA
      --mt_ref [path]               Mitochondrion reference FASTA

    Organelle assembly options:
      --organelle_assembler [str]   Organelle assembler: 'flye' or 'oatk'           [default: flye]
                                    'oatk' uses the embryophyta OatkDB (v20230921).
      --filter_organelles [bool]    Apply reference-based filter post-assembly      [default: true]
      --organelle_min_qcov [float]  Min fraction of contig covered by ref alignment [default: 0.5]
      --organelle_min_ident [float] Min identity (matches / alignment length)       [default: 0.7]

    Nuclear assembly options:
      --genome_size [str]           Estimated nuclear genome size                   [default: 800m]
      --nuclear_ref [path]          Nuclear reference FASTA (for scaffolding)       [default: none]
      --run_hapdup [bool]           Run HapDup phasing step                         [default: false]

    Polishing & QC options:
      --medaka_model [str]          Medaka model (null = auto-detect from reads)    [default: null]
      --busco_lineage [str]         BUSCO lineage dataset                           [default: poales_odb10]

    General options:
      --outdir [path]               Output directory                                [default: results]
      --help                        Show this message and exit
    ===================================================================
    """.stripIndent()
}

// ---- Module imports ----
// Imports are grouped by source file (one include per module), and within each
// group processes are listed in the order they appear in the workflow below.
// Module order itself follows the pipeline stages: QC → mapping → assembly → polishing.

// QC: read stats, filtering, BUSCO completeness, MultiQC aggregation
include {NANOPLOT; FILTER_READS; FILTER_ORGANELLE_CONTIGS; BUSCO_NUCLEAR; MULTIQC; BANDAGE_IMAGE; QUAST_ORGANELLE; QUAST_NUCLEAR; ALIGN_FOR_QC; SORT_FOR_QC; QUALIMAP_BAMQC; BLOBTOOLS_COVERAGE} from './modules/qc.nf'
include {FETCH_OATKDB} from './modules/dbs.nf'
include {ALIGN_TO_ORGANELLES; SORT_INDEX_BAM; EXTRACT_RAW_READSETS; DEDUP_ORGANELLE_READS; READSET_STATS} from './modules/mapping.nf'
include {ASSEMBLE_CP_FLYE; ASSEMBLE_MT_FLYE; ASSEMBLE_ORGANELLES_OATK; ASSEMBLE_NUCLEAR} from './modules/assembly.nf'
include {POLISH_MEDAKA; PURGE_DUPS; SKIP_PURGE_MARKER; ALIGN_FOR_HAPDUP; SORT_FOR_HAPDUP; HAPDUP} from './modules/polishing.nf'
include {RAGTAG_SCAFFOLD; RAGTAG_SCAFFOLD as RAGTAG_PREPURGE} from './modules/scaffolding.nf'
include {FINAL_SUMMARY; TOOLS_REPORT} from './modules/reports.nf'

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
        reads dir           : ${params.reads}
        cp reference        : ${params.cp_ref}
        mt reference        : ${params.mt_ref}
        nuclear ref         : ${params.nuclear_ref ?: '(none — scaffolding off)'}
        organelle assembler : ${params.organelle_assembler}${params.organelle_assembler == 'oatk' ? "  (OatkDB ${params.oatkdb_version})" : ''}
        filter organelles   : ${params.filter_organelles}
        output dir          : ${params.outdir}
        genome size         : ${params.genome_size}
        busco lineage       : ${params.busco_lineage}
        medaka model        : ${params.medaka_model}
        run HapDup          : ${params.run_hapdup}
        skip purge_dups     : ${params.skip_purge}
        """.stripIndent()

    // ---- Reference channels ----
    cp_ref_ch = Channel.value(file(params.cp_ref, checkIfExists: true))
    mt_ref_ch = Channel.value(file(params.mt_ref, checkIfExists: true))
    // Nuclear reference is optional — only needed for RagTag scaffolding.
    nuclear_ref_ch = params.nuclear_ref \
        ? Channel.value(file(params.nuclear_ref, checkIfExists: true))
        : Channel.empty()

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

   // 4. Organelle assemblies — branch on assembler choice
    if (params.organelle_assembler == "oatk") {
        // Always use the pre-built embryophyta OatkDB.
        // Building custom HMMs from whole-genome FASTAs (hmmbuild on 140-468 kb
        // sequences) creates profiles too large for HMMER's DP matrix — integer
        // overflow at scan time. The standard DB covers sorghum organelle genes.
        FETCH_OATKDB()
        oatk_mito_db = FETCH_OATKDB.out.mito
        oatk_pltd_db = FETCH_OATKDB.out.pltd
        ASSEMBLE_ORGANELLES_OATK(
            FILTER_READS.out.reads,
            oatk_mito_db,
            oatk_pltd_db
        )
        cp_raw_ch = ASSEMBLE_ORGANELLES_OATK.out.cp_assembly
        mt_raw_ch = ASSEMBLE_ORGANELLES_OATK.out.mt_assembly

        // Visualize assembly graphs with Bandage (published next to filtered assemblies)
        bandage_in = ASSEMBLE_ORGANELLES_OATK.out.cp_gfa
            .map { id, gfa -> tuple(id, 'chloroplast', gfa) }
            .mix(
                ASSEMBLE_ORGANELLES_OATK.out.mt_gfa
                    .map { id, gfa -> tuple(id, 'mitochondria', gfa) }
            )
        BANDAGE_IMAGE(bandage_in)
    } else {
        ASSEMBLE_CP_FLYE(DEDUP_ORGANELLE_READS.out.cp_reads)
        ASSEMBLE_MT_FLYE(DEDUP_ORGANELLE_READS.out.mt_reads)
        cp_raw_ch = ASSEMBLE_CP_FLYE.out.assembly
        mt_raw_ch = ASSEMBLE_MT_FLYE.out.assembly
    }

    // 4b. Reference-based filtering — remove nuclear contamination (NUPTs/NUMTs)
    if (params.filter_organelles) {
        cp_filter_in = cp_raw_ch.map { id, fa -> tuple(id, 'chloroplast', fa) }
                                .combine(cp_ref_ch)
        mt_filter_in = mt_raw_ch.map { id, fa -> tuple(id, 'mitochondria', fa) }
                                .combine(mt_ref_ch)

        all_organelle_in = cp_filter_in.mix(mt_filter_in)
        FILTER_ORGANELLE_CONTIGS(all_organelle_in)

        // Final outputs (filtered)
        cp_final = FILTER_ORGANELLE_CONTIGS.out.filtered
            .filter { id, comp, fa -> comp == 'chloroplast' }
            .map { id, comp, fa -> tuple(id, fa) }
        mt_final = FILTER_ORGANELLE_CONTIGS.out.filtered
            .filter { id, comp, fa -> comp == 'mitochondria' }
            .map { id, comp, fa -> tuple(id, fa) }
    } else {
        cp_final = cp_raw_ch
        mt_final = mt_raw_ch
    }

    // 4c. QUAST on filtered organelle assemblies (uses cp/mt references)
    quast_organelle_in = cp_final
        .map { id, fa -> tuple(id, 'chloroplast', fa) }
        .combine(cp_ref_ch)
        .mix(
            mt_final
                .map { id, fa -> tuple(id, 'mitochondria', fa) }
                .combine(mt_ref_ch)
        )
    QUAST_ORGANELLE(quast_organelle_in)

    // 5. Nuclear assembly (no dedup needed — unmapped reads only appear once)
    ASSEMBLE_NUCLEAR(EXTRACT_RAW_READSETS.out.nuclear_reads)

    // 6. Polish nuclear assembly with Medaka
    polish_in = ASSEMBLE_NUCLEAR.out.assembly
        .join(EXTRACT_RAW_READSETS.out.nuclear_reads)
    POLISH_MEDAKA(polish_in)

    // 7. Purge haplotigs (or bypass with --skip_purge to go medaka → scaffold directly)
    if (params.skip_purge) {
        purge_final      = POLISH_MEDAKA.out.assembly
        SKIP_PURGE_MARKER(POLISH_MEDAKA.out.assembly)
        purge_cutoffs_ch = SKIP_PURGE_MARKER.out.cutoffs
        purge_calcuts_ch = SKIP_PURGE_MARKER.out.calcuts_log
    } else {
        purge_in = POLISH_MEDAKA.out.assembly
            .join(EXTRACT_RAW_READSETS.out.nuclear_reads)
        PURGE_DUPS(purge_in)
        purge_final      = PURGE_DUPS.out.assembly
        purge_cutoffs_ch = PURGE_DUPS.out.cutoffs
        purge_calcuts_ch = PURGE_DUPS.out.calcuts_log
    }

    // 8. Optional: HapDup phasing
    if (params.run_hapdup) {
        // Align reads to purge-stage assembly
        align_in = purge_final
            .join(EXTRACT_RAW_READSETS.out.nuclear_reads)
        ALIGN_FOR_HAPDUP(align_in)

        // Sort and index BAM
        SORT_FOR_HAPDUP(ALIGN_FOR_HAPDUP.out.sam)

        // Run HapDup with the prepared BAM
        HAPDUP(SORT_FOR_HAPDUP.out.bam)
    }

    // 8b. Chromosome scaffolding with RagTag (correct + scaffold), if a nuclear
    //     reference is provided. QC then runs on the scaffolded assembly.
    if (params.nuclear_ref) {
        // Canonical (purge-stage) scaffold — the pipeline's final nuclear assembly.
        RAGTAG_SCAFFOLD(purge_final.map { id, fa -> tuple(id, 'purge', fa) }, nuclear_ref_ch)
        nuclear_final = RAGTAG_SCAFFOLD.out.scaffold

        // When purge_dups runs, also scaffold the pre-purge (Medaka) assembly so the
        // QUAST table shows medaka_scaf vs purge scaffold side-by-side. This exposes
        // over-purging: if purge_dups discards unique sequence, medaka_scaf retains a
        // much higher genome fraction. Skipped under --skip_purge (no purge to compare).
        if (!params.skip_purge) {
            RAGTAG_PREPURGE(POLISH_MEDAKA.out.assembly.map { id, fa -> tuple(id, 'medaka', fa) }, nuclear_ref_ch)
            medaka_scaffold_ch  = RAGTAG_PREPURGE.out.scaffold
            has_medaka_scaffold = true
        } else {
            medaka_scaffold_ch  = RAGTAG_SCAFFOLD.out.scaffold   // placeholder, excluded via flag
            has_medaka_scaffold = false
        }

        // All assembly stages + reference for direct comparison
        quast_nuclear_in = ASSEMBLE_NUCLEAR.out.assembly
            .join(POLISH_MEDAKA.out.assembly)
            .join(medaka_scaffold_ch)
            .join(purge_final)
            .join(RAGTAG_SCAFFOLD.out.scaffold)
            .map { id, flye, medaka, mscaf, purge, scaffold ->
                tuple(id, flye, medaka, mscaf, purge, scaffold)
            }
        QUAST_NUCLEAR(quast_nuclear_in, Channel.value(true), Channel.value(has_medaka_scaffold), nuclear_ref_ch, Channel.value(params.skip_purge))
    } else {
        nuclear_final = purge_final

        // Three stages without scaffold or reference (medaka_scaf / scaffold are placeholders)
        quast_nuclear_in = ASSEMBLE_NUCLEAR.out.assembly
            .join(POLISH_MEDAKA.out.assembly)
            .join(purge_final)
            .map { id, flye, medaka, purge ->
                tuple(id, flye, medaka, purge, purge, purge)  // asm_medaka_scaffold + asm_scaffold unused
            }
        QUAST_NUCLEAR(quast_nuclear_in, Channel.value(false), Channel.value(false), Channel.value([]), Channel.value(params.skip_purge))
    }

    // 9b. BUSCO — nuclear only
    BUSCO_NUCLEAR(nuclear_final)

    // 9c. Optional: BAM-level QC (Qualimap) and coverage blob plots (BlobTools)
    //     One sorted BAM per assembly stage is computed once and shared between both tools.
    if (params.run_qualimap || params.run_blobtools) {
        qc_stages = purge_final
            .map { sid, fa -> tuple(sid, "purge", fa) }

        if (params.nuclear_ref) {
            qc_stages = qc_stages.mix(
                RAGTAG_SCAFFOLD.out.scaffold
                    .map { sid, fa -> tuple(sid, "scaffold", fa) }
            )
        }

        qc_align_in = qc_stages.combine(EXTRACT_RAW_READSETS.out.nuclear_reads, by: 0)
        ALIGN_FOR_QC(qc_align_in)
        SORT_FOR_QC(ALIGN_FOR_QC.out.sam)

        if (params.run_qualimap)  QUALIMAP_BAMQC(SORT_FOR_QC.out.bam)
        if (params.run_blobtools) BLOBTOOLS_COVERAGE(SORT_FOR_QC.out.bam)
    }

    // 10. MultiQC aggregation
    // QUAST reports are staged as their whole output directory (uniquely named
    // per compartment) rather than the bare report.tsv — three files all named
    // report.tsv collide when staged flat. MultiQC recurses into each dir.
    multiqc_in = Channel.empty()
        .mix( NANOPLOT.out.report.map                        { sample, f -> f } )
        .mix( READSET_STATS.out.mqc.map                      { sample, f -> f } )
        .mix( BUSCO_NUCLEAR.out.summary.map                  { sample, f -> f } )
        .mix( QUAST_ORGANELLE.out.report.map                 { sample, comp, d -> d } )
        .mix( QUAST_NUCLEAR.out.report.map                   { sample, d -> d } )
    if (!params.skip_purge) {
        multiqc_in = multiqc_in
            .mix( PURGE_DUPS.out.pbstat.map  { sample, f -> f } )
            .mix( PURGE_DUPS.out.cutoffs.map { sample, f -> f } )
    }

    if (params.run_qualimap) {
        multiqc_in = multiqc_in.mix(
            QUALIMAP_BAMQC.out.report.map { sid, stage, d -> d }
        )
    }

    MULTIQC(multiqc_in.collect())

    // 11. Human-readable per-sample summary
    nano_stats_ch = NANOPLOT.out.report
        .map { id, files ->
            def list = files instanceof List ? files : [files]
            tuple(id, list.find { it.name == 'NanoStats.txt' })
        }
    summary_in = nano_stats_ch
        .join(QUAST_NUCLEAR.out.report)
        .join(BUSCO_NUCLEAR.out.summary)
        .join(purge_cutoffs_ch)
        .join(purge_calcuts_ch)
        .join(RAGTAG_SCAFFOLD.out.stats)
    FINAL_SUMMARY(summary_in)
    TOOLS_REPORT()

    // Capture outdir into a local — `params` resolves to null inside the
    // onComplete closure when it fires, so reference the captured value.
    def outdir = params.outdir ?: 'NA'
    workflow.onComplete {
        log.info "Pipeline completed | Output: ${outdir}"
    }
}
