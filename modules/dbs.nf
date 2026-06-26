// modules/db.nf
// ============================================================
// Reference database fetchers.
// Processes here download external reference data needed by the pipeline.
// They use `storeDir` so each version is cached persistently per
// outdir — re-running the pipeline reuses cached databases instead of
// re-downloading. Versions are pinned via params for reproducibility.
// ============================================================


process DOWNLOAD_OATKDB {
    label    'fetch'
    storeDir "${params.outdir}/databases/oatkdb/${params.oatkdb_version}/raw"
    // no container — uses host curl

    output:
    path "embryophyta_mito.fam", emit: mito
    path "embryophyta_pltd.fam", emit: pltd

    script:
    def base = "https://raw.githubusercontent.com/c-zhou/OatkDB/${params.oatkdb_commit}/${params.oatkdb_version}"
    """
    curl -fsSL -o embryophyta_mito.fam ${base}/embryophyta_mito.fam
    curl -fsSL -o embryophyta_pltd.fam ${base}/embryophyta_pltd.fam
    """
}

process PRESS_OATKDB {
    label     'fetch'
    storeDir  "${params.outdir}/databases/oatkdb/${params.oatkdb_version}/pressed"
    container 'docker://assteindorff/oatk:1.0'

    input:
    path mito_fam
    path pltd_fam

    output:
    tuple path("embryophyta_mito.fam"), path("embryophyta_mito.fam.h3*"), emit: mito
    tuple path("embryophyta_pltd.fam"), path("embryophyta_pltd.fam.h3*"), emit: pltd

    script:
    """
    hmmpress ${mito_fam}
    hmmpress ${pltd_fam}
    """
}

// Wrapper sub-workflow so main.nf still calls one thing
workflow FETCH_OATKDB {
    main:
    DOWNLOAD_OATKDB()
    PRESS_OATKDB(DOWNLOAD_OATKDB.out.mito, DOWNLOAD_OATKDB.out.pltd)

    emit:
    mito = PRESS_OATKDB.out.mito
    pltd = PRESS_OATKDB.out.pltd
}

// ============================================================
// Build OATK HMM databases from user-supplied FASTA sequences
// (chloroplast and mitochondrion complete genome assemblies).
// Emits the same channel format as FETCH_OATKDB so the rest of
// the pipeline is unchanged.
// ============================================================

process BUILD_OATKDB_FROM_FASTA {
    label     'fetch'
    storeDir  "${params.outdir}/databases/custom_oatkdb"
    container 'docker://assteindorff/oatk:1.0'

    input:
    tuple val(prefix), path(fasta)

    output:
    tuple val(prefix), path("${prefix}.fam"), path("${prefix}.fam.h3*")

    script:
    """
    hmmbuild --dna ${prefix}.fam ${fasta}
    hmmpress ${prefix}.fam
    """
}

workflow BUILD_OATKDB {
    take:
    cp_fasta   // Channel<Path>  — chloroplast reference FASTA
    mt_fasta   // Channel<Path>  — mitochondrion reference FASTA

    main:
    mito_in = mt_fasta.map { f -> tuple("custom_mito", f) }
    pltd_in = cp_fasta.map { f -> tuple("custom_pltd", f) }

    BUILD_OATKDB_FROM_FASTA(mito_in.mix(pltd_in))

    emit:
    mito = BUILD_OATKDB_FROM_FASTA.out
               .filter { label, fam, idx -> label == "custom_mito" }
               .map    { label, fam, idx -> tuple(fam, idx) }
    pltd = BUILD_OATKDB_FROM_FASTA.out
               .filter { label, fam, idx -> label == "custom_pltd" }
               .map    { label, fam, idx -> tuple(fam, idx) }
}