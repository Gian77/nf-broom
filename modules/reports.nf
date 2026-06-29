// Human-readable final summary.
// Pulls the key metrics from NanoPlot (read set), QUAST (nuclear assembly), and
// BUSCO (completeness) into one markdown table per sample. Reporting only —
// errorStrategy 'ignore' so a missing/odd input never fails the pipeline.
// Reuses the already-pulled multiqc image (has bash/grep/awk/python).

process FINAL_SUMMARY {
    tag           { sample_id }
    label         'qc'
    errorStrategy 'ignore'
    publishDir    { "${params.outdir}/reports/${sample_id}" }, mode: 'symlink'
    container     'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'

    input:
    tuple val(sample_id), path(nanostats), path(quast_nuclear), path(busco_summary),
          path(purge_cutoffs), path(purge_calcuts_log), path(ragtag_stats)

    output:
    tuple val(sample_id), path("${sample_id}_assembly_summary.md"), emit: summary

    script:
    def busco_lin = params.busco_lineage
    def sample    = sample_id
    def gsize     = params.genome_size
    """
    #!/usr/bin/env bash
    set -uo pipefail

    NANO="${nanostats}"
    QUAST="${quast_nuclear}/report.tsv"
    BUSCO=\$(ls ${busco_summary} 2>/dev/null | head -1)
    CUTOFFS="${purge_cutoffs}"
    CALCUTS="${purge_calcuts_log}"
    RAGTAG="${ragtag_stats}"

    # ── helper: NanoPlot field lookup ───────────────────────────────────────
    nano() { awk -F'\\t' -v k="\$1" 'tolower(\$1) ~ tolower(k) {print \$NF; exit}' "\$NANO" 2>/dev/null || echo NA; }

    # ── helper: QUAST row as markdown ───────────────────────────────────────
    qrow() {
        awk -F'\\t' -v k="\$1" '
            NR==1 { ncols=NF }
            \$1==k { printf "| **%s** |", \$1; for(i=2;i<=ncols;i++) printf " %s |", \$i; print ""; exit }
        ' "\$QUAST" 2>/dev/null
    }

    qheader() {
        awk -F'\\t' '
            NR==1 { printf "| Metric |"; for(i=2;i<=NF;i++) printf " **%s** |", \$i; print "";
                    printf "|---|";       for(i=2;i<=NF;i++) printf "---|";           print "" }
        ' "\$QUAST" 2>/dev/null
    }

    # ── helper: extract single QUAST cell by metric + column-header pattern ─
    qval() {
        local metric="\$1" colpat="\$2"
        awk -F'\\t' -v m="\$metric" -v p="\$colpat" '
            NR==1 { for(i=2;i<=NF;i++) if(tolower(\$i) ~ tolower(p)) { col=i; break } }
            \$1==m && col { print \$col; exit }
        ' "\$QUAST" 2>/dev/null || echo 0
    }

    # ── extract metrics ──────────────────────────────────────────────────────
    TOTAL_BASES=\$(nano 'number.of.bases')
    READ_N50=\$(nano 'n50')
    MEAN_QUAL=\$(nano 'mean.qual')

    GFRAC_FLYE=\$(qval "Genome fraction (%)" "flye")
    GFRAC_PURGE=\$(qval "Genome fraction (%)" "purge")
    GFRAC_SCAFFOLD=\$(qval "Genome fraction (%)" "scaffold")
    # Pre-purge (Medaka) RagTag scaffold, when present — column label "medaka_scaf"
    # deliberately lacks the substring "scaffold" so the lookup above is unaffected.
    GFRAC_MEDAKA_SCAF=\$(qval "Genome fraction (%)" "medaka_scaf")
    N50_FLYE=\$(qval "N50" "flye")
    N50_SCAFFOLD=\$(qval "N50" "scaffold")
    LEN_FLYE=\$(qval "Total length" "flye")
    LEN_PURGE=\$(qval "Total length" "purge")
    NCONTIGS_FLYE=\$(qval "# contigs" "flye")
    MISASM_FLYE=\$(qval "# misassemblies" "flye")
    MISASM_SCAFFOLD=\$(qval "# misassemblies" "scaffold")

    COVERAGE=\$(echo "${gsize}" | awk '{
        s=tolower(\$0)
        if (sub(/g\$/, "", s)) s=s*1e9
        else if (sub(/m\$/, "", s)) s=s*1e6
        else if (sub(/k\$/, "", s)) s=s*1e3
        print s
    }' | awk -v tb="\$TOTAL_BASES" 'BEGIN{tb+=0} {if(\$1>0) printf "%.0fx", tb/\$1; else print "NAx"}')

    PURGE_RETAINED=\$(echo "\$LEN_PURGE \$LEN_FLYE" | awk '{if(\$2>0) printf "%.1f", \$1/\$2*100; else print "NA"}')
    PURGE_REMOVED=\$(echo "\$LEN_PURGE \$LEN_FLYE"  | awk '{if(\$2>0) printf "%.1f", (1-\$1/\$2)*100; else print "NA"}')
    C_JUNK=\$(awk '{print \$1}' "\$CUTOFFS" 2>/dev/null || echo "?")
    C_HAP_LOW=\$(awk '{print \$2}' "\$CUTOFFS" 2>/dev/null || echo "?")
    C_HAP_HIGH=\$(awk '{print \$3}' "\$CUTOFFS" 2>/dev/null || echo "?")
    C_DIP_LOW=\$(awk '{print \$4}' "\$CUTOFFS" 2>/dev/null || echo "?")
    C_DIP_HIGH=\$(awk '{print \$5}' "\$CUTOFFS" 2>/dev/null || echo "?")
    C_REPEAT=\$(awk '{print \$6}' "\$CUTOFFS" 2>/dev/null || echo "?")

    PLACED_MB=\$(awk -F'\\t' 'NR==2 {printf "%.0f", \$2/1e6}' "\$RAGTAG" 2>/dev/null || echo 0)
    UNPLACED_MB=\$(awk -F'\\t' 'NR==2 {printf "%.0f", \$4/1e6}' "\$RAGTAG" 2>/dev/null || echo 0)
    PLACED_PCT=\$(echo "\$PLACED_MB \$UNPLACED_MB" | awk '{t=\$1+\$2; if(t>0) printf "%.0f", \$1/t*100; else print "NA"}')

    BUSCO_C=\$(grep 'C:[0-9]' "\$BUSCO" 2>/dev/null | grep -o 'C:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)
    BUSCO_S=\$(grep 'S:[0-9]' "\$BUSCO" 2>/dev/null | grep -o 'S:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)
    BUSCO_D=\$(grep 'D:[0-9]' "\$BUSCO" 2>/dev/null | grep -o 'D:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)
    BUSCO_F=\$(grep 'F:[0-9]' "\$BUSCO" 2>/dev/null | grep -o 'F:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)
    BUSCO_M=\$(grep 'M:[0-9]' "\$BUSCO" 2>/dev/null | grep -o 'M:[0-9][0-9.]*' | cut -d: -f2 | head -1 || echo 0)

    PEAK=\$(grep "autotune" "\$CALCUTS" 2>/dev/null | grep -o 'Peak: [0-9]*x' | grep -o '[0-9]*' || echo "unknown")
    SKIP_PURGE=0
    grep -q "skip_purge" "\$CALCUTS" 2>/dev/null && SKIP_PURGE=1 || true

    QUALITY_TIER=\$(awk -v c="\$BUSCO_C" -v gf="\$GFRAC_SCAFFOLD" -v n50="\$N50_SCAFFOLD" 'BEGIN {
        c=c+0; gf=gf+0; n50=n50+0; score=0
        if (c  >= 95) score+=2; else if (c  >= 90) score++
        if (gf >= 80) score+=2; else if (gf >= 60) score++
        if (n50 >= 20000000) score+=2; else if (n50 >= 5000000) score++
        if (score >= 5) print "high"
        else if (score >= 3) print "moderate"
        else print "low"
    }')

    # ── interpretation helpers ───────────────────────────────────────────────
    read_interp() {
        awk -v q="\$MEAN_QUAL" -v n50="\$READ_N50" -v cov="\$COVERAGE" 'BEGIN {
            q=q+0; n50=n50+0; cov=cov+0
            if (q >= 20)      qmsg="excellent (Q>=20)"
            else if (q >= 15) qmsg="good (Q15-Q20)"
            else if (q >= 12) qmsg="acceptable (Q12-Q15)"
            else              qmsg="low (Q<12) -- polishing essential"
            if (n50 >= 30000)      nmsg="excellent (>=30 kb)"
            else if (n50 >= 15000) nmsg="good (15-30 kb)"
            else if (n50 >= 8000)  nmsg="acceptable (8-15 kb)"
            else                   nmsg="short (<8 kb) -- may limit contiguity"
            if (cov >= 50)      cmsg="deep -- well-suited for polishing and purging"
            else if (cov >= 25) cmsg="adequate for ONT assembly"
            else if (cov >= 15) cmsg="moderate -- assembly may be fragmented"
            else                cmsg="low -- additional sequencing recommended"
            printf "Mean quality: %s. Read N50: %s. Coverage ~%dx (%s).", qmsg, nmsg, cov, cmsg
        }'
    }

    assembly_interp() {
        awk -v gf="\$GFRAC_SCAFFOLD" -v n50="\$N50_SCAFFOLD" -v misasm="\$MISASM_SCAFFOLD" -v mf="\$MISASM_FLYE" 'BEGIN {
            gf=gf+0; n50=n50+0; misasm=misasm+0; mf=mf+0
            if (gf >= 85)      gfmsg="high -- assembly covers most of the reference"
            else if (gf >= 70) gfmsg="moderate -- some regions absent or divergent"
            else if (gf >= 50) gfmsg="low -- significant sequence missing or highly divergent"
            else               gfmsg="very low -- highly divergent from reference, or assembly fragmented"
            if (n50 >= 50000000)      n50msg="chromosome-scale (>=50 Mb)"
            else if (n50 >= 10000000) n50msg="near-chromosome-scale (10-50 Mb)"
            else if (n50 >= 1000000)  n50msg="Mb-scale (1-10 Mb)"
            else                      n50msg="sub-Mb -- scaffolding may be limited by contig fragmentation"
            printf "Scaffold N50 is %.1f Mb. Reference genome fraction: %.1f%% (%s). Misassemblies: %d (Flye) -> %d (scaffold).", n50/1e6, gf, gfmsg, mf, misasm
        }'
    }

    purge_interp() {
        awk -v ret="\$PURGE_RETAINED" -v gf_pre="\$GFRAC_FLYE" -v gf_post="\$GFRAC_PURGE" -v peak="\$PEAK" \
            -v junk="\$C_JUNK" -v hl="\$C_HAP_LOW" -v hh="\$C_HAP_HIGH" \
            -v dl="\$C_DIP_LOW" -v dh="\$C_DIP_HIGH" -v rep="\$C_REPEAT" 'BEGIN {
            ret=ret+0; gf_pre=gf_pre+0; gf_post=gf_post+0; drop=gf_pre-gf_post
            removed=sprintf("%.1f", 100-ret)
            if (peak != "unknown") {
                pmsg = "Autotune detected haploid coverage peak at " peak "x. "
                zmsg = "Coverage zones: junk (<" junk "x), haploid-kept (" hl "-" hh "x, centered on " peak "x), haplotig-purged (" dl "-" dh "x), repeat (>" rep "x). "
            } else {
                pmsg = "Coverage peak auto-detection unavailable. "
                zmsg = ""
            }
            mech = "Removed " removed "% of assembled sequence (genome fraction: " sprintf("%.1f",gf_pre) "% -> " sprintf("%.1f",gf_post) "%, drop " sprintf("%.1f",drop) "%). "
            mech = mech "Purging is driven primarily by the self-alignment overlap step: contigs that map redundantly against the primary assembly are classified as haplotigs and discarded; the coverage cutoffs assign each contig to a depth class before overlap testing. "
            if (drop > 15)
                verdict = "Genome fraction fell " sprintf("%.1f",drop) " points -- purge_dups removed UNIQUE reference sequence, not just haplotigs (over-purging). A true haplotig purge leaves genome fraction roughly flat. Strongly consider --skip_purge (or manual --calcuts_args), especially at low coverage. See the pre-purge comparison below."
            else if (ret < 55)
                verdict = "Large removal -- possible over-purging for a homozygous sample. Consider --calcuts_args to set manual cutoffs, or skip purge_dups."
            else if (ret < 75)
                verdict = "Moderate removal, typical for a heterozygous plant genome where alternate-haplotype contigs are assembled separately at similar depth to the primary contigs."
            else
                verdict = "Conservative removal -- sample is likely largely homozygous with minimal haplotig duplication in the assembly."
            print pmsg zmsg mech verdict
        }'
    }

    scaffold_interp() {
        awk -v pct="\$PLACED_PCT" -v n50c="\$N50_FLYE" -v n50s="\$N50_SCAFFOLD" 'BEGIN {
            pct=pct+0; n50c=n50c+0; n50s=n50s+0
            if (pct >= 85)      anch="excellent anchoring"
            else if (pct >= 70) anch="good anchoring"
            else if (pct >= 50) anch="moderate anchoring"
            else                anch="low anchoring -- many contigs unplaced, possibly from highly divergent or novel sequence"
            fold=int(n50s/n50c)
            printf "%s (%d%% of bases placed on chromosomes). N50 improved %d-fold: %.1f kb (contigs) -> %.1f Mb (scaffolds).", anch, pct, fold, n50c/1e3, n50s/1e6
        }'
    }

    busco_interp() {
        awk -v c="\$BUSCO_C" -v s="\$BUSCO_S" -v d="\$BUSCO_D" -v f="\$BUSCO_F" -v m="\$BUSCO_M" 'BEGIN {
            c=c+0; s=s+0; d=d+0; f=f+0; m=m+0
            if (c >= 95 && d <= 3)      qual="Excellent gene-space recovery"
            else if (c >= 90 && d <= 5) qual="Good gene-space recovery"
            else if (c >= 80)           qual="Moderate gene-space recovery"
            else                        qual="Poor gene-space recovery -- significant gene loss detected"
            notes=""
            if (d > 5)  notes=notes " Duplication (" d "%) elevated -- residual haplotigs may remain."
            if (m > 10) notes=notes " Missing BUSCOs (" m "%) -- check assembly completeness."
            if (f > 5)  notes=notes " Fragmented BUSCOs (" f "%) indicate assembly fragmentation."
            printf "%s (complete: %s%%, single: %s%%, duplicated: %s%%, fragmented: %s%%, missing: %s%%).", qual, c, s, d, f, m
            if (notes != "") printf " Notes:%s", notes
            printf "\\n"
        }'
    }

    {
      echo "# Assembly summary — ${sample}"
      echo
      echo "_Generated by nf-broom on \$(date -u '+%Y-%m-%d %H:%M UTC')_"
      echo

      # ── 1. Read set ───────────────────────────────────────────────────────
      echo "## 1. Read set (NanoPlot)"
      echo
      echo "| Metric | Value |"
      echo "|---|---|"
      echo "| Number of reads     | \$(nano 'number.of.reads') |"
      echo "| Total bases         | \$TOTAL_BASES |"
      echo "| Estimated coverage  | \$COVERAGE |"
      echo "| Read N50            | \$READ_N50 bp |"
      echo "| Median read length  | \$(nano 'median.read.length') bp |"
      echo "| Mean read quality   | Q\$MEAN_QUAL |"
      echo
      echo "> **Interpretation:** \$(read_interp)"
      echo

      # ── 2. Nuclear assembly progression ──────────────────────────────────
      echo "## 2. Nuclear assembly progression (QUAST)"
      echo
      qheader
      qrow "# contigs"
      qrow "Total length"
      qrow "Largest contig"
      qrow "N50"
      qrow "GC (%)"
      qrow "Genome fraction (%)"
      qrow "# misassemblies"
      qrow "# N's per 100 kbp"
      echo
      echo "> **Interpretation:** \$(assembly_interp)"
      echo

      # ── 3. Purge_dups ─────────────────────────────────────────────────────
      echo "## 3. Haplotig removal (Purge_dups)"
      echo
      if [ "\$SKIP_PURGE" -gt 0 ]; then
        echo '> **Note:** Purge_dups was skipped (`--skip_purge true`). Assembly proceeded directly from Medaka polishing to scaffolding without haplotig removal. To compare, rerun without this flag.'
      else
        echo "Coverage cutoffs (junk / hap-low / hap-high / dip-low / dip-high / repeat):"
        echo
        echo '```'
        cat "\$CUTOFFS" 2>/dev/null || echo "NA"
        echo '```'
        echo
        WARN=\$(grep -i "warn" "\$CALCUTS" 2>/dev/null | head -3 || true)
        if [ -n "\$WARN" ]; then
          echo "> **Warning:** \$WARN"
          echo
        fi
        echo "> **Interpretation:** \$(purge_interp)"
        # Pre-purge comparison: RagTag on the Medaka assembly (no purge_dups).
        if [ -n "\$GFRAC_MEDAKA_SCAF" ] && [ "\$GFRAC_MEDAKA_SCAF" != "0" ]; then
          echo
          echo "> **Pre-purge comparison:** scaffolding the Medaka assembly *without* purge_dups (the medaka_scaf column) retains genome fraction \${GFRAC_MEDAKA_SCAF}% vs \${GFRAC_SCAFFOLD}% after purge_dups. A large gap means purge_dups discarded unique sequence (over-purging) rather than haplotigs -- consider running with --skip_purge."
        fi
      fi
      echo

      # ── 4. Scaffolding ────────────────────────────────────────────────────
      echo "## 4. Chromosomal scaffolding (RagTag)"
      echo
      echo "| | Sequences | Bases |"
      echo "|---|---|---|"
      awk -F'\\t' 'NR==2 {
          placed_mb   = sprintf("%.1f Mb", \$2/1e6);
          unplaced_mb = sprintf("%.1f Mb", \$4/1e6);
          gap_kb      = sprintf("%.1f kb", \$5/1e3);
          printf "| Placed on chromosomes | %s | %s |\\n", \$1, placed_mb;
          printf "| Unplaced              | %s | %s |\\n", \$3, unplaced_mb;
          printf "| N-gaps introduced     | %s | %s |\\n", \$6, gap_kb
      }' "\$RAGTAG" 2>/dev/null || echo "| NA | NA | NA |"
      echo
      echo "> **Interpretation:** \$(scaffold_interp)"
      echo

      # ── 5. BUSCO ──────────────────────────────────────────────────────────
      echo "## 5. Gene space completeness (BUSCO — ${busco_lin})"
      echo
      echo '```'
      grep -E "C:|Complete|Fragmented|Missing|Total BUSCO" "\$BUSCO" 2>/dev/null | sed 's/^[[:space:]]*//' || echo "NA"
      echo '```'
      echo
      echo "> **Interpretation:** \$(busco_interp)"
      echo

      # ── Conclusion ────────────────────────────────────────────────────────
      echo "---"
      echo
      echo "## Conclusion"
      echo
      echo "The **${sample}** nuclear genome assembly is of **\${QUALITY_TIER} overall quality** based on BUSCO completeness (\${BUSCO_C}%), scaffold N50 (\${N50_SCAFFOLD} bp), and reference genome fraction (\${GFRAC_SCAFFOLD}%)."
      echo
      if [ "\$SKIP_PURGE" -gt 0 ]; then
        HAPLO_SENTENCE="**Haplotig removal:** Purge_dups was skipped (--skip_purge true) -- no haplotigs removed."
      else
        PURGE_MB=\$(echo \$LEN_PURGE | awk '{printf "%.0f Mb", \$1/1e6}')
        HAPLO_SENTENCE="**Haplotig removal:** Purge_dups removed \${PURGE_REMOVED}% of the assembly (autotune peak: \${PEAK}x), retaining \${PURGE_MB} of primary sequence."
      fi
      echo "**Sequencing:** \$COVERAGE ONT reads at mean quality Q\${MEAN_QUAL} (read N50 \${READ_N50} bp). **Assembly:** Flye produced \$NCONTIGS_FLYE contigs totalling \$(echo \$LEN_FLYE | awk '{printf "%.0f Mb", \$1/1e6}') from a \$(echo \$TOTAL_BASES | awk '{printf "%.1f Gb", \$1/1e9}') read set. **Polishing:** Medaka corrected base-level errors without structural changes. \${HAPLO_SENTENCE} **Scaffolding:** RagTag anchored \${PLACED_PCT}% of assembled bases to chromosomal positions using the reference, raising scaffold N50 to \$(echo \$N50_SCAFFOLD | awk '{printf "%.1f Mb", \$1/1e6}'). **Gene space:** BUSCO completeness \${BUSCO_C}% (${busco_lin}), with \${BUSCO_D}% duplication and \${BUSCO_M}% missing genes."
      echo
      awk -v c="\$BUSCO_C" -v gf="\$GFRAC_SCAFFOLD" -v ret="\$PURGE_RETAINED" -v dup="\$BUSCO_D" -v m="\$BUSCO_M" 'BEGIN {
          c=c+0; gf=gf+0; ret=ret+0; dup=dup+0; m=m+0
          rec=""
          if (ret < 55)
              rec=rec "- Purge_dups removed " sprintf("%.1f",100-ret) "% of the assembly. If the sample is largely homozygous, try --calcuts_args to set manual cutoffs, or consider skipping purge_dups.\\n"
          if (dup > 5)
              rec=rec "- BUSCO duplication (" dup "%) is elevated -- residual haplotigs may remain. Consider stricter purging parameters.\\n"
          if (gf < 60)
              rec=rec "- Reference genome fraction (" gf "%) is low. If BUSCO is high, this likely reflects genuine divergence from reference (not poor assembly quality).\\n"
          if (c < 90)
              rec=rec "- BUSCO completeness (" c "%) is below 90%. Consider deeper sequencing, additional polishing, or verifying the busco lineage param.\\n"
          if (m > 10)
              rec=rec "- " m "% of BUSCO genes are missing. Check for contamination, coverage gaps, or lineage mismatch.\\n"
          if (rec != "")
              printf "**Suggestions:**\\n\\n%s\\n", rec
      }'
      echo "_No contamination screening was performed. For publishable assemblies, consider running BlobTools or Kraken2 on the final scaffold FASTA._"
    } > ${sample}_assembly_summary.md
    """
}

process TOOLS_REPORT {
    label         'qc'
    errorStrategy 'ignore'
    publishDir    "${params.outdir}/reports", mode: 'copy'
    container     'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'

    output:
    path "pipeline_tools.md"

    script:
    """
    cat > pipeline_tools.md << 'TOOLSEOF'
# Pipeline tools -- nf-broom

Describes the role and rationale of each tool in the nf-broom ONT plant genome assembly pipeline.

---

## Read QC and filtering

### NanoPlot
Generates statistics and plots for the raw ONT read set: read length histogram, N50, total
base count, and per-read quality scores (Phred scale). Run before filtering so the raw
distribution is visible in MultiQC alongside filtered results.

### Filtlong
Filters reads by minimum length (>=1 kb) and minimum mean quality. Removing short, low-quality
reads reduces assembly fragmentation and lowers the chance of incorporating sequencing errors
into the draft assembly. The filtered reads are used for all downstream steps.

---

## Organelle read separation

### minimap2 (organelle alignment)
Aligns filtered reads to the combined chloroplast + mitochondria reference using the map-ont
preset (tuned for noisy long reads against a short reference). Reads that align are extracted
as cp/mt read sets; unmapped reads form the nuclear read set.

### samtools
Sorts and indexes the organelle alignment BAM so reads can be extracted by genomic region.

---

## Organelle assembly

### Flye (organelle mode)
De novo assembler for long reads. Assembled separately on deduplicated cp and mt read sets.
Produces FASTA + assembly graph (GFA format). Used when --organelle_assembler flye (default).

### OATK (--organelle_assembler oatk)
HMM-based organelle assembler. Identifies organelle reads directly from the full filtered
read set using pHMM profiles from the OatkDB (embryophyta plant gene database), then assembles
them. More sensitive than minimap2-based extraction for highly divergent species. Produces GFA
files for both chloroplast and mitochondria simultaneously.

### Bandage
Renders the GFA assembly graph as a PNG image. A complete circular chloroplast typically
appears as a single circular node or a two-bubble structure (reflecting the inverted repeat).
Useful for a quick sanity check on organelle topology.

### FILTER_ORGANELLE_CONTIGS
After assembly, aligns organelle contigs back to the reference and retains only those with
>=50% query coverage and >=70% identity. Removes nuclear-derived contigs (NUPTs/NUMTs) that
Flye or OATK occasionally incorporates into the organelle assembly. Controlled by
--filter_organelles (default: true), --organelle_min_qcov, and --organelle_min_ident.

---

## Nuclear assembly

### Flye (nuclear mode)
Same assembler as for organelles, but run on the nuclear read set with the full estimated
genome size (--genome_size). Uses a repeat graph approach to handle the high repeat content
of plant genomes.

### Medaka
Neural-network-based polisher trained on ONT signal data. Re-aligns the original reads to
the Flye draft and calls a consensus sequence to correct base-level errors. Improves raw
accuracy from approximately Q20 (Flye output) toward Q30+. The model must match the flowcell
chemistry and basecaller used; set via --medaka_model.

---

## Haplotig removal (--skip_purge to bypass)

### purge_dups
Removes redundant haplotig contigs from diploid assemblies. In a heterozygous organism, Flye
may assemble both haplotypes of a locus as separate contigs. purge_dups identifies redundant
contigs by two criteria:

1. Coverage depth: reads map to a haplotig at roughly half the expected haploid depth because
   reads from both haplotypes align to the primary contig, leaving the haplotig undercovered.
   Coverage cutoffs define depth zones: junk / haploid-kept / haplotig-purged / repeat.

2. Self-alignment overlap: contigs that align redundantly against each other at >50% overlap
   are flagged as haplotigs regardless of depth class. This step drives most of the removal.

The pipeline uses autotune mode by default: the haploid coverage peak is detected from the
depth distribution (PB.stat) and thresholds are set automatically. Use --calcuts_args to
override, or --skip_purge to bypass this step entirely and pass Medaka output directly to
scaffolding.

### HapDup (--run_hapdup)
Phases the purged primary assembly into haplotype-resolved contigs using long-read signal.
Aligns reads to the purged assembly then separates them by haplotype using heterozygous variant
sites. Produces two haplotype FASTA files. More compute-intensive than purge_dups.

---

## Scaffolding (requires --nuclear_ref)

### RagTag
Orders and orients purged contigs into chromosome-scale scaffolds by alignment to a reference
genome. Two steps are run in sequence:

- ragtag.py correct (enabled by default via --ragtag_correct): breaks contigs at positions
  where the assembly disagrees with the reference, reducing misassemblies before joining.
- ragtag.py scaffold: orders corrected contigs into scaffolds, inserting N-gaps at joins.

The reference does not need to be from the same species -- a related genome provides useful
synteny information even at moderate sequence identity.

---

## Assembly QC

### QUAST
Computes standard assembly statistics (N50, total length, contig count, GC content) and,
when a reference is provided, alignment-based metrics: genome fraction covered,
misassemblies, and structural variants. Run simultaneously on all nuclear assembly stages
(Flye, Medaka, purge_dups/medaka_nopurge, scaffold) to track how each step affects quality.

### BUSCO
Searches the final scaffolded assembly for a lineage-specific set of near-universal
single-copy orthologs. Completeness categories:
- Complete single-copy (S): present exactly once -- ideal
- Complete duplicated (D): present in multiple copies -- may indicate residual haplotigs
- Fragmented (F): partially recovered -- may indicate assembly fragmentation
- Missing (M): absent -- may indicate genuinely missing sequence or lineage mismatch

Use --busco_lineage to select the appropriate database for your organism.

### MultiQC
Aggregates QC reports from NanoPlot, BUSCO, QUAST, and (optionally) Qualimap into a single
interactive HTML report. Enables comparison across samples when multiple samples are processed
in the same batch.

---

## Optional BAM-level QC

### Qualimap bamqc (--run_qualimap)
Generates per-base coverage statistics, GC bias plots, and insert-size distributions from
the read-to-assembly BAM. Run at both the purge_dups and scaffold stages to compare coverage
uniformity before and after scaffolding. Results are integrated into the MultiQC report
automatically via native Qualimap support.

### BlobTools (--run_blobtools)
Generates GC-vs-coverage blob plots from the read-to-assembly BAM and assembly FASTA.
Coverage-only mode (no taxonomy database required) visualises whether all contigs cluster at
the expected sequencing depth and GC content. Unexpected clusters may indicate contamination
or organelle sequence carry-over into the nuclear assembly. Results are published as PNG plots
and TSV tables in the output directory.

---

## References and links

Citations for each tool, in order of appearance, with project homepages.

- **NanoPlot / NanoPack** -- De Coster W, et al. NanoPack: visualizing and processing long-read sequencing data. Bioinformatics. 2018;34(15):2666-2669. <https://github.com/wdecoster/NanoPlot>
- **Filtlong** -- Wick RR. Filtlong (no associated publication). <https://github.com/rrwick/Filtlong>
- **minimap2** -- Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094-3100. <https://github.com/lh3/minimap2>
- **SAMtools** -- Danecek P, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008. <https://www.htslib.org>
- **Flye** -- Kolmogorov M, et al. Assembly of long, error-prone reads using repeat graphs. Nat Biotechnol. 2019;37:540-546. <https://github.com/fenderglass/Flye>
- **OATK** -- Zhou C. Oatk: organelle assembly toolkit (no associated publication). <https://github.com/c-zhou/oatk>
- **Bandage** -- Wick RR, et al. Bandage: interactive visualization of de novo genome assemblies. Bioinformatics. 2015;31(20):3350-3352. <https://rrwick.github.io/Bandage>
- **Medaka** -- Oxford Nanopore Technologies. Medaka (no associated publication). <https://github.com/nanoporetech/medaka>
- **purge_dups** -- Guan D, et al. Identifying and removing haplotypic duplication in primary genome assemblies. Bioinformatics. 2020;36(9):2896-2898. <https://github.com/dfguan/purge_dups>
- **HapDup** -- Kolmogorov M, et al. Scalable nanopore sequencing of human genomes provides a comprehensive view of haplotype-resolved variation and methylation. Nat Methods. 2023;20:1483-1492. <https://github.com/KolmogorovLab/hapdup>
- **RagTag** -- Alonge M, et al. Automated assembly scaffolding using RagTag elevates a new tomato system for high-throughput genome editing. Genome Biol. 2022;23:258. <https://github.com/malonge/RagTag>
- **QUAST** -- Gurevich A, et al. QUAST: quality assessment tool for genome assemblies. Bioinformatics. 2013;29(8):1072-1075. <https://github.com/ablab/quast>
- **BUSCO** -- Manni M, et al. BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage. Mol Biol Evol. 2021;38(10):4647-4654. <https://busco.ezlab.org>
- **MultiQC** -- Ewels P, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047-3048. <https://multiqc.info>
- **Qualimap** -- Okonechnikov K, et al. Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data. Bioinformatics. 2016;32(2):292-294. <http://qualimap.conf.es>
- **BlobTools** -- Laetsch DR, Blaxter ML. BlobTools: Interrogation of genome assemblies. F1000Research. 2017;6:1287. <https://github.com/DRL/blobtools>
TOOLSEOF
    """
}
