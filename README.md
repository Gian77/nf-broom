# nf-broom 

<!-- badges: start -->
[![Visits](https://img.shields.io/badge/visits-count-blue)](https://github.com/OWNER/REPO)
[![Latest Release](https://img.shields.io/github/v/release/OWNER/REPO)](https://github.com/OWNER/REPO/releases)
[![Total Downloads](https://img.shields.io/github/downloads/OWNER/REPO/total)](https://github.com/OWNER/REPO/releases)
[![Open Issues](https://img.shields.io/github/issues/OWNER/REPO)](https://github.com/OWNER/REPO/issues)
<!-- badges: end -->

## a nextflow pipeline for plnt genome assembly

`nf-broom` is a [nextflow](https://www.nextflow.io/) pipeline for plnt genome assembly. 
This pipeline is currently under development. At the moment this is ONT-only de novo 
assembly developed for sorghum genomes with explicit organelle separation. The name 
**nf-broom** comes from an older common name for sorhgum i.e. broom, as it was used to 
make brooms.

## What nf-broom does (for now)

For each sample:

1. **QC** — NanoPlot read-level QC
2. **Filter** — filtlong (min length 1 kb, min mean quality 80)
3. **Map to organelles** — minimap2 against combined cp + mt reference
4. **Partition reads** — split into chloroplast / mitochondrial / nuclear sets
5. **Assemble organelles** — Flye per compartment (default), or Oatk HMM-based read
   detection on all filtered reads (`--organelle_assembler oatk`)
6. **Filter organelle contigs** — reference-based filter removes nuclear-inserted
   organelle sequence (NUPTs/NUMTs)
7. **Assemble nuclear genome** — Flye
8. **Polish nuclear** — Medaka
9. **Purge duplicates** — purge_dups always runs; both the purged and unpurged
   (Medaka) genomes are carried forward for comparison
10. **Phase (optional)** — HapDup for diploid output (`--run_hapdup`, off by default)
11. **Scaffold (optional)** — RagTag correct + scaffold against `--nuclear_ref`, run on
    both candidate genomes
12. **Compare & select final** — QUAST + BUSCO score both candidates side-by-side (this
    exposes purge_dups over-purging); `--final_assembly medaka|purge` (default `medaka`)
    picks which one is published as the final genome
13. **Contaminant screening (optional)** — Kraken2 classification + BlobTools
    coverage/taxonomy blob plots on the final genome (`--run_kraken2`, `--run_blobtools`,
    `--run_qualimap`)
14. **Reports** — MultiQC aggregation, a human-readable per-sample assembly summary
    (`FINAL_SUMMARY`, requires `--nuclear_ref`), and a tool-citations report

## Directory structure expected

```
sorghum-assembly/
├── main.nf
├── nextflow.config
├── modules/
│   ├── qc.nf
│   ├── mapping.nf
│   ├── assembly.nf
│   └── polishing.nf
└── reads/                       <-- you provide
    ├── sample_A/
    │   └── *.fastq.gz
    ├── sample_B/
    │   └── *.fastq.gz
    └── ...
```

## Running

### Local test (small data, validate the wiring)

This may fail depending on how big your subsetted data is. For example, if too small it is hard to assemble the genome and all the coverage parameters will be off.

```bash
nextflow run main.nf \
    --reads /path/to/reads \
    --cp_ref /path/to/sorghum_chloroplast.fasta \
    --mt_ref /path/to/sorghum_mitochondrion.fasta \
    --outdir results_local
```

### HTCondor (production)

```bash
nextflow run main.nf \
    -profile condor \
    --reads /path/to/reads \
    --cp_ref /path/to/sorghum_chloroplast.fasta \
    --mt_ref /path/to/sorghum_mitochondrion.fasta \
    --outdir results \
    -w /scratch/$USER/nf-work-sorghum \
    -resume
```

### Enable HapDup phasing

Add `--run_hapdup` to either command above.

## Key parameters

| Parameter               | Default                                | Notes                                                          |
|------------------------|-----------------------------------------|-----------------------------------------------------------------|
| `--reads`               | `./reads`                              | Dir containing `<sample_id>/*.fastq.gz`                        |
| `--cp_ref`              | (required)                             | Chloroplast reference FASTA                                    |
| `--mt_ref`              | (required)                             | Mitochondrion reference FASTA                                  |
| `--nuclear_ref`         | none                                    | Nuclear reference FASTA; enables RagTag scaffolding, QUAST genome-fraction, and `FINAL_SUMMARY` |
| `--organelle_assembler` | `flye`                                 | `flye` (per-compartment) or `oatk` (HMM-based read detection)  |
| `--genome_size`         | `720m`                                 | Estimated nuclear genome size for Flye                         |
| `--medaka_model`        | `r1041_e82_400bps_sup_v5.0.0`          | Match your basecaller + chemistry                               |
| `--busco_lineage`       | `poales_odb10`                         | Plant lineage; downloaded auto by BUSCO                        |
| `--run_hapdup`          | `false`                                | Enable for diploid phasing                                     |
| `--calcuts_args`        | `""` (autotune)                        | Manual purge_dups cutoffs, e.g. `"-l 5 -m 22 -u 120"`           |
| `--final_assembly`      | `medaka`                               | `medaka` (unpurged) or `purge` — selects the published final nuclear genome; purge_dups always runs and both are compared |
| `--run_qualimap`        | `false`                                | BAM-level coverage QC on the final genome                      |
| `--run_blobtools`       | `false`                                | Coverage-vs-GC blob plot on the final genome                   |
| `--run_kraken2`         | `false`                                | Taxonomic contaminant screening on the final genome (needs `--kraken2_db` / `--taxdump_dir`) |
| `--kraken2_db`          | PlusPF DB path                         | Kraken2 database; small enough to load fully into RAM          |
| `--taxdump_dir`         | PlusPF DB path                         | NCBI taxdump (nodes.dmp/names.dmp) BlobTools uses to resolve Kraken2 taxids |
| `--outdir`              | `results`                              | Output directory                                                |

> `--skip_purge` is retired — purge_dups always runs now; use `--final_assembly medaka` (equivalent to the old skip behavior) or `--final_assembly purge`.

## Resource classes (configured in nextflow.config)

| Label             | CPUs | RAM     | Time   | Used by                                         |
|-------------------|------|---------|--------|--------------------------------------------------|
| `qc`              | 16   | 128 GB  | 12h    | NanoPlot, filtlong, MultiQC                      |
| `qc_heavy`        | 16   | 128 GB  | 24h    | BUSCO, QUAST, Qualimap, BlobTools, Kraken2       |
| `map`             | 24   | 128 GB  | 12h    | minimap2 + samtools                              |
| `assemble_small`  | 16   | 128 GB  | 6h     | Flye on cp / mt                                  |
| `assemble_heavy`  | 32   | 448 GB  | 120h   | Flye nuclear, HapDup                             |
| `polish`          | 24   | 256 GB  | 48h    | Medaka, purge_dups                               |

`KRAKEN2_CLASSIFY` overrides `qc_heavy` to 160 GB (`withName` in nextflow.config) so the
104 GB PlusPF database loads fully into RAM instead of relying on `--memory-mapping`.

Adjust to your cluster's queue limits and node capacity.

## HPC notes (GLBRC / HTCondor)

- **Run Nextflow from the submit node** — it dispatches jobs, doesn't compute. Use `tmux` or `screen` to keep it alive.
- **`work/` directory MUST be on scratch** — assembly intermediates are 100s of GB. Use `-w /scratch/$USER/nf-work-sorghum`.
- **Apptainer install required** — install via `conda install -c conda-forge apptainer -y` in the Nextflow env.
- **First run downloads ~6 container images** — takes 10–20 min, cached in `~/.apptainer_cache`.

## What's NOT yet included (Phase 3+)

- Merqury QV estimation
- Organelle annotation (GeSeq / PGA)
- Illumina polishing (Pilon)
- Variant calling against reference

## Sanity checks before running

```bash
# Verify config parses
nextflow config -profile condor

# Show the DAG without running anything
nextflow run main.nf -preview --cp_ref ... --mt_ref ...

# Test on one sample first
nextflow run main.nf -profile condor --reads ./reads_single_sample ...
```

## Use the `test_channel.nf` to test
```bash
nextflow run test_channel.nf --reads reads/
```

## Test run
# Pick one sample, subsample to ~100,000 reads (again, this may fail becasue to small of a subset)

```
mkdir -p reads_test/test_sample
seqkit sample -n 50000 reads/genotype_1/*.fastq.gz -o reads_test/test_sample/test.fastq.gz

# Run
nextflow run main.nf \
    --reads $PWD/reads_test \
    -profile condor \
    --cp_ref refs/sorghum/sorghum_cp_NC008602.fasta \
    --mt_ref refs/sorghum/sorghum_mt_NC008360.fasta \
    --outdir results_test \
    -w /scratch/$USER/nf-work-test
```

## For containers 
### 1) search the hash and then use the link to fidn the quay

For example for both samtools and minimap

```
curl -sO https://raw.githubusercontent.com/BioContainers/multi-package-containers/master/combinations/hash.tsv
grep -E "^minimap2=.*samtools=|^samtools=.*minimap2=" hash.tsv
```

Then paste the packages in:
https://midnighter.github.io/mulled

### 2) Use the quay_tools.sh 

For single-tool images, for example for `flye` just check `quay.io` directly at:
https://quay.io/repository/biocontainers/flye?tab=tags

Or you can find the tag directly, with an API call

```
conda activate base
curl -s 'https://quay.io/api/v1/repository/biocontainers/flye/tag/?limit=20&onlyActiveTags=true'  | \
python -m json.tool | grep '"name"' | head -10

# Then pick the most recent one matching the version you want, then build the image reference:
quay.io/biocontainers/flye:<paste-tag-here>

# That turns into:
container 'quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1'
```

#### Alternative option

Source the `quay_tools.sh` in `helper_functions/` to search throuhg the quay, verify an image exist and is pullable, and test a command for the tool you are looking into using. Remember, this works for images that have one single tool. If you need/want more tool in one image you need to rely on the mulled images, see above. 

To use the script, you need to have a conda nextflow environment with `apptainer` installed.

For example:
```
[benucci@scarcity-ap-1 ~]$ source helper-functions/quay_tools.sh 
[benucci@scarcity-ap-1 ~]$ conda activate nextflow
(nextflow) [benucci@scarcity-ap-1 ~]$ quay_tags flye
quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1
quay.io/biocontainers/flye:2.9.6--py312h734f728_1
quay.io/biocontainers/flye:2.9.6--py311h93bbee8_1
quay.io/biocontainers/flye:2.9.6--py310h5850263_1
quay.io/biocontainers/flye:2.9.6--py310h275bdba_0
quay.io/biocontainers/flye:2.9.6--py39h475c85d_0
quay.io/biocontainers/flye:2.9.6--py311h2de2dd3_0
quay.io/biocontainers/flye:2.9.5--py310h275bdba_2
quay.io/biocontainers/flye:2.9.5--py39h475c85d_2
quay.io/biocontainers/flye:2.9.5--py312h5e9d817_2
...

(nextflow) [benucci@scarcity-ap-1 ~]$ verify_image quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1
[OK pull] quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1  (no command check)
(nextflow) [benucci@scarcity-ap-1 ~]$ verify_image quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1 "flye --version"
[OK] quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1  (flye --version works)
(nextflow) [benucci@scarcity-ap-1 ~]$ 
```

# Test the pipeline

See the use of resume, to resume previous chashed data.

```
nextflow run main.nf \
-profile condor \
--reads $PWD/tests \
--cp_ref refs/sorghum/sorghum_cp_NC008602.fasta \
--mt_ref refs/sorghum/sorghum_mt_NC008360.fasta \
--outdir results_test \
-w nf-work-test \
-resume
```

Or simply:

```
nextflow run main.nf -profile condor --reads $PWD/B11077_test/ --cp_ref $PWD/refs/sorghum/sorghum_cp_NC008602.fasta --mt_ref $PWD/refs/sorghum/sorghum_mt_NC008360.fasta --outdir $PWD/results_B11077 -w $PWD/nf-work-test -resume
```

Using a specific nextflow log:

```
nextflow run main.nf -profile condor     --reads /home/glbrc.org/benucci/nf-broom/B11077_test/     --cp_ref /home/glbrc.org/benucci/nf-broom/refs/sorghum/sorghum_cp_NC008602.fasta     --mt_ref /home/glbrc.org/benucci/nf-broom/refs/sorghum/sorghum_mt_NC008360.fasta     --outdir /home/glbrc.org/benucci/nf-broom/results_B11077     -w /home/glbrc.org/benucci/nf-broom/nf-work-test  --run_hapdup true   -resume 22d67d52-0eb7-48c7-b602-567d87268e30
```


# Start interactive sessions in scarcity

## Method 1: tmux (recommended)

```
# Start a new session named "nf-broom"
tmux new -s nf-broom

# Inside tmux, launch the pipeline as usual
cd /home/glbrc.org/benucci/nf-broom
conda activate nextflow
nextflow run main.nf -profile condor \
    --reads $PWD/tests \
    --cp_ref refs/sorghum/sorghum_cp_NC008602.fasta \
    --mt_ref refs/sorghum/sorghum_mt_NC008360.fasta \
    --outdir results_test \
    -w nf-work-test \
    -resume
```

<div style="padding: 15px; border: 1px solid #007bcc; background-color: #f0f8ff; border-radius: 5px;"> 
    <strong>More about tmux use:</strong> To Detach from tmux (leaves it running), press <code>Ctrl-B</code>, then <code>D</code>. Now you can close your laptop, log out, whatever. If you want to copy output on the terminal you can press <code>Ctrl + b</code> then release. Press the <code>[</code> key (this enters "copy mode"). Use your <code>Up/Down</code> arrow keys or <code>Page Up/Page Down</code> to scroll through your output. Press <code>q</code> to exit scroll mode and return to typing. 
</div>

To reconnect later from anywhere:
```
ssh scarcity-ap-1.glbrc.org
tmux attach -t nf-broom
```

## Method 2: nohup
Simpler but less interactive — no live progress bars to look at:

```
nohup nextflow run main.nf -profile condor \
    --reads $PWD/tests \
    --cp_ref refs/sorghum/sorghum_cp_NC008602.fasta \
    --mt_ref refs/sorghum/sorghum_mt_NC008360.fasta \
    --outdir results_test \
    -w nf-work-test \
    -resume \
    > nf.log 2>&1 &

# Note the PID
echo $! > nf.pid

# Check progress later:
tail -f nf.log
ps -p $(cat nf.pid)
kill $(cat nf.pid)
```

# To clean up an start a complete new session

```
cd /home/glbrc.org/benucci/nf-broom

# The work directory (cached task outputs — this is the big one)
rm -rf /home/glbrc.org/benucci/nf-broom/nf-work-test/

# The .nextflow hidden directory (history, cache metadata, session info)
rm -rf .nextflow/

# The published results from previous runs (you backed this up above)
rm -rf results_B11077/

# Any leftover log files
rm -f .nextflow.log* nextflow_report*.html timeline*.html trace*.txt
```

And to cancel a `tmux` session, after closing it

```
tmux kill-session -t nf-broom
tmux ls    # should now say "no server running"
```


# Additional cleanups

## Check and clean the apptainer cache
```
apptainer cache list -v
apptainer cache clean
```




