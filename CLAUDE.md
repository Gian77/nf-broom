# nf-broom — Claude Code guidance

## Project overview
ONT-based plant genome assembly pipeline (Nextflow DSL2, HTCondor cluster).
Assembles chloroplast, mitochondria, and nuclear genomes from long reads.

## Repository layout
- `main.nf` — top-level workflow
- `modules/` — one file per stage: `assembly.nf`, `dbs.nf`, `mapping.nf`, `polishing.nf`, `qc.nf`
- `nextflow.config` — all params and Condor label profiles
- `refs/sorghum/` — reference FASTAs: `sorghum_cp_NC008602.fasta`, `sorghum_mt_NC008360.fasta`, `Sbicolor_730_v5.0.fa`
- `reads/` — per-sample subdirectories of `.fastq.gz` files

## Container policy — always verify before changing an image tag

Before editing any `container` directive in a module, verify the image works on the cluster nodes using the helper script:

```bash
source ~/helper-functions/quay_tools.sh

# Fast API check (no download — confirms tag exists):
image_exists quay.io/biocontainers/flye:2.9.5--py39h6935b12_1

# Full pull + binary check (slow — confirms it actually runs):
verify_image quay.io/biocontainers/flye:2.9.4--py310h2b6aa90_0 "flye --version"
```

`verify_image` pulls into a temp SIF and discards it; set `VERIFY_VERBOSE=1` to see pull errors.
Use `quay_tags <tool>` to list available tags when looking for alternatives.

**Known broken image**: `quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1` — dies with
`SIGILL` (illegal CPU instruction) on the Condor compute nodes. Use `2.9.4--py310h2b6aa90_0` instead.

## Running the pipeline

Typical test run (HTCondor, with resume):
```bash
nextflow run main.nf -profile condor \
  --reads B11077_test/ \
  --cp_ref refs/sorghum/sorghum_cp_NC008602.fasta \
  --mt_ref refs/sorghum/sorghum_mt_NC008360.fasta \
  --outdir results_B11077 \
  -w nf-work-test \
  --organelle_assembler oatk \
  -resume
```

Always use `-resume` when re-running after a failure — Nextflow will reuse all cached tasks
and only re-run what changed or failed.

## Organelle assembly

Two assembler modes (set with `--organelle_assembler`):
- `flye` (default) — runs `ASSEMBLE_CP_FLYE` and `ASSEMBLE_MT_FLYE` on deduplicated organelle reads
- `oatk` — runs `ASSEMBLE_ORGANELLES_OATK` on all filtered reads using HMM-based read detection

When `--organelle_assembler oatk` AND `--cp_ref`/`--mt_ref` are provided, the pipeline
automatically builds custom HMM databases from those FASTAs (`BUILD_OATKDB` in `modules/dbs.nf`)
instead of downloading the generic embryophyta OatkDB. The sorghum-specific profiles give
better read classification than the broad plant database.

## Cluster notes
- Executor: HTCondor (`-profile condor`)
- Node `scarcity-14.glbrc.org` is excluded in config (known bad node)
- Node `scarcity-9.glbrc.org` has an older CPU: the medaka torch/ONNX binary
  SIGILLs (exit 132) there. It is excluded from the `polish` label only (it is
  still allowlisted for `assemble_heavy` nuclear assembly, which it handles fine).
  Retrying alone didn't help — Condor kept rematching the same node.
- Retry on SIGILL/OOM/timeout exit codes: 132, 137, 139, 143, 247
