# New features to implement

# Filter & Map to Organelles

- Running `filtlong` before filtering chloroplast (cp) and mitochondrial (mt) reads wastes CPU cycles. Swap the order, map raw reads directly to your organelle references using `minimap2` first, partition them, and then run `filtlong` exclusively on the isolated nuclear reads.

- When running `minimap2` to catch organelle sequences, with `-ax map-ont` preset, maybe we should append `-K 20M` to process blocks in parallel and lower runtime memory spikes.

- consider adding options for using `Hifiasm` (in its -nt mode for Nanopore) and/or `Verkko` instead of `Flye` for the nuclear compartment. They handle structural haplotypes natively, saving you a massive amount of cleanup time during the purge_dups step.

# Chromosome Scaffolding

- insert `RagTag` directly after the duplicate purging/phasing step.

* Misassembly Correction (`ragtag.py correction`): Uses the reference to identify and fix chimeric contigs in your Flye output.

* Scaffolding (`ragtag.py scaffold`): Orders and orients your contigs into 10 distinct chromosome-scale groups based on the reference coordinates.

* Patching: It inserts a standard string of Ns (usually 100 bp) between contigs to bridge them into continuous chromosome sequences.

# Human-Readable Reports & Stats Generation

- Introduce `QUAST` (Quality Assessment Tool for Genome Assemblies) into your pipeline right after chromosome assignment. You can supply your reference sorghum genome, and `QUAST` will output beautiful `HTML`, `PDF`, and `Markdown` tracking tables showing structural accuracy, misassemblies, and genomic synteny.

# Compiling Custom Final Reports via Nextflow

You can use Nextflow's internal `.collect()` operator to pull metadata text blocks from `BUSCO`, `QUAST`, and `NanoPlot`, and pipe them directly into a neat, human-readable markdown/text summary inside a final process.


# Linear Pipeline Flow
[ Raw ONT Reads ] ──► [ Map to Organelles ] ──► [ Partition Reads ]
                                                       │
                                        ┌──────────────┴──────────────┐
                                        ▼                             ▼
                                 [ Organelle Sets ]            [ Nuclear Reads ]
                                        │                             │
                                        ▼                             ▼ 
                              (Flye / oatk Assembly)             (Filtlong QC)
                                        |                             │
                                        ▼                             ▼
                                (QUAST oragenelle)              (Nuclear Flye)
                                        |                             │
                                        ▼                             ▼
                                (Bandage organelle)            (Medaka Polish)
                                                                      │
                                                                      ▼
                                                            (Purge Dups / HapDup)
                                                                      │
                                                                      ▼
                                                           [ Sorghum Reference ]
                                                                      │
                                                                      ▼
                                                             (RagTag Scaffolder)
                                                                      │
                                                                      ▼
                                                        (QUAST & BUSCO QC Evaluators)
                                                                      │
                                                                      ▼
                                                          (Final Summary & MultiQC)


