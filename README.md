# bioinfo

Small collection of bioinformatics pipeline scripts and C++ utilities for variant-calling and gene-position analyses.

## Repository layout

- `Scripts/` — Bash workflows for project-specific alignment, variant calling, and post-processing tasks.
  - `Liu_prj_pipeline.sh` and `ayesha_prj_pipeline.sh` set project paths, references, FASTQ discovery patterns, read-group construction, and then source the shared GATK/Picard pipeline steps from `step2-21.sh`.
  - `step2-21.sh` contains the shared processing stages: SAM sorting, BAM merging, alignment metrics, duplicate marking, indexing, realignment, variant calling, filtering, recalibration, annotation, and coverage generation.
  - `tomato.sh` contains a tomato resequencing workflow using BWA/SAMtools/Picard/BCFtools.
- `gene_match/` — C++ command-line tools and headers.
  - `genematch.cpp` matches chromosome positions to gene intervals from a reference gene table.
  - `Perm-test.cpp` runs a simple permutation-test utility.
  - `alphanum_comp.h` provides alphanumeric sorting helpers used by `genematch.cpp`.

## Requirements

The scripts assume a Linux shell environment with common bioinformatics command-line tools installed and available on `PATH` or configured in the scripts:

- GNU Bash, GNU Find, GNU Parallel, AWK, and Java
- BWA
- SAMtools
- BCFtools
- BEDTools
- Picard
- GATK 3.x
- snpEff
- A C++ compiler such as `g++` for the utilities in `gene_match/`

> Note: several scripts contain absolute local paths for project directories, reference genomes, Picard, GATK, and snpEff. Update those variables before running the workflows on a new system.

## Building the C++ utilities

From the repository root:

```bash
g++ -std=c++11 -O2 -o gene_match/gene_match gene_match/genematch.cpp
g++ -std=c++11 -O2 -o gene_match/permtest gene_match/Perm-test.cpp
```

## Usage examples

### Match variants or positions to genes

Prepare a tab-delimited reference gene file with columns like:

```text
Gene stable ID	Chromosome	Gene start	Gene end	Gene name
ENSG00000283891	15	55372940	55373034	MIR628
```

Then run:

```bash
gene_match/gene_match -cf Gene_location.txt mutation.txt > out_matched.txt
```

Input position files should have chromosome and position in the first two tab-delimited columns, for example:

```text
chr17	51391264
```

### Run a project pipeline

Edit the variables at the top of the desired script, especially `PRJDIR`, `GE_REF`, `RESULT`, `LOGFILE`, and any tool paths. Then run the script from a configured environment:

```bash
bash Scripts/ayesha_prj_pipeline.sh
```

The project pipeline scripts use `STARTSTEP` and Bash `case` fall-through to resume from selected stages.

## License

This repository is licensed under the MIT License. See `LICENSE` for details.
