# scOUTseq Pipeline

Portable workspace version of the `scOUTseq` analysis pipeline for processing 10X and PARSE-based scOUTseq runs.

This repository contains:

- a cleaned-up shell pipeline entrypoint
- helper scripts grouped by function
- an example config
- example FASTQ inputs
- a portable genome directory convention

## Repository Layout

```text
.
├── .gitignore
├── README.md
├── example.config
├── hdr_filter_example.ini
├── requirements.txt
├── data/
│   ├── scOUT_seq_example_R1.fastq.gz
│   └── scOUT_seq_example_R2.fastq.gz
├── genome/
│   └── <reference fasta and bwa index files>
└── scripts/
    ├── required-core/
    ├── optional-parse/
    └── optional-hdr-downstream/
```

## What The Pipeline Does

At a high level, the pipeline:

1. loads a run config
2. links the input FASTQs into a sample output directory
3. optionally filters reads by sequence patterns
4. runs `umi_tools whitelist` and `umi_tools extract`
5. optionally annotates PARSE barcodes
6. optionally separates HDR-barcode-containing reads
7. runs `CRISPResso`
8. optionally performs downstream HDR/off-target/translocation/editing outcome analysis

## Entry Point

The main entrypoint is:

- [scripts/required-core/scOUT_pipeline.sh](/Users/moritzschlapansky/Documents/Playground/scripts/required-core/scOUT_pipeline.sh)

Run it from the repository root with:

```bash
bash ./scripts/required-core/scOUT_pipeline.sh ./example.config
```

## Requirements

The pipeline expects these tools to be available in your environment:

- `bash`
- `conda`
- `python`
- `seqtk` `1.3`
- `pigz` `2.6`
- `umi_tools` `1.1.4`
- `CRISPResso` `2.3.1`
- `bwa` `0.7.19`
- `samtools` `1.6`
- `htslib` `1.6`

Some helper scripts also depend on Python packages such as:

- `pandas` `2.2.3`
- `numpy` `2.1.3`
- `matplotlib` `3.10.0`
- `seaborn` `0.13.2`
- `openpyxl` `3.1.0`
- `pyfastx` `2.1.0`

Install the Python dependencies with:

```bash
pip install -r requirements.txt
```

The repository-level [requirements.txt](/Users/moritzschlapansky/Documents/Playground/requirements.txt) also documents the required command-line tools as comments.

## Conda Environments

The pipeline attempts to activate:

- `scout`
- `crispresso2_env`

If either environment is missing, the script continues with the current shell environment. In practice, the required tools still need to be available on `PATH`.

## Required Input Files

For a basic run, you need:

- paired FASTQ files
- a config file
- a reference genome FASTA in `./genome/`
- BWA index files for that FASTA

Example minimal input layout:

```text
.
├── example.config
├── data/
│   ├── sample_R1.fastq.gz
│   └── sample_R2.fastq.gz
└── genome/
    ├── hg38.fa
    ├── hg38.fa.amb
    ├── hg38.fa.ann
    ├── hg38.fa.bwt
    ├── hg38.fa.pac
    └── hg38.fa.sa
```

If you are using PARSE mode or HDR downstream analysis, more files are needed as described below.

## Genome Files

The pipeline now assumes genome files live under:

- [genome/](/Users/moritzschlapansky/Documents/Playground/genome)

The helper scripts look for common filenames such as:

- `hg38.fa`
- `hg38.fasta`
- `human.fa`
- `genome.fa`
- `mm10.fa`
- `mouse.fa`
- `AAV*.fa`

For `bwa mem`, the reference must already be indexed:

```bash
bwa index ./genome/hg38.fa
samtools faidx ./genome/hg38.fa
```

## Configuration

An example config is provided at:

- [example.config](/Users/moritzschlapansky/Documents/Playground/example.config)

Optional HDR barcode settings live in:

- [hdr_filter_example.ini](/Users/moritzschlapansky/Documents/Playground/hdr_filter_example.ini)

Paths in the config may be:

- absolute paths
- paths relative to the config file

Important config variables:

- `R1`, `R2`: input FASTQ files
- `SAMPLENAME`: output folder name
- `LIBTYPE`: `10X` or `PARSE`
- `EXPECTEDCELLS` or `SETCELLNUMBER`: whitelist behavior
- `BCPATTERN`: barcode pattern for 10X mode
- `TARGET`, `AMPLICONSEQ`, `GUIDE`, `CRISPRESSOWINDOW`: CRISPResso settings
- `SEQFILTER`: optional pre-filter pattern file
- `FILTERHDRREADSCONFIG`: enables HDR barcode extraction when set
- `CBCPATH`: metadata directory for PARSE or downstream cell barcode annotation
- `GENOME_DIR`: genome directory, defaults to `./genome`

## Script Groups

The scripts are separated into three groups.

### `required-core`

Core scripts needed for routine runs:

- pipeline entrypoint
- generic sequence filtering
- GAPDH pseudogene filtering
- genome path resolution
- translocation + editing outcome integration

See:

- [scripts/required-core/](/Users/moritzschlapansky/Documents/Playground/scripts/required-core)

### `optional-parse`

Only needed for `LIBTYPE="PARSE"`:

- barcode-to-ID annotation
- PARSE barcode lookup tables
- PARSE translocation extraction
- PARSE editing outcome extraction

See:

- [scripts/optional-parse/](/Users/moritzschlapansky/Documents/Playground/scripts/optional-parse)

### `optional-hdr-downstream`

Only needed when HDR barcode filtering and downstream remapping are enabled:

- HDR read extraction
- CRISPResso FASTQ-to-table conversion
- off-target remapping
- translocation remapping
- 10X downstream extractors

See:

- [scripts/optional-hdr-downstream/](/Users/moritzschlapansky/Documents/Playground/scripts/optional-hdr-downstream)

## Typical Run Modes

### 1. Minimal 10X smoke test

Use this when you only want to verify the core pipeline runs on test data.

Recommended config changes:

- `LIBTYPE="10X"`
- `FILTERHDRREADSCONFIG=""`
- `CBCPATH=""`
- keep or disable `SEQFILTER` depending on whether you want that pre-filter step

Run:

```bash
bash ./scripts/required-core/scOUT_pipeline.sh ./example.config
```

### 2. Full 10X run with HDR downstream analysis

Use:

- `LIBTYPE="10X"`
- `FILTERHDRREADSCONFIG="./hdr_filter_example.ini"` or another real HDR config
- valid `CBCPATH`
- valid genome FASTA and BWA index files

This mode uses scripts from:

- `required-core`
- `optional-hdr-downstream`

### 3. Full PARSE run

Use:

- `LIBTYPE="PARSE"`
- valid `CBCPATH`
- optional `PARSEKIT="mini"` when appropriate

This mode uses scripts from:

- `required-core`
- `optional-parse`
- `optional-hdr-downstream` if HDR downstream analysis is enabled

## Logging

Each run now writes a per-run log file automatically.

Log files are created under:

```text
<SAMPLENAME>/logs/
```

Filename format:

```text
<SAMPLENAME>_YYYYMMDD_HHMMSS.log
```

Pipeline log lines are timestamped, for example:

```text
2026-03-18 14:22:31 [scOUT] Activating conda environment 'scout'.
```

Tool output from programs like `bwa` and `CRISPResso` is also captured in the same log file.

Additional tool-specific logs are written by some steps, for example:

- `fastq_corrected_bc/<sample>_whitelist.log`
- `fastq_corrected_bc/<sample>.extracted.log`
- `editing_outcomes/EditingOutcomeAssignment.log`

## Outputs

The pipeline creates a sample output directory in the repository root:

```text
<SAMPLENAME>/
```

Common subdirectories include:

- `fastq/`
- `fastq_corrected_bc/`
- `<SAMPLENAME>_crispresso_out/`
- `logs/`

Optional downstream analysis may also create:

- translocation tables
- off-target tables
- `editing_outcomes/`

Within `editing_outcomes/`, the pipeline can also generate:

- `EditingOutcomeFrequencies.csv`
- `EditingOutcomesWithTranslocations.csv`
- `FilteredEditingOutcomesWithTranslocations.csv`

`EditingOutcomesWithTranslocations.csv` is the merged table of top editing outcomes plus translocation summaries.

`FilteredEditingOutcomesWithTranslocations.csv` is the downstream reporting table. It also includes:

- `GroupedRepairOutcome_1`
- `RepairClass_1`
- `GroupedRepairOutcome_2`
- `RepairClass_2`

`EditingOutcomeFrequencies.csv` is a tidy frequency table with these outcome families:

- `RepairOutcome`
- `GroupedRepairOutcome`
- `RepairClass`

## Example Test Run

With the bundled example FASTQs:

```bash
bash ./scripts/required-core/scOUT_pipeline.sh ./example.config
```

Before running, make sure:

- the FASTQs exist in [data/](/Users/moritzschlapansky/Documents/Playground/data)
- a reference FASTA exists in [genome/](/Users/moritzschlapansky/Documents/Playground/genome)
- the reference has BWA index files
- the required tools are installed

If you only want a lightweight test, disable HDR downstream analysis in the config by setting:

```bash
FILTERHDRREADSCONFIG=""
CBCPATH=""
```

## Common Errors

### `fail to locate the index files`

This means `bwa` found the FASTA but not the corresponding BWA index files.

Fix by indexing the FASTA:

```bash
bwa index ./genome/hg38.fa
```

### Missing conda environment

The pipeline prints a message and continues, but the required binaries still need to be available in the current environment.

### Missing metadata directory

If `CBCPATH` is set for PARSE or downstream annotation, the referenced directory must exist and contain the expected metadata CSV files.

## Notes

- The pipeline does not bundle genome FASTA files.
- The example FASTQs are included for testing, but not all downstream metadata required for a full production-style run is bundled here.
- Some downstream helper scripts are still domain-specific and expect particular input table formats produced by upstream steps.
- Runtime outputs are ignored by git via [.gitignore](/Users/moritzschlapansky/Documents/Playground/.gitignore).

## License

This project is MIT licensed.
