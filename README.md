# scOUTseq Pipeline

Analysis pipeline for processing 10X and Parse-based scOUT-seq data to map Cas9 editing outcomes at single cell and single nucleotide resolution.

## What The Pipeline Does

At a high level, the pipeline:

1. Salvages and assigns cell barcodes and UMIs.
2. Optionally filters reads by sequence patterns.
3. Optionally assigns HDR-barcode-containing reads.
4. Uses CRISPResso to align reads against target to assign editing outcomes and categorizes them at several levels of granularity.
5. Optionally performs downstream HDR/large deletion/translocation/editing outcome analysis.


## Entry Point

The main entrypoint is:

- [scripts/required-core/scOUT_pipeline.sh]

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

The requirements text file [requirements.txt] also documents the required command-line tools as comments.

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

The pipeline now assumes genome files (e.g. hg38 or mm10) live under:

- [`genome/`](./genome/)

The genome files are not included in the pipeline, so you will have to download and index them yourself.
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

When using AAV for editing, the pipeline can also check against the AAV genome (e.g. to check for AAV integrations). Therefore, include also your indexed AAV genome.

## Configuration

An example config is provided at:

- [example.config]

Optional HDR barcode settings live in:

- [hdr_filter_example.ini]

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
- `CBCPATH`: metadata directory for Parse or downstream cell barcode annotation (to match with transcriptome data)
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

- [scripts/required-core/]

### `optional-parse`

Only needed for `LIBTYPE="PARSE"`:

- barcode-to-ID annotation
- PARSE barcode lookup tables
- PARSE translocation extraction
- PARSE editing outcome extraction

See:

- [scripts/optional-parse/]

### `optional-hdr-downstream`

Only needed when HDR barcode filtering and downstream remapping are enabled:

- HDR read extraction
- CRISPResso FASTQ-to-table conversion
- off-target remapping
- translocation remapping
- 10X downstream extractors

See:

- [scripts/optional-hdr-downstream/]

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

### 3. Full Parse run

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

Within `editing_outcomes/`, the pipeline can also generate among other file:

- `EditingOutcomeFrequencies.csv`
- `FilteredEditingOutcomesWithTranslocations.csv`

I included an example of the different editing categories as a picture under [outcome_categories.pdf].

## Example Test Run

With the bundled example FASTQs:

```bash
bash ./scripts/required-core/scOUT_pipeline.sh ./example.config
```

Before running, make sure:

- the FASTQs exist in [data/]
- a reference FASTA exists in [genome/]
- the reference has BWA index files
- the required tools are installed

If you only want a lightweight test, disable HDR downstream analysis in the config by setting:

```bash
FILTERHDRREADSCONFIG=""
CBCPATH=""
```

## Notes

- The pipeline does not bundle genome FASTA files.
- The example FASTQs are included for testing, but not all downstream metadata required for a full production-style run is bundled here.
- Some downstream helper scripts are still domain-specific and expect particular input table formats produced by upstream steps.

## License

This project is MIT licensed.
