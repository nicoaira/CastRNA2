# CastRNA

A Nextflow pipeline for RNA secondary structure inference by projecting Rfam consensus structures onto query sequences.

## Overview

CastRNA derives RNA secondary structures by projecting pre-computed consensus structures from the [Rfam database](https://rfam.org/) onto query sequences via covariance model alignment. Rather than predicting structures from scratch using thermodynamic models, CastRNA leverages the evolutionary information encoded in Rfam family alignments to produce biologically meaningful structure annotations.

**Key concept:** Given RNA sequences, the pipeline identifies which Rfam families they belong to (or uses provided annotations) and "casts" the consensus secondary structures from those families onto the query sequences through sequence alignment.

## Features

- **Multiple input modes**: Pure FASTA (inference mode), FASTA with mapping file, or structured CSV samplesheet
- **Automatic Rfam database management**: Downloads and indexes Rfam covariance models automatically
- **Dual alignment methods**:
  - `cmalign`: Computes new alignments using Infernal (default)
  - `seed`: Uses pre-computed Rfam seed alignments for reproducibility with bpRNA
- **Flexible handling of edge cases**:
  - Non-canonical base pair strategies (keep, remove, or flag)
  - Pseudoknot preservation (letter notation)
  - Configurable inference thresholds
- **MSA output**: Optionally generates multiple sequence alignments grouped by Rfam family
- **Container support**: Docker, Singularity, and Conda environments
- **Comprehensive reporting**: JSON reports, TSV summaries, and HTML execution visualizations

## Requirements

### Software

- [Nextflow](https://www.nextflow.io/) >= 23.04.0
- One of the following execution environments:
  - [Docker](https://www.docker.com/)
  - [Singularity](https://sylabs.io/singularity/)
  - [Conda](https://docs.conda.io/)
  - Local installation of [Infernal](http://eddylab.org/infernal/) 1.1.5

### Hardware

Default resource limits (configurable):
- Memory: up to 128 GB
- CPUs: up to 16 cores
- Time: up to 240 hours

## Installation

1. **Install Nextflow**:
   ```bash
   curl -s https://get.nextflow.io | bash
   chmod +x nextflow
   # Optionally move to a directory in your PATH
   mv nextflow ~/bin/
   ```

2. **Clone the repository**:
   ```bash
   git clone https://github.com/your-username/CastRNA.git
   cd CastRNA
   ```

3. **Choose an execution profile** (see [Usage](#usage) below)

## Usage

### Basic Examples

**Mode 1: Pure FASTA (inference mode)**

All sequences are scanned against the Rfam database to identify their family:
```bash
nextflow run main.nf --input sequences.fasta -profile docker
```

**Mode 2: FASTA with known Rfam IDs**

Provide a mapping file to skip inference for sequences with known families:
```bash
nextflow run main.nf --input sequences.fasta --mapping mapping.csv -profile docker
```

**Mode 3: CSV samplesheet**

Use a structured samplesheet for complex inputs:
```bash
nextflow run main.nf --input samplesheet.csv -profile docker
```

**Mode 4: Reproduce bpRNA methodology**

Use seed alignments for exact reproducibility:
```bash
nextflow run main.nf --input sequences.fasta --alignment_method seed --rfam_release 12.2 -profile docker
```

**Mode 5: Generate MSA outputs**

Group sequences by family and generate multiple sequence alignments:
```bash
nextflow run main.nf --input sequences.fasta --mapping mapping.csv --output_msa -profile docker
```

### Execution Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers |
| `singularity` | Run with Singularity containers |
| `conda` | Run with Conda environments |
| `test` | Run with minimal test dataset |
| `debug` | Enable debug mode with verbose output |

### Running the Test Dataset

```bash
nextflow run main.nf -profile test,docker
```

## Input Formats

### Option A: FASTA file

Standard FASTA format with one or more sequences:
```fasta
>sequence_1
ACGUACGUACGUACGUACGU
>sequence_2
GCUAGCUAGCUAGCUAGCUA
```

### Option B: FASTA + Mapping CSV

A CSV file mapping sequence IDs to Rfam families. Sequences without an `rfam_id` will be inferred via cmscan:
```csv
seq_id,rfam_id
sequence_1,RF00067
sequence_2,RF00655
sequence_3,
```

### Option C: CSV Samplesheet

A CSV file with paths to individual FASTA files:
```csv
seq_id,fasta,rfam_id
seq1,/path/to/seq1.fasta,RF00001
seq2,/path/to/seq2.fasta,
seq3,/path/to/seq3.fasta,RF00003
```

## Output

Results are written to the `results/` directory (configurable with `--outdir`):

```
results/
├── structures/
│   ├── *.dbn              # Dot-bracket notation files
│   └── *.json             # Detailed projection reports
├── alignments/
│   └── *.sto              # Stockholm format alignments
├── msa/                   # (if --output_msa)
│   ├── {rfam_id}.msa.sto  # Multi-sequence Stockholm
│   ├── {rfam_id}.aligned.fa
│   └── {rfam_id}.structures.tsv
├── summary.tsv            # Aggregated results
└── pipeline_info/
    ├── execution_report_*.html
    ├── execution_timeline_*.html
    ├── execution_trace_*.txt
    └── pipeline_dag_*.svg
```

### Output Formats

**Dot-bracket notation (`.dbn`)**:
```
>sequence_id
ACGUACGUACGU
((((....))))
```

**JSON report**: Contains detailed statistics including sequence length, number of base pairs, unpaired bases, insertions, and non-canonical pairs.

**Summary TSV**: Tab-separated file with columns:
- `seq_id`, `rfam_id`, `status`, `sequence_length`
- `num_base_pairs`, `num_unpaired`, `num_insertions`, `num_noncanonical`

## Parameters

### Input/Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | *required* | Path to FASTA file or CSV samplesheet |
| `--mapping` | - | Path to CSV mapping file (seq_id, rfam_id) |
| `--outdir` | `results` | Output directory |
| `--publish_dir_mode` | `copy` | How to publish results (copy, symlink, etc.) |

### Rfam Database

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--rfam_cm` | - | Path to Rfam.cm or Rfam.cm.gz (downloads if not provided) |
| `--rfam_seed` | - | Path to Rfam.seed (required for seed method) |
| `--rfam_clanin` | - | Path to Rfam.clanin (optional) |
| `--rfam_release` | `CURRENT` | Rfam version (e.g., `14.9`, `12.2`) |
| `--rfam_cache_dir` | `${launchDir}/rfam_db` | Cache directory for downloaded files |

### Alignment & Inference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--alignment_method` | `cmalign` | Alignment method: `cmalign` or `seed` |
| `--cmscan_evalue` | `1e-5` | E-value threshold for cmscan |
| `--cmscan_score` | `20.0` | Bit score threshold for cmscan |

### Structure Processing

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--noncanonical_pairs` | `keep` | How to handle non-canonical pairs: `keep`, `remove`, or `flag` |
| `--keep_all_bracket_types` | `false` | Preserve `<>`, `[]`, `{}` brackets (vs. converting to `()`) |
| `--output_msa` | `false` | Generate MSA outputs grouped by Rfam family |

### Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_memory` | `128.GB` | Maximum memory limit |
| `--max_cpus` | `16` | Maximum CPU limit |
| `--max_time` | `240.h` | Maximum time limit |

## Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│                         CastRNA Pipeline                         │
└─────────────────────────────────────────────────────────────────┘
                                │
                                ▼
                    ┌───────────────────────┐
                    │   1. PREPARE_RFAM     │
                    │  Download/index Rfam  │
                    └───────────────────────┘
                                │
                                ▼
                    ┌───────────────────────┐
                    │   2. Parse Input      │
                    │  FASTA/CSV/Samplesheet│
                    └───────────────────────┘
                                │
                                ▼
                    ┌───────────────────────┐
                    │   3. Branch Sequences │
                    └───────────────────────┘
                          │           │
                Known IDs │           │ Unknown IDs
                          ▼           ▼
                          │   ┌───────────────────┐
                          │   │ 4. INFERNAL_CMSCAN│
                          │   │  Identify families│
                          │   └───────────────────┘
                          │           │
                          └─────┬─────┘
                                │
                                ▼
                    ┌───────────────────────┐
                    │   5. Alignment        │
                    │  cmalign OR seed      │
                    └───────────────────────┘
                                │
                                ▼
                    ┌───────────────────────┐
                    │ 6. CONSENSUS_PROJECTION│
                    │  Project structures   │
                    └───────────────────────┘
                                │
                                ▼
                    ┌───────────────────────┐
                    │ 7. GENERATE_SUMMARY   │
                    │  Aggregate results    │
                    └───────────────────────┘
```

## Alignment Methods

### cmalign (default)

Uses Infernal's `cmalign` to compute new alignments between query sequences and their identified Rfam covariance models. This method:
- Works for any sequence that matches an Rfam family
- Computes optimal alignments using the full covariance model
- Is the recommended method for general use

### seed

Uses pre-computed seed alignments from Rfam. This method:
- Mimics the bpRNA database methodology
- Only works for sequences that are already in Rfam seed alignments
- Provides exact reproducibility with bpRNA annotations
- Requires the `--rfam_seed` file (downloaded automatically if not provided)

## Tips

1. **Large datasets**: For many sequences, consider increasing `--max_cpus` and using a cluster executor (SLURM, PBS, etc.)

2. **Specific Rfam version**: Use `--rfam_release` to ensure reproducibility across runs

3. **Known families**: If you already know the Rfam families for your sequences, provide a mapping file to skip the computationally expensive cmscan step

4. **Pseudoknots**: CastRNA preserves pseudoknot annotations using letter notation (A...a, B...b, etc.)

5. **Non-canonical pairs**: Use `--noncanonical_pairs remove` if you need strictly canonical (AU, GC, GU) base pairs only

## Citation

If you use CastRNA in your research, please cite:

- The Rfam database: [Rfam 14.0: updates to the RNA families database](https://doi.org/10.1093/nar/gkaa1047)
- Infernal: [Infernal 1.1: 100-fold faster RNA homology searches](https://doi.org/10.1093/bioinformatics/btt509)

## License

[Add license information here]

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.
