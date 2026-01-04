# rust-featurecounts

A fast, memory-efficient feature counting tool for **prokaryotic RNA-seq analysis**, written in Rust. This is a clean-room reimplementation inspired by [featureCounts](https://subread.sourceforge.net/featureCounts.html) from the Subread package.

> **Note**: This tool is designed specifically for bacterial and archaeal genomes. It has not been tested with eukaryotic genomes and may not handle complex features such as alternative splicing.

> **Vibe Coding**: This project was developed through AI-assisted programming (vibe coding) using Claude.

## Features

- **High Performance**: Multi-threaded processing with automatic load-based thread allocation
- **Low Memory Footprint**: Efficient memory usage through Rust's zero-cost abstractions
- **featureCounts Compatible**: Output format fully compatible with downstream tools (DESeq2, edgeR, etc.)
- **Multiple Annotation Formats**: Supports GFF, GTF, and GenBank (.gbk, .gb, .genbank, .gbff) formats
- **Strand-Specific Counting**: Unstranded, forward, and reverse strand options
- **Flexible Counting Modes**: Intersection and union overlap methods
- **Paired-End Support**: Fragment-based counting for paired-end data
- **Quality Filtering**: MAPQ-based read filtering

## Performance

### Benchmarks

Tested on *E. coli* K-12 MG1655 RNA-seq data (913,011 fragments, 4,506 genes):

| Metric | rust-featurecounts | Notes |
|--------|-------------------|-------|
| **Processing Speed** | ~300ms (4 threads) | Region-based parallel processing |
| **Processing Speed** | ~140ms (30 threads) | Auto-scaling with `-T 0` |
| **Throughput** | ~3 million fragments/sec | With multi-threading |
| **Memory Usage** | ~50 MB | For typical bacterial genome |
| **Binary Size** | ~4 MB | Statically linked, no runtime dependencies |

### Performance Features

- **Region-Based Parallelism**: BAM files are split into genomic regions and processed concurrently
- **Interval Tree Indexing**: O(log n + k) overlap queries using the bio crate's interval tree
- **Zero-Copy Parsing**: Efficient BAM parsing with noodles (pure Rust, no htslib dependency)
- **Smart Thread Allocation**: With `-T 0`, automatically determines optimal thread count based on current system load
- **FxHash**: Uses rustc-hash for faster hash map operations compared to standard HashMap
- **SmallVec Optimisation**: Reduces heap allocations for small collections

### Memory Efficiency

rust-featurecounts is designed to minimise memory usage:

- Annotation features are loaded once and shared across threads
- Read data is processed in a streaming fashion
- No intermediate files are created
- Rust's ownership model ensures no memory leaks

## Installation

### From Source

Ensure you have Rust installed (https://rustup.rs/), then:

```bash
git clone https://github.com/necoli1822/rust-featurecounts.git
cd rust-featurecounts
cargo build --release
```

The binary will be available at `target/release/rust-featurecounts`.

### From crates.io

```bash
cargo install rust-featurecounts
```

## Usage

```bash
rust-featurecounts [options] -a <annotation> -o <output> input.bam [input2.bam ...]
```

### Required Arguments

| Option | Description |
|--------|-------------|
| `-a <file>` | Annotation file (GFF/GTF/GenBank format) |
| `-o <file>` | Output file for count matrix |

### Optional Arguments

| Option | Description | Default |
|--------|-------------|---------|
| `-t <string>` | Feature type to count (exon, gene, CDS, auto) | exon |
| `-g <string>` | Attribute for feature ID (gene_id, locus_tag, etc.) | gene_id |
| `-s <int>` | Strandedness: 0=unstranded, 1=forward, 2=reverse | 0 |
| `-Q <int>` | Minimum mapping quality | 0 |
| `-p` | Paired-end mode (counts fragments) | false |
| `-T <int>` | Number of threads (0=auto based on system load) | 1 |
| `-m <string>` | Counting mode: intersection, union | intersection |
| `-h, --help` | Show help message | |
| `-v, --version` | Show version | |

### Examples

Basic usage with GFF annotation:

```bash
rust-featurecounts -a annotation.gff -o counts.tsv aligned.bam
```

Strand-specific counting with multiple samples:

```bash
rust-featurecounts \
    -a annotation.gff \
    -o counts.tsv \
    -t gene \
    -g locus_tag \
    -s 2 \
    -p \
    -T 0 \
    sample1.bam sample2.bam sample3.bam
```

Using GenBank annotation (automatically detected):

```bash
rust-featurecounts -a genome.gbk -o counts.tsv -t gene -g locus_tag aligned.bam
```

## Output Format

The output is a tab-separated file compatible with featureCounts and downstream analysis tools:

```
Geneid    Chr           Start   End     Strand  Length  sample1.bam  sample2.bam
gene001   NC_000913.3   190     255     +       66      22           18
gene002   NC_000913.3   337     2799    +       2463    453          512
```

A summary file (`.summary`) is also generated with statistics:

```
rust-featurecounts Summary
==========================

Overall Statistics:
  Total reads processed: 919574
  Total fragments: 913011
  Assigned: 707488

Unassigned reads:
  No feature: 33544
  Ambiguous: 171979
  Low quality: 0
```

## Supported Annotation Formats

### GFF/GTF
Standard GFF3 and GTF formats are supported. Features are extracted based on the `-t` option (default: exon).

### GenBank
GenBank flat file format (.gbk, .gb, .genbank, .gbff) is automatically detected and parsed. The VERSION field is used for chromosome/sequence matching with BAM references.

Supported qualifiers for gene identification:
- `locus_tag` (recommended for bacterial genomes)
- `gene`
- `protein_id`

## Differences from Original featureCounts

This is a clean-room Rust reimplementation, not a port of the original C code. Key differences:

1. **Prokaryote Focus**: Optimised for bacterial and archaeal genomes; not tested with eukaryotes
2. **GenBank Support**: Native support for GenBank annotation files
3. **Auto Thread Detection**: Intelligent thread allocation based on system load
4. **Pure Rust**: No HTS library dependency; uses noodles for BAM parsing
5. **Simplified Options**: Focused on bacterial RNA-seq workflows

## Licence

This project is licensed under the MIT Licence. See the [LICENCE](LICENCE) file for details.

Note: This is an independent reimplementation and is not affiliated with or derived from the original [Subread/featureCounts](https://subread.sourceforge.net/) project (GPL-3.0).

## Dependencies

- [noodles](https://github.com/zaeleus/noodles) - BAM/SAM file handling (MIT)
- [rayon](https://github.com/rayon-rs/rayon) - Parallel processing (MIT/Apache-2.0)
- [bio](https://github.com/rust-bio/rust-bio) - Bioinformatics algorithms (MIT)
- [rustc-hash](https://github.com/rust-lang/rustc-hash) - Fast hashing (MIT/Apache-2.0)
- [smallvec](https://github.com/servo/rust-smallvec) - Small vector optimisation (MIT/Apache-2.0)

## Citation

If you use rust-featurecounts in your research, please cite:

> Kim S. (2025). rust-featurecounts: A fast feature counting tool for prokaryotic RNA-seq analysis. https://github.com/necoli1822/rust-featurecounts

For the original featureCounts algorithm, please also cite:

> Liao Y, Smyth GK and Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30, 2014.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request at [GitHub](https://github.com/necoli1822/rust-featurecounts).
