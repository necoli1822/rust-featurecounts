# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2025-01-04

### Added

- Initial release of rust-featurecounts
- Core feature counting functionality compatible with featureCounts output format
- Support for GFF, GTF, and GenBank annotation formats
- Automatic annotation format detection based on file extension
- GenBank parsing with VERSION field for BAM reference matching
- Pseudo gene filtering for GenBank files to match GFF output
- Strand-specific counting modes (unstranded, forward, reverse)
- Paired-end fragment counting with `-p` option
- MAPQ-based quality filtering with `-Q` option
- Multi-threaded processing with `-T` option
- Smart thread allocation based on system load (`-T 0`)
- Intersection and union counting modes
- Automatic feature type detection with `-t auto`
- Summary statistics output file
- featureCounts-compatible output format with Geneid, Chr, Start, End, Strand, Length columns

### Technical Details

- Pure Rust implementation using noodles for BAM/SAM parsing
- Interval tree-based overlap detection using bio crate
- Parallel region-based BAM processing with rayon
- Memory-efficient design with rustc-hash and smallvec

[0.1.0]: https://github.com/necoli1822/rust-featurecounts/releases/tag/v0.1.0
