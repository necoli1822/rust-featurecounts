//! BAM file reading module using noodles (pure Rust)

use crate::error::{Context, Result, AppError};
use noodles::bam;
use noodles::bgzf;
use noodles::sam::alignment::Record;
use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::BufReader;
use std::path::Path;

/// Region information for parallel processing
#[derive(Debug, Clone)]
pub struct BamRegion {
    pub name: String,
    pub start: u64,
    pub end: u64,
}

/// Represents a read or fragment from BAM
#[derive(Debug, Clone)]
pub struct AlignedRead {
    pub qname: String,
    pub tid: i32,
    pub start: i64,
    pub end: i64,
    pub mapq: u8,
    pub is_reverse: bool,
    pub is_read1: bool,
    pub is_paired: bool,
    pub is_unmapped: bool,
}

/// Compact read struct - uses hash instead of string for qname (saves ~40 bytes per read)
#[derive(Debug, Clone, Copy)]
pub struct CompactRead {
    pub qname_hash: u64,
    pub start: i64,
    pub end: i64,
    pub mapq: u8,
    pub flags: u8, // packed: is_reverse(1), is_read1(2), is_paired(4), is_unmapped(8)
}

impl CompactRead {
    #[inline]
    pub fn is_reverse(&self) -> bool {
        self.flags & 1 != 0
    }

    #[inline]
    pub fn is_read1(&self) -> bool {
        self.flags & 2 != 0
    }

    #[inline]
    pub fn is_paired(&self) -> bool {
        self.flags & 4 != 0
    }

    #[inline]
    pub fn is_unmapped(&self) -> bool {
        self.flags & 8 != 0
    }
}

/// Fast hash function for qname string
#[inline]
fn hash_qname_str(qname: &str) -> u64 {
    let mut hasher = DefaultHasher::new();
    qname.hash(&mut hasher);
    hasher.finish()
}

/// Simple BAM reader
pub struct BamReader {
    reader: bam::io::Reader<bgzf::Reader<BufReader<File>>>,
    header: noodles::sam::Header,
    reference_names: Vec<String>,
}

impl BamReader {
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(&path).context("Failed to open BAM file")?;
        let mut reader = bam::io::reader::Builder::default()
            .build_from_reader(BufReader::new(file));

        let header = reader.read_header().context("Failed to read BAM header")?;

        // Extract reference names
        let reference_names: Vec<String> = header
            .reference_sequences()
            .keys()
            .map(|name| name.to_string())
            .collect();

        Ok(Self {
            reader,
            header,
            reference_names,
        })
    }

    /// Get reference name by tid
    pub fn tid_to_name(&self, tid: i32) -> Option<String> {
        if tid < 0 {
            return None;
        }
        self.reference_names.get(tid as usize).cloned()
    }

    /// Iterate over all reads
    pub fn records(&mut self) -> BamRecordIterator<'_> {
        BamRecordIterator {
            reader: &mut self.reader,
            header: &self.header,
            record: bam::Record::default(),
        }
    }
}

pub struct BamRecordIterator<'a> {
    reader: &'a mut bam::io::Reader<bgzf::Reader<BufReader<File>>>,
    #[allow(dead_code)]
    header: &'a noodles::sam::Header,
    record: bam::Record,
}

impl<'a> Iterator for BamRecordIterator<'a> {
    type Item = Result<AlignedRead>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record(&mut self.record) {
            Ok(0) => None, // EOF
            Ok(_) => Some(self.convert_record()),
            Err(e) => Some(Err(AppError::new(format!("Error reading BAM record: {}", e)))),
        }
    }
}

impl<'a> BamRecordIterator<'a> {
    fn convert_record(&self) -> Result<AlignedRead> {
        let flags = self.record.flags();

        // Get read name
        let qname = self
            .record
            .name()
            .map(|n| n.to_string())
            .unwrap_or_else(|| "*".to_string());

        // Get reference sequence id (tid)
        let tid = self
            .record
            .reference_sequence_id()
            .transpose()
            .context("Invalid reference sequence ID")?
            .map(|id| id as i32)
            .unwrap_or(-1);

        // Get alignment position (0-based)
        let start = self
            .record
            .alignment_start()
            .transpose()
            .context("Invalid alignment start")?
            .map(|pos| (usize::from(pos) - 1) as i64)  // Convert to 0-based
            .unwrap_or(0);

        // Calculate end position from CIGAR
        let end = self
            .record
            .alignment_end()
            .transpose()
            .context("Invalid alignment end")?
            .map(|pos| usize::from(pos) as i64)
            .unwrap_or(start);

        // Get mapping quality
        let mapq = self
            .record
            .mapping_quality()
            .map(|q| q.get())
            .unwrap_or(255);

        Ok(AlignedRead {
            qname,
            tid,
            start,
            end,
            mapq,
            is_reverse: flags.is_reverse_complemented(),
            is_read1: flags.is_first_segment(),
            is_paired: flags.is_segmented(),
            is_unmapped: flags.is_unmapped(),
        })
    }
}

/// Read BAM region into Vec of CompactRead (memory efficient)
pub fn read_region_compact<P: AsRef<Path>>(path: P, region: &BamRegion) -> Result<Vec<CompactRead>> {
    use noodles::bam::bai;
    use noodles::core::Region;

    let path = path.as_ref();

    let mut reader = File::open(path)
        .map(BufReader::new)
        .map(bam::io::Reader::new)
        .context("Failed to open BAM file")?;

    let header = reader.read_header().context("Failed to read BAM header")?;

    let index_path = path.with_extension("bam.bai");
    let index = bai::read(&index_path)
        .or_else(|_| {
            let alt_path = format!("{}.bai", path.display());
            bai::read(&alt_path)
        })
        .context("Failed to open BAM index")?;

    let region_str = format!("{}:{}-{}", region.name, region.start + 1, region.end);
    let query_region: Region = region_str.parse().context("Failed to parse region")?;

    let query = reader
        .query(&header, &index, &query_region)
        .context("Failed to query BAM region")?;

    let mut reads = Vec::new();

    for record_result in query {
        let record = match record_result {
            Ok(r) => r,
            Err(_) => continue,
        };

        let flags = record.flags();

        let qname_hash = record
            .name()
            .map(|n| hash_qname_str(&n.to_string()))
            .unwrap_or(0);

        let start = record
            .alignment_start()
            .transpose()
            .ok()
            .flatten()
            .map(|pos| (usize::from(pos) - 1) as i64)
            .unwrap_or(0);

        let end = record
            .alignment_end()
            .transpose()
            .ok()
            .flatten()
            .map(|pos| usize::from(pos) as i64)
            .unwrap_or(start);

        let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(255);

        let packed_flags = (flags.is_reverse_complemented() as u8)
            | ((flags.is_first_segment() as u8) << 1)
            | ((flags.is_segmented() as u8) << 2)
            | ((flags.is_unmapped() as u8) << 3);

        reads.push(CompactRead {
            qname_hash,
            start,
            end,
            mapq,
            flags: packed_flags,
        });
    }

    Ok(reads)
}

/// Get BAM regions for parallel processing
pub fn get_bam_regions<P: AsRef<Path>>(path: P, num_chunks: usize) -> Result<Vec<BamRegion>> {
    let file = File::open(&path).context("Failed to open BAM file")?;
    let mut reader = bam::io::Reader::new(BufReader::new(file));

    let header = reader.read_header().context("Failed to read BAM header")?;

    let mut regions = Vec::new();

    for (_tid, (name, ref_seq)) in header.reference_sequences().iter().enumerate() {
        let length = ref_seq.length().get() as u64;

        if length == 0 {
            continue;
        }

        let chunk_size = (length / num_chunks as u64).max(1);

        for i in 0..num_chunks {
            let start = i as u64 * chunk_size;
            let end = if i == num_chunks - 1 {
                length
            } else {
                (i as u64 + 1) * chunk_size
            };

            if start < length {
                regions.push(BamRegion {
                    name: name.to_string(),
                    start,
                    end,
                });
            }
        }
    }

    Ok(regions)
}
