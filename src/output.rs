//! Output writing module

use crate::counter::{CountStats, FeatureMeta};
use crate::error::{Context, Result};
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Count matrix writer
pub struct CountMatrixWriter {
    path: std::path::PathBuf,
}

impl CountMatrixWriter {
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(Self {
            path: path.as_ref().to_path_buf(),
        })
    }

    /// Write count matrix to file (featureCounts compatible format)
    /// Format: Geneid Chr Start End Strand Length sample1.bam sample2.bam ...
    pub fn write(
        &self,
        counts: &[(String, FxHashMap<String, u64>)],
        feature_ids: &[String],
        feature_meta: &FxHashMap<String, FeatureMeta>,
    ) -> Result<()> {
        let file = File::create(&self.path)
            .context(format!("Failed to create output file: {:?}", self.path))?;
        let mut writer = BufWriter::new(file);

        // Write header (featureCounts format)
        write!(writer, "Geneid\tChr\tStart\tEnd\tStrand\tLength")?;
        for (sample_name, _) in counts {
            write!(writer, "\t{}.bam", sample_name)?;
        }
        writeln!(writer)?;

        // Write counts for each feature
        for feature_id in feature_ids {
            write!(writer, "{}", feature_id)?;

            // Write metadata columns
            if let Some(meta) = feature_meta.get(feature_id) {
                write!(
                    writer,
                    "\t{}\t{}\t{}\t{}\t{}",
                    meta.chr, meta.start, meta.end, meta.strand, meta.length
                )?;
            } else {
                write!(writer, "\tNA\tNA\tNA\t.\t0")?;
            }

            // Write counts
            for (_, sample_counts) in counts {
                let count = sample_counts.get(feature_id).copied().unwrap_or(0);
                write!(writer, "\t{}", count)?;
            }
            writeln!(writer)?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Write summary statistics
    pub fn write_summary<P: AsRef<Path>>(
        &self,
        path: P,
        counts: &[(String, FxHashMap<String, u64>)],
        stats: &CountStats,
    ) -> Result<()> {
        let file = File::create(path.as_ref())?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "rust-featurecounts Summary")?;
        writeln!(writer, "==========================")?;
        writeln!(writer)?;

        writeln!(writer, "Overall Statistics:")?;
        writeln!(writer, "  Total reads processed: {}", stats.total_reads)?;
        writeln!(writer, "  Total fragments: {}", stats.total_fragments)?;
        writeln!(writer, "  Assigned: {}", stats.assigned)?;
        writeln!(writer)?;

        writeln!(writer, "Unassigned reads:")?;
        writeln!(writer, "  No feature: {}", stats.unassigned_no_feature)?;
        writeln!(writer, "  Ambiguous: {}", stats.unassigned_ambiguous)?;
        writeln!(writer, "  Low quality: {}", stats.unassigned_low_quality)?;
        writeln!(writer, "  Unmapped: {}", stats.unassigned_unmapped)?;
        writeln!(writer)?;

        let assign_rate = if stats.total_fragments > 0 {
            (stats.assigned as f64 / stats.total_fragments as f64) * 100.0
        } else {
            0.0
        };

        writeln!(
            writer,
            "Assignment rate: {:.2}% ({}/{})",
            assign_rate, stats.assigned, stats.total_fragments
        )?;
        writeln!(writer)?;

        writeln!(writer, "Per-sample counts:")?;
        for (sample_name, sample_counts) in counts {
            let total: u64 = sample_counts.values().sum();
            let non_zero = sample_counts.values().filter(|&&v| v > 0).count();
            writeln!(
                writer,
                "  {}: {} total counts, {} features with counts",
                sample_name, total, non_zero
            )?;
        }

        writer.flush()?;
        Ok(())
    }
}
