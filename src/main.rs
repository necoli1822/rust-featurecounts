//! rust-featurecounts: A fast feature counting tool compatible with featureCounts

mod bam;
mod cli;
mod counter;
mod error;
mod gff;
mod output;

use cli::Args;
use counter::{CountingMode, FeatureCounter, Strandedness};
use error::Result;
use gff::GffReader;
use output::CountMatrixWriter;
use std::time::Instant;

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

/// Smart thread count based on system load
fn smart_thread_count() -> usize {
    let cores = std::thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(1);

    // Read 1-minute load average from /proc/loadavg
    let load = std::fs::read_to_string("/proc/loadavg")
        .ok()
        .and_then(|s| s.split_whitespace().next()?.parse::<f32>().ok())
        .unwrap_or(0.0);

    // Calculate available threads: cores - current_load, leave 1 for system
    let available = ((cores as f32 - load).max(1.0) as usize).min(cores.saturating_sub(1)).max(1);

    eprintln!("  System: {} cores, load avg: {:.2}, using {} threads", cores, load, available);
    available
}

fn run() -> Result<()> {
    let args = Args::parse().map_err(|e| error::AppError::new(e))?;

    // Determine number of threads
    let num_threads = if args.threads <= 0 {
        smart_thread_count()
    } else {
        args.threads as usize
    };

    // Set thread pool size
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .map_err(|e| error::AppError::new(format!("Failed to set thread pool: {}", e)))?;

    eprintln!("rust-featurecounts v0.1.0");
    eprintln!("Reading annotation: {:?}", args.annotation);

    // Parse strandedness (featureCounts compatible: 0, 1, 2)
    let strandedness = match args.stranded.as_str() {
        "0" | "unstranded" => Strandedness::Unstranded,
        "1" | "forward" => Strandedness::Forward,
        "2" | "reverse" => Strandedness::Reverse,
        _ => {
            eprintln!("Warning: Unknown strandedness '{}', using unstranded", args.stranded);
            Strandedness::Unstranded
        }
    };

    // Read GFF annotation
    let start = Instant::now();
    let gff_reader = GffReader::new(&args.annotation)?;

    // Handle feature type (support "auto" for backward compatibility)
    let feature_type = if args.feature_type == "auto" {
        gff_reader.detect_best_feature_type()?
    } else {
        args.feature_type.clone()
    };

    eprintln!("Using feature type: {}", feature_type);

    let features = gff_reader.parse_features(&feature_type, &args.attribute)?;
    eprintln!("Loaded {} features in {:?}", features.len(), start.elapsed());

    if features.is_empty() {
        return Err(error::AppError::new(format!(
            "No features found with type '{}' and attribute '{}'",
            feature_type, args.attribute
        )));
    }

    // Parse counting mode
    let counting_mode = match args.counting_mode.to_lowercase().as_str() {
        "intersection" | "intersect" => CountingMode::Intersection,
        "union" => CountingMode::Union,
        _ => {
            eprintln!("Warning: Unknown counting mode '{}', using intersection", args.counting_mode);
            CountingMode::Intersection
        }
    };

    // Initialize counter
    let mut counter = FeatureCounter::new(
        features,
        strandedness,
        args.min_mapq,
        args.paired,
        counting_mode,
    );
    counter.set_threads(num_threads);

    if num_threads > 1 {
        eprintln!("Using {} threads", num_threads);
    }

    // Process each BAM file
    let mut all_counts: Vec<(String, rustc_hash::FxHashMap<String, u64>)> = Vec::new();

    for bam_path in &args.bam_files {
        eprintln!("Processing: {:?}", bam_path);
        let start = Instant::now();

        let sample_name = bam_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        let counts = counter.count_reads(bam_path)?;

        let total_assigned: u64 = counts.values().sum();
        let stats = counter.get_stats();

        eprintln!(
            "{}: {} assigned, {} unassigned (no_feature: {}, ambiguous: {}, low_quality: {}) in {:?}",
            sample_name,
            total_assigned,
            stats.unassigned_no_feature + stats.unassigned_ambiguous + stats.unassigned_low_quality,
            stats.unassigned_no_feature,
            stats.unassigned_ambiguous,
            stats.unassigned_low_quality,
            start.elapsed()
        );

        all_counts.push((sample_name, counts));
    }

    // Write output
    eprintln!("Writing: {:?}", args.output);
    let writer = CountMatrixWriter::new(&args.output)?;
    writer.write(&all_counts, counter.get_feature_ids(), counter.get_feature_meta())?;

    // Write summary
    let summary_path = args.output.with_extension("summary");
    writer.write_summary(&summary_path, &all_counts, counter.get_stats())?;

    eprintln!("Done!");
    Ok(())
}
