//! Feature counting module with parallel processing and memory optimization

use crate::bam::{get_bam_regions, read_region_compact, BamReader};
use crate::gff::{Feature, Strand};
use crate::error::Result;
use bio::data_structures::interval_tree::IntervalTree;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use smallvec::SmallVec;
use std::ops::Range;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

/// Strandedness mode
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strandedness {
    Unstranded,
    Forward,
    Reverse,
}

/// Counting mode for paired-end reads
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CountingMode {
    Intersection,
    Union,
}

/// Statistics for counting
#[derive(Debug, Clone, Default)]
pub struct CountStats {
    pub total_reads: u64,
    pub total_fragments: u64,
    pub assigned: u64,
    pub unassigned_no_feature: u64,
    pub unassigned_ambiguous: u64,
    pub unassigned_low_quality: u64,
    pub unassigned_unmapped: u64,
}

/// Feature metadata for output (featureCounts compatible)
#[derive(Debug, Clone)]
pub struct FeatureMeta {
    pub chr: String,
    pub start: String,
    pub end: String,
    pub strand: char,
    pub length: u64,
}

/// Feature entry in interval tree
#[derive(Debug, Clone)]
struct FeatureEntry {
    feature_id: String,
    strand: Strand,
}

/// Shared data for parallel processing
struct SharedData {
    trees: FxHashMap<String, IntervalTree<u64, FeatureEntry>>,
    strandedness: Strandedness,
    feature_to_idx: FxHashMap<String, u16>,
}

/// Type alias for memory-efficient feature list (most reads have 0-2 features)
type FeatureIndices = SmallVec<[u16; 2]>;

impl SharedData {
    /// Find overlapping features and return their indices (memory efficient)
    /// Uses SmallVec to avoid heap allocation for common case of 0-2 features
    #[inline]
    fn find_overlapping_feature_indices(
        &self,
        ref_name: &str,
        start: u64,
        end: u64,
        is_reverse: bool,
        is_read1: bool,
    ) -> FeatureIndices {
        let tree = match self.trees.get(ref_name) {
            Some(t) => t,
            None => return FeatureIndices::new(),
        };

        let overlaps = tree.find(start..end);
        let mut matching_indices = FeatureIndices::new();

        for entry in overlaps {
            let feature = entry.data();

            if self.strandedness != Strandedness::Unstranded {
                let read_strand = get_read_strand(self.strandedness, is_reverse, is_read1);
                if !strand_matches(read_strand, feature.strand) {
                    continue;
                }
            }

            if let Some(&idx) = self.feature_to_idx.get(&feature.feature_id) {
                // Avoid duplicates (use linear search for small vec)
                if !matching_indices.contains(&idx) {
                    matching_indices.push(idx);
                }
            }
        }

        matching_indices
    }
}

/// Feature counter using interval tree for efficient overlap detection
pub struct FeatureCounter {
    trees: FxHashMap<String, IntervalTree<u64, FeatureEntry>>,
    feature_ids: Vec<String>,
    feature_meta: FxHashMap<String, FeatureMeta>,
    strandedness: Strandedness,
    min_mapq: u8,
    paired_mode: bool,
    counting_mode: CountingMode,
    stats: CountStats,
    num_threads: usize,
}

impl FeatureCounter {
    pub fn new(
        features: Vec<Feature>,
        strandedness: Strandedness,
        min_mapq: u8,
        paired_mode: bool,
        counting_mode: CountingMode,
    ) -> Self {
        let mut trees: FxHashMap<String, Vec<(Range<u64>, FeatureEntry)>> = FxHashMap::default();
        let mut feature_ids = Vec::new();
        let mut feature_meta: FxHashMap<String, FeatureMeta> = FxHashMap::default();
        let mut seen_ids: FxHashSet<String> = FxHashSet::default();

        for feature in &features {
            let start = feature.start.saturating_sub(1);
            let end = feature.end;

            let entry = FeatureEntry {
                feature_id: feature.feature_id.clone(),
                strand: feature.strand,
            };

            trees
                .entry(feature.seqname.clone())
                .or_default()
                .push((start..end, entry));

            // Build feature metadata (handle multiple exons per gene)
            let length = feature.end.saturating_sub(feature.start) + 1;
            let strand_char = match feature.strand {
                Strand::Plus => '+',
                Strand::Minus => '-',
                Strand::Unknown => '.',
            };

            if let Some(meta) = feature_meta.get_mut(&feature.feature_id) {
                // Append to existing metadata (multiple exons)
                meta.chr.push(';');
                meta.chr.push_str(&feature.seqname);
                meta.start.push(';');
                meta.start.push_str(&feature.start.to_string());
                meta.end.push(';');
                meta.end.push_str(&feature.end.to_string());
                meta.length += length;
            } else {
                feature_meta.insert(
                    feature.feature_id.clone(),
                    FeatureMeta {
                        chr: feature.seqname.clone(),
                        start: feature.start.to_string(),
                        end: feature.end.to_string(),
                        strand: strand_char,
                        length,
                    },
                );
            }

            if !seen_ids.contains(&feature.feature_id) {
                feature_ids.push(feature.feature_id.clone());
                seen_ids.insert(feature.feature_id.clone());
            }
        }

        // Build interval trees
        let trees: FxHashMap<String, IntervalTree<u64, FeatureEntry>> = trees
            .into_par_iter()
            .map(|(seqname, intervals)| {
                let tree: IntervalTree<u64, FeatureEntry> = intervals.into_iter().collect();
                (seqname, tree)
            })
            .collect();

        Self {
            trees,
            feature_ids,
            feature_meta,
            strandedness,
            min_mapq,
            paired_mode,
            counting_mode,
            stats: CountStats::default(),
            num_threads: 1,
        }
    }

    /// Set number of threads for parallel processing
    pub fn set_threads(&mut self, num_threads: usize) {
        self.num_threads = num_threads;
    }

    pub fn get_feature_ids(&self) -> &[String] {
        &self.feature_ids
    }

    pub fn get_feature_meta(&self) -> &FxHashMap<String, FeatureMeta> {
        &self.feature_meta
    }

    pub fn get_stats(&self) -> &CountStats {
        &self.stats
    }

    /// Count reads from a BAM file
    pub fn count_reads<P: AsRef<Path>>(
        &mut self,
        bam_path: P,
    ) -> Result<FxHashMap<String, u64>> {
        if self.num_threads > 1 {
            self.count_reads_parallel(bam_path)
        } else {
            self.count_reads_sequential(bam_path)
        }
    }

    /// Sequential counting
    fn count_reads_sequential<P: AsRef<Path>>(
        &mut self,
        bam_path: P,
    ) -> Result<FxHashMap<String, u64>> {
        let mut reader = BamReader::new(&bam_path)?;

        let ref_names: FxHashMap<i32, String> = (0..100)
            .filter_map(|tid| reader.tid_to_name(tid).map(|name| (tid, name)))
            .collect();

        let mut counts: FxHashMap<String, u64> = self
            .feature_ids
            .iter()
            .map(|id| (id.clone(), 0))
            .collect();

        self.stats = CountStats::default();
        let mut processed_fragments: FxHashMap<String, FxHashSet<String>> = FxHashMap::default();
        let mut read_count = 0u64;

        for result in reader.records() {
            let read = result?;
            read_count += 1;
            if read_count % 500_000 == 0 {
                eprintln!("  Processed {} reads...", read_count);
            }

            self.stats.total_reads += 1;

            if read.is_unmapped {
                self.stats.unassigned_unmapped += 1;
                continue;
            }

            if read.mapq < self.min_mapq {
                self.stats.unassigned_low_quality += 1;
                continue;
            }

            let ref_name = match ref_names.get(&read.tid) {
                Some(name) => name,
                None => continue,
            };

            let overlapping = self.find_overlapping_features(
                ref_name,
                read.start as u64,
                read.end as u64,
                read.is_reverse,
                read.is_read1,
            );

            if self.paired_mode && read.is_paired {
                let fragment_key = read.qname.clone();

                if let Some(existing_features) = processed_fragments.remove(&fragment_key) {
                    let overlapping_set: FxHashSet<String> = overlapping.into_iter().collect();
                    let combined = self.combine_features(existing_features, overlapping_set);
                    self.count_fragment(&combined, &mut counts);
                    self.stats.total_fragments += 1;
                } else {
                    let feature_set: FxHashSet<String> = overlapping.into_iter().collect();
                    processed_fragments.insert(fragment_key, feature_set);
                }
            } else {
                self.stats.total_fragments += 1;
                let feature_set: FxHashSet<String> = overlapping.into_iter().collect();
                self.count_fragment(&feature_set, &mut counts);
            }
        }

        // Process orphans
        if self.paired_mode {
            for (_qname, features) in processed_fragments {
                self.stats.total_fragments += 1;
                self.count_fragment(&features, &mut counts);
            }
        }

        Ok(counts)
    }

    /// Parallel counting using indexed BAM reading
    fn count_reads_parallel<P: AsRef<Path>>(
        &mut self,
        bam_path: P,
    ) -> Result<FxHashMap<String, u64>> {
        let bam_path = bam_path.as_ref();

        // Get regions to process in parallel
        let regions = get_bam_regions(bam_path, self.num_threads)?;
        eprintln!("  Splitting BAM into {} regions for parallel processing", regions.len());

        // Build feature_id -> index mapping for memory efficiency
        let feature_to_idx: FxHashMap<String, u16> = self
            .feature_ids
            .iter()
            .enumerate()
            .map(|(i, id)| (id.clone(), i as u16))
            .collect();

        // Create shared data
        let shared = Arc::new(SharedData {
            trees: std::mem::take(&mut self.trees),
            strandedness: self.strandedness,
            feature_to_idx,
        });

        self.stats = CountStats::default();
        let counting_mode = self.counting_mode;
        let paired_mode = self.paired_mode;
        let min_mapq = self.min_mapq;

        // Atomic counters for stats
        let total_reads = AtomicU64::new(0);
        let unmapped = AtomicU64::new(0);
        let low_quality = AtomicU64::new(0);

        eprintln!("  Processing {} regions in parallel...", regions.len());

        // Process regions in parallel
        type PendingReads = FxHashMap<u64, FeatureIndices>;
        type LocalCounts = FxHashMap<u16, u64>;

        struct RegionResult {
            pending: PendingReads,
            counts: LocalCounts,
            fragments: u64,
            assigned: u64,
            no_feature: u64,
            ambiguous: u64,
        }

        let region_results: Vec<RegionResult> = regions
            .par_iter()
            .map(|region| {
                let shared = Arc::clone(&shared);
                let mut pending: PendingReads = FxHashMap::default();
                let mut counts: LocalCounts = FxHashMap::default();
                let mut fragments = 0u64;
                let mut assigned = 0u64;
                let mut no_feature = 0u64;
                let mut ambiguous = 0u64;

                // Use CompactRead with pre-computed qname_hash (memory optimized)
                if let Ok(reads) = read_region_compact(bam_path, region) {
                    pending.reserve(reads.len() / 2);

                    for read in reads {
                        let read_start = read.start as u64;
                        if read_start < region.start || read_start >= region.end {
                            continue;
                        }

                        total_reads.fetch_add(1, Ordering::Relaxed);

                        if read.is_unmapped() {
                            unmapped.fetch_add(1, Ordering::Relaxed);
                            continue;
                        }

                        if read.mapq < min_mapq {
                            low_quality.fetch_add(1, Ordering::Relaxed);
                            continue;
                        }

                        let features = shared.find_overlapping_feature_indices(
                            &region.name,
                            read_start,
                            read.end as u64,
                            read.is_reverse(),
                            read.is_read1(),
                        );

                        let qname_hash = read.qname_hash;

                        if paired_mode && read.is_paired() {
                            if let Some(mate_features) = pending.remove(&qname_hash) {
                                let combined = combine_indices_smallvec(
                                    mate_features,
                                    features,
                                    counting_mode,
                                );
                                fragments += 1;
                                count_fragment_indices(
                                    &combined,
                                    counting_mode,
                                    &mut counts,
                                    &mut assigned,
                                    &mut no_feature,
                                    &mut ambiguous,
                                );
                            } else {
                                pending.insert(qname_hash, features);
                            }
                        } else {
                            fragments += 1;
                            count_fragment_indices(
                                &features,
                                counting_mode,
                                &mut counts,
                                &mut assigned,
                                &mut no_feature,
                                &mut ambiguous,
                            );
                        }
                    }
                }

                RegionResult {
                    pending,
                    counts,
                    fragments,
                    assigned,
                    no_feature,
                    ambiguous,
                }
            })
            .collect();

        self.stats.total_reads = total_reads.load(Ordering::Relaxed);
        self.stats.unassigned_unmapped = unmapped.load(Ordering::Relaxed);
        self.stats.unassigned_low_quality = low_quality.load(Ordering::Relaxed);

        // Phase 2: Merge pending reads from all regions

        let mut total_fragments = 0u64;
        let mut total_assigned = 0u64;
        let mut total_no_feature = 0u64;
        let mut total_ambiguous = 0u64;
        let mut aggregated_counts: FxHashMap<u16, u64> = FxHashMap::default();
        let mut global_pending: FxHashMap<u64, FeatureIndices> = FxHashMap::default();

        for result in region_results {
            for (idx, count) in result.counts {
                *aggregated_counts.entry(idx).or_insert(0) += count;
            }
            total_fragments += result.fragments;
            total_assigned += result.assigned;
            total_no_feature += result.no_feature;
            total_ambiguous += result.ambiguous;

            for (qname_hash, features) in result.pending {
                if let Some(mate_features) = global_pending.remove(&qname_hash) {
                    let combined = combine_indices_smallvec(mate_features, features, counting_mode);
                    total_fragments += 1;
                    count_fragment_indices(
                        &combined,
                        counting_mode,
                        &mut aggregated_counts,
                        &mut total_assigned,
                        &mut total_no_feature,
                        &mut total_ambiguous,
                    );
                } else {
                    global_pending.insert(qname_hash, features);
                }
            }
        }

        // Process remaining orphan reads
        for (_qname_hash, features) in global_pending {
            total_fragments += 1;
            count_fragment_indices(
                &features,
                counting_mode,
                &mut aggregated_counts,
                &mut total_assigned,
                &mut total_no_feature,
                &mut total_ambiguous,
            );
        }

        eprintln!("  Processed {} fragments", total_fragments);

        // Convert indices back to feature IDs

        let mut final_counts: FxHashMap<String, u64> = self
            .feature_ids
            .iter()
            .map(|id| (id.clone(), 0))
            .collect();

        for (idx, count) in aggregated_counts {
            if let Some(feature_id) = self.feature_ids.get(idx as usize) {
                *final_counts.entry(feature_id.clone()).or_insert(0) += count;
            }
        }

        // Restore trees
        self.trees = Arc::try_unwrap(shared)
            .map(|s| s.trees)
            .unwrap_or_default();

        // Update stats
        self.stats.total_fragments = total_fragments;
        self.stats.assigned = total_assigned;
        self.stats.unassigned_no_feature = total_no_feature;
        self.stats.unassigned_ambiguous = total_ambiguous;

        Ok(final_counts)
    }

    /// Find overlapping features for a read
    fn find_overlapping_features(
        &self,
        ref_name: &str,
        start: u64,
        end: u64,
        is_reverse: bool,
        is_read1: bool,
    ) -> Vec<String> {
        let tree = match self.trees.get(ref_name) {
            Some(t) => t,
            None => return Vec::new(),
        };

        let overlaps = tree.find(start..end);
        let mut matching_features: FxHashSet<String> = FxHashSet::default();

        for entry in overlaps {
            let feature = entry.data();

            if self.strandedness != Strandedness::Unstranded {
                let read_strand = get_read_strand(self.strandedness, is_reverse, is_read1);
                if !strand_matches(read_strand, feature.strand) {
                    continue;
                }
            }

            matching_features.insert(feature.feature_id.clone());
        }

        matching_features.into_iter().collect()
    }

    /// Combine features from two reads
    fn combine_features(
        &self,
        features1: FxHashSet<String>,
        features2: FxHashSet<String>,
    ) -> FxHashSet<String> {
        combine_features_static(features1, features2, self.counting_mode)
    }

    /// Count a single fragment
    fn count_fragment(&mut self, features: &FxHashSet<String>, counts: &mut FxHashMap<String, u64>) {
        match features.len() {
            0 => {
                self.stats.unassigned_no_feature += 1;
            }
            1 => {
                let feature_id = features.iter().next().unwrap();
                *counts.entry(feature_id.clone()).or_insert(0) += 1;
                self.stats.assigned += 1;
            }
            _ => {
                self.stats.unassigned_ambiguous += 1;
            }
        }
    }
}

/// Get read strand based on strandedness mode
fn get_read_strand(strandedness: Strandedness, is_reverse: bool, is_read1: bool) -> Strand {
    match strandedness {
        Strandedness::Unstranded => Strand::Unknown,
        Strandedness::Forward => {
            if is_read1 {
                if is_reverse {
                    Strand::Minus
                } else {
                    Strand::Plus
                }
            } else {
                if is_reverse {
                    Strand::Plus
                } else {
                    Strand::Minus
                }
            }
        }
        Strandedness::Reverse => {
            if is_read1 {
                if is_reverse {
                    Strand::Plus
                } else {
                    Strand::Minus
                }
            } else {
                if is_reverse {
                    Strand::Minus
                } else {
                    Strand::Plus
                }
            }
        }
    }
}

/// Check if read strand matches feature strand
fn strand_matches(read_strand: Strand, feature_strand: Strand) -> bool {
    match (read_strand, feature_strand) {
        (Strand::Unknown, _) | (_, Strand::Unknown) => true,
        (Strand::Plus, Strand::Plus) | (Strand::Minus, Strand::Minus) => true,
        _ => false,
    }
}

/// Static version of combine_features for use in parallel context
fn combine_features_static(
    features1: FxHashSet<String>,
    features2: FxHashSet<String>,
    counting_mode: CountingMode,
) -> FxHashSet<String> {
    match counting_mode {
        CountingMode::Intersection => {
            if features1.is_empty() {
                features2
            } else if features2.is_empty() {
                features1
            } else {
                features1.intersection(&features2).cloned().collect()
            }
        }
        CountingMode::Union => {
            let mut union = features1;
            union.extend(features2);
            union
        }
    }
}

/// Combine feature indices using SmallVec (avoids heap allocation for common cases)
#[inline]
fn combine_indices_smallvec(
    indices1: FeatureIndices,
    indices2: FeatureIndices,
    counting_mode: CountingMode,
) -> FeatureIndices {
    match counting_mode {
        CountingMode::Intersection => {
            if indices1.is_empty() {
                indices2
            } else if indices2.is_empty() {
                indices1
            } else {
                // Intersection: keep only elements in both
                let mut result = FeatureIndices::new();
                for idx in &indices1 {
                    if indices2.contains(idx) {
                        result.push(*idx);
                    }
                }
                result
            }
        }
        CountingMode::Union => {
            // Union: combine both, avoiding duplicates
            let mut result = indices1;
            for idx in indices2 {
                if !result.contains(&idx) {
                    result.push(idx);
                }
            }
            result
        }
    }
}

/// Count a fragment using feature indices (inline for performance)
#[inline]
fn count_fragment_indices(
    features: &FeatureIndices,
    _counting_mode: CountingMode,
    counts: &mut FxHashMap<u16, u64>,
    assigned: &mut u64,
    no_feature: &mut u64,
    ambiguous: &mut u64,
) {
    match features.len() {
        0 => {
            *no_feature += 1;
        }
        1 => {
            *counts.entry(features[0]).or_insert(0) += 1;
            *assigned += 1;
        }
        _ => {
            *ambiguous += 1;
        }
    }
}
