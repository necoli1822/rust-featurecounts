//! Annotation file parsing module (GFF/GTF/GenBank)

use crate::error::{Context, Result, AppError};
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Annotation file format
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotationFormat {
    Gff,
    GenBank,
}

impl AnnotationFormat {
    /// Detect format from file extension
    pub fn from_path<P: AsRef<Path>>(path: P) -> Self {
        let ext = path.as_ref()
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.to_lowercase())
            .unwrap_or_default();

        match ext.as_str() {
            "gbk" | "gb" | "genbank" | "gbff" => AnnotationFormat::GenBank,
            _ => AnnotationFormat::Gff,
        }
    }
}

/// Represents a genomic feature from GFF/GTF
#[derive(Debug, Clone)]
pub struct Feature {
    pub seqname: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub feature_id: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

impl Strand {
    pub fn from_char(c: char) -> Self {
        match c {
            '+' => Strand::Plus,
            '-' => Strand::Minus,
            _ => Strand::Unknown,
        }
    }
}

/// Annotation file reader (GFF/GTF/GenBank)
pub struct GffReader {
    path: std::path::PathBuf,
    format: AnnotationFormat,
    feature_counts: FxHashMap<String, usize>,
    seqname: String, // For GenBank: the LOCUS name
}

impl GffReader {
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        if !path.exists() {
            return Err(AppError::new(format!("Annotation file not found: {:?}", path)));
        }

        let format = AnnotationFormat::from_path(&path);

        let mut reader = Self {
            path,
            format,
            feature_counts: FxHashMap::default(),
            seqname: String::new(),
        };

        match format {
            AnnotationFormat::Gff => reader.scan_feature_types_gff()?,
            AnnotationFormat::GenBank => reader.scan_feature_types_genbank()?,
        }

        Ok(reader)
    }

    /// Scan GFF/GTF file to count feature types
    fn scan_feature_types_gff(&mut self) -> Result<()> {
        let file = File::open(&self.path)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 9 {
                continue;
            }

            let feature_type = fields[2];
            *self.feature_counts.entry(feature_type.to_string()).or_insert(0) += 1;
        }

        Ok(())
    }

    /// Scan GenBank file to count feature types
    fn scan_feature_types_genbank(&mut self) -> Result<()> {
        let file = File::open(&self.path)?;
        let reader = BufReader::new(file);

        let mut in_features = false;
        let mut locus_name = String::new();

        for line in reader.lines() {
            let line = line?;

            // Get LOCUS name (fallback)
            if line.starts_with("LOCUS") {
                if let Some(name) = line.split_whitespace().nth(1) {
                    locus_name = name.to_string();
                }
                continue;
            }

            // Get VERSION (preferred - includes version number like NC_000913.3)
            if line.starts_with("VERSION") {
                if let Some(version) = line.split_whitespace().nth(1) {
                    self.seqname = version.to_string();
                }
                continue;
            }

            // Start of FEATURES section
            if line.starts_with("FEATURES") {
                in_features = true;
                continue;
            }

            // End of FEATURES section
            if in_features && !line.starts_with(' ') && !line.is_empty() {
                break;
            }

            if !in_features {
                continue;
            }

            // Feature line starts with 5 spaces then feature type
            if line.len() > 5 && line.chars().take(5).all(|c| c == ' ') {
                let rest = line[5..].trim();
                if let Some(feature_type) = rest.split_whitespace().next() {
                    if !feature_type.starts_with('/') {
                        *self.feature_counts.entry(feature_type.to_string()).or_insert(0) += 1;
                    }
                }
            }
        }

        // Use LOCUS if VERSION not found
        if self.seqname.is_empty() {
            self.seqname = locus_name;
        }

        eprintln!("GenBank format detected, seqname: {}", self.seqname);
        Ok(())
    }

    /// Detect the best feature type to use (gene or CDS, whichever has more entries)
    pub fn detect_best_feature_type(&self) -> Result<String> {
        let gene_count = self.feature_counts.get("gene").copied().unwrap_or(0);
        let cds_count = self.feature_counts.get("CDS").copied().unwrap_or(0);

        eprintln!("Feature counts - gene: {}, CDS: {}", gene_count, cds_count);

        if gene_count == 0 && cds_count == 0 {
            // Try other feature types
            if let Some((feature_type, count)) = self
                .feature_counts
                .iter()
                .filter(|(k, _)| *k != "region" && *k != "exon")
                .max_by_key(|(_, v)| *v)
            {
                eprintln!(
                    "Neither gene nor CDS found. Using '{}' ({} features)",
                    feature_type, count
                );
                return Ok(feature_type.clone());
            }
            return Err(AppError::new("No suitable feature type found in annotation file"));
        }

        // Use whichever has more entries
        let best = if gene_count >= cds_count {
            "gene"
        } else {
            "CDS"
        };

        eprintln!(
            "Auto-selected feature type: {} ({} features)",
            best,
            if best == "gene" { gene_count } else { cds_count }
        );

        Ok(best.to_string())
    }

    /// Parse features of a specific type
    pub fn parse_features(&self, feature_type: &str, id_attribute: &str) -> Result<Vec<Feature>> {
        match self.format {
            AnnotationFormat::Gff => self.parse_features_gff(feature_type, id_attribute),
            AnnotationFormat::GenBank => self.parse_features_genbank(feature_type, id_attribute),
        }
    }

    /// Parse features from GFF/GTF file
    fn parse_features_gff(&self, feature_type: &str, id_attribute: &str) -> Result<Vec<Feature>> {
        let file = File::open(&self.path).context("Failed to open annotation file")?;
        let reader = BufReader::new(file);

        let mut features = Vec::new();
        let id_attributes: Vec<&str> = id_attribute.split(',').collect();

        for (line_num, line) in reader.lines().enumerate() {
            let line = line?;
            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 9 {
                continue;
            }

            if fields[2] != feature_type {
                continue;
            }

            let seqname = fields[0].to_string();
            let start: u64 = fields[3].parse().context(format!(
                "Invalid start position at line {}",
                line_num + 1
            ))?;
            let end: u64 = fields[4].parse().context(format!(
                "Invalid end position at line {}",
                line_num + 1
            ))?;
            let strand = Strand::from_char(fields[6].chars().next().unwrap_or('.'));

            // Parse attributes
            let attributes = parse_attributes(fields[8]);

            // Try to find feature ID using the specified attribute(s)
            let feature_id = id_attributes
                .iter()
                .filter_map(|attr| attributes.get(*attr))
                .next()
                .cloned();

            let feature_id = match feature_id {
                Some(id) => id,
                None => {
                    // Fallback: try common attributes
                    let fallback_attrs = ["locus_tag", "gene_id", "ID", "Name", "gene"];
                    fallback_attrs
                        .iter()
                        .filter_map(|attr| attributes.get(*attr))
                        .next()
                        .cloned()
                        .unwrap_or_else(|| format!("unknown_{}", line_num + 1))
                }
            };

            features.push(Feature {
                seqname,
                start,
                end,
                strand,
                feature_id,
            });
        }

        eprintln!(
            "Parsed {} {} features from GFF/GTF",
            features.len(),
            feature_type
        );

        Ok(features)
    }

    /// Parse features from GenBank file
    fn parse_features_genbank(&self, feature_type: &str, id_attribute: &str) -> Result<Vec<Feature>> {
        let file = File::open(&self.path).context("Failed to open annotation file")?;
        let reader = BufReader::new(file);

        let mut features = Vec::new();
        let id_attributes: Vec<&str> = id_attribute.split(',').collect();

        let mut in_features = false;
        let mut current_feature_type = String::new();
        let mut current_location = String::new();
        let mut current_qualifiers: FxHashMap<String, String> = FxHashMap::default();
        let mut seqname = self.seqname.clone();
        let mut feature_count = 0;

        let lines: Vec<String> = reader.lines().filter_map(|l| l.ok()).collect();
        let mut i = 0;

        while i < lines.len() {
            let line = &lines[i];

            // Get LOCUS name (fallback)
            if line.starts_with("LOCUS") {
                if seqname.is_empty() {
                    if let Some(name) = line.split_whitespace().nth(1) {
                        seqname = name.to_string();
                    }
                }
                i += 1;
                continue;
            }

            // Get VERSION (preferred)
            if line.starts_with("VERSION") {
                if let Some(version) = line.split_whitespace().nth(1) {
                    seqname = version.to_string();
                }
                i += 1;
                continue;
            }

            // Start of FEATURES section
            if line.starts_with("FEATURES") {
                in_features = true;
                i += 1;
                continue;
            }

            // End of FEATURES section (ORIGIN or other top-level keyword)
            if in_features && !line.starts_with(' ') && !line.is_empty() {
                // Save last feature if any
                if !current_feature_type.is_empty() && current_feature_type == feature_type {
                    if let Some(f) = parse_genbank_feature(
                        &seqname,
                        &current_location,
                        &current_qualifiers,
                        &id_attributes,
                        feature_count,
                    ) {
                        features.push(f);
                    }
                }
                break;
            }

            if !in_features {
                i += 1;
                continue;
            }

            // Feature line starts with 5 spaces then feature type
            if line.len() > 5 && &line[..5] == "     " {
                let rest = &line[5..];

                // New feature (type starts at column 5, location after spaces)
                if !rest.starts_with(' ') && !rest.starts_with('/') {
                    // Save previous feature
                    if !current_feature_type.is_empty() && current_feature_type == feature_type {
                        if let Some(f) = parse_genbank_feature(
                            &seqname,
                            &current_location,
                            &current_qualifiers,
                            &id_attributes,
                            feature_count,
                        ) {
                            features.push(f);
                        }
                    }

                    // Parse new feature type and location
                    let parts: Vec<&str> = rest.split_whitespace().collect();
                    if parts.len() >= 2 {
                        current_feature_type = parts[0].to_string();
                        current_location = parts[1..].join("");
                        current_qualifiers.clear();
                        feature_count += 1;
                    }
                }
                // Qualifier line (starts with /)
                else if rest.trim().starts_with('/') {
                    let qual_line = rest.trim();
                    parse_genbank_qualifier(qual_line, &mut current_qualifiers);
                }
                // Continuation of location or qualifier
                else {
                    let trimmed = rest.trim();
                    // Check if it's a qualifier continuation (has quotes or is inside a qualifier value)
                    if trimmed.starts_with('/') {
                        parse_genbank_qualifier(trimmed, &mut current_qualifiers);
                    } else if current_location.contains('(') && !current_location.contains(')') {
                        // Location continuation
                        current_location.push_str(trimmed);
                    }
                }
            }

            i += 1;
        }

        eprintln!(
            "Parsed {} {} features from GenBank",
            features.len(),
            feature_type
        );

        Ok(features)
    }
}

/// Parse GFF/GTF attributes field
fn parse_attributes(attr_string: &str) -> FxHashMap<String, String> {
    let mut attrs = FxHashMap::default();

    // GFF3 format: key=value;key=value
    // GTF format: key "value"; key "value"

    if attr_string.contains('=') {
        // GFF3 format
        for part in attr_string.split(';') {
            let part = part.trim();
            if let Some((key, value)) = part.split_once('=') {
                let value = value.trim_matches('"');
                attrs.insert(key.to_string(), value.to_string());
            }
        }
    } else {
        // GTF format
        for part in attr_string.split(';') {
            let part = part.trim();
            if part.is_empty() {
                continue;
            }
            let mut tokens = part.splitn(2, ' ');
            if let (Some(key), Some(value)) = (tokens.next(), tokens.next()) {
                let value = value.trim().trim_matches('"');
                attrs.insert(key.to_string(), value.to_string());
            }
        }
    }

    attrs
}

/// Parse GenBank qualifier line (e.g., /locus_tag="b0001")
fn parse_genbank_qualifier(line: &str, qualifiers: &mut FxHashMap<String, String>) {
    let line = line.trim_start_matches('/');

    if let Some((key, value)) = line.split_once('=') {
        let value = value.trim_matches('"');
        qualifiers.insert(key.to_string(), value.to_string());
    } else {
        // Flag-style qualifier (e.g., /pseudo)
        qualifiers.insert(line.to_string(), "true".to_string());
    }
}

/// Parse GenBank location and create a Feature
fn parse_genbank_feature(
    seqname: &str,
    location: &str,
    qualifiers: &FxHashMap<String, String>,
    id_attributes: &[&str],
    feature_count: usize,
) -> Option<Feature> {
    // Skip pseudo genes (these are separate "pseudogene" type in GFF)
    if qualifiers.contains_key("pseudo") || qualifiers.contains_key("pseudogene") {
        return None;
    }

    // Parse location: can be "190..255", "complement(190..255)", "join(190..255,300..400)"
    let (start, end, strand) = parse_genbank_location(location)?;

    // Find feature ID from qualifiers
    let feature_id = id_attributes
        .iter()
        .filter_map(|attr| qualifiers.get(*attr))
        .next()
        .cloned()
        .or_else(|| {
            // Fallback: try common attributes
            let fallback_attrs = ["locus_tag", "gene", "protein_id", "db_xref"];
            fallback_attrs
                .iter()
                .filter_map(|attr| qualifiers.get(*attr))
                .next()
                .cloned()
        })
        .unwrap_or_else(|| format!("feature_{}", feature_count));

    Some(Feature {
        seqname: seqname.to_string(),
        start,
        end,
        strand,
        feature_id,
    })
}

/// Parse GenBank location string
/// Returns (start, end, strand)
fn parse_genbank_location(location: &str) -> Option<(u64, u64, Strand)> {
    let location = location.trim();

    // Check for complement
    let (inner, is_complement) = if location.starts_with("complement(") && location.ends_with(')') {
        (&location[11..location.len() - 1], true)
    } else {
        (location, false)
    };

    // Handle join() - for simplicity, use the overall range
    let range_str = if inner.starts_with("join(") && inner.ends_with(')') {
        &inner[5..inner.len() - 1]
    } else {
        inner
    };

    // Parse the range(s) - get min start and max end
    let mut min_start = u64::MAX;
    let mut max_end = 0u64;

    for part in range_str.split(',') {
        let part = part.trim();
        // Handle "190..255" or "<190..255" or "190..>255"
        let part = part.trim_start_matches('<').trim_start_matches('>');

        if let Some((start_str, end_str)) = part.split_once("..") {
            let start_str = start_str.trim().trim_start_matches('<').trim_start_matches('>');
            let end_str = end_str.trim().trim_start_matches('<').trim_start_matches('>');

            if let (Ok(start), Ok(end)) = (start_str.parse::<u64>(), end_str.parse::<u64>()) {
                min_start = min_start.min(start);
                max_end = max_end.max(end);
            }
        } else if let Ok(pos) = part.parse::<u64>() {
            // Single position
            min_start = min_start.min(pos);
            max_end = max_end.max(pos);
        }
    }

    if min_start == u64::MAX || max_end == 0 {
        return None;
    }

    let strand = if is_complement {
        Strand::Minus
    } else {
        Strand::Plus
    };

    Some((min_start, max_end, strand))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_gff3_attributes() {
        let attr = "ID=gene-b0001;Name=thrL;locus_tag=b0001;gene_biotype=protein_coding";
        let parsed = parse_attributes(attr);
        assert_eq!(parsed.get("ID"), Some(&"gene-b0001".to_string()));
        assert_eq!(parsed.get("locus_tag"), Some(&"b0001".to_string()));
    }

    #[test]
    fn test_parse_gtf_attributes() {
        let attr = r#"gene_id "b0001"; gene_name "thrL"; locus_tag "b0001";"#;
        let parsed = parse_attributes(attr);
        assert_eq!(parsed.get("gene_id"), Some(&"b0001".to_string()));
        assert_eq!(parsed.get("locus_tag"), Some(&"b0001".to_string()));
    }
}
