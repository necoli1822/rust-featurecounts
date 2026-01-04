//! Simple CLI argument parsing without clap

use std::env;
use std::path::PathBuf;

#[derive(Debug)]
pub struct Args {
    pub bam_files: Vec<PathBuf>,
    pub annotation: PathBuf,
    pub output: PathBuf,
    pub feature_type: String,
    pub attribute: String,
    pub stranded: String,
    pub min_mapq: u8,
    pub paired: bool,
    pub threads: i32,
    pub counting_mode: String,
}

impl Args {
    pub fn parse() -> Result<Self, String> {
        let args: Vec<String> = env::args().collect();

        let mut annotation: Option<PathBuf> = None;
        let mut output: Option<PathBuf> = None;
        let mut feature_type = "exon".to_string();
        let mut attribute = "gene_id".to_string();
        let mut stranded = "0".to_string();
        let mut min_mapq: u8 = 0;
        let mut paired = false;
        let mut threads: i32 = 1;
        let mut counting_mode = "intersection".to_string();
        let mut bam_files: Vec<PathBuf> = Vec::new();

        let mut i = 1;
        while i < args.len() {
            match args[i].as_str() {
                "-a" => {
                    i += 1;
                    if i < args.len() {
                        annotation = Some(PathBuf::from(&args[i]));
                    }
                }
                "-o" => {
                    i += 1;
                    if i < args.len() {
                        output = Some(PathBuf::from(&args[i]));
                    }
                }
                "-t" => {
                    i += 1;
                    if i < args.len() {
                        feature_type = args[i].clone();
                    }
                }
                "-g" => {
                    i += 1;
                    if i < args.len() {
                        attribute = args[i].clone();
                    }
                }
                "-s" => {
                    i += 1;
                    if i < args.len() {
                        stranded = args[i].clone();
                    }
                }
                "-Q" => {
                    i += 1;
                    if i < args.len() {
                        min_mapq = args[i].parse().unwrap_or(0);
                    }
                }
                "-p" => {
                    paired = true;
                }
                "-T" => {
                    i += 1;
                    if i < args.len() {
                        threads = args[i].parse().unwrap_or(1);
                    }
                }
                "-m" | "--mode" => {
                    i += 1;
                    if i < args.len() {
                        counting_mode = args[i].clone();
                    }
                }
                "-h" | "--help" => {
                    print_help();
                    std::process::exit(0);
                }
                "-v" | "--version" => {
                    println!("rust-featurecounts v0.1.0");
                    std::process::exit(0);
                }
                arg if !arg.starts_with('-') => {
                    bam_files.push(PathBuf::from(arg));
                }
                _ => {
                    // Ignore unknown options for compatibility
                }
            }
            i += 1;
        }

        let annotation = annotation.ok_or("Missing required argument: -a <annotation_file>")?;
        let output = output.ok_or("Missing required argument: -o <output_file>")?;

        if bam_files.is_empty() {
            return Err("No input BAM files provided".to_string());
        }

        Ok(Self {
            bam_files,
            annotation,
            output,
            feature_type,
            attribute,
            stranded,
            min_mapq,
            paired,
            threads,
            counting_mode,
        })
    }
}

fn print_help() {
    println!(r#"rust-featurecounts v0.1.0 - A fast feature counting tool

Usage: rust-featurecounts [options] -a <annotation> -o <output> input.bam [input2.bam ...]

Required arguments:
  -a <file>       Annotation file (GFF/GTF/GenBank format)
                  Supported: .gff, .gff3, .gtf, .gbk, .gb, .genbank, .gbff
  -o <file>       Output file for count matrix (featureCounts compatible)

Optional arguments:
  -t <string>     Feature type to count [default: exon]
                  Options: exon, gene, CDS, auto, or any feature type
  -g <string>     Attribute for feature ID [default: gene_id]
                  GenBank: locus_tag, gene, protein_id
  -s <int>        Strandedness: 0=unstranded, 1=stranded, 2=reversely stranded [default: 0]
  -Q <int>        Minimum mapping quality [default: 0]
  -p              Paired-end mode (counts fragments)
  -T <int>        Number of threads [default: 1, 0=auto based on load]
  -m <string>     Counting mode [default: intersection]
                  Options: intersection, union
  -h, --help      Show this help message
  -v, --version   Show version
"#);
}
