#[cfg(not(debug_assertions))]
use human_panic::setup_panic;

#[cfg(debug_assertions)]
extern crate better_panic;

use anyhow::{anyhow, Context, Result};
use clap::Parser;
use debruijn::dna_string::*;
use debruijn::{kmer::Kmer16, Mer, Vmer};
use simple_logger::SimpleLogger;
use std::cmp::Ordering;
use std::{
    collections::hash_map::Entry, collections::HashMap, collections::HashSet, fs::File,
    io::BufReader, path::Path, path::PathBuf,
};

#[derive(Parser)]
#[clap(author, version, about)]
struct Cli {
    /// File containing reads to examine (or /dev/stdin)
    #[clap(short, long, parse(from_os_str), value_name = "FILE", required = true)]
    reads: PathBuf,

    /// Files containing primer sets to check against
    #[clap(
        short,
        long,
        parse(from_os_str),
        value_name = "FILE",
        multiple_values = true,
        required = true
    )]
    primer_sets: Vec<PathBuf>,

    /// Increase logging verbosity with -d or -dd
    #[clap(short, long, parse(from_occurrences))]
    debug: usize,
}

struct PrimerSet {
    name: String,
    primer_counter: HashMap<Kmer16, i64>,
    num_consistent_reads: i64,
    num_inconsistent_reads: i64,
    frac_consistent: f32,
}

const EXPECTED_NON_MATCHING_RATIO: f32 = 0.005;
const DEFAULT_PRIMER_SET: &str = "unknown";
const MER_SIZE: usize = 16;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Human Panic. Only enabled when *not* debugging.
    #[cfg(not(debug_assertions))]
    {
        setup_panic!();
    }

    // Better Panic. Only enabled *when* debugging.
    #[cfg(debug_assertions)]
    {
        better_panic::Settings::debug()
            .most_recent_first(false)
            .lineno_suffix(true)
            .verbosity(better_panic::Verbosity::Full)
            .install();
    }

    let args = Cli::parse();
    //TODO: add reading directly from STDIN

    let default_log_level = if args.debug >= 2 {
        log::LevelFilter::Trace
    } else if args.debug >= 1 {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Warn
    };
    // Setup Logging, UTC timestamps to avoid platform-specific problems
    SimpleLogger::new()
        .with_level(default_log_level)
        .env()
        .with_utc_timestamps()
        .init()
        .unwrap();

    check_inputs(&args)?;

    let mut primer_set_counters = import_primer_sets(&args.primer_sets)?;

    classify_reads(&args.reads, &mut primer_set_counters)?;

    let ps_detected = identify_primer_set(&primer_set_counters);

    println!("{:?}", ps_detected);

    Ok(())
}

/// checks the passed input structure for reasonableness, printing error if necessary.
fn check_inputs(args: &Cli) -> Result<(), anyhow::Error> {
    let mut error_messages = Vec::new();

    if args.reads.exists() && args.primer_sets.iter().all(|ps| ps.exists()) {
        log::info!(
            "Searching for primers from {:?} in reads from: {:?}",
            args.primer_sets,
            args.reads
        );
    } else {
        for ps in &args.primer_sets {
            if !ps.exists() {
                error_messages.push(format!("Could not find primer sets at {:?}", ps.as_path()));
            }
        }
        if !args.reads.exists() {
            error_messages.push(format!(
                "Could not find reads at {:?}",
                args.reads.as_path()
            ));
        }
    }
    if !error_messages.is_empty() {
        log::error!("{}", error_messages.join("\n"));
        //Err(anyhow!("Invalid input:{:?}", error_messages.join("\n")));
        Err(anyhow!("Invalid input"))
    } else {
        Ok(())
    }
}

/// imports primer sets, creating primer_set_counter objects containing k-mers to search for
fn import_primer_sets(primer_set_paths: &[PathBuf]) -> Result<Vec<PrimerSet>, anyhow::Error> {
    let primer_set_counters: Vec<PrimerSet> = primer_set_paths
        .iter()
        .map(|ps_filename| {
            let primer_reader = File::open(ps_filename.as_path())
                .map(BufReader::new)
                .map(noodles::fasta::Reader::new)
                .with_context(|| {
                    anyhow!(
                        "Failed to open primer_sets file: {:?}",
                        ps_filename.as_path()
                    )
                })
                .unwrap();
            //TODO: consider an array of size 65536 instead and just index into that array
            let primer_counts: HashMap<Kmer16, i64> = primer_counts_from_file(primer_reader)
                .with_context(|| anyhow!("Failed to read records"))
                .unwrap();
            PrimerSet {
                name: ps_filename
                    .file_stem()
                    .unwrap()
                    .to_string_lossy()
                    .into_owned(),
                primer_counter: primer_counts,
                num_consistent_reads: 0,
                num_inconsistent_reads: 0,
                frac_consistent: 0.0,
            }
        })
        .collect();

    //TODO: compare all primer sets, removing ambiguous primers
    Ok(primer_set_counters)
}

/// Adds an 8-bit (16 nt) representation of the primer to the counter hash.
fn primer_counts_from_file(
    mut primer_reader: noodles::fasta::Reader<BufReader<File>>,
) -> Result<HashMap<Kmer16, i64>, anyhow::Error> {
    let mut primer_counts: HashMap<Kmer16, i64> = HashMap::new();
    for result in primer_reader.records() {
        let record = result?;
        let primer_seq = DnaString::from_acgt_bytes(record.sequence().as_ref());
        let record_name = record.name().to_lowercase();
        let key: Kmer16 = if record_name.contains("left") {
            assert!(MER_SIZE <= record.sequence().len());
            primer_seq.slice(0, MER_SIZE).get_kmer(0)
        } else if record_name.contains("right") {
            let offset = record.sequence().len() - MER_SIZE;
            assert!(offset > 0);
            primer_seq
                .slice(offset, record.sequence().len())
                .get_kmer(0)
        } else {
            return Err(anyhow!(
                "Unable to identify left/right primer from {}",
                record.name()
            ));
        };

        for key in [key, key.rc()] {
            if primer_counts.insert(key, 0).is_some() {
                log::info!("Ambiguous primer: {:?}", primer_seq);
            }
        }
    }
    Ok(primer_counts)
}

/// populates counts of primers observed in reads
fn classify_reads(
    reads_file: &Path,
    primer_set_counters: &mut Vec<PrimerSet>,
) -> Result<(), anyhow::Error> {
    let mut fastq_reader = File::open(reads_file)
        .map(BufReader::new)
        .map(noodles::fastq::Reader::new)
        .with_context(|| anyhow!("Failed to open reads file: {:?}", reads_file))?;

    for result in fastq_reader.records() {
        let record = result?;
        let read_seq = DnaString::from_acgt_bytes(record.sequence());
        if read_seq.len() < MER_SIZE {
            log::warn!("skipping short read {:?}", read_seq);
            break;
        }
        // extracts beginning and ending bases of each read
        let left_key: Kmer16 = read_seq.slice(0, MER_SIZE).get_kmer(0);
        let right_key: Kmer16 = read_seq
            .slice(read_seq.len() - MER_SIZE, read_seq.len())
            .get_kmer(0);

        //primer_set_counters.par_iter_mut().for_each(|psc| {
        for psc in &mut *primer_set_counters {
            for key in [left_key, right_key] {
                match psc.primer_counter.entry(key) {
                    Entry::Occupied(val) => {
                        *val.into_mut() += 1;
                        psc.num_consistent_reads += 1;
                    }
                    Entry::Vacant(_) => psc.num_inconsistent_reads += 1,
                }
            }
            psc.frac_consistent = psc.num_consistent_reads as f32
                / (psc.num_consistent_reads + psc.num_inconsistent_reads) as f32
        }
    }
    Ok(())
}

/// summarizes primer set observations deciding which primer set was used
fn identify_primer_set(primer_set_counters: &[PrimerSet]) -> (String, f32) {
    if primer_set_counters.is_empty() {
        return (DEFAULT_PRIMER_SET.to_string(), 0.0);
    };

    //TODO: add requirement that X fraction of primers have been observed allowing for some dropouts
    // a primer set with fewer reads classified but with all primers represented is more confident
    // than one with more raw counts on fewer primers

    //TODO: infer fragmentation or full length

    //TODO: check that ratio of left and right are similar for each primer pair

    let mut ps_fracs = primer_set_counters
        .iter()
        .map(|psc| {
            log::debug!(
                "{} consistent reads: {:?}",
                psc.name,
                psc.primer_counter
                    .iter()
                    .filter(|&(_, &count)| count > 0)
                    .collect::<HashMap<&Kmer16, &i64>>()
            );
            log::info!(
                "{} con/inconsistent reads: {}/{}",
                psc.name,
                psc.num_consistent_reads,
                psc.num_inconsistent_reads
            );
            (String::from(&psc.name), psc.frac_consistent)
        })
        .collect::<Vec<(String, f32)>>();
    //TODO add a bootstrapping confidence calculation
    ps_fracs.sort_unstable_by_key(|ps_frac| (ps_frac.1 * -1000.0) as i32);
    log::debug!("matching fractions: {:?}", ps_fracs);
    match ps_fracs.len().cmp(&1) {
        Ordering::Equal => {
            // can't use a ratio here - only one being checked
            // typical ratio for non-matching library is 0.001
            let confidence = ps_fracs[0].1 / EXPECTED_NON_MATCHING_RATIO;
            if ps_fracs[0].1 >= EXPECTED_NON_MATCHING_RATIO {
                (String::from(&ps_fracs[0].0), confidence)
            } else {
                (String::from(DEFAULT_PRIMER_SET), 0.0)
            }
        }
        Ordering::Greater => {
            let confidence = ps_fracs[0].1 / ps_fracs[1].1;
            if ps_fracs[0].1 / ps_fracs[1].1 > EXPECTED_NON_MATCHING_RATIO * 1000.0 {
                (String::from(&ps_fracs[0].0), confidence)
            } else {
                log::debug!("Resolving related primer sets");
                let (primer_set, confidence) = compare_only_unique_primers(primer_set_counters);
                (primer_set, confidence)
            }
        }
        _ => (String::from(DEFAULT_PRIMER_SET), 0.0),
    }
}

/// compares a ratio of only the unique keys to see if we can differentiate between two similar sets
fn compare_only_unique_primers(primer_set_counters: &[PrimerSet]) -> (String, f32) {
    if primer_set_counters.len() < 2 {
        (String::from(DEFAULT_PRIMER_SET), 0.0)
    } else {
        let mut top = &primer_set_counters[0];
        let mut second = &primer_set_counters[0];
        let mut max_consistent_reads = 0;
        let mut second_consistent_reads = 0;
        for psc in primer_set_counters {
            if psc.num_consistent_reads > max_consistent_reads {
                top = psc;
                max_consistent_reads = top.num_consistent_reads;
            }
        }
        for psc in primer_set_counters {
            if psc.name == top.name {
                continue;
            }
            if psc.num_consistent_reads > second_consistent_reads {
                second = psc;
                second_consistent_reads = psc.num_consistent_reads;
            }
        }

        log::debug!(
            "top primer_set {:?}({:?}), second primer_set: {:?}({:?})",
            top.name,
            top.num_consistent_reads,
            second.name,
            second.num_consistent_reads
        );

        let top_ps_keys: HashSet<Kmer16> = top.primer_counter.keys().cloned().collect();
        let second_ps_keys: HashSet<Kmer16> = second.primer_counter.keys().cloned().collect();
        let top_unique_keys = &top_ps_keys - &second_ps_keys;
        let second_unique_keys = &second_ps_keys - &top_ps_keys;
        let mut uniq_top_count = 0;
        let mut uniq_second_count = 0;
        for uniq_key in top_unique_keys {
            let top_count = top.primer_counter.get(&uniq_key).unwrap_or(&0);
            let second_count = second.primer_counter.get(&uniq_key).unwrap_or(&0);
            uniq_top_count += top_count;
            uniq_second_count += second_count;
            log::debug!(
                "top uniq_key: {:?} top/second {:?}/{:?}",
                uniq_key,
                top_count,
                second_count
            );
        }
        for uniq_key in second_unique_keys {
            let top_count = top.primer_counter.get(&uniq_key).unwrap_or(&0);
            let second_count = second.primer_counter.get(&uniq_key).unwrap_or(&0);
            uniq_top_count += top_count;
            uniq_second_count += second_count;
            log::debug!(
                "second uniq_key: {:?} top/second {:?}/{:?}",
                uniq_key,
                top_count,
                second_count
            );
        }

        let count_ratio = uniq_top_count as f32 / uniq_second_count as f32;
        log::debug!(
            "{}/{}= {}/{}",
            top.name,
            second.name,
            uniq_top_count,
            uniq_second_count
        );
        //TODO: figure out probability threshold to choose which one...
        if count_ratio > EXPECTED_NON_MATCHING_RATIO * 100.0 {
            //first set's uniq counts
            (String::from(&top.name), 0.0)
        } else if 1.0 / count_ratio > EXPECTED_NON_MATCHING_RATIO * 100.0 {
            //second set
            (String::from(&second.name), 0.0)
        } else {
            (String::from(DEFAULT_PRIMER_SET), 0.0)
        }
    }
}
