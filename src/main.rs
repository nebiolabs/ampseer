#[cfg(not(debug_assertions))]
use human_panic::setup_panic;

#[cfg(debug_assertions)]
extern crate better_panic;

// use std::io::{stdin, stdout, BufWriter}; TODO add support for reads on stdin
use anyhow::{anyhow, Context, Result};
use clap::Parser;
use debruijn::dna_string::*;
use debruijn::{kmer::Kmer16, Mer, Vmer};
use simple_logger::SimpleLogger;
use std::{
    collections::hash_map::Entry, collections::HashMap, fs::File, io::BufReader, path::PathBuf,
};

#[derive(Parser)]
#[clap(author, version, about)]
struct Cli {
    /// file containing reads to examine
    #[clap(short, long, parse(from_os_str), value_name = "FILE", required = true)]
    reads: PathBuf,

    /// files containing primer sets to check against
    #[clap(
        short,
        long,
        parse(from_os_str),
        value_name = "FILE",
        multiple_values = true,
        required = true
    )]
    primer_sets: Vec<PathBuf>,

    /// Turn debugging information on
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

    // Setup Logging, UTC timestamps to avoid platform-specific problems
    SimpleLogger::new().with_utc_timestamps().init().unwrap();

    let args = Cli::parse();
    check_inputs(&args)?;

    let mut primer_set_counters = import_primer_sets(&args.primer_sets)?;

    classify_reads(&args.reads, &mut primer_set_counters)?;

    let ps_detected = identify_primer_set(primer_set_counters);

    println!("{}", ps_detected);

    Ok(())
}

/// populates counts of primers observed in reads
fn classify_reads(
    reads_file: &PathBuf,
    primer_set_counters: &mut Vec<PrimerSet>,
) -> Result<(), anyhow::Error> {
    let mut fastq_reader = File::open(reads_file.as_path())
        .map(BufReader::new)
        .map(noodles::fastq::Reader::new)
        .with_context(|| {
            anyhow!(
                "Failed to open reads file: {:?}",
                reads_file.as_path()
            )
        })?;
    Ok(for result in fastq_reader.records() {
        let record = result?;
        let read_seq = DnaString::from_acgt_bytes(record.sequence().as_ref());
        if read_seq.len() < 16 {
            log::warn!("skipping short read {:?}", read_seq);
            break;
        }
        // extracts beginning and ending bases of each read
        let left_key: Kmer16 = read_seq.slice(0, 16).get_kmer(0);
        let right_key: Kmer16 = read_seq
            .slice(read_seq.len() - 16, read_seq.len())
            .get_kmer(0);

        for psc in &mut *primer_set_counters {
            for key in [left_key, right_key] {
                match psc.primer_counter.entry(key) {
                    Entry::Occupied(val) => {
                        *val.into_mut() += 1;
                        psc.num_consistent_reads += 1
                    }
                    Entry::Vacant(_) => psc.num_inconsistent_reads += 1,
                }
                psc.frac_consistent = psc.num_consistent_reads as f32
                    / (psc.num_consistent_reads + psc.num_inconsistent_reads) as f32
            }
        }
    })
}

/// imports primer sets, creating primer_set_counter objects containing k-mers to search for
fn import_primer_sets(primer_set_paths: &Vec<PathBuf>) -> Result<Vec<PrimerSet>, anyhow::Error> {
    let mut primer_set_counters: Vec<PrimerSet> = Vec::new();
    for ps_filename in primer_set_paths {
        let mut primer_reader = File::open(ps_filename.as_path())
            .map(BufReader::new)
            .map(noodles::fasta::Reader::new)
            .with_context(|| {
                anyhow!(
                    "Failed to open primer_sets file: {:?}",
                    ps_filename.as_path()
                )
            })?;

        //TODO: consider an array of size 65536 instead and just index into that array
        let mut primer_counts: HashMap<Kmer16, i64> = HashMap::new();
        for result in primer_reader.records() {
            let record = result?;
            populate_primer_count_hash(&mut primer_counts, record)?;
        }
        primer_set_counters.push(PrimerSet {
            name: ps_filename
                .file_stem()
                .unwrap()
                .to_string_lossy()
                .into_owned(),
            primer_counter: primer_counts,
            num_consistent_reads: 0,
            num_inconsistent_reads: 0,
            frac_consistent: 0.0,
        });
    }
    //TODO: compare all primer sets, removing ambiguous primers
    Ok(primer_set_counters)
}

/// summarizes primer set observations deciding which primer set was used
fn identify_primer_set(primer_set_counters: Vec<PrimerSet>) -> String {
    //TODO: add requirement that X fraction of primers have been observed allowing for some dropouts
    // a primer set with fewer reads classified but with all primers represented is more confident than one with more raw counts.

    //TODO: infer fragmentation or full length

    //TODO: check that ratio of left and right are similar for each primer pair

    let mut ps_fracs: Vec<(String, f32)> = Vec::new();
    for psc in primer_set_counters {
        log::debug!(
            "{} consistent reads: {:?}",
            psc.name,
            psc.primer_counter
                .iter()
                .filter(|&(_, &count)| count > 0)
                .collect::<HashMap<&Kmer16, &i64>>()
        );
        log::info!(
            "{} inconsistent reads: {}",
            psc.name,
            psc.num_inconsistent_reads
        );
        ps_fracs.push((psc.name, psc.frac_consistent));
    }

    ps_fracs.sort_unstable_by_key(|ps_frac| (ps_frac.1 * -1000.0) as i32);
    let ps_detected = if ps_fracs[0].1 > 0.0 {
        if ps_fracs.len() == 1 || ps_fracs[1].1 <= 0.0 {
            &ps_fracs[0].0
        } else if ps_fracs.len() == 1 || ps_fracs[0].1 / ps_fracs[1].1 > 5.0 {
            &ps_fracs[0].0
        } else {
            "unknown"
        }
    } else {
        "unknown"
    };
    ps_detected.to_string()
}

/// Adds an 8-bit (16 nt) representation of the primer to the counter hash.
fn populate_primer_count_hash(
    primer_counts: &mut HashMap<debruijn::kmer::IntKmer<u32>, i64>,
    record: noodles::fasta::Record,
) -> Result<(), anyhow::Error> {
    let key_length: usize = 16;
    let key: Kmer16;
    let primer_seq = DnaString::from_acgt_bytes(record.sequence().as_ref());
    if record.name().to_lowercase().contains("left") {
        assert!(key_length <= record.sequence().len());
        key = primer_seq.slice(0, key_length).get_kmer(0);
        log::debug!(
            "Adding left key: {:?} from {:?}:{:?}",
            key,
            record.name(),
            primer_seq
        );
    } else if record.name().to_lowercase().contains("right") {
        let offset = record.sequence().len() - key_length;
        assert!(offset > 0);
        key = primer_seq
            .slice(offset, record.sequence().len())
            .get_kmer(0);
        log::debug!(
            "Adding right key: {:?} from {:?}:{:?}",
            key,
            record.name(),
            primer_seq
        );
    } else {
        return Err(anyhow!(
            "Unable to identify left/right primer from {}",
            record.name()
        ));
    }
    for key in [key, key.rc()] {
        if primer_counts.contains_key(&key) {
            log::warn!("Ambiguous primer: {:?}", primer_seq);
        } else {
            primer_counts.insert(key, 0);
        }
    }
    Ok(())
}
