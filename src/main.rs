#[cfg(not(debug_assertions))]
use human_panic::setup_panic;

#[cfg(debug_assertions)]
extern crate better_panic;

// use std::io::{stdin, stdout, BufWriter}; TODO add support for reads on stdin
use clap::Parser;
use simple_logger::SimpleLogger;
use anyhow::{anyhow, Context, Result};
use std::{
    path::PathBuf,
    fs::File,
    io::BufReader,
    collections::HashMap,
};
use debruijn::dna_string::*;
use debruijn::kmer::Kmer16;
use debruijn::Vmer;

#[derive(Parser)]
#[clap(author, version, about)]
struct Cli {
    /// file containing reads to examine
    #[clap(short, long, parse(from_os_str), value_name = "FILE")]
    reads: PathBuf,

    /// file containing primer sets to check against
    #[clap(short, long, parse(from_os_str), value_name = "FILE")]
    primer_sets: PathBuf,

    /// Turn debugging information on
    #[clap(short, long, parse(from_occurrences))]
    debug: usize,
}

/// checks the passed input structure for reasonableness, printing error if necessary.
fn check_inputs(args: &Cli) -> Result<(), anyhow::Error> {
    let mut error_messages = Vec::new();

    if args.reads.exists() && args.primer_sets.exists() {
        log::info!(
            "Searching for primers from {:?} in reads from: {:?}",
            args.primer_sets,
            args.reads
        );
    } else {
        if !args.primer_sets.exists() {
            error_messages.push(format!("Could not find primer sets at {:?}", args.primer_sets.as_path()));
        }
        if !args.reads.exists() {
            error_messages.push(format!("Could not find reads at {:?}", args.reads.as_path()));
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

    //TODO: add a list of multiple primer sets so we can differentiate between them
    let mut reader = File::open(args.primer_sets.as_path())
        .map(BufReader::new)
        .map(noodles::fasta::Reader::new)
        .with_context(|| anyhow!("Failed to open primer_sets file: {:?}", args.primer_sets.as_path()))?;

    let mut primer_counts:HashMap<Kmer16,i32> = HashMap::new();
    for result in reader.records() {
        let record = result?;
        popuplate_primer_count_hash(&mut primer_counts, record)?;
    }
    //TODO: compare all primer sets, remove ambiguous primers 
    //TODO: count up all other start/ends reporting the top N - compare to generate confidence estimate
    //TODO: infer fragmentation or full length
    //TODO: check that ratio of left and right are similar for each primer pair

    // extract beginning bases of R1 and R2

    //count in batches of N reads
    // let lines = io::stdin().lines();
    // for line in lines {
    //     println!("got a line: {}", line.unwrap());
    // }
    Ok(())
}

/// Adds an 8-bit representation of the primer to the counter hash.
fn popuplate_primer_count_hash(primer_counts: &mut HashMap<debruijn::kmer::IntKmer<u32>, i32>, record: noodles::fasta::Record,) -> Result<(), anyhow::Error> {
    let key_length:usize = 16;
    let key:Kmer16;
    let primer_seq = DnaString::from_acgt_bytes(record.sequence().as_ref());
    if record.name().to_lowercase().contains("left") {
        assert!(key_length <= record.sequence().len());
        key = primer_seq.slice(0, key_length).get_kmer(0); 
        log::debug!("Adding left key: {:?} from {:?}:{:?}", key, record.name(), primer_seq);
    } else if record.name().to_lowercase().contains("right") {
        let offset = record.sequence().len() - key_length;
        assert!(offset > 0);
        key = primer_seq.slice(offset,record.sequence().len()).get_kmer(0); 
        log::debug!("Adding right key: {:?} from {:?}:{:?}", key, record.name(), primer_seq);
    } else {
        return Err(anyhow!("Unable to identify left/right primer from {}",record.name()));
    }
    if primer_counts.contains_key(&key) {
        log::warn!("Ambiguous primer: {:?}", primer_seq);
    } else {
        primer_counts.insert(key,0);
    }
    Ok(())
}