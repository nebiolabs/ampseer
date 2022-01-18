#[cfg(not(debug_assertions))]
use human_panic::setup_panic;

#[cfg(debug_assertions)]
extern crate better_panic;

// use std::io::{stdin, stdout, BufWriter}; TODO add support for reads on stdin
use anyhow::{anyhow, Context, Result};
use clap::Parser;
use debruijn::dna_string::*;
use debruijn::kmer::Kmer16;
use debruijn::Vmer;
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
    reads_classified: i64,
    reads_unclassified: i64,
}

/// checks the passed input structure for reasonableness, printing error if necessary.
fn check_inputs(args: &Cli) -> Result<(), anyhow::Error> {
    let mut error_messages = Vec::new();

    if args.reads.exists()
        && args.primer_sets.iter().all(|ps| ps.exists())
    {
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

    let mut primer_set_counters: Vec<PrimerSet> = Vec::new();
    for ps_filename in args.primer_sets {
        let mut primer_reader = File::open(ps_filename.as_path())
            .map(BufReader::new)
            .map(noodles::fasta::Reader::new)
            .with_context(|| {
                anyhow!(
                    "Failed to open primer_sets file: {:?}",
                    ps_filename.as_path()
                )
            })?;

        //TODO: maybe consider an array of size 65536 instead and just index into that array
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
            reads_classified: 0,
            reads_unclassified: 0,
        });
    }
    //TODO: compare all primer sets, remove ambiguous primers
    //TODO: count up all other start/ends reporting the top N - compare to generate confidence estimate
    //TODO: infer fragmentation or full length
    //TODO: check that ratio of left and right are similar for each primer pair

    // extract beginning bases of R1 and R2
    let mut fastq_reader = File::open(args.reads.as_path())
        .map(BufReader::new)
        .map(noodles::fastq::Reader::new)
        .with_context(|| {
            anyhow!(
                "Failed to open primer_sets file: {:?}",
                args.reads.as_path()
            )
        })?;
    for result in fastq_reader.records() {
        let record = result?;
        let read_seq = DnaString::from_acgt_bytes(record.sequence().as_ref());
        let key: Kmer16 = read_seq.slice(0, 16).get_kmer(0);
        //TODO consider last 16 bases too - for full-length inserts
        for psc in &mut primer_set_counters {
            match psc.primer_counter.entry(key) {
                Entry::Occupied(val) => {
                    *val.into_mut() += 1;
                    psc.reads_classified += 1
                }
                Entry::Vacant(_) => psc.reads_unclassified += 1,
            }
        }
    }
    for psc in primer_set_counters {
        log::info!(
            "{} classified reads: {:?}",
            psc.name,
            psc.primer_counter
                .iter()
                .filter(|&(_, &count)| count > 0)
                .collect::<HashMap<&Kmer16, &i64>>()
        );
        log::info!(
            "{} unclassified reads: {}",
            psc.name,
            psc.reads_unclassified
        );
    }

    //count in batches of N reads
    // let lines = io::stdin().lines();
    // for line in lines {
    //     println!("got a line: {}", line.unwrap());
    // }
    Ok(())
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
    if primer_counts.contains_key(&key) {
        log::warn!("Ambiguous primer: {:?}", primer_seq);
    } else {
        primer_counts.insert(key, 0);
    }
    Ok(())
}
