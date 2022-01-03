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
};

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

    let mut reader = File::open(args.primer_sets.as_path())
        .map(BufReader::new)
        .map(noodles::fasta::Reader::new)
        .with_context(|| anyhow!("Failed to open primer_sets file: {:?}", args.primer_sets.as_path()))?;

    for result in reader.records() {
        let record = result?;
        println!("{}\t{}", record.name(), record.sequence().len());
    }
    // sniff the input file format
    //import primer schemes to hash table
    // extract beginning bases of R1 and R2

    //count in batches of 1000 reads
    // let lines = io::stdin().lines();
    // for line in lines {
    //     println!("got a line: {}", line.unwrap());
    // }
    Ok(())
}