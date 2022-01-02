#[cfg(not(debug_assertions))]
use human_panic::setup_panic;

#[cfg(debug_assertions)]
extern crate better_panic;


use clap::Parser;
use simple_logger::SimpleLogger;
use anyhow::{anyhow, Context, Result};
use std::{
    path::PathBuf,
    fs::File,
    io::BufReader,
};
use noodles;



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

    if args.reads.exists() && args.primer_sets.exists() {
        log::info!(
            "Searching for primers from {:?} in reads from: {:?}",
            args.primer_sets,
            args.reads
        );
    } else {
        if !args.primer_sets.exists() {
            log::error!("Could not find primer sets at {:?}", args.primer_sets);
            //return Err();
        }
        if !args.reads.exists() {
            log::error!("Could not find reads at {:?}", args.reads);
        }
    }
    if let Some(primer_sets_filename) = args.primer_sets.to_str() {
        let mut reader = File::open(primer_sets_filename)
            .map(BufReader::new)
            .map(noodles::fasta::Reader::new)
            .with_context(|| anyhow!("Failed to open primer_sets file: {:?}", primer_sets_filename))?;

        for result in reader.records() {
            let record = result?;
            println!("{}\t{}", record.name(), record.sequence().len());
        }
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

// use std::io::{stdin, stdout, BufWriter};
// use structopt::StructOpt;

// /// This program compares fastq read files (on stdin or from first argument) with expected
// /// mulitiplex PCR amplicon start and end sequences to identify which of
// /// the specified PCR schemes was used during library preparation


// fn main() {
//     let args = Cli::from_args();

//     println!(
//         "AmpSeer comparing: {reads} to {primers}",
//         reads = args.reads.display(),
//         primers = args.primer_path.display()
//     );

// }
