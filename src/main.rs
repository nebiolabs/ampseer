use structopt::StructOpt;
use std::io::{stdout, stdin, BufWriter};


/// This program compares fastq read files (on stdin or from first argument) with expected 
/// mulitiplex PCR amplicon start and end sequences to identify which of 
/// the specified PCR schemes was used during library preparation
#[derive(StructOpt)]
struct Cli {
    /// The path to the file to read
    #[structopt(parse(from_os_str))]
    reads: std::path::PathBuf,
    #[structopt(parse(from_os_str), short = "p", long = "primers")]
    primer_path: std::path::PathBuf,
}

fn main() {
    let args = Cli::from_args();

    println!("AmpSeer comparing: {reads} to {primers}", 
            reads=args.reads.display(), 
            primers=args.primer_path.display()
        );

    // sniff the input file format 
    //import primer schemes to hash table
    // extract beginning bases of R1 and R2 

    //count in batches of 1000 reads
    // let lines = io::stdin().lines();
    // for line in lines {
    //     println!("got a line: {}", line.unwrap());
    // }

}
