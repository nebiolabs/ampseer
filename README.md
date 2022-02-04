# ampseer
Ampseer examines reads in fastq format and identifies which multiplex PCR primer set was used to generate the SARS-CoV-2 sequencing library they are read from. 
It is intended to differentiate between ARTIC v3, ARTIC v4, ARTIC v4.1, VarSkip 1a, VarSkip 2a, Midnight, and VarSkip Long primer sets.
## This program is not yet fully tested, it's shared now to enable commentary from the scientific community.

Pull requests and issues are welcome.

When compiled with --release optimizations, Ampseer processes reads at the same speed as samtools fastq ( less than 4s for a 155M bam file on 2019 Macbook Pro)

```sh
time samtools fastq tests/fixtures/vss2_large.bam \
| target/release/ampseer --reads /dev/stdin \
 --primer-sets primer_sets/*.fasta
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 4348571 reads
("neb_vss2a", 0.0)
samtools fastq tests/fixtures/vss2_large.bam  3.55s user 0.11s system 99% cpu 3.662 total
target/release/ampseer --reads /dev/stdin --primer-sets primer_sets/*.fasta  3.40s user 0.09s system 95% cpu 3.661 total
```

Note: Ampseer will produce "unknown" unless one primer set can be clearly separated from other candidates. It will not be able to identify differences between related sets unless both candidate sets are included. For example, ampseer will identify a ARTIC v4.1 library as ARTIC v4 unless both primer sets are included as candidates.

## Example Commands:
This tool does not yet have any binary releases. To try it, you will need to [install rustup](https://www.rust-lang.org/tools/install), or `rustup update` if you are using an older rust installation.

### run the program: 
```sh
cargo build --release
samtools fastq tests/fixtures/vss2_small.bam \
| target/release/ampseer --reads /dev/stdin \
          --primer-sets primer_sets/*.fasta
```
### view ampseer help:
```sh
cargo build --release # may take some time to compile the first time
target/release/ampseer -h
```

### run the tests:
```sh
cargo test # may take some time to compile the first time
```

### make a flamegraph (--root needed on MacOS):
```sh
samtools fastq tests/fixtures/vss2_small.bam \
| CARGO_PROFILE_RELEASE_DEBUG=true cargo flamegraph --root \
  -- --reads /dev/stdin \
     --primer-sets primer_sets/*.fasta
```

## How does it work?
Ampseer examines the ends of reads in 16 bp chunks and compares these with the expected sequences specified in the primer sets to identify which primer set is most consistent with the observed reads.
