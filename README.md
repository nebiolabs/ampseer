# ampseer
Ampseer examines reads in fastq format and identifies which multiplex PCR primer set was used to generate the SARS-CoV-2 sequencing library they are read from. 
It is intended to differentiate between ARTIC v3, ARTIC v4, ARTIC v4.1, VarSkip 1a, VarSkip 2, Midnight and VarSkip Long primer sets.
When compiled with --release optimizations, Ampseer processes reads at the same speed as samtools fastq ( < 4s for a 155M bam file on 2019 Macbook Pro)

## This program is not yet fully tested, it's shared now to enable commentary from the scientific community.

```
time samtools fastq tests/fixtures/vss2_large.bam \
| target/release/ampseer --reads /dev/stdin \
          --primer-sets primer_sets/*.fasta
...
samtools fastq tests/fixtures/vss2_large.bam  3.69s user 0.12s system 99% cpu 3.839 total
target/release/ampseer --reads /dev/stdin --primer-sets primer_sets/*.fasta  3.44s user 0.11s system 92% cpu 3.841 total
```

Note: Ampseer will produce "unknown" unless one primer set can be clearly separated from other candidates. It will not be able to identify differences between related sets unless both candidate sets are included. For example, ampseer will identify a ARTIC v4.1 library as ARTIC v4 unless both primer sets are included as candidates. 

Ampseer can be called like this:

## Example Commands:
### run the program: 
```sh
cargo build --release;
samtools fastq tests/fixtures/vss2_small.bam
| target/release/ampseer --reads /dev/stdin \
          --primer-sets primer_sets/*.fasta
```

### run the tests
```
cargo test
```

### make a flamegraph (--root needed on MacOS)
```sh
samtools fastq tests/fixtures/vss2_small.bam
| CARGO_PROFILE_RELEASE_DEBUG=true cargo flamegraph --root \
  -- --reads /dev/stdin \
     --primer-sets primer_sets/*.fasta
```
