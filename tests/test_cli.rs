use assert_cmd::Command;
#[cfg(test)]
use predicates::prelude::*;
use std::path::{Path, PathBuf};

fn path_to_fixtures() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures"))
}
fn set_cwd_to_fixtures() {
    std::env::set_current_dir(Path::new(&path_to_fixtures())).unwrap();
}
fn path_to_ampseer() -> PathBuf {
    let mut bin_path = std::path::PathBuf::from(
        std::env::var("CARGO_TARGET_DIR")
            .unwrap_or_else(|_| concat!(env!("CARGO_MANIFEST_DIR"), "/target").to_string()),
    );
    #[cfg(debug_assertions)]
    bin_path.push("debug");
    #[cfg(not(debug_assertions))]
    bin_path.push("release");
    bin_path.push("ampseer");
    bin_path
}

#[test]
fn test_cli() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    cmd.assert().failure();
}

#[test]
fn test_version() {
    let expected_version = format!("ampseer {}\n", env!("CARGO_PKG_VERSION"));
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    cmd.arg("--version").assert().stdout(expected_version);
}

#[test]
fn test_help_exit_code() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    cmd.arg("--help").assert().code(0);
}

#[test]
fn test_insufficient_arguments() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    cmd.arg("--reads")
        .arg("x")
        .assert()
        .stderr(predicate::str::contains("required argument"));
}

#[test]
fn missing_fastq_file() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    cmd.arg("--primer-sets")
        .arg("primer_sets/Midnight_1200.fasta");
    cmd.arg("--reads").arg("missing.fastq");
    cmd.assert()
        .stderr(predicate::str::contains("Could not find reads"));
}

#[test]
fn test_classify_reads_from_stdin() {
    let mut cmd = Command::cargo_bin("ampseer").unwrap();
    set_cwd_to_fixtures();

    cmd.arg("--primer-sets").arg("primer_sets/neb_vss1a.fasta");

    cmd.pipe_stdin("vss.fastq")
        .unwrap()
        .assert()
        .stdout(predicate::str::contains("neb_vss1a"));
}

#[test]
fn non_matching_primer_sets() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    set_cwd_to_fixtures();

    cmd.arg("--primer-sets").arg("primer_sets/ARTIC_v3.fasta");
    cmd.arg("--reads").arg("vss.fastq");
    cmd.assert().stdout(predicate::str::contains("unknown"));
}

#[test]
fn ont_amps_find_both_orientations() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    set_cwd_to_fixtures();

    cmd.arg("--primer-sets").arg("vss_18_28.fasta");
    cmd.arg("--reads")
        .arg("ont_vss_full_length_amp18rev_amp28for.fastq");
    cmd.assert().stdout(predicate::str::contains("vss_18_28"));
}

#[test]
fn differentiate_vss_from_artic_v3() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    set_cwd_to_fixtures();

    cmd.arg("--primer-sets")
        .arg("primer_sets/ARTIC_v3.fasta")
        .arg("primer_sets/neb_vss1a.fasta");
    cmd.arg("--reads").arg("vss1a.fastq");
    cmd.assert().stdout(predicate::str::contains("neb_vss1a"));
}
#[test]
fn differentiate_artic_v3_from_vss() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    set_cwd_to_fixtures();

    cmd.arg("--primer-sets")
        .arg("primer_sets/ARTIC_v3.fasta")
        .arg("primer_sets/neb_vss2a.fasta")
        .arg("primer_sets/ARTIC_v4.fasta")
        .arg("primer_sets/Midnight_1200.fasta")
        .arg("primer_sets/neb_vsl1a.fasta")
        .arg("primer_sets/neb_vss1a.fasta");
    cmd.arg("--reads").arg("artic_v3.fastq");
    cmd.assert().stdout(predicate::str::contains("ARTIC_v3"));
}

#[test]
fn differentiate_vss2_from_vss1a() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    set_cwd_to_fixtures();

    cmd.arg("--primer-sets")
        .arg("primer_sets/ARTIC_v3.fasta")
        .arg("primer_sets/neb_vss1a.fasta")
        .arg("primer_sets/neb_vss2a.fasta")
        .arg("primer_sets/ARTIC_v4.fasta");
    cmd.arg("--reads").arg("vss2.fastq");
    cmd.assert().stdout(predicate::str::contains("neb_vss2a"));
}

#[test]
fn vss2_within_common_primer_sets_2023() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    set_cwd_to_fixtures();

    cmd.arg("--primer-sets")
        .arg("primer_sets/ARTIC_v4.fasta")
        .arg("primer_sets/neb_vss1a.fasta")
        .arg("primer_sets/neb_vss2a.fasta");
    cmd.arg("--reads").arg("vss2.fastq");
    cmd.assert().stdout(predicate::str::contains("neb_vss2a"));
}
#[test]
fn vss1_within_common_primer_sets_2023() {
    set_cwd_to_fixtures();
    let decompress = duct::cmd!("zstd", "-d", "-c", "broad_vss1a.fastq.zstd");

    let ampseer_cmd = duct::cmd!(
        path_to_ampseer(),
        "--primer-sets",
        "primer_sets/ARTIC_v3.fasta",
        "primer_sets/ARTIC_v4.fasta",
        "primer_sets/Midnight_1200.fasta",
        "primer_sets/neb_vsl1a.fasta",
        "primer_sets/neb_vss1a.fasta",
        "primer_sets/neb_vss2a.fasta"
    );

    let pipeline = decompress.pipe(ampseer_cmd);

    let output = pipeline.read().expect("Failed to execute pipeline");
    assert!(output.contains("neb_vss1a"));
}

//TODO: add ARTICv4.1, v5, v5.1, and Vss2b
