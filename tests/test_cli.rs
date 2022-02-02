#[cfg(test)]
use assert_cmd;
use predicates;

use assert_cmd::prelude::*;
use predicates::prelude::*;

use std::{path::Path, process::Command};

fn path_to_fixtures() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures"))
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
    let expected_output_predicate = || {
        predicate::function(|x: &str| {
            return x.contains("required arguments");
        })
    };
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    cmd.arg("--reads")
        .arg("y")
        .assert()
        .stderr(expected_output_predicate());
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    cmd.arg("--primer-sets")
        .arg("x")
        .assert()
        .stderr(expected_output_predicate());
}
#[test]
fn missing_reads() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    cmd.arg("--primer-sets").arg("midnight.bed");
    cmd.arg("--reads").arg("missing.fastq");
    cmd.assert().stderr(predicate::function(|x: &str| {
        x.contains("Could not find reads")
    }));
}

#[test]
fn non_matching_primer_sets() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");

    cmd.arg("--primer-sets")
        .arg(path_to_fixtures().join("primer_sets/ARTIC_v3.fasta"));
    cmd.arg("--reads").arg(path_to_fixtures().join("vss.fastq"));
    cmd.assert()
        .stdout(predicate::function(|x: &str| x.contains("unknown")));
}

#[test]
fn ont_amps_find_both_orientations() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");

    cmd.arg("--primer-sets")
        .arg(path_to_fixtures().join("vss_18_28.fasta"));
    cmd.arg("--reads")
        .arg(path_to_fixtures().join("ont_vss_full_length_amp18rev_amp28for.fastq"));
    cmd.assert()
        .stdout(predicate::function(|x: &str| x.contains("vss_18_28")));
}

#[test]
fn differentiate_vss_from_artic_v3() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");

    cmd.arg("--primer-sets")
        .arg(path_to_fixtures().join("primer_sets/ARTIC_v3.fasta"))
        .arg(path_to_fixtures().join("primer_sets/neb_vss1a.fasta"));
    cmd.arg("--reads").arg(path_to_fixtures().join("vss1a.fastq"));
    cmd.assert()
        .stdout(predicate::function(|x: &str| x.contains("neb_vss1a")));
}
#[test]
fn differentiate_artic_v3_from_vss() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");

    cmd.arg("--primer-sets")
        .arg(path_to_fixtures().join("primer_sets/ARTIC_v3.fasta"))
        .arg(path_to_fixtures().join("primer_sets/neb_vss1a.fasta"));
    cmd.arg("--reads").arg(path_to_fixtures().join("artic_v3.fastq"));
    cmd.assert()
        .stdout(predicate::function(|x: &str| x.contains("ARTIC_v3")));
}

#[test]
fn differentiate_vss2_from_vss1a() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");

    cmd.arg("--primer-sets")
        .arg(path_to_fixtures().join("primer_sets/ARTIC_v3.fasta"))
        .arg(path_to_fixtures().join("primer_sets/neb_vss1a.fasta"))
        .arg(path_to_fixtures().join("primer_sets/neb_vss2a.fasta"));
    cmd.arg("--reads").arg(path_to_fixtures().join("vss2.fastq"));
    cmd.assert()
        .stdout(predicate::function(|x: &str| x.contains("neb_vss2a")));
}

//TODO: add ARTICv4.1
