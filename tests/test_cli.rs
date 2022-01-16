#[cfg(test)]
use assert_cmd;
use predicates;

use assert_cmd::prelude::*;
use predicates::{prelude::*, path};

use std::{process::Command, path::{Path}};

#[test]
fn test_cli() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    cmd.assert().failure();
}

#[test]
fn test_version() {

    let expected_version = format!("ampseer {}\n",env!("CARGO_PKG_VERSION"));
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
    let expected_output_predicate = || predicate::function(|x: &str| {
        return x.contains("required arguments");
    });
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    cmd.arg("--reads").arg("y").assert().stderr(expected_output_predicate());
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    cmd.arg("--primer-sets").arg("x").assert().stderr(expected_output_predicate());
}
#[test]
fn missing_reads() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    cmd.arg("--primer-sets").arg("midnight.bed");
    cmd.arg("--reads").arg("missing.fastq");
    cmd.assert().stderr( predicate::function(|x: &str| x.contains("Could not find reads")));
}

fn path_to_fixtures() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures"))
}

#[test]
fn non_matching_primer_sets() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");

    cmd.arg("--primer-sets").arg(path_to_fixtures().join("ARTIC_v3.bed.fasta"));
    cmd.arg("--reads").arg(path_to_fixtures().join("vss.fastq"));
    cmd.assert().stdout( predicate::function(|x: &str| x.contains("unknown")));
}
