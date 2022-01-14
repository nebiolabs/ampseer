#[cfg(test)]
extern crate assert_cmd;
extern crate predicates;

use assert_cmd::prelude::*;
use predicates::prelude::*;

use std::process::Command;

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
fn test_hazard_exit_code() {
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    cmd.arg("--help").assert().code(0);
}

#[test]
fn test_insufficient_arguments() {
    let expected_output_predicate = predicate::function(|x: &str| {
        return x.contains("required arguments");
    });
    let mut cmd = Command::cargo_bin("ampseer").expect("Calling binary failed");
    cmd.arg("--reads").arg("y").assert().stderr(expected_output_predicate);
    //cmd.arg("--primer-sets").arg("x").assert().stdout(expected_output_predicate);
}
