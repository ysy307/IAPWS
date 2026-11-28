use std::{
    ffi::OsString,
    fs,
    io::{self, Write},
    path::{Path, PathBuf},
    process::{Command, ExitCode},
};

use anyhow::{Context, Result, ensure};
use clap::{Args, Parser, Subcommand, ValueEnum};
use fortran_check::{
    config::{Config, RULES},
    diagnostic::{Diagnostic, Location, sarif},
    input::{Backend, InputOptions},
    output::write_json,
    pipeline,
};
use serde::{Deserialize, Serialize};

#[derive(Parser)]
#[command(version, about)]
struct Cli {
    #[command(subcommand)]
    command: Subcommands,
}

#[derive(Subcommand)]
enum Subcommands {
    /// Check sources together as one Fortran project.
    Check(Box<Check>),
    /// Print the effective TOML policy or its JSON representation.
    Config {
        #[command(flatten)]
        policy: Policy,
        #[arg(long, value_enum, default_value = "toml")]
        format: ConfigFormat,
    },
    /// List native rules, their phase and default severity.
    Rules {
        #[arg(long, value_enum, default_value = "text")]
        format: ListFormat,
    },
    /// Print strict compiler options for a CMake toolchain.
    CompilerFlags {
        #[command(flatten)]
        policy: Policy,
        #[arg(long)]
        compiler_id: Backend,
        #[arg(long)]
        windows: bool,
        #[arg(long, value_enum, default_value = "text")]
        format: ListFormat,
    },
    /// Invoke a compiler after a successful project check.
    Launch {
        #[arg(long)]
        gate: PathBuf,
        #[arg(last = true, required = true, num_args = 1.., value_name = "COMPILER_AND_ARGS")]
        command: Vec<OsString>,
    },
}

#[derive(Clone, Copy, ValueEnum)]
enum Profile {
    Generic,
    Ftcms,
}

#[derive(Clone, Copy, ValueEnum)]
enum ConfigFormat {
    Toml,
    Json,
}

#[derive(Clone, Copy, ValueEnum)]
enum ListFormat {
    Text,
    Json,
}

#[derive(Clone, Copy, ValueEnum)]
enum DiagnosticFormat {
    Gcc,
    Json,
}

#[derive(Args)]
struct Policy {
    /// Project root for sources, policy and path patterns.
    #[arg(long, default_value = ".")]
    root: PathBuf,
    /// Built-in settings to apply before the local policy.
    #[arg(long, value_enum, default_value = "generic")]
    profile: Profile,
    /// Local TOML file containing policy overrides.
    #[arg(long, visible_alias = "policy")]
    config: Option<PathBuf>,
}

impl Policy {
    fn load(&self, root: &Path) -> Result<Config> {
        let path = self.config.as_ref().map(|path| root.join(path));
        let profile = match self.profile {
            Profile::Generic => "generic",
            Profile::Ftcms => "ftcms",
        };
        Config::load(profile, path.as_deref())
    }
}

#[derive(Args)]
struct Check {
    #[command(flatten)]
    policy: Policy,
    /// Files and directories to analyze. With no paths, scan the project root.
    #[arg(value_name = "SOURCE", conflicts_with = "compile_db")]
    paths: Vec<PathBuf>,
    /// Use native preprocessing from a CMake compilation database.
    #[arg(long, requires = "compiler_id")]
    compile_db: Option<PathBuf>,
    /// GNU (gfortran) or IntelLLVM (ifx).
    #[arg(long)]
    compiler_id: Option<Backend>,
    /// Use this native preprocessor for directly specified files or Fypp output.
    #[arg(long, requires = "compiler_id", conflicts_with = "compile_db")]
    cpp_executable: Option<PathBuf>,
    /// Supply a native preprocessor option; repeat for multiple options.
    #[arg(long, requires = "cpp_executable", allow_hyphen_values = true)]
    cpp_arg: Vec<String>,
    /// Write SARIF, metrics JSON and dependency DOT here.
    #[arg(long)]
    reports_dir: Option<PathBuf>,
    /// Publish the compiler launch status here.
    #[arg(long)]
    gate: Option<PathBuf>,
    /// Use this directory for the incremental cache.
    #[arg(long, conflicts_with = "no_cache")]
    cache_dir: Option<PathBuf>,
    #[arg(long)]
    no_cache: bool,
    #[arg(long, value_parser = clap::value_parser!(u32).range(1..))]
    jobs: Option<u32>,
    /// Enable Fortitude with the executable at this path or on PATH.
    #[arg(long, conflicts_with = "no_fortitude")]
    fortitude: Option<PathBuf>,
    #[arg(long)]
    no_fortitude: bool,
    #[arg(long, value_enum, default_value = "gcc")]
    format: DiagnosticFormat,
}

#[derive(Serialize, Deserialize)]
struct Gate {
    protocol: u8,
    passed: bool,
    sources: usize,
}

fn json(value: &impl Serialize) -> Result<()> {
    let mut stdout = io::stdout().lock();
    serde_json::to_writer_pretty(&mut stdout, value)?;
    writeln!(stdout)?;
    Ok(())
}

fn print_error(diagnostic: &Diagnostic, format: DiagnosticFormat) -> Result<()> {
    match format {
        DiagnosticFormat::Gcc => eprintln!("{}", diagnostic.gcc()),
        DiagnosticFormat::Json => json(&serde_json::json!({
            "complete": false,
            "diagnostics": [diagnostic],
        }))?,
    }
    Ok(())
}

fn check(args: Check) -> Result<u8> {
    let gate = args.gate.as_ref().map(|path| args.policy.root.join(path));
    if let Some(path) = &gate {
        write_json(
            path,
            &Gate {
                protocol: 1,
                passed: false,
                sources: 0,
            },
        )?;
    }
    let root = fs::canonicalize(&args.policy.root)
        .with_context(|| format!("cannot open project root {}", args.policy.root.display()))?;
    ensure!(root.is_dir(), "project root must be a directory");
    let reports = args.reports_dir.as_ref().map(|path| root.join(path));
    let mut write_error_sarif = true;
    let result = (|| -> Result<u8> {
        let mut config = args.policy.load(&root)?;
        write_error_sarif = config.reports.include_sarif;
        if let Some(jobs) = args.jobs {
            config.runtime.jobs = jobs as usize;
        }
        if args.no_cache {
            config.runtime.cache = false;
        }
        if let Some(executable) = args.fortitude {
            config.fortitude.enabled = true;
            config.fortitude.executable = executable;
        }
        if args.no_fortitude {
            config.fortitude.enabled = false;
        }
        let cache = args
            .cache_dir
            .or_else(|| config.runtime.cache_directory.clone())
            .map(|path| root.join(path))
            .or_else(|| reports.as_ref().map(|path| path.join("cache")));
        let inputs = InputOptions {
            paths: args.paths,
            compile_db: args.compile_db,
            backend: args.compiler_id,
            cpp_executable: args.cpp_executable,
            cpp_args: args.cpp_arg,
        };
        let report = pipeline::run(&root, &inputs, &config, cache.as_deref())?;
        if let Some(directory) = &reports {
            report.write(directory, &root, &config)?;
        }
        let code = if !report.complete {
            2
        } else if report.errors() > 0 {
            1
        } else {
            0
        };
        match args.format {
            DiagnosticFormat::Gcc => {
                for diagnostic in &report.diagnostics {
                    eprintln!("{}", diagnostic.gcc());
                }
                println!(
                    "FortranCheck: {} sources, {} errors, {} warnings{}{}",
                    report.sources,
                    report.errors(),
                    report.warnings(),
                    if report.cached { " (cached)" } else { "" },
                    if report.complete {
                        ""
                    } else {
                        " (incomplete analysis)"
                    },
                );
            }
            DiagnosticFormat::Json => json(&report)?,
        }
        if let Some(path) = &gate {
            write_json(
                path,
                &Gate {
                    protocol: 1,
                    passed: code == 0,
                    sources: report.sources,
                },
            )?;
        }
        Ok(code)
    })();
    match result {
        Ok(code) => Ok(code),
        Err(error) => {
            let source = args
                .policy
                .config
                .as_ref()
                .map_or_else(|| root.clone(), |path| root.join(path));
            let diagnostic =
                Diagnostic::new("TOOL001", format!("{error:#}"), Location::new(source, 1, 1));
            if let Some(directory) = reports
                && write_error_sarif
            {
                write_json(
                    &directory.join("diagnostics.sarif"),
                    &sarif(std::slice::from_ref(&diagnostic), &root),
                )
                .context("could not write failure diagnostics")?;
            }
            print_error(&diagnostic, args.format)?;
            Ok(2)
        }
    }
}

fn launch(gate: &Path, command: &[OsString]) -> Result<u8> {
    let gate: Gate = serde_json::from_slice(
        &fs::read(gate).with_context(|| format!("no successful check at {}", gate.display()))?,
    )
    .with_context(|| format!("invalid check status at {}", gate.display()))?;
    ensure!(
        gate.protocol == 1 && gate.passed,
        "Fortran compilation requires a successful project check"
    );
    let status = Command::new(&command[0])
        .args(&command[1..])
        .status()
        .context("could not start the compiler")?;
    Ok(status
        .code()
        .and_then(|code| u8::try_from(code).ok())
        .unwrap_or(2))
}

fn run(command: Subcommands) -> Result<u8> {
    match command {
        Subcommands::Check(args) => check(*args),
        Subcommands::Config { policy, format } => {
            let config = policy.load(&fs::canonicalize(&policy.root)?)?;
            match format {
                ConfigFormat::Toml => print!("{}", toml::to_string_pretty(&config)?),
                ConfigFormat::Json => json(&config)?,
            }
            Ok(0)
        }
        Subcommands::Rules { format } => {
            match format {
                ListFormat::Json => json(&RULES)?,
                ListFormat::Text => {
                    for rule in RULES {
                        println!(
                            "{:<8} {:>2} {:<12} {}",
                            rule.id, rule.phase, rule.default_severity, rule.name
                        );
                    }
                }
            }
            Ok(0)
        }
        Subcommands::CompilerFlags {
            policy,
            compiler_id,
            windows,
            format,
        } => {
            let config = policy.load(&fs::canonicalize(&policy.root)?)?;
            let flags = compiler_id.flags(&config, windows || cfg!(windows));
            match format {
                ListFormat::Json => json(&flags)?,
                ListFormat::Text => {
                    for flag in flags {
                        println!("{flag}");
                    }
                }
            }
            Ok(0)
        }
        Subcommands::Launch { gate, command } => launch(&gate, &command),
    }
}

fn main() -> ExitCode {
    match run(Cli::parse().command) {
        Ok(code) => ExitCode::from(code),
        Err(error) => {
            eprintln!("fortran-check: error [TOOL001]: {error:#}");
            ExitCode::from(2)
        }
    }
}
