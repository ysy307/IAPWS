use std::{
    collections::{BTreeMap, BTreeSet},
    ffi::OsStr,
    fs,
    path::{Path, PathBuf},
    process::Command,
    str::FromStr,
    sync::OnceLock,
};

use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
use regex::Regex;
use serde::Deserialize;
use sha2::{Digest, Sha256};
use walkdir::WalkDir;

use crate::{
    config::{self, Config},
    diagnostic::slash_path,
    source::Source,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend {
    Gnu,
    Ifx,
}

impl FromStr for Backend {
    type Err = anyhow::Error;

    fn from_str(name: &str) -> Result<Self> {
        match name.to_ascii_lowercase().as_str() {
            "gnu" | "gfortran" => Ok(Self::Gnu),
            "intel" | "intelllvm" | "ifx" => Ok(Self::Ifx),
            _ => bail!("unsupported Fortran compiler '{name}', expected GNU or IntelLLVM"),
        }
    }
}

impl Backend {
    pub fn flags(self, config: &Config, windows: bool) -> Vec<String> {
        let flags = match (self, windows) {
            (Self::Gnu, _) => &config.compiler.gnu_flags,
            (Self::Ifx, false) => &config.compiler.ifx_flags,
            (Self::Ifx, true) => &config.compiler.ifx_windows_flags,
        };
        let standard = &config.language.standard;
        let intel = match standard.as_str() {
            "f2003" => "f03",
            "f2008" => "f08",
            "f2018" => "f18",
            "f2023" => "f23",
            value => value,
        };
        flags
            .iter()
            .map(|flag| {
                flag.replace("{standard}", standard)
                    .replace("{intel_standard}", intel)
            })
            .collect()
    }

    fn cpp(self, config: &Config, windows: bool) -> &[String] {
        match (self, windows) {
            (Self::Gnu, _) => &config.compiler.gnu_cpp,
            (Self::Ifx, false) => &config.compiler.ifx_cpp,
            (Self::Ifx, true) => &config.compiler.ifx_windows_cpp,
        }
    }
}

#[derive(Debug, Default)]
pub struct InputOptions {
    pub paths: Vec<PathBuf>,
    pub compile_db: Option<PathBuf>,
    pub backend: Option<Backend>,
    pub cpp_executable: Option<PathBuf>,
    pub cpp_args: Vec<String>,
}

#[derive(Deserialize)]
struct CompileEntry {
    directory: PathBuf,
    file: PathBuf,
    arguments: Option<Vec<String>>,
    command: Option<String>,
}

#[derive(Debug)]
struct CompileJob {
    directory: PathBuf,
    file: PathBuf,
    arguments: Vec<String>,
}

fn extension(path: &Path) -> String {
    path.extension()
        .and_then(OsStr::to_str)
        .unwrap_or_default()
        .to_ascii_lowercase()
}

fn supported_extension(path: &Path, config: &Config) -> bool {
    let extension = extension(path);
    config
        .inputs
        .extensions
        .iter()
        .any(|ext| ext.trim_start_matches('.').eq_ignore_ascii_case(&extension))
}

fn source_included(root: &Path, path: &Path, config: &Config) -> bool {
    let relative = slash_path(path.strip_prefix(root).unwrap_or(path));
    config::included(&config.inputs.include, &relative)
        && !config::any_match(&config.inputs.exclude, &relative)
}

fn sources_from_paths(
    root: &Path,
    selected: &[PathBuf],
    config: &Config,
) -> Result<BTreeSet<PathBuf>> {
    let defaults = [root.to_path_buf()];
    let mut sources = BTreeSet::new();
    for entry in if selected.is_empty() {
        &defaults
    } else {
        selected
    } {
        let path = fs::canonicalize(root.join(entry))
            .with_context(|| format!("cannot open input {}", entry.display()))?;
        if path.is_file() {
            ensure!(
                supported_extension(&path, config) && source_included(root, &path, config),
                "explicit input {} was excluded or has an unsupported extension",
                path.display()
            );
            sources.insert(path);
            continue;
        }
        for item in WalkDir::new(&path)
            .follow_links(config.inputs.follow_links)
            .into_iter()
            .filter_entry(|item| {
                if item.depth() == 0 || !item.file_type().is_dir() {
                    return true;
                }
                let path = item.path();
                let relative = slash_path(path.strip_prefix(root).unwrap_or(path));
                let directory = format!("{relative}/");
                !config::any_match(&config.inputs.exclude, &directory)
                    && !config::any_match(&config.inputs.discovery_exclude, &directory)
            })
        {
            let item = item.with_context(|| format!("cannot traverse {}", path.display()))?;
            if item.file_type().is_file()
                && supported_extension(item.path(), config)
                && source_included(root, item.path(), config)
            {
                sources.insert(item.into_path());
            }
        }
    }
    Ok(sources)
}

fn read_database(
    path: &Path,
    root: &Path,
    config: &Config,
    windows: bool,
) -> Result<Vec<CompileJob>> {
    let text =
        fs::read_to_string(path).with_context(|| format!("cannot read {}", path.display()))?;
    let entries: Vec<CompileEntry> = serde_json::from_str(&text)
        .with_context(|| format!("invalid compilation database {}", path.display()))?;
    let mut jobs = Vec::new();
    for entry in entries {
        let directory = path.parent().unwrap_or(root).join(entry.directory);
        let file = directory.join(entry.file);
        let fixed = ["f", "for", "f77", "ftn"].contains(&extension(&file).as_str());
        ensure!(
            !fixed || !source_included(root, &file, config),
            "fixed-form Fortran input is unsupported: {}",
            file.display()
        );
        if !supported_extension(&file, config) {
            continue;
        }
        let file = fs::canonicalize(&file)
            .with_context(|| format!("missing generated or source file {}", file.display()))?;
        let args = match entry.arguments {
            Some(arguments) => arguments,
            None => split_command(
                &entry
                    .command
                    .context("compile command has neither arguments nor command")?,
                windows,
            )?,
        };
        ensure!(
            !args.is_empty(),
            "empty compile command for {}",
            file.display()
        );
        jobs.push(CompileJob {
            directory,
            file,
            arguments: args,
        });
    }
    Ok(jobs)
}

fn split_command(command: &str, windows: bool) -> Result<Vec<String>> {
    if !windows {
        return Ok(shell_words::split(command)?);
    }
    let mut input = command.chars().peekable();
    let mut args = Vec::new();
    while input.peek().is_some() {
        while input.peek().is_some_and(|chr| chr.is_ascii_whitespace()) {
            input.next();
        }
        if input.peek().is_none() {
            break;
        }
        let mut argument = String::new();
        let mut quoted = false;
        while let Some(chr) = input.peek().copied() {
            if chr.is_ascii_whitespace() && !quoted {
                break;
            }
            let mut slashes = 0;
            while input.peek() == Some(&'\\') {
                slashes += 1;
                input.next();
            }
            if input.peek() == Some(&'"') {
                input.next();
                argument.extend(std::iter::repeat_n('\\', slashes / 2));
                if slashes % 2 == 1 {
                    argument.push('"');
                } else if quoted && input.peek() == Some(&'"') {
                    input.next();
                    argument.push('"');
                } else {
                    quoted = !quoted;
                }
            } else {
                argument.extend(std::iter::repeat_n('\\', slashes));
                if let Some(chr) = input.peek().copied() {
                    if chr.is_ascii_whitespace() && !quoted {
                        break;
                    }
                    input.next();
                    argument.push(chr);
                }
            }
        }
        ensure!(!quoted, "unmatched quote in compiler command");
        args.push(argument);
    }
    Ok(args)
}

fn expand_response(
    args: &[String],
    directory: &Path,
    windows_quoting: bool,
    active: &mut BTreeSet<PathBuf>,
) -> Result<Vec<String>> {
    let mut expanded = Vec::new();
    for argument in args {
        if let Some(name) = argument.strip_prefix('@') {
            let path = fs::canonicalize(directory.join(name))
                .with_context(|| format!("cannot open response file {name}"))?;
            ensure!(
                active.insert(path.clone()),
                "recursive response file {}",
                path.display()
            );
            let bytes = fs::read(&path)?;
            let text = if bytes.starts_with(&[0xff, 0xfe]) || bytes.starts_with(&[0xfe, 0xff]) {
                ensure!(
                    bytes.len() % 2 == 0,
                    "truncated UTF-16 response file {}",
                    path.display()
                );
                let little = bytes[0] == 0xff;
                let words: Vec<_> = bytes[2..]
                    .chunks_exact(2)
                    .map(|pair| {
                        if little {
                            u16::from_le_bytes([pair[0], pair[1]])
                        } else {
                            u16::from_be_bytes([pair[0], pair[1]])
                        }
                    })
                    .collect();
                String::from_utf16(&words)?
            } else {
                std::str::from_utf8(&bytes)?
                    .trim_start_matches('\u{feff}')
                    .into()
            };
            let nested = split_command(&text, windows_quoting)?;
            expanded.extend(expand_response(
                &nested,
                directory,
                windows_quoting,
                active,
            )?);
            active.remove(&path);
        } else {
            expanded.push(argument.clone());
        }
    }
    Ok(expanded)
}

fn cpp_command(
    job: &CompileJob,
    config: &Config,
    backend: Backend,
    windows: bool,
) -> Result<Vec<String>> {
    let intel_windows = windows && backend == Backend::Ifx;
    let expanded = expand_response(
        &job.arguments[1..],
        &job.directory,
        intel_windows,
        &mut BTreeSet::new(),
    )?;
    let mut arguments = vec![job.arguments[0].clone()];
    arguments.extend_from_slice(backend.cpp(config, windows));
    let mut options = expanded.into_iter();
    while let Some(option) = options.next() {
        if matches!(
            option.as_str(),
            "-c" | "-E" | "-S" | "-MD" | "-MMD" | "-MP" | "-fpreprocessed"
        ) {
            continue;
        }
        if matches!(option.as_str(), "-o" | "-MF" | "-MT" | "-MQ") {
            ensure!(options.next().is_some(), "missing argument for {option}");
            continue;
        }
        if ["-MF", "-MT", "-MQ"]
            .iter()
            .any(|prefix| option.starts_with(prefix))
            || option.starts_with("-o") && option.ends_with(".o")
        {
            continue;
        }
        let lower = option.to_ascii_lowercase();
        if intel_windows
            && (matches!(lower.as_str(), "/c" | "/e")
                || lower.starts_with("/fo")
                || lower.starts_with("/object:"))
        {
            continue;
        }
        ensure!(
            !(matches!(
                option.as_str(),
                "-nocpp" | "-nofpp" | "-P" | "-M" | "-MM" | "-preprocess-only" | "-EP"
            ) || intel_windows
                && matches!(lower.as_str(), "/nofpp" | "/p" | "/ep" | "/preprocess-only")),
            "option {option} prevents native preprocessing with line markers"
        );
        arguments.push(option);
    }
    Ok(arguments)
}

fn marker_regex() -> &'static Regex {
    static LINE: OnceLock<Regex> = OnceLock::new();
    LINE.get_or_init(|| {
        Regex::new(r#"^\s*#\s*(?:line\s+)?(\d+)(?:\s+"((?:[^"\\]|\\.)*)")?(?:\s+\d+)*\s*$"#)
            .unwrap()
    })
}

fn marker_path(name: &str) -> String {
    let mut chars = name.chars().peekable();
    let mut decoded = String::new();
    while let Some(chr) = chars.next() {
        if chr == '\\' && chars.peek().is_some_and(|next| matches!(next, '\\' | '"')) {
            decoded.push(chars.next().unwrap());
        } else {
            decoded.push(chr);
        }
    }
    decoded
}

fn unmap_markers(
    file: &Path,
    directory: &Path,
    text: &str,
    require_markers: bool,
    original: &str,
) -> Result<Source> {
    let mut aliases = BTreeMap::new();
    for line in original.lines() {
        if let Some(marker) = marker_regex()
            .captures(line)
            .and_then(|marker| marker.get(2))
        {
            let path = PathBuf::from(marker_path(marker.as_str()));
            if path.is_absolute() {
                let parent = file.parent().unwrap_or(directory).display();
                aliases.insert(format!("{parent}/{}", path.display()), path);
            }
        }
    }
    let mut source = Source::new(file, String::new());
    let mut number = 1;
    let mut filename = Some(file.to_path_buf());
    let mut markers = false;
    for line in text.lines() {
        if let Some(marker) = marker_regex().captures(line) {
            markers = true;
            number = marker[1].parse::<usize>()?;
            if let Some(name) = marker.get(2) {
                let name = marker_path(name.as_str());
                filename = if name.starts_with('<') && name.ends_with('>') {
                    None
                } else {
                    Some(
                        aliases
                            .get(&name)
                            .cloned()
                            .unwrap_or_else(|| directory.join(name)),
                    )
                };
            }
        } else {
            source.text.push_str(line);
            source.text.push('\n');
            source.line_map.push(if number > 0 {
                filename.clone().map(|path| (path, number))
            } else {
                None
            });
            number = number
                .checked_add(1)
                .context("source line counter overflow")?;
        }
    }
    ensure!(
        !require_markers || markers,
        "{}: native preprocessor emitted no source line markers",
        file.display()
    );
    if !markers {
        source.line_map.clear();
    }
    Ok(source)
}

fn output(command: &mut Command, label: &str) -> Result<String> {
    let result = command
        .output()
        .with_context(|| format!("could not execute {label} ({command:?})"))?;
    ensure!(
        result.status.success(),
        "{label} failed ({}): {}",
        result.status,
        String::from_utf8_lossy(&result.stderr).trim()
    );
    String::from_utf8(result.stdout).with_context(|| format!("{label} did not produce UTF-8"))
}

pub fn tool_path(root: &Path, executable: &Path) -> PathBuf {
    if executable.is_relative() && executable.components().count() > 1 {
        root.join(executable)
    } else {
        executable.into()
    }
}

fn cpp_cache_key(job: &CompileJob, argv: &[String], input: &str) -> Option<String> {
    let volatile = ["__DATE__", "__TIME__", "__TIMESTAMP__"];
    let unknown_options = [
        "-include",
        "-imacros",
        "-Xpreprocessor",
        "-Wp,",
        "-Qoption",
        "/Qoption",
        "/FI",
        "-specs",
        "--specs",
        "-B",
    ];
    if input
        .lines()
        .any(|line| line.contains('#') && !marker_regex().is_match(line))
        || argv.iter().skip(1).any(|argument| {
            unknown_options
                .iter()
                .any(|prefix| argument.starts_with(prefix))
        })
        || volatile
            .iter()
            .any(|name| input.contains(name) || argv.iter().any(|item| item.contains(name)))
    {
        return None;
    }
    let compiler = PathBuf::from(&argv[0]);
    if !compiler.is_absolute() {
        return None;
    }
    let metadata = fs::metadata(compiler).ok()?;
    let version = metadata
        .modified()
        .ok()?
        .duration_since(std::time::UNIX_EPOCH)
        .ok()?;
    let vars: BTreeMap<_, _> = [
        "PATH",
        "CPATH",
        "FPP",
        "IFXCFG",
        "INCLUDE",
        "GCC_EXEC_PREFIX",
        "COMPILER_PATH",
        "SOURCE_DATE_EPOCH",
    ]
    .into_iter()
    .map(|name| (name, std::env::var(name).ok()))
    .collect();
    if vars["FPP"].is_some() || vars["IFXCFG"].is_some() {
        return None;
    }
    let key = serde_json::to_vec(&(
        job.directory.as_path(),
        job.file.as_path(),
        argv,
        input,
        vars,
        version.as_nanos(),
        metadata.len(),
    ))
    .ok()?;
    Some(format!("{:x}", Sha256::digest(key)))
}

fn native_source(
    job: CompileJob,
    backend: Backend,
    config: &Config,
    cache: Option<&Path>,
) -> Result<Source> {
    let raw = fs::read_to_string(&job.file)
        .with_context(|| format!("cannot read {}", job.file.display()))?;
    let argv = cpp_command(&job, config, backend, cfg!(windows))?;
    let cached = cache.and_then(|base| {
        cpp_cache_key(&job, &argv, &raw).map(|key| base.join("cpp").join(format!("{key}.json")))
    });
    if let Some(source) = cached
        .as_ref()
        .and_then(|path| fs::read(path).ok())
        .and_then(|data| serde_json::from_slice(&data).ok())
    {
        return Ok(source);
    }
    let expanded = output(
        Command::new(&argv[0])
            .current_dir(&job.directory)
            .args(&argv[1..]),
        &format!("native preprocessing of {}", job.file.display()),
    )?;
    let source = unmap_markers(&job.file, &job.directory, &expanded, true, &raw)?;
    if let Some(path) = cached {
        crate::output::atomic_write(&path, &serde_json::to_vec(&source)?)?;
    }
    Ok(source)
}

fn direct_source(
    path: &Path,
    root: &Path,
    options: &InputOptions,
    config: &Config,
) -> Result<Source> {
    let raw = if extension(path) == "fypp" {
        let relative = slash_path(path.strip_prefix(root).unwrap_or(path));
        let args = config
            .fypp
            .args
            .iter()
            .chain(
                config
                    .fypp
                    .overrides
                    .iter()
                    .filter(|entry| config::any_match(&entry.paths, &relative))
                    .flat_map(|entry| &entry.args),
            )
            .map(|arg| {
                arg.replace("{line_length}", &config.language.line_length.to_string())
                    .replace("{root}", &slash_path(root))
            });
        output(
            Command::new(tool_path(root, &config.fypp.executable))
                .current_dir(root)
                .args(args)
                .arg(path),
            &format!("Fypp expansion of {}", path.display()),
        )?
    } else {
        fs::read_to_string(path).with_context(|| format!("cannot read {}", path.display()))?
    };
    if let Some(executable) = &options.cpp_executable {
        let backend = options
            .backend
            .context("--cpp-executable also requires --compiler-id")?;
        let staged = tempfile::Builder::new().suffix(".F90").tempfile()?;
        let file = if extension(path) == "fypp" {
            fs::write(staged.path(), &raw)?;
            staged.path()
        } else {
            path
        };
        let mut command = Command::new(tool_path(root, executable));
        let include = if backend == Backend::Ifx && cfg!(windows) {
            "/I"
        } else {
            "-I"
        };
        command
            .current_dir(root)
            .args(backend.cpp(config, cfg!(windows)))
            .args(&options.cpp_args)
            .arg(format!(
                "{include}{}",
                path.parent().unwrap_or(root).display()
            ))
            .arg(file);
        let expanded = output(
            &mut command,
            &format!("native preprocessing of {}", path.display()),
        )?;
        unmap_markers(file, root, &expanded, true, &raw).map(|mut source| {
            source.path = path.into();
            source
        })
    } else {
        unmap_markers(path, root, &raw, extension(path) == "fypp", &raw)
    }
}

pub fn collect(
    root: &Path,
    options: &InputOptions,
    config: &Config,
    cache: Option<&Path>,
) -> Result<Vec<Source>> {
    let results: Result<Vec<_>> = if let Some(path) = &options.compile_db {
        let backend = options
            .backend
            .context("--compile-db requires --compiler-id")?;
        let db = root.join(path);
        read_database(&db, root, config, cfg!(windows))?
            .into_par_iter()
            .map(|job| native_source(job, backend, config, cache))
            .collect()
    } else {
        sources_from_paths(root, &options.paths, config)?
            .into_par_iter()
            .map(|path| direct_source(&path, root, options, config))
            .collect()
    };
    let mut by_path = BTreeMap::new();
    for source in results? {
        if options.compile_db.is_some() {
            let original = source
                .text
                .lines()
                .zip(&source.line_map)
                .find_map(|(line, location)| {
                    let line = line.trim_start();
                    let (path, _) = location.as_ref()?;
                    (!line.is_empty()
                        && !line.starts_with('!')
                        && supported_extension(path, config))
                    .then_some(path)
                })
                .unwrap_or(&source.path);
            if !source_included(root, original, config) {
                continue;
            }
        }
        if let Some(previous) = by_path.get(&source.path) {
            ensure!(
                previous == &source,
                "{} has different expanded variants in the same project",
                source.path.display()
            );
        } else {
            by_path.insert(source.path.clone(), source);
        }
    }
    ensure!(
        config.inputs.allow_empty || !by_path.is_empty(),
        "no Fortran sources matched the input policy"
    );
    Ok(by_path.into_values().collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn windows_quoted_arguments_keep_first_option_and_paths() {
        let result = split_command(
            r#" /I"C:\Program Files\Headers"  /DTEST -c "a b.F90" "" "#,
            true,
        )
        .unwrap();
        assert_eq!(
            result,
            [r"/IC:\Program Files\Headers", "/DTEST", "-c", "a b.F90", ""]
        );
    }

    #[test]
    fn native_and_fypp_markers_map_locations() {
        let input = "# 1 \"<built-in>\"\n# 71 \"src/file.fypp\"\nuse other\n#line 18 \"src/header.h\"\ninvalid\n";
        let source = unmap_markers(
            Path::new("/project/output.f90"),
            Path::new("/project"),
            input,
            true,
            "",
        )
        .unwrap();
        assert_eq!(
            source.location(1, 1).path,
            Path::new("/project/src/file.fypp")
        );
        assert_eq!(source.location(1, 1).line, 71);
        assert_eq!(source.location(2, 3).line, 18);
        assert!(unmap_markers(Path::new("out.f90"), Path::new("."), "end\n", true, "").is_err());
    }
}
