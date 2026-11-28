use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
};

use anyhow::{Context, Result, ensure};
use serde::Deserialize;
use sha2::{Digest, Sha256};

use crate::{
    config::Config, diagnostic::Diagnostic, input::tool_path, output::write_json, source::Source,
};

pub struct Fortitude {
    executable: PathBuf,
    user_config: Option<PathBuf>,
    version: String,
}

impl Fortitude {
    pub fn discover(config: &Config, root: &Path) -> Result<Option<Self>> {
        if !config.fortitude.enabled {
            return Ok(None);
        }
        let executable = tool_path(root, &config.fortitude.executable);
        let output = run(Command::new(&executable).arg("--version"))?;
        let version = String::from_utf8(output.stdout)?.trim().to_owned();
        ensure!(
            output.status.success() && version == format!("fortitude {}", config.fortitude.version),
            "expected Fortitude {}, got '{}' ({})",
            config.fortitude.version,
            version,
            output.status
        );
        Ok(Some(Self {
            executable,
            version,
            user_config: config.fortitude.config.as_ref().map(|path| root.join(path)),
        }))
    }

    pub fn fingerprint(&self) -> Result<Vec<u8>> {
        let config = self
            .user_config
            .as_ref()
            .map(|path| fs::read(path).with_context(|| format!("cannot read {}", path.display())))
            .transpose()?;
        let metadata = fs::metadata(&self.executable).ok().and_then(|file| {
            let time = file
                .modified()
                .ok()?
                .duration_since(std::time::UNIX_EPOCH)
                .ok()?;
            Some((time.as_nanos(), file.len()))
        });
        Ok(serde_json::to_vec(&(
            &self.executable,
            &self.version,
            config,
            metadata,
            std::env::var("PATH").ok(),
        ))?)
    }

    pub fn check(
        &self,
        sources: &[Source],
        config: &Config,
        cache: Option<&Path>,
    ) -> Result<Vec<Diagnostic>> {
        let key = cache
            .map(|_| -> Result<_> {
                Ok(serde_json::to_vec(&(
                    env!("CARGO_PKG_VERSION"),
                    self.fingerprint()?,
                    &config.fortitude,
                    &config.language.standard,
                    config.language.line_length,
                ))?)
            })
            .transpose()?;
        let mut pending = Vec::new();
        let mut diagnostics = Vec::new();
        for source in sources {
            let path = cache
                .zip(key.as_ref())
                .map(|(directory, key)| -> Result<_> {
                    let mut hash = Sha256::new();
                    hash.update(key);
                    hash.update(serde_json::to_vec(source)?);
                    Ok(directory
                        .join("fortitude")
                        .join(format!("{:x}.json", hash.finalize())))
                })
                .transpose()?;
            if let Some(cached) = path
                .as_ref()
                .and_then(|file| fs::read(file).ok())
                .and_then(|bytes| serde_json::from_slice::<Vec<Diagnostic>>(&bytes).ok())
            {
                diagnostics.extend(cached);
            } else {
                pending.push((source, path));
            }
        }
        if pending.is_empty() {
            return Ok(diagnostics);
        }
        let new = self.run_check(
            &pending
                .iter()
                .map(|(source, _)| *source)
                .collect::<Vec<_>>(),
            config,
        )?;
        for ((_, path), found) in pending.into_iter().zip(new) {
            if let Some(path) = path {
                write_json(&path, &found)?;
            }
            diagnostics.extend(found);
        }
        Ok(diagnostics)
    }

    fn run_check(&self, sources: &[&Source], config: &Config) -> Result<Vec<Vec<Diagnostic>>> {
        let stage = tempfile::Builder::new()
            .prefix("fortran-check-fortitude-")
            .tempdir()?;
        let policy_path = match &self.user_config {
            Some(path) => path.clone(),
            None => {
                let path = stage.path().join("fortitude.toml");
                fs::write(&path, "[check]\nfix = false\nunsafe-fixes = false\n")?;
                path
            }
        };
        let files = sources
            .iter()
            .enumerate()
            .map(|(index, source)| -> Result<_> {
                let name = format!("{index:06}.f90");
                fs::write(stage.path().join(&name), &source.text)?;
                Ok((name, (index, source)))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;
        let file_names: Vec<_> = files.keys().collect();
        let mut diagnostics = vec![Vec::new(); sources.len()];
        for names in file_names.chunks(config.fortitude.batch_size) {
            let mut command = Command::new(&self.executable);
            command
                .current_dir(stage.path())
                .arg("check")
                .arg("--config")
                .arg(&policy_path)
                .args(["--config", "check.force-exclude = false", "--config"])
                .arg(format!(
                    "check.line-length = {}",
                    config.language.line_length
                ))
                .arg("--target-std")
                .arg(&config.language.standard)
                .arg("--select")
                .arg(config.fortitude.select.join(","))
                .args([
                    "--no-respect-gitignore",
                    "--no-fix",
                    "--no-fix-only",
                    "--output-format",
                    "json",
                    "--quiet",
                ])
                .arg(if config.fortitude.preview {
                    "--preview"
                } else {
                    "--no-preview"
                });
            let mut ignored = config.fortitude.ignore.clone();
            if !ignored.iter().any(|rule| rule.eq_ignore_ascii_case("C121")) {
                ignored.push("C121".into());
            }
            command.arg("--ignore").arg(ignored.join(","));
            if !config.fortitude.allow_comments {
                command.arg("--ignore-allow-comments");
            }
            command.arg("--").args(names);
            let output = run(&mut command)?;
            ensure!(
                matches!(output.status.code(), Some(0 | 1)),
                "Fortitude exited with {}: {}",
                output.status,
                String::from_utf8_lossy(&output.stderr)
            );
            let found: Vec<FortitudeDiagnostic> =
                serde_json::from_slice(&output.stdout).context("invalid Fortitude diagnostics")?;
            ensure!(
                output.status.success() || !found.is_empty(),
                "Fortitude returned a failed check without diagnostics"
            );
            for item in found {
                let name = item
                    .filename
                    .file_name()
                    .and_then(|name| name.to_str())
                    .context("missing Fortitude source filename")?;
                let (index, source) = files
                    .get(name)
                    .with_context(|| format!("Fortitude returned an unstaged source '{name}'"))?;
                ensure!(
                    item.location.row > 0 && item.location.column > 0 && !item.code.is_empty(),
                    "invalid Fortitude diagnostic location or ID"
                );
                let mut diagnostic = Diagnostic::new(
                    item.code,
                    item.message,
                    source.location(item.location.row, item.location.column),
                );
                diagnostic.severity = config.fortitude.severity;
                diagnostics[*index].push(diagnostic);
            }
        }
        Ok(diagnostics)
    }
}

fn run(command: &mut Command) -> Result<Output> {
    for (name, _) in std::env::vars_os() {
        if name
            .to_string_lossy()
            .to_ascii_uppercase()
            .starts_with("FORTITUDE_")
        {
            command.env_remove(name);
        }
    }
    command
        .output()
        .with_context(|| format!("could not run Fortitude ({command:?})"))
}

#[derive(Deserialize)]
struct FortitudeDiagnostic {
    code: String,
    message: String,
    filename: PathBuf,
    location: Position,
}

#[derive(Deserialize)]
struct Position {
    row: usize,
    column: usize,
}
