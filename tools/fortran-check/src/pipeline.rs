use std::{
    fs,
    path::{Path, PathBuf},
};

use anyhow::Result;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use sha2::{Digest, Sha256};

use crate::{
    config::Config,
    diagnostic::{Diagnostic, Severity, sarif},
    fortitude::Fortitude,
    input::{self, InputOptions},
    output::{atomic_write, write_json},
    semantic::Project,
    source::Unit,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Report {
    pub sources: usize,
    pub complete: bool,
    pub diagnostics: Vec<Diagnostic>,
    pub metrics: Value,
    pub dependencies_dot: String,
    #[serde(default)]
    pub cached: bool,
}

impl Report {
    pub fn errors(&self) -> usize {
        self.diagnostics
            .iter()
            .filter(|item| item.severity == Severity::Error)
            .count()
    }

    pub fn warnings(&self) -> usize {
        self.diagnostics
            .iter()
            .filter(|item| item.severity == Severity::Warning)
            .count()
    }

    pub fn write(&self, directory: &Path, root: &Path, config: &Config) -> Result<()> {
        if config.reports.include_sarif {
            write_json(
                &directory.join("diagnostics.sarif"),
                &sarif(&self.diagnostics, root),
            )?;
        }
        if config.reports.include_metrics {
            write_json(&directory.join("metrics.json"), &self.metrics)?;
        }
        if config.reports.include_graph {
            atomic_write(
                &directory.join("dependencies.dot"),
                self.dependencies_dot.as_bytes(),
            )?;
        }
        Ok(())
    }
}

fn cache_file(
    cache: &Path,
    root: &Path,
    config: &Config,
    sources: &[crate::source::Source],
    fortitude: &[u8],
) -> Result<Option<PathBuf>> {
    let stamp = std::env::current_exe().ok().and_then(|executable| {
        let file = fs::metadata(&executable).ok()?;
        let modified = file
            .modified()
            .ok()?
            .duration_since(std::time::UNIX_EPOCH)
            .ok()?
            .as_nanos();
        Some((executable, modified, file.len()))
    });
    let Some(stamp) = stamp else { return Ok(None) };
    let mut digest = Sha256::new();
    digest.update(env!("CARGO_PKG_VERSION").as_bytes());
    digest.update(serde_json::to_vec(&(root, config, stamp, fortitude))?);
    for source in sources {
        digest.update(serde_json::to_vec(&source)?);
    }
    Ok(Some(
        cache
            .join("report")
            .join(format!("{:x}.json", digest.finalize())),
    ))
}

pub fn run(
    root: &Path,
    options: &InputOptions,
    config: &Config,
    cache: Option<&Path>,
) -> Result<Report> {
    let executor = rayon::ThreadPoolBuilder::new()
        .num_threads(config.runtime.jobs)
        .build()?;
    executor.install(|| {
        let fortitude = Fortitude::discover(config, root)?;
        let cache = cache.filter(|_| config.runtime.cache);
        let sources = input::collect(root, options, config, cache)?;
        let file = cache
            .map(|cache| {
                let extra = fortitude
                    .as_ref()
                    .map(Fortitude::fingerprint)
                    .transpose()?
                    .unwrap_or_default();
                cache_file(cache, root, config, &sources, &extra)
            })
            .transpose()?
            .flatten();
        if let Some(mut report) = file
            .as_ref()
            .and_then(|file| fs::read(file).ok())
            .and_then(|bytes| serde_json::from_slice::<Report>(&bytes).ok())
        {
            report.cached = true;
            return Ok(report);
        }
        let (units, fortitude_results) = rayon::join(
            || {
                sources
                    .par_iter()
                    .map(|source| Unit::parse(source.clone()))
                    .collect::<Result<Vec<_>>>()
            },
            || {
                fortitude
                    .as_ref()
                    .map(|tool| tool.check(&sources, config, cache))
                    .transpose()
            },
        );
        let mut diagnostics = fortitude_results?.unwrap_or_default();
        let mut valid = Vec::new();
        for unit in units? {
            let problems = unit.diagnostics();
            if problems.is_empty() {
                valid.push(unit);
            } else {
                diagnostics.extend(problems);
            }
        }
        let complete = valid.len() == sources.len();
        let project = Project::new(&valid, config);
        diagnostics.extend(project.diagnostics.iter().cloned());
        diagnostics.extend(crate::language::check(&project, config));
        let design = crate::design::analyze(&project, config, root);
        diagnostics.extend(design.diagnostics);
        diagnostics.extend(crate::memory::check(&project, config));
        diagnostics.extend(crate::ownership::check(&project, config));
        diagnostics.extend(crate::numerical::check(&project, config));
        diagnostics.extend(crate::parallel::check(&project, config));
        let mut diagnostics: Vec<_> = diagnostics
            .into_iter()
            .filter_map(|item| config.filter(item, root))
            .collect();
        diagnostics.sort();
        diagnostics.dedup();
        let report = Report {
            sources: sources.len(),
            complete,
            diagnostics,
            metrics: design.metrics,
            dependencies_dot: design.graph,
            cached: false,
        };
        if let Some(path) = file {
            write_json(&path, &report)?;
        }
        Ok(report)
    })
}
