use std::{fs, io::Write, path::Path};

use anyhow::{Context, Result};
use serde::Serialize;

pub(crate) fn atomic_write(path: &Path, contents: &[u8]) -> Result<()> {
    let directory = path
        .parent()
        .filter(|base| !base.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(directory)
        .with_context(|| format!("cannot create {}", directory.display()))?;
    let mut temporary = tempfile::NamedTempFile::new_in(directory)?;
    temporary.as_file_mut().write_all(contents)?;
    temporary
        .persist(path)
        .with_context(|| format!("cannot publish {}", path.display()))?;
    Ok(())
}

pub fn write_json(path: &Path, value: &impl Serialize) -> Result<()> {
    let mut text = serde_json::to_vec_pretty(value)?;
    text.push(b'\n');
    atomic_write(path, &text)
}
