"""Thin wrappers around the external bioinformatics tools easy45 depends on.

Design contract (driven by the bioconda packaging goal):
    The Python package contains *no* bundled binaries. Every heavy tool is a
    conda run-dependency declared in recipe/meta.yaml. This module is the only
    place that shells out, and it fails *fast and clearly* at startup if a
    required tool is missing from PATH — never mid-pipeline.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass


@dataclass(frozen=True)
class Tool:
    name: str            # executable name on PATH
    version_args: tuple  # args that print version info (for the startup check)
    optional: bool = False


# Tools used by the pipeline. `optional` ones are only needed for certain modes.
REQUIRED_TOOLS = (
    Tool("minimap2", ("--version",)),
    Tool("seqkit", ("version",)),
    Tool("vsearch", ("--version",)),
    Tool("ITSx", ("-h",)),
    Tool("barrnap", ("--version",)),
    Tool("abpoa", ("-v",)),
    Tool("cmsearch", ("-h",)),   # infernal; CM-based mature-boundary trimming (S5)
)


class DependencyError(RuntimeError):
    """Raised when a required external tool is missing from PATH."""


def check_dependencies(include_optional: bool = False) -> dict[str, str | None]:
    """Verify required tools are on PATH. Returns {tool: resolved_path_or_None}.

    Raises DependencyError listing every missing *required* tool at once.
    """
    found: dict[str, str | None] = {}
    missing: list[str] = []
    for tool in REQUIRED_TOOLS:
        path = shutil.which(tool.name)
        found[tool.name] = path
        if path is None and not tool.optional and not include_optional:
            missing.append(tool.name)
        elif path is None and tool.optional and include_optional:
            missing.append(tool.name)
    if missing:
        raise DependencyError(
            "Missing required tool(s): "
            + ", ".join(missing)
            + "\nInstall everything with: conda env create -f environment.yml"
        )
    return found


def run(cmd: list[str], **kwargs) -> subprocess.CompletedProcess:
    """Run an external command, raising CalledProcessError on non-zero exit.

    Output is captured by default so callers can inspect / log it.
    """
    kwargs.setdefault("check", True)
    kwargs.setdefault("text", True)
    kwargs.setdefault("capture_output", True)
    return subprocess.run([str(c) for c in cmd], **kwargs)
