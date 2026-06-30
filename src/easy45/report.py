"""Run reporting: collects warnings and writes the human-readable summary.

Warnings are accumulated across stages (organelle depletion skipped, hybrid
ribotype signal, IGS incomplete, ...) and surfaced together at the end so the
user can judge how much to trust each output — essential given the tool's
ribotype-diversity positioning.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Report:
    warnings: list[str] = field(default_factory=list)
    stats: dict = field(default_factory=dict)

    def warn(self, message: str) -> None:
        self.warnings.append(message)

    def write(self, path) -> None:
        """Write the summary report (TSV/markdown). Phase 2."""
        raise NotImplementedError
