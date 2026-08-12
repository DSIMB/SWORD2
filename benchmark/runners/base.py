from __future__ import annotations

import math
import os
import re
import subprocess
import time
from dataclasses import dataclass, replace
from pathlib import Path


@dataclass(frozen=True)
class ParsedPrediction:
    chopping: str
    n_domains: int | None = None
    raw: dict[str, str] | None = None


@dataclass(frozen=True)
class PartitionPrediction:
    tool: str
    name: str
    variant: str
    chopping: str
    n_domains: int
    raw_path: Path | None = None
    runtime_s: float | None = None


@dataclass(frozen=True)
class ToolRunResult:
    command: list[str]
    cwd: Path
    returncode: int
    runtime_s: float
    peak_rss_mb: float | None
    stdout: str
    stderr: str
    peak_rss_kb: int | None = None


_RSS_RE = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")
_CUDA_CHECK_CODE = """
import sys
import torch

available = torch.cuda.is_available()
print(
    f"torch={torch.__version__} "
    f"cuda_available={available} "
    f"cuda_version={torch.version.cuda} "
    f"device_count={torch.cuda.device_count()}"
)
raise SystemExit(0 if available else 1)
""".strip()


def parse_peak_rss_mb(time_stderr: str) -> float | None:
    match = _RSS_RE.search(time_stderr)
    if not match:
        return None
    return int(match.group(1)) / 1024.0


def parse_locked_peak_rss_kb(time_stderr: str) -> int:
    labels = [
        line
        for line in time_stderr.splitlines()
        if "Maximum resident set size" in line
    ]
    if len(labels) != 1:
        raise ValueError("locked GNU time evidence must contain exactly one RSS line")
    match = re.fullmatch(r"\s*Maximum resident set size \(kbytes\):\s*(\d+)\s*", labels[0])
    if match is None:
        raise ValueError("locked GNU time RSS line has the wrong label or unit")
    value = int(match.group(1))
    if value <= 0:
        raise ValueError("locked GNU time RSS must be positive")
    return value


def ensure_cuda_available(python: Path, cwd: Path) -> None:
    command = [str(python.absolute()), "-c", _CUDA_CHECK_CODE]
    completed = subprocess.run(
        command,
        cwd=cwd.resolve(),
        text=True,
        capture_output=True,
        check=False,
        timeout=30,
    )
    if completed.returncode == 0:
        return

    details = "\n".join(
        part for part in [completed.stdout.strip(), completed.stderr.strip()] if part
    )
    if not details:
        details = f"CUDA check exited with status {completed.returncode}"
    raise RuntimeError(
        "CUDA is required for this benchmark run, but the tool Python could not initialize it.\n"
        f"Python: {python}\n"
        f"CWD: {cwd}\n"
        f"{details}\n"
        "Fix NVIDIA driver/CUDA visibility, or pass --allow-dl-cpu to permit CPU fallback."
    )


def run_timed(
    command: list[str],
    cwd: Path,
    env: dict[str, str] | None = None,
    timeout_s: int | None = None,
) -> ToolRunResult:
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)

    timed_command = ["/usr/bin/time", "-v", *command]
    start = time.perf_counter()
    completed = subprocess.run(
        timed_command,
        cwd=cwd.resolve(),
        env=merged_env,
        timeout=timeout_s,
        text=True,
        capture_output=True,
        check=False,
    )
    runtime_s = time.perf_counter() - start
    peak_rss_mb = parse_peak_rss_mb(completed.stderr)
    return ToolRunResult(
        command=command,
        cwd=cwd.resolve(),
        returncode=completed.returncode,
        runtime_s=runtime_s,
        peak_rss_mb=peak_rss_mb,
        stdout=completed.stdout,
        stderr=completed.stderr,
        peak_rss_kb=(
            int(round(peak_rss_mb * 1024))
            if peak_rss_mb is not None
            else None
        ),
    )


def run_timed_locked(
    command: list[str],
    cwd: Path,
    env: dict[str, str],
) -> ToolRunResult:
    result = run_timed(command, cwd=cwd, env=env, timeout_s=None)
    if not math.isfinite(result.runtime_s) or result.runtime_s <= 0:
        raise ValueError("locked wall duration must be finite and positive")
    peak_rss_kb = parse_locked_peak_rss_kb(result.stderr)
    return replace(
        result,
        peak_rss_kb=peak_rss_kb,
        peak_rss_mb=peak_rss_kb / 1024.0,
    )
