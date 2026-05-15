#!/usr/bin/env python3
"""
Inspect Cursor terminal snapshot files to see whether shells look alive or stuck.

Typical use (from repo root):
  python scripts/terminal_probe.py --auto
  python scripts/terminal_probe.py --dir "$HOME/.cursor/projects/<slug>/terminals"

Heuristics:
  - Parse the leading --- metadata block for pid, cwd, command, running_for_ms.
  - If pid is present, ask the OS whether that pid still exists (Unix: signal 0).
  - Scan the file tail for an exit_code footer (agent/shell completed).
  - Show mtime + last lines of transcript for "is output still moving?" spot checks.
"""

from __future__ import annotations

import argparse
import errno
import os
import subprocess
import sys
import time
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Probe Cursor terminal snapshot files for liveness.")
    p.add_argument(
        "--dir",
        type=str,
        default="",
        help="Path to a Cursor terminals/ folder containing *.txt snapshots.",
    )
    p.add_argument(
        "--auto",
        action="store_true",
        help="Locate this Git repo root and scan ~/.cursor/projects/*/terminals/*.txt for matching cwd.",
    )
    p.add_argument("--tail", type=int, default=12, help="Lines of each transcript tail to print.")
    return p.parse_args()


def git_toplevel(start: Path) -> Path | None:
    try:
        out = subprocess.check_output(
            ["git", "-C", str(start), "rev-parse", "--show-toplevel"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        return Path(out).resolve()
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return None


def default_cursor_terminals_dir() -> Path | None:
    env = os.environ.get("CURSOR_TERMINALS_DIR", "").strip()
    if env:
        return Path(env).expanduser().resolve()
    return None


def iter_terminal_files(term_dir: Path) -> list[Path]:
    if not term_dir.is_dir():
        return []
    return sorted(term_dir.glob("*.txt"), key=lambda p: p.stat().st_mtime, reverse=True)


def parse_meta_and_body(text: str) -> tuple[dict[str, str], str]:
    lines = text.splitlines()
    meta: dict[str, str] = {}
    if not lines or lines[0].strip() != "---":
        return meta, text
    i = 1
    while i < len(lines) and lines[i].strip() != "---":
        line = lines[i]
        if ":" in line:
            key, _, val = line.partition(":")
            meta[key.strip()] = val.strip()
        i += 1
    body = "\n".join(lines[i + 1 :])
    return meta, body


def file_has_exit_footer(text: str) -> bool:
    tail = "\n".join(text.splitlines()[-40:])
    return "exit_code:" in tail


def pid_exists(pid: int) -> bool | str:
    if pid <= 0:
        return "invalid"
    try:
        os.kill(pid, 0)
    except OSError as exc:
        if exc.errno == errno.ESRCH:
            return False
        if exc.errno == errno.EPERM:
            return "unknown(no-perm)"
        return f"unknown({exc.errno})"
    return True


def discover_auto_term_dirs(repo_root: Path) -> list[Path]:
    projects = Path.home() / ".cursor" / "projects"
    if not projects.is_dir():
        return []
    hits: list[Path] = []
    for term_dir in projects.glob("*/terminals"):
        if not term_dir.is_dir():
            continue
        for snap in term_dir.glob("*.txt"):
            try:
                text = snap.read_text(errors="replace")
            except OSError:
                continue
            meta, _ = parse_meta_and_body(text)
            cwd = meta.get("cwd", "").strip().strip('"')
            if not cwd:
                continue
            try:
                cwdp = Path(cwd).resolve()
            except OSError:
                continue
            try:
                cwdp.relative_to(repo_root)
                hits.append(term_dir)
                break
            except ValueError:
                if repo_root == cwdp or str(cwdp).startswith(str(repo_root) + os.sep):
                    hits.append(term_dir)
                    break
    seen: set[str] = set()
    out: list[Path] = []
    for d in hits:
        key = str(d.resolve())
        if key not in seen:
            seen.add(key)
            out.append(d)
    return out


def probe_file(path: Path, tail_n: int) -> None:
    try:
        st = path.stat()
    except OSError as exc:
        print(f"{path}: cannot stat ({exc})")
        return
    text = path.read_text(errors="replace")
    meta, body = parse_meta_and_body(text)
    pid_s = meta.get("pid", "").strip()
    try:
        pid = int(pid_s) if pid_s else 0
    except ValueError:
        pid = 0

    alive = pid_exists(pid) if pid else "n/a"
    exited = file_has_exit_footer(text)
    print(f"\n== {path} ==")
    print(f"mtime_age_s={max(0, int(time.time() - st.st_mtime))} size_bytes={st.st_size}")
    if meta:
        for k in ("pid", "cwd", "last_command", "command", "started_at", "running_for_ms", "last_exit_code"):
            if k in meta:
                print(f"{k}: {meta[k]}")
    print(f"pid_alive={alive} footer_exit_hint={exited}")
    tail_lines = body.splitlines()[-tail_n:]
    if tail_lines:
        print("--- tail ---")
        for ln in tail_lines:
            print(ln)


def main() -> int:
    args = parse_args()
    term_dirs: list[Path] = []

    if args.dir:
        term_dirs.append(Path(args.dir).expanduser().resolve())
    elif args.auto:
        here = Path.cwd().resolve()
        root = git_toplevel(here)
        if root is None:
            print("error: --auto requires a Git checkout (git rev-parse failed).", file=sys.stderr)
            return 2
        term_dirs = discover_auto_term_dirs(root)
        if not term_dirs:
            print(
                "warning: no Cursor terminals dirs matched this repo under ~/.cursor/projects.\n"
                "Set CURSOR_TERMINALS_DIR or pass --dir explicitly; see docs/LONG_RUNNING_JOBS.md",
                file=sys.stderr,
            )
    else:
        d = default_cursor_terminals_dir()
        if d:
            term_dirs.append(d)

    if not term_dirs:
        print(
            "usage: terminal_probe.py (--dir PATH | --auto | set CURSOR_TERMINALS_DIR)\n"
            "See docs/LONG_RUNNING_JOBS.md",
            file=sys.stderr,
        )
        return 2

    any_files = False
    for d in term_dirs:
        files = iter_terminal_files(d)
        print(f"# terminals_dir={d} files={len(files)}")
        if not files:
            continue
        any_files = True
        for f in files:
            probe_file(f, args.tail)

    if not any_files:
        print("no *.txt snapshots found in given dir(s).", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
