#!/usr/bin/env bash
# Emit a wall-clock line to stderr every INTERVAL seconds while COMMAND runs.
# Usage: with_heartbeat.sh <interval_s> <command...>
set -euo pipefail
interval="${1:?interval seconds required}"
shift
(
  while true; do
    echo "[heartbeat] $(date -Iseconds)" >&2
    sleep "$interval"
  done
) &
hb_pid=$!
trap 'kill "$hb_pid" 2>/dev/null || true' EXIT
"$@"
status=$?
kill "$hb_pid" 2>/dev/null || true
wait "$hb_pid" 2>/dev/null || true
exit "$status"
