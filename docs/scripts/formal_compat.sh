#!/usr/bin/env bash
# Bring the call graphs of a FORMAL release older than the class-based highlight to the form FORMAL now writes
# (szaghi/formal, generator.py format_call_diagram): left to right, the procedure marked with the class `focus`, no
# colours. Before it, FORMAL wrote `flowchart TD` (a procedure with many callers came out too wide to read) and an
# inline `style X fill:#3e63dd,...` that no stylesheet can override and that fails contrast in both themes.
# Idempotent: on the output of a FORMAL that writes the new form it changes nothing. Drop it once the FORMAL used by
# the docs build (the toolchain image's) writes the new form.
#
# Usage: bash scripts/formal_compat.sh <api dir>      (npm run docs:api calls it after formal generate)
set -euo pipefail
api="${1:?usage: formal_compat.sh <api dir>}"
find "$api" -name '*.md' -print0 | xargs -0 -r sed -i \
  -e 's/^flowchart TD$/flowchart LR/' \
  -e 's/^  style \([A-Za-z0-9_]*\) fill:#3e63dd,stroke:#99b,stroke-width:2px$/  classDef focus stroke-width:3px,font-weight:bold\n  class \1 focus/'
