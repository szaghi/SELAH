#!/usr/bin/env python3
"""Check that docs/guide/api.md lists every name the facade module src/lib/selah.F90 exports (`public ::` lines).

Why: the API page is hand-written (the generated reference cannot show the facade's re-exports), so a name added to the
facade and forgotten on the page would be invisible to users. Exit status 1, with the missing names, if any.

Also checks that README.md and docs/index.md list every exported SCHEME_* constant (the landing-page table).

Usage: python3 scripts/check_api_doc.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
names = []
for line in (root / 'src/lib/selah.F90').read_text().splitlines():
    m = re.match(r'\s*public\s*::\s*(.*)', line)
    if m:
        names += [n.strip() for n in m.group(1).split(',') if n.strip()]
doc = (root / 'docs/guide/api.md').read_text()
missing = [n for n in names if not re.search(r'`' + re.escape(n) + r'[`(\[ ]', doc)]
print(f'{len(names)} names exported by selah.F90; missing from docs/guide/api.md: {missing or "none"}')
status = 1 if missing else 0
# The landing pages list the time integrators: every exported SCHEME_* constant must appear in both.
schemes = [n for n in names if n.startswith('SCHEME_')]
for page in ('README.md', 'docs/index.md'):
    text = (root / page).read_text()
    absent = [n for n in schemes if f'`{n}`' not in text]
    print(f'{len(schemes)} schemes exported; missing from {page}: {absent or "none"}')
    status |= 1 if absent else 0
sys.exit(status)
