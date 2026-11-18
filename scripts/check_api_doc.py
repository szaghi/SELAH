#!/usr/bin/env python3
"""Check that the generated API reference lists every name the facade module src/lib/selah.F90 exports.

Why: `use selah` is the whole API, and its generated page (docs/api/src/lib/selah.md, written by formal: run
`npm run docs:api` in docs/ first) must list each exported name under Re-exports; a formal that dropped one would leave
it invisible to users. Exit status 1, with the missing names, if any.

Also checks that README.md, docs/index.md and docs/ford.md (the API overview) list every exported SCHEME_* constant
(the landing-page tables).

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
page = root / 'docs/api/src/lib/selah.md'
if not page.is_file():
    print(f'{page.relative_to(root)} not found: run `npm run docs:api` in docs/ first')
    sys.exit(1)
doc = page.read_text()
reexports = doc[doc.index('## Re-exports'):] if '## Re-exports' in doc else ''
missing = [n for n in names if not re.search(r'`' + re.escape(n) + r'`', reexports, re.IGNORECASE)]
print(f'{len(names)} names exported by selah.F90; missing from the Re-exports of {page.relative_to(root)}: '
      f'{missing or "none"}')
status = 1 if missing else 0
# The landing pages list the time integrators: every exported SCHEME_* constant must appear in both.
schemes = [n for n in names if n.startswith('SCHEME_')]
for page in ('README.md', 'docs/index.md', 'docs/ford.md'):
    text = (root / page).read_text()
    absent = [n for n in schemes if f'`{n}`' not in text]
    print(f'{len(schemes)} schemes exported; missing from {page}: {absent or "none"}')
    status |= 1 if absent else 0
sys.exit(status)
