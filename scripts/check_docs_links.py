#!/usr/bin/env python3
"""Check the #anchors of the documentation's internal links against the ids of the built pages.

Why: VitePress fails the build on a link to a missing page, but not on a link to a missing anchor of an existing page
(a renamed heading). Run after `npm run docs:build`; exit status 1, with the broken links, if any.

Usage: python3 scripts/check_docs_links.py [--docs docs]
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--docs', type=Path, default=Path(__file__).resolve().parents[1] / 'docs')
    args = ap.parse_args()
    dist = args.docs / '.vitepress' / 'dist'
    if not dist.is_dir():
        print('no built site: run npm run docs:build first')
        return 1
    ids: dict[str, set[str]] = {}
    for f in dist.rglob('*.html'):
        page = '/' + str(f.relative_to(dist)).replace('index.html', '').replace('.html', '')
        ids[page.rstrip('/') or '/'] = set(re.findall(r'id="([^"]+)"', f.read_text(errors='ignore')))
    bad = []
    for f in args.docs.rglob('*.md'):
        rel = f.relative_to(args.docs)
        if rel.parts[0] in ('node_modules', 'api', 'examples') or '.vitepress' in rel.parts:
            continue
        for m in re.finditer(r'\]\(([^)#\s]*)#([^)\s]+)\)', f.read_text()):
            path, anchor = m.group(1), m.group(2)
            if path.startswith(('http://', 'https://')):
                continue
            if path.startswith('/'):
                target = path
            elif path:
                target = os.path.normpath('/' + str(rel.parent) + '/' + path)
            else:
                target = '/' + str(rel.with_suffix(''))
            target = target.replace('.md', '').rstrip('/')
            target = target[:-6] if target.endswith('/index') else target
            target = target or '/'
            if target not in ids:
                bad.append(f'{rel}: {path}#{anchor} (no page {target})')
            elif anchor not in ids[target]:
                bad.append(f'{rel}: {path}#{anchor} (no anchor on {target})')
    for b in bad:
        print(b)
    print(f'{len(bad)} broken anchors')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
