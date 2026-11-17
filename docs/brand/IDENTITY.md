# SELAH · new identity: Hamon 波紋 in Tsuyukusa 露草

The mark is one drop falling on still water and the four rings it sends out: one evaluation, and the
stages spreading from it. Each ring is a single brush stroke, pressed in, left open at a different
point and paler as it widens, with the spacing growing outward as on real water. The drop falls along
a faint broken trace. It carries the Sumi family seal (the refined quadrant) and the brushed
three-level AMR ground, painted one stroke per cell edge and refined only under the rings. Ink:
**Tsuyukusa**, dayflower blue `#6fb8ff` on silver ink (the dark theme); on washi (the light theme) it
deepens to `#1f6fc0` for the drop and the seal, and `#1b62aa` for links and the primary button.

The bundle targets a VitePress `docs/` with `base: '/SELAH/'`, as the published site uses.

## What is in it

```
MIGRATION.md                       this file
apply-identity.py                  edits to two tracked files (config.mts, index.md); see below
tokens.json                        Tsuyukusa values for both themes, with contrast on ground
source/                            snapshot of the generator (reference only)
docs/.vitepress/theme/index.ts     new: default theme + custom.css
docs/.vitepress/theme/custom.css   new: Tsuyukusa tokens mapped to VitePress variables
docs/public/logo.svg               hero mark, washi (sumi AMR ground + seal)
docs/public/logo-dark.svg          hero mark, silver ink
docs/public/logo-nav.svg           navbar mark, washi (no ground, no seal, strokes thickened)
docs/public/logo-nav-dark.svg      navbar mark, silver ink
docs/public/favicon.svg            favicon tile, washi: the drop and two thick rings only
docs/public/favicon-dark.svg       favicon tile, silver ink
docs/public/favicon-16.png         PNG fallbacks rendered from favicon.svg
docs/public/favicon-32.png
docs/public/apple-touch-icon.png   180 × 180, opaque (iOS rounds it)
docs/public/icons/*.svg            nine icons, light and -dark
```

## Applying it on trunk

From the repository root, on `trunk` (or whatever your default branch is), with a clean tree:

```sh
unzip selah-identity.zip -d /tmp/selah-identity
cp -r /tmp/selah-identity/docs/. docs/
python3 /tmp/selah-identity/apply-identity.py --check   # dry run
python3 /tmp/selah-identity/apply-identity.py
npx vitepress dev docs                                  # look at it
git add docs
git commit -m "docs: new identity, Hamon in Tsuyukusa (Sumi family)"
git push
```

One commit, straight on trunk, like the other projects of the family. `tokens.json`, `source/` and
this file are for reference and do not need to be committed (put them under `docs/brand/` if you
want them in the repository).

### Why a script and not a patch

The SELAH repository was not reachable from the session that built this bundle, so there was no
trunk to diff against. Instead of a `git diff`, `apply-identity.py` makes the same edits by anchors
and is idempotent (a second run changes nothing; `--check` is a dry run; if an anchor is missing it
stops, names it and writes nothing):

- `docs/.vitepress/config.mts`: adds `head` after `base: '/SELAH/'` (SVG favicons switched by
  `prefers-color-scheme`, PNG fallbacks, apple-touch icon, `theme-color`), or prepends those entries
  to an existing `head`; adds `themeConfig.logo` and `siteTitle: 'SELAH'`. VitePress does not add the
  base to `head` hrefs, so they carry `/SELAH/` explicitly; logo and image paths get it automatically.
- `docs/index.md`: adds the hero image (light/dark) after the tagline (a folded `>-` tagline is
  handled) and sets the brushed SVG icons on the three features, replacing an emoji icon where there
  is one. No text is changed.

| Feature | Icon |
|---|---|
| Verified tableaux | `tableaux` (strictly lower-triangular A, then the weights b) |
| Kernels over your arrays | `kernels` (a block, its halo left untouched) |
| A stepper you drive | `stepper` (begin_stage, end_stage, finish_step) |

The other icons (`ssp`, `low-storage`, `splitting`, `cpu-gpu`, `order-conditions`, and `selah`, the
mark at icon size) are not used by the home page. They are there for the guide: the tableaux
catalogue (SSP, low-storage, splitting), the backends page (CPU and GPU), the verification of order
conditions. Every icon draws its stages with the same small ring.

The script and the files were tested on a stand-in `docs/` with SELAH's home page (hero, the three
features above, one with an emoji icon) and built with VitePress 1: logos, favicons and icons all
resolve under `/SELAH/`.

## Colours

| | Washi | Silver ink |
|---|---|---|
| Drop and seal (fills only) | `#1f6fc0` | `#6fb8ff` |
| Links, active nav | `#1b62aa` | `#8cc6ff` |
| Primary button | white on `#1b62aa` (hover `#164e88`) | `#06121f` on `#6fb8ff` (hover darker) |

Contrast on the page ground (all in `tokens.json`):

- Washi: ink 14.7, muted 5.9, link 5.3 (4.9 on `ground-alt`), button text 6.2 (8.5 on hover),
  drop and seal 4.4, control borders 4.6.
- Silver ink: ink 15.9, muted 8.7, link 10.5, button text 9.0 (6.7 on hover), drop and seal 9.0,
  control borders 6.2.

## Rules, as for the rest of the family

- Do not hand-edit the SVG paths. They are generated geometry; change the drawing only by
  regenerating it (`source/`).
- No gradients, glows or shadows behind the mark: it is flat ink.
  `--vp-home-hero-image-background-image` is set to `none`.
- Use the full logo (with ground and seal) from 48 px up. Below 48 px use `logo-nav*.svg` or the
  favicons. The favicon is a simplified drawing (the drop and two thick rings), because the four thin
  rings do not survive 16 px.
- Only the drop and the seal take the Tsuyukusa ink; the rings stay sumi.
- `source/` is a snapshot of the generator.
