import os, json, shutil, subprocess, math
import svg
from svg import mark_svg, icon_svg
from marks5 import MARKS5
from icons import pack
from pal import cr, shade
from pal6 import PAL6
B = "/home/claude/bundle/selah-identity"
shutil.rmtree(B, ignore_errors=True)
pub = f"{B}/docs/public"; os.makedirs(f"{pub}/icons", exist_ok=True); os.makedirs(f"{B}/docs/.vitepress/theme", exist_ok=True)
p = PAL6["tsuyukusa"]
WASHI = dict(ground="#f3ede1", raised="#f7f2e8", line="#ddd4c3", ink="#1d1b18", muted="#5f5a52", mesh="#6f6a62")
SILVER = dict(ground="#0e1116", raised="#141820", line="#222831", ink="#e8ecf0", muted="#a9b1bb", mesh="#8b949e")
fn = MARKS5["hamon"]
cw = dict(mesh=WASHI["mesh"], ink=WASHI["ink"], ground=WASHI["ground"], acc=p["acc"])
cn = dict(mesh=SILVER["mesh"], ink=SILVER["ink"], ground=SILVER["ground"], acc=p["accn"])
def clean(s):
    return s.replace(' width="200" height="200"', '')
def write(path, s):
    open(path, "w").write(s + "\n")
# logos
svg.SUMI_DEFAULT = True
write(f"{pub}/logo.svg", mark_svg(fn(), cw, 256, label="SELAH: one drop, four rings on still water"))
write(f"{pub}/logo-dark.svg", mark_svg(fn(), cn, 256, label="SELAH: one drop, four rings on still water"))
write(f"{pub}/logo-nav.svg", mark_svg(fn(1.6), cw, 64, mesh_levels=0, seal=False, label="SELAH", vb="2 2 58 58"))
write(f"{pub}/logo-nav-dark.svg", mark_svg(fn(1.6), cn, 64, mesh_levels=0, seal=False, label="SELAH", vb="2 2 58 58"))
# favicons: the bare mark, rings thickened, on the family tile
def fav(col, T, ws=2.1):
    inner = mark_svg(fn(ws), col, 64, mesh_levels=0, seal=False, label="SELAH", vb="2 2 58 58")
    body = inner[inner.index(">") + 1:-len("</svg>")]
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64" width="64" height="64" role="img" aria-label="SELAH">'
            f'<rect width="64" height="64" rx="12" fill="{T["raised"]}"/><g transform="translate(-0.4 0) scale(1.1) translate(-2.9 -2.9)">{body}</g></svg>')
write(f"{pub}/favicon.svg", fav(cw, WASHI))
write(f"{pub}/favicon-dark.svg", fav(cn, SILVER))
for sz in (16, 32):
    subprocess.run(["rsvg-convert", "-w", str(sz), "-h", str(sz), f"{pub}/favicon.svg", "-o", f"{pub}/favicon-{sz}.png"], check=True)
at = fav(cw, WASHI).replace('rx="12" ', '')
open("/tmp/at.svg", "w").write(at)
subprocess.run(["rsvg-convert", "-w", "180", "-h", "180", "-b", WASHI["raised"], "/tmp/at.svg", "-o", f"{pub}/apple-touch-icon.png"], check=True)
# icons
names = {"Butcher tableaux": "tableaux", "Strong stability (SSP)": "ssp", "Low storage": "low-storage", "Operator splitting": "splitting",
         "Interior-only kernels": "kernels", "CPU and GPU": "cpu-gpu", "Stage contract": "stepper", "Order conditions": "order-conditions"}
icons = pack("ring")
for lab, ink in icons:
    write(f"{pub}/icons/{names[lab]}.svg", icon_svg(ink, dict(ink=WASHI["ink"], acc=p["acc"]), 24, lab))
    write(f"{pub}/icons/{names[lab]}-dark.svg", icon_svg(ink, dict(ink=SILVER["ink"], acc=p["accn"]), 24, lab))
write(f"{pub}/icons/selah.svg", mark_svg(fn(1.6), cw, 24, mesh_levels=0, seal=False, label="SELAH", vb="2 2 58 58"))
write(f"{pub}/icons/selah-dark.svg", mark_svg(fn(1.6), cn, 24, mesh_levels=0, seal=False, label="SELAH", vb="2 2 58 58"))
# colours
brand_w, drop_w = p["link"], p["acc"]
press_w = shade(brand_w, 0.8)
brand_n, drop_n = p["linkn"], p["accn"]
press_n = shade(drop_n, 0.85)
on_n = "#06121f"
C = {
 "washi": {"drop": drop_w, "link": brand_w, "button": brand_w, "button-pressed": press_w, "on-button": "#ffffff"},
 "silver-ink": {"drop": drop_n, "link": brand_n, "button": drop_n, "button-pressed": press_n, "on-button": on_n},
}
contrast = {
 "washi": {"ink": cr(WASHI["ink"], WASHI["ground"]), "muted": cr(WASHI["muted"], WASHI["ground"]), "link": cr(brand_w, WASHI["ground"]),
           "link on ground-alt": cr(brand_w, "#ebe3d3"), "button text": cr("#ffffff", brand_w), "button text pressed": cr("#ffffff", press_w),
           "drop and seal": cr(drop_w, WASHI["ground"]), "control borders": cr(WASHI["mesh"], WASHI["ground"])},
 "silver-ink": {"ink": cr(SILVER["ink"], SILVER["ground"]), "muted": cr(SILVER["muted"], SILVER["ground"]), "link": cr(brand_n, SILVER["ground"]),
           "link on ground-alt": cr(brand_n, "#0a0d11"), "button text": cr(on_n, drop_n), "button text pressed": cr(on_n, press_n),
           "drop and seal": cr(drop_n, SILVER["ground"]), "control borders": cr(SILVER["mesh"], SILVER["ground"])},
}
tokens = {"name": "SELAH · Hamon in Tsuyukusa", "family": "Sumi family", "ink": {"name": "Tsuyukusa", "kanji": "露草", "gloss": p["gloss"]},
          "themes": {"washi": {**WASHI, **C["washi"]}, "silver-ink": {**SILVER, **C["silver-ink"]}},
          "contrast": {t: {k: round(v, 2) for k, v in d.items()} for t, d in contrast.items()}}
json.dump(tokens, open(f"{B}/tokens.json", "w"), ensure_ascii=False, indent=2)
print(json.dumps(tokens["contrast"], indent=1))
def rgba(h, a):
    return f"rgba({int(h[1:3],16)}, {int(h[3:5],16)}, {int(h[5:7],16)}, {a})"
css = f"""/* SELAH · Hamon 波紋 in Tsuyukusa 露草 — Sumi family tokens for VitePress.
   Washi is :root, silver ink is .dark. Every text pair is at least 4.5:1 on its ground.
   Dayflower blue {drop_n} is the night ink; on washi it deepens to {drop_w} (drop, seal)
   and {brand_w} (links, primary button). */
@import url('https://fonts.googleapis.com/css2?family=Shippori+Mincho:wght@800&family=Inter:wght@400;500;600&family=JetBrains+Mono:wght@400;500&display=swap');

:root {{
  --vp-font-family-base: Inter, ui-sans-serif, system-ui, -apple-system, 'Segoe UI', sans-serif;
  --vp-font-family-mono: 'JetBrains Mono', ui-monospace, 'SF Mono', Menlo, monospace;
  --selah-font-display: 'Shippori Mincho', 'Hiragino Mincho ProN', serif;
  /* ground */
  --vp-c-bg: #f3ede1;
  --vp-c-bg-alt: #ebe3d3;
  --vp-c-bg-elv: #f7f2e8;
  --vp-c-bg-soft: #ebe3d3;
  --vp-c-divider: #ddd4c3;
  --vp-c-gutter: #ddd4c3;
  --vp-c-border: #6f6a62;
  /* ink */
  --vp-c-text-1: #1d1b18;
  --vp-c-text-2: #5f5a52;
  --vp-c-text-3: #5f5a52;
  /* tsuyukusa: links, active nav, primary button */
  --vp-c-brand-1: {brand_w};
  --vp-c-brand-2: {press_w};
  --vp-c-brand-3: {brand_w};
  --vp-c-brand-soft: {rgba(drop_w, 0.14)};
  --vp-button-brand-bg: {brand_w};
  --vp-button-brand-border: {brand_w};
  --vp-button-brand-text: #ffffff;
  --vp-button-brand-hover-bg: {press_w};
  --vp-button-brand-hover-border: {press_w};
  --vp-button-brand-hover-text: #ffffff;
  --vp-button-brand-active-bg: {press_w};
  --vp-button-brand-active-border: {press_w};
  --vp-button-brand-active-text: #ffffff;
  --vp-button-alt-border: #6f6a62;
  /* the drop and the seal: fills only, never text */
  --selah-drop: {drop_w};
  --selah-seal: {drop_w};
  --selah-focus: #1d1b18;
  --vp-home-hero-name-color: #1d1b18;
  --vp-code-block-bg: #ebe3d3;
  /* no glow behind the mark: the ink carries it */
  --vp-home-hero-image-background-image: none;
  --vp-home-hero-image-filter: none;
}}

.dark {{
  /* ground */
  --vp-c-bg: #0e1116;
  --vp-c-bg-alt: #0a0d11;
  --vp-c-bg-elv: #141820;
  --vp-c-bg-soft: #0a0d11;
  --vp-c-divider: #222831;
  --vp-c-gutter: #222831;
  --vp-c-border: #8b949e;
  /* ink */
  --vp-c-text-1: #e8ecf0;
  --vp-c-text-2: #a9b1bb;
  --vp-c-text-3: #a9b1bb;
  /* tsuyukusa: links, active nav, primary button */
  --vp-c-brand-1: {brand_n};
  --vp-c-brand-2: {drop_n};
  --vp-c-brand-3: {drop_n};
  --vp-c-brand-soft: {rgba(drop_n, 0.14)};
  --vp-button-brand-bg: {drop_n};
  --vp-button-brand-border: {drop_n};
  --vp-button-brand-text: {on_n};
  --vp-button-brand-hover-bg: {press_n};
  --vp-button-brand-hover-border: {press_n};
  --vp-button-brand-hover-text: {on_n};
  --vp-button-brand-active-bg: {press_n};
  --vp-button-brand-active-border: {press_n};
  --vp-button-brand-active-text: {on_n};
  --vp-button-alt-border: #8b949e;
  /* the drop and the seal: fills only */
  --selah-drop: {drop_n};
  --selah-seal: {drop_n};
  --selah-focus: {brand_n};
  --vp-home-hero-name-color: #e8ecf0;
  --vp-code-block-bg: #0a0d11;
}}

/* names in the display face, tracked wide */
.VPNavBarTitle .title,
.VPHero .name,
.VPHero .name .clip {{
  font-family: var(--selah-font-display);
  font-weight: 800;
  letter-spacing: 0.04em;
}}
.VPHero .name .clip {{
  background: none;
  -webkit-text-fill-color: var(--vp-home-hero-name-color);
  color: var(--vp-home-hero-name-color);
}}

/* the mark in the hero: square, no shadow */
.VPHero .image-src {{
  max-width: 320px;
  max-height: 320px;
}}

/* feature icons are drawn at 24 units: show them at 44 px, no tile */
.VPFeature .icon {{
  width: 52px;
  height: 52px;
  background-color: transparent;
}}

/* keyboard focus: solid 2px ring with 2px of ground */
:focus-visible {{
  outline: 2px solid var(--selah-focus);
  outline-offset: 2px;
}}
"""
write(f"{B}/docs/.vitepress/theme/custom.css", css.rstrip())
write(f"{B}/docs/.vitepress/theme/index.ts", "// SELAH docs theme: VitePress default theme + the Sumi family tokens in the Tsuyukusa ink.\nimport DefaultTheme from 'vitepress/theme'\nimport './custom.css'\n\nexport default DefaultTheme")
# source snapshot
os.makedirs(f"{B}/source", exist_ok=True)
for f in ["brush.py", "svg.py", "icons.py", "marks2.py", "marks4.py", "marks5.py", "pal.py", "pal6.py", "bundle.py"]:
    shutil.copy(f, f"{B}/source/{f}")
print(sorted(os.listdir(pub)), len(os.listdir(f"{pub}/icons")))
