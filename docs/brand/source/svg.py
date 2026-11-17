from brush import amr, seal_d
SUMI_DEFAULT = False

def mark_svg(ink, col, size=200, mesh_levels=3, seal=True, label="", vb="0 0 64 64", sumi=None):
    sumi = SUMI_DEFAULT if sumi is None else sumi
    from brush import Ink
    parts = []
    if seal:
        # the seal is ink for the AMR as well
        for x in range(53, 62, 2):
            for y in range(53, 62, 2):
                ink.pts.append((x, y, 0.5))
    if mesh_levels:
        from brush import amr_sumi
        layers = amr_sumi(ink)[:mesh_levels] if sumi else amr(ink, levels=mesh_levels)
        for d, op in layers:
            parts.append(f'<path fill="{col["mesh"]}" opacity="{op}" d="{d}"/>')
    order = ["accb", "ink55", "ink", "ground", "acc"]
    for role in order:
        for (r, eo), ds in ink.layers.items():
            if r != role:
                continue
            fill = {"accb": col["acc"], "ink55": col["ink"], "ink": col["ink"], "ground": col["ground"], "acc": col["acc"]}[role]
            op = ' opacity="0.55"' if role == "ink55" else ""
            parts.append(f'<path fill="{fill}"{op} d="{" ".join(ds)}"/>')
    if seal:
        parts.append(f'<path fill="{col["acc"]}" fill-rule="evenodd" d="{seal_d()}"/>')
    aria = f' role="img" aria-label="{label}"' if label else ' aria-hidden="true"'
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vb}" width="{size}" height="{size}"{aria}>' + "".join(parts) + "</svg>"

def icon_svg(ink, col, size=40, label=""):
    parts = []
    for role in ["ink55", "ink", "acc"]:
        for (r, eo), ds in ink.layers.items():
            if r != role:
                continue
            fill = col["acc"] if role == "acc" else col["ink"]
            op = ' opacity="0.55"' if role == "ink55" else ""
            parts.append(f'<path fill="{fill}"{op} d="{" ".join(ds)}"/>')
    aria = f' role="img" aria-label="{label}"' if label else ' aria-hidden="true"'
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" width="{size}" height="{size}"{aria}>' + "".join(parts) + "</svg>"
