"""The colours of SELAH's documentation figures and animations: the identity (Hamon in Tsuyukusa) applied to data.

Why this module exists: every figure script (make_doc_figures, make_convention_figures, make_performance_figures) and
the animations (make_doc_animations) take their colours from here, once.

The grounds and inks are the Sumi family's, as in docs/brand/tokens.json (washi light, silver ink dark). The accent is
the identity's ink, Tsuyukusa (dayflower blue): #1f6fc0 on washi, #6fb8ff on silver ink, for fills only, as the drop
and the seal of the mark. The categorical data lines are those validated for SETH's docs on these same grounds with the
dataviz skill's validate_palette.js (lightness band, chroma floor, colour-blind separation of adjacent pairs,
contrast); their first hue is a blue next to Tsuyukusa, so the SSP family, the core of the catalogue, wears the
identity's colour.

The static figures are drawn twice by their scripts, with the light colours and, inside figure_theme.theme("dark")
(ADAM's harness, unchanged), with the dark ones (STYLE_DARK, LINES_DARK); the animations draw each theme with its own.
"""

from __future__ import annotations

# grounds and inks (custom.css)
WASHI, WASHI_ALT, INK_LIGHT, MESH_LIGHT = '#f3ede1', '#ebe3d3', '#1d1b18', '#6f6a62'
SILVER, SILVER_ALT, INK_DARK, MESH_DARK = '#0e1116', '#0a0d11', '#e8ecf0', '#8b949e'
ACCENT_LIGHT, ACCENT_DARK = '#1f6fc0', '#6fb8ff'  # Tsuyukusa (custom.css --selah-drop)

# categorical data colours, in fixed order (never cycled): blue, green, plum, umber
BLUE, GREEN, PLUM, UMBER = '#1f63a8', '#4a7a06', '#8a4f9e', '#a8661c'
BLUE_D, GREEN_D, PLUM_D, UMBER_D = '#3f8fd6', '#6aa41c', '#b06ac4', '#c27a2a'
LINES_LIGHT = (BLUE, GREEN, PLUM, UMBER)
LINES_DARK = (BLUE_D, GREEN_D, PLUM_D, UMBER_D)

# one colour and marker per scheme: the family is the hue (the SSP Runge-Kutta schemes blue, the classical rk4 umber,
# the low-storage 2N schemes green, the splittings plum), the order a darker shade, the marker the family again
STYLE = {
    'euler': ('0.55', '.'), 'ssp-22': ('#8db6e0', 'o'), 'ssp-33': ('#4f8fcf', 'o'), 'ssp-54': (BLUE, 'o'),
    'ssp-104': ('#123e6b', 'o'), 'rk4': (UMBER, 'D'), 'lsrk-33': ('#76a52a', 's'), 'lsrk-54': (GREEN, 's'),
    'strang': ('#c49bd2', '^'), 'yoshida-4': (PLUM, '^'), 'blanes-moan-s6': ('#6b3a7d', 'v'),
    'blanes-moan-s10': ('#45214f', 'v'),
}
STYLE_DARK = {
    'euler': ('0.62', '.'), 'ssp-22': ('#9cc6f0', 'o'), 'ssp-33': ('#6aa8e6', 'o'), 'ssp-54': (BLUE_D, 'o'),
    'ssp-104': ('#2f77c0', 'o'), 'rk4': (UMBER_D, 'D'), 'lsrk-33': ('#8fc23e', 's'), 'lsrk-54': (GREEN_D, 's'),
    'strang': ('#dcb6e8', '^'), 'yoshida-4': (PLUM_D, '^'), 'blanes-moan-s6': ('#c98ad9', 'v'),
    'blanes-moan-s10': ('#a35fb6', 'v'),
}

# the storage forms (register maps, convention figures): one hue each
FORM_COLOUR = {'butcher': BLUE, 'shu-osher-2r': GREEN, 'williamson-2n': PLUM, 'ketcheson-2r': UMBER, 'splitting': '0.45'}
FORM_COLOUR_DARK = {'butcher': BLUE_D, 'shu-osher-2r': GREEN_D, 'williamson-2n': PLUM_D, 'ketcheson-2r': UMBER_D,
                    'splitting': '0.7'}
