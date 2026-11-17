// docs/.vitepress/theme/index.ts
// The redirects of moved pages (redirects.ts, client side: GitHub Pages has no server redirects), and
// the SELAH identity (custom.css: Hamon in Tsuyukusa, the Sumi family tokens of the identity bundle, docs/brand/) plus the dark figures of the ADAM docs
// harness (both copied from ADAM's docs/.vitepress/theme/index.ts): in the dark theme a page figure <img src=".../x.png"> shows x-dark.png, the silver-ink twin
// drawn by docs/scripts/figure_theme.py, when it exists (dark-figures.data.ts), and stays a washi card otherwise
// (custom.css). The verification catalog, cards and casts are styled in custom.css. Wide Mermaid diagrams shrink to the
// column down to 0.75 of their size, and scroll past it (floorDiagrams).
import DefaultTheme from 'vitepress/theme'
import { h, nextTick, watch } from 'vue'
import { onContentUpdated, useData, useRoute, withBase } from 'vitepress'
import { pageRedirects } from './redirects'
import { data as darkFigures } from './dark-figures.data'
import SelahIcon from './SelahIcon.vue'
import './custom.css'

// The twins are swapped on the client only, after the content is mounted (the server-rendered page always carries the
// light figure, so hydration never mismatches): on every content update (first mount, route change) and theme change.
// The light source is kept in data-light-src; a twin that fails to load falls back to the washi card.
const twins = new Set<string>(darkFigures)
const TWIN = 'selah-dark-twin'

function swapFigures(dark: boolean, base: string) {
  if (typeof document === 'undefined') return
  document.querySelectorAll<HTMLImageElement>('.vp-doc img').forEach((img) => {
    const light = img.dataset.lightSrc ?? img.getAttribute('src') ?? ''
    if (!img.dataset.lightSrc) {
      const path = new URL(light, window.location.href).pathname
      if (!path.startsWith(base) || !twins.has(path.slice(base.length))) return
      img.dataset.lightSrc = light
    }
    if (dark && !img.classList.contains(TWIN)) {
      img.addEventListener('error', () => {
        img.classList.remove(TWIN)
        img.setAttribute('src', light)
      }, { once: true })
      img.classList.add(TWIN)
      img.setAttribute('src', light.replace(/\.png$/, '-dark.png'))
    } else if (!dark && img.classList.contains(TWIN)) {
      img.classList.remove(TWIN)
      img.setAttribute('src', light)
    }
  })
}

// '/SELAH/guide/verification.html' -> '/guide/verification'
function sitePath(pathname: string, base: string): string {
  let p = pathname.startsWith(base) ? '/' + pathname.slice(base.length) : pathname
  p = p.replace(/\.html$/, '').replace(/\/index$/, '/')
  if (p.length > 1 && p.endsWith('/')) p = p.slice(0, -1)
  return p || '/'
}

// A moved page renders the 404 page; replace it by its new location (redirects.ts), keeping the anchor.
function redirect(isNotFound: boolean, base: string) {
  if (typeof window === 'undefined' || !isNotFound) return
  const target = pageRedirects[sitePath(window.location.pathname, base)]
  if (target) window.location.replace(withBase(target) + (window.location.hash || ''))
}

// Mermaid draws each diagram with width 100% and its natural width as max-width, so a wide one shrinks to the column
// and its labels with it (a 1800 px call graph came out with 6 px labels). A floor on the shrink: min-width at
// MIN_SCALE of the natural width, past which the diagram keeps its size and its container scrolls (custom.css).
// The diagrams are rendered on the client after mount and again on a theme change, hence the observer.
const MIN_SCALE = 0.75

function floorDiagrams(root: ParentNode) {
  root.querySelectorAll<SVGSVGElement>('.mermaid svg').forEach((svg) => {
    const natural = svg.viewBox?.baseVal?.width
    if (natural) svg.style.minWidth = `${Math.round(MIN_SCALE * natural)}px`
  })
}

let diagramObserver: MutationObserver | undefined
function watchDiagrams() {
  if (typeof document === 'undefined' || diagramObserver) return
  floorDiagrams(document)
  diagramObserver = new MutationObserver((records) => {
    for (const r of records) r.addedNodes.forEach((n) => { if (n instanceof Element) floorDiagrams(n.parentElement ?? n) })
  })
  diagramObserver.observe(document.body, { childList: true, subtree: true })
}

const Layout = {
  setup() {
    const { page, site, isDark } = useData()
    const route = useRoute()
    watch(() => [route.path, page.value.isNotFound],
          () => nextTick(() => redirect(!!page.value.isNotFound, site.value.base)),
          { immediate: true })
    onContentUpdated(() => { swapFigures(isDark.value, site.value.base); watchDiagrams() })
    watch(isDark, (dark) => nextTick(() => swapFigures(dark, site.value.base)))
    // a page with `icon: <name>` in its frontmatter shows that icon of the identity above its title
    const { frontmatter } = useData()
    return () => h(DefaultTheme.Layout, null, {
      'doc-before': () => frontmatter.value.icon
        ? h('div', { class: 'selah-page-icon' }, [h(SelahIcon, { name: frontmatter.value.icon, size: 56 })])
        : null,
    })
  },
}

export default {
  extends: DefaultTheme,
  Layout,
  enhanceApp({ app }) {
    app.component('SelahIcon', SelahIcon)
  },
}
