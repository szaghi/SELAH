// docs/.vitepress/theme/index.ts
// The redirects of moved pages (redirects.ts, client side: GitHub Pages has no server redirects), and
// the SELAH identity (custom.css: Hamon in Tsuyukusa, the Sumi family tokens of the identity bundle, docs/brand/) plus the dark figures of the ADAM docs
// harness (both copied from ADAM's docs/.vitepress/theme/index.ts): in the dark theme a page figure <img src=".../x.png"> shows x-dark.png, the silver-ink twin
// drawn by docs/scripts/figure_theme.py, when it exists (dark-figures.data.ts), and stays a washi card otherwise
// (custom.css). The verification catalog, cards and casts are styled in custom.css.
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

const Layout = {
  setup() {
    const { page, site, isDark } = useData()
    const route = useRoute()
    watch(() => [route.path, page.value.isNotFound],
          () => nextTick(() => redirect(!!page.value.isNotFound, site.value.base)),
          { immediate: true })
    onContentUpdated(() => swapFigures(isDark.value, site.value.base))
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
