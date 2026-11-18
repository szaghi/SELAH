import { existsSync, readFileSync } from 'node:fs'
import { resolve } from 'node:path'
import { withMermaid } from 'vitepress-plugin-mermaid'

// formal's sidebar of the API pages (docs/api/_sidebar.json, written by `npm run docs:api`), when the API has been
// generated; its one group (the source directory, src/lib) is shown open, as "Modules"
const apiSidebarFile = resolve(__dirname, '../api/_sidebar.json')
const apiModules = existsSync(apiSidebarFile)
  ? JSON.parse(readFileSync(apiSidebarFile, 'utf-8')).flatMap((group: { items: unknown[] }) => group.items)
  : []

export default withMermaid({
  title: 'SELAH',
  description: 'Stage-wise Evolution Library for Accelerators and HPC: time-evolution operators for block-structured field solvers',
  base: '/SELAH/',
  // ford.md is the API generator's project file; _tables/ are the generated fragments the pages include
  srcExclude: ['ford.md', '**/_tables/**', 'examples/**', 'brand/**'],
  // dead links fail the build; SELAH_DOCS_DRAFT=1 tolerates them while a refactor is half done (never in CI)
  ignoreDeadLinks: process.env.SELAH_DOCS_DRAFT === '1',
  // VitePress 1.6 prefetches the page of every link in view through link.pathname, which the SVG links of the API's
  // Mermaid diagrams lack (a TypeError in production; fixed in VitePress 2, vuejs/vitepress#5442)
  router: { prefetchLinks: false },

  // the identity (Hamon in Tsuyukusa, docs/brand/): favicons by colour scheme with PNG fallbacks, the touch icon, the
  // theme colour; and the theme-aware videos of the verification pages (<video class="light-only"> / "dark-only").
  // VitePress does not add the base to head hrefs: they carry /SELAH/ explicitly.
  head: [
    ['link', { rel: 'icon', type: 'image/svg+xml', href: '/SELAH/favicon.svg', media: '(prefers-color-scheme: light)' }],
    ['link', { rel: 'icon', type: 'image/svg+xml', href: '/SELAH/favicon-dark.svg', media: '(prefers-color-scheme: dark)' }],
    ['link', { rel: 'icon', type: 'image/png', sizes: '32x32', href: '/SELAH/favicon-32.png' }],
    ['link', { rel: 'icon', type: 'image/png', sizes: '16x16', href: '/SELAH/favicon-16.png' }],
    ['link', { rel: 'apple-touch-icon', sizes: '180x180', href: '/SELAH/apple-touch-icon.png' }],
    ['meta', { name: 'theme-color', content: '#f3ede1', media: '(prefers-color-scheme: light)' }],
    ['meta', { name: 'theme-color', content: '#0e1116', media: '(prefers-color-scheme: dark)' }],
    ['style', {}, '.dark-only{display:none}.dark .dark-only{display:block}.dark .light-only{display:none}'],
  ],

  markdown: {
    math: true,
    languages: ['fortran-free-form', 'fortran-fixed-form'],
    languageAlias: {
      fortran: 'fortran-free-form',
      f90: 'fortran-free-form',
      f03: 'fortran-free-form',
      f08: 'fortran-free-form',
    },
  },

  themeConfig: {
    logo: { light: '/logo-nav.svg', dark: '/logo-nav-dark.svg', alt: 'SELAH' },
    siteTitle: 'SELAH',
    nav: [
      { text: 'Home',         link: '/' },
      { text: 'Guide',        link: '/guide/' },
      { text: 'Tutorial',     link: '/tutorial/' },
      { text: 'Cookbook',     link: '/cookbook/' },
      { text: 'Verification', link: '/verification/' },
      { text: 'API',          link: '/api/' },
      { text: 'GitHub',       link: 'https://github.com/szaghi/SELAH' },
    ],

    sidebar: {
      '/guide/': [
        {
          text: 'Introduction',
          items: [
            { text: 'About',        link: '/guide/' },
            { text: 'Installation', link: '/guide/installation' },
            { text: 'Quick start',  link: '/guide/quickstart' },
          ],
        },
        {
          text: 'Reference',
          items: [
            { text: 'Concepts',  link: '/guide/concepts' },
            { text: 'Schemes',   link: '/guide/schemes' },
            { text: 'Kernels',   link: '/guide/kernels' },
            { text: 'Stepper',   link: '/guide/stepper' },
            { text: 'Precision', link: '/guide/precision' },
            { text: 'API',       link: '/api/' },
          ],
        },
        {
          text: 'Quality',
          items: [
            { text: 'Verification', link: '/verification/' },
            { text: 'Performance',  link: '/guide/performance' },
            { text: 'Sources',      link: '/guide/sources' },
          ],
        },
        {
          text: 'Project',
          items: [
            { text: 'Contributing', link: '/guide/contributing' },
            { text: 'Changelog',    link: '/guide/changelog' },
          ],
        },
      ],
      '/tutorial/': [
        {
          text: 'Tutorial',
          items: [
            { text: 'Overview',                  link: '/tutorial/' },
            { text: '1. A step by hand',         link: '/tutorial/01-by-hand' },
            { text: '2. The catalogue',          link: '/tutorial/02-catalogue' },
            { text: '3. A field and a box',      link: '/tutorial/03-box' },
            { text: '4. The stepper',            link: '/tutorial/04-stepper' },
            { text: '5. Low-storage forms',      link: '/tutorial/05-low-storage' },
            { text: '6. Splittings',             link: '/tutorial/06-splittings' },
            { text: '7. Masks',                  link: '/tutorial/07-masks' },
            { text: '8. The flux register',      link: '/tutorial/08-flux-register' },
            { text: '9. On the GPU',             link: '/tutorial/09-gpu' },
            { text: '10. Production',            link: '/tutorial/10-production' },
          ],
        },
      ],
      '/cookbook/': [
        {
          text: 'Cookbook',
          items: [
            { text: 'Overview',     link: '/cookbook/' },
            { text: 'Schemes',      link: '/cookbook/schemes' },
            { text: 'Storage',      link: '/cookbook/storage' },
            { text: 'AMR',          link: '/cookbook/amr' },
            { text: 'Diagnostics',  link: '/cookbook/diagnostics' },
            { text: 'Porting from ADAM', link: '/cookbook/porting' },
            { text: 'Performance',  link: '/cookbook/performance' },
          ],
        },
      ],
      '/api/': [
        { text: 'API', items: [
          { text: 'Overview', link: '/api/' },
        ] },
        { text: 'Modules', items: apiModules },
      ],
      '/verification/': [
        { text: 'Verification', items: [
          { text: 'Catalogue',     link: '/verification/' },
          { text: 'Test programs', link: '/verification/unit-tests' },
          { text: 'Regression',    link: '/verification/regression' },
        ] },
        { text: 'Accuracy', collapsed: false, items: [
          { text: 'Overview',                       link: '/verification/accuracy/' },
          { text: 'A1 Order conditions',            link: '/verification/accuracy/order-conditions' },
          { text: 'A2 Orders on scalar ODEs',       link: '/verification/accuracy/scalar-orders' },
          { text: 'A3 Orders through the stepper',  link: '/verification/accuracy/stepper-orders' },
          { text: 'A4 Splittings',                  link: '/verification/accuracy/splittings' },
          { text: 'A5 Linear stability',            link: '/verification/accuracy/stability' },
        ] },
        { text: 'Properties', collapsed: false, items: [
          { text: 'Overview',                       link: '/verification/properties/' },
          { text: 'P1 Strong stability (TVD)',      link: '/verification/properties/ssp-tvd' },
          { text: 'P2 SSP coefficients',            link: '/verification/properties/ssp-coefficients' },
          { text: 'P3 Splitting energy',            link: '/verification/properties/splitting-energy' },
          { text: 'P4 Conservation',                link: '/verification/properties/conservation' },
        ] },
        { text: 'Contracts', collapsed: false, items: [
          { text: 'Overview',                       link: '/verification/contracts/' },
          { text: 'C1 Equivalence with ADAM',       link: '/verification/contracts/adam-equivalence' },
          { text: 'C2 Storage forms',               link: '/verification/contracts/forms' },
          { text: 'C3 Stage contract',              link: '/verification/contracts/stage-contract' },
          { text: 'C4 Golden regression',           link: '/verification/contracts/golden' },
          { text: 'C5 Device',                      link: '/verification/contracts/device' },
          { text: 'C6 Mutation testing',            link: '/verification/contracts/mutation' },
          { text: 'C7 Precision',                   link: '/verification/contracts/precision' },
        ] },
      ],
    },

    socialLinks: [
      { icon: 'github', link: 'https://github.com/szaghi/SELAH' },
    ],

    footer: {
      message: 'GPL-3.0-or-later, BSD-2-Clause, BSD-3-Clause or MIT, at your choice',
      copyright: '© 2026 Stefano Zaghi',
    },

    search: { provider: 'local' },
  },

  mermaid: {},

  vite: {
    // Build with an explicit modern JS target so the docs compile regardless of
    // which mermaid/vitepress/esbuild versions npm resolves. Vite's default
    // es2020 target forces esbuild to down-level modern syntax (e.g. the
    // destructuring mermaid 11.16+ emits), which it refuses to do and the build
    // dies. es2022 needs no lowering and is within VitePress's browser floor.
    build: {
      target: 'es2022',
    },
    // The dev server pre-bundles the dependencies with its own esbuild target (es2020), ignoring build.target:
    // without this, `vitepress dev` dies on the destructuring of mark.js and @vueuse. Bundling mermaid whole also
    // covers the CJS dependencies vitepress-plugin-mermaid does not list (fastdom since mermaid 11.17), which leave
    // a white dev page otherwise.
    optimizeDeps: {
      include: ['mermaid'],
      esbuildOptions: {
        target: 'es2022',
      },
    },
  },
})
