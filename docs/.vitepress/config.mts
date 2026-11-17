import { withMermaid } from 'vitepress-plugin-mermaid'

export default withMermaid({
  title: 'SELAH',
  description: 'Stage-wise Evolution Library for Accelerators and HPC: time-evolution operators for block-structured field solvers',
  base: '/SELAH/',
  srcExclude: ['ford.md'],

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
    siteTitle: 'SELAH',
    nav: [
      { text: 'Home',   link: '/' },
      { text: 'Guide',  link: '/guide/' },
      { text: 'API',    link: '/api/' },
      { text: 'GitHub', link: 'https://github.com/szaghi/SELAH' },
    ],

    sidebar: {
      '/guide/': [
        {
          text: 'Introduction',
          items: [
            { text: 'About',        link: '/guide/' },
            { text: 'Installation', link: '/guide/installation' },
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
