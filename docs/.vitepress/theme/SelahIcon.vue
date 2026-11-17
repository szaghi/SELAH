<!-- docs/.vitepress/theme/SelahIcon.vue: one icon of the identity (docs/public/icons/<name>.svg and <name>-dark.svg,
     Hamon in Tsuyukusa), drawn in the theme of the page. Used inline in Markdown (<SelahIcon name="ssp" />, registered
     globally by index.ts; alone on a line before a paragraph with `float`, the paragraph wraps around it) and as the
     page icon of the guide pages (frontmatter `icon: <name>`, the doc-before slot). -->
<script setup lang="ts">
import { withBase } from 'vitepress'

const props = withDefaults(defineProps<{ name: string; size?: number; label?: string; float?: boolean }>(),
                           { size: 28, label: '', float: false })
</script>

<template>
  <span :class="['selah-icon', { 'selah-icon-float': props.float }]" :style="{ width: props.size + 'px', height: props.size + 'px' }" :role="props.label ? 'img' : undefined"
        :aria-label="props.label || undefined" :aria-hidden="props.label ? undefined : 'true'">
    <img class="selah-icon-light" :src="withBase(`/icons/${props.name}.svg`)" alt="" :width="props.size" :height="props.size" />
    <img class="selah-icon-dark" :src="withBase(`/icons/${props.name}-dark.svg`)" alt="" :width="props.size" :height="props.size" />
  </span>
</template>
