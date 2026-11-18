---
title: Cookbook
---

# Cookbook

Short recipes for specific tasks. The schemes, storage, AMR and diagnostics recipes are each one program in
`docs/examples/src/`, with the recipe chosen on the command line. The code and output on each page are that program
and what it printed when the documentation was built. The porting and performance pages are references, with code
taken from the tests and the benchmarks.

| Page | Recipes |
|---|---|
| [Schemes](./schemes) | choosing a scheme and its storage form; the largest strong-stability-preserving time step |
| [Storage](./storage) | registers from your own pool; the working precision of a build |
| [Any rank](./ranks) | a 1-D, 2-D or 3-D solver: the box and the stepper at the field's rank; a periodic field with no ghost cells |
| [AMR](./amr) | refluxing with `stage_weight`; where the ghost cells of a stage must come from |
| [Diagnostics](./diagnostics) | the step residual as a convergence monitor; contract checks; what happens on a request the stepper cannot honour |
| [Porting from ADAM](./porting) | ADAM's RK calls and their SELAH equivalents, call by call |
| [Performance](./performance) | benchmark discipline on a GPU; the paths that save memory traffic |
