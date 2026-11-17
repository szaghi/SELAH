! the flux register accumulates F_right - F_left of stage k, with the stage's effective weight
flux_w = flux_w + st%stage_weight(k) * boundary_flux(u)
flux_u = flux_u + boundary_flux(u) / st%stages()             ! the wrong weights: uniform, for comparison
