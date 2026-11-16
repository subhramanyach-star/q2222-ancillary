# Rigorous single-round min-entropy intervals (bits) replacing the NPA level-3 numbers of prl_v13

Each cell: exact value lies in [H_lo, H_hi]; a single number means both ends agree to 6 decimals (interval width < 10⁻⁶).
H_lo: exact-rational NPA certificate (level in brackets). H_hi: explicit qubit attack. Paper value (NPA-3 float) below each row.

### θ = π/8: H_min(A₁|X=1,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| tilted CHSH g_θ | 0.776634 | 0.383682 | 0.185275 | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.777 | 0.384 | 0.185 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| paper "best in cone" (its weights) | [0.832439, 0.832440] | 0.438340 | 0.206499 | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.832 | 0.438 | 0.206 | 0.000 | 0.000 |
| *(NPA level used)* | 2 | 2 | 2 | 1+AB | 1+AB |
| **sup over the whole hexagon C_θ** | 0.832443 | 0.438352 | 0.206530 | 0 (exact) | 0 (exact) |
| *(NPA level used)* | 2,3 | 2 | 2 | – | – |
| CHSH | 0.226659 | 0.210974 | 0.194352 | 0.148883 | 0.084669 |
| *(paper NPA-3)* | 0.227 | 0.211 | 0.194 | 0.149 | 0.085 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| full statistics | [0.835087, 0.835088] | 0.514289 | 0.376273 | 0.212696 | 0.103973 |
| *(paper NPA-3)* | 0.835 | 0.514 | 0.376 | 0.213 | 0.104 |
| *(NPA level used)* | 3 | 1+AB | 1+AB | 2 | 2 |

### θ = 27.76° (Liu et al.): H_min(A₁|X=1,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| tilted CHSH g_θ | 0.842291 | 0.548525 | 0.393120 | 0.118590 | 0 (exact) |
| *(paper NPA-3)* | 0.842 | 0.549 | 0.393 | 0.119 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| paper "best in cone" (its weights) | 0.871352 | 0.606007 | 0.439178 | 0.126896 | 0 (exact) |
| *(paper NPA-3)* | 0.871 | 0.606 | 0.439 | 0.127 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 2,3 | 1+AB |
| **sup over the whole hexagon C_θ** | 0.871353 | 0.606026 | 0.439191 | 0.126918 | 0 (exact) |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 3 | – |
| CHSH | 0.350085 | 0.326529 | 0.302144 | 0.237970 | 0.152205 |
| *(paper NPA-3)* | 0.350 | 0.327 | 0.302 | 0.238 | 0.152 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| full statistics | 0.871943 | [0.630677, 0.630678] | 0.505821 | 0.322627 | 0.186098 |
| *(paper NPA-3)* | 0.872 | 0.631 | 0.506 | 0.323 | 0.186 |
| *(NPA level used)* | 1+AB | 2 | 1+AB | 2 | 2 |

### θ = 3π/16 (ancillary): H_min(A₁|X=1,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| tilted CHSH g_θ | 0.884187 | 0.660479 | 0.538610 | 0.318483 | 0.094460 |
| *(paper NPA-3)* | 0.884 | 0.660 | 0.539 | 0.318 | 0.094 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| paper "best in cone" (its weights) | 0.895283 | 0.692396 | 0.578470 | 0.346602 | 0.098757 |
| *(paper NPA-3)* | 0.895 | 0.692 | 0.578 | 0.347 | 0.099 |
| *(NPA level used)* | 1+AB | 1+AB | 2 | 1+AB | 1+AB |
| **sup over the whole hexagon C_θ** | 0.895283 | 0.692398 | 0.578473 | 0.346607 | 0.098841 |
| *(NPA level used)* | 1+AB | 1+AB | 2 | 1+AB | 1+AB |
| CHSH | 0.527518 | 0.485829 | 0.445285 | 0.347412 | 0.229271 |
| *(paper NPA-3)* | 0.528 | 0.486 | 0.445 | 0.347 | 0.229 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| full statistics | [0.895491, 0.895492] | 0.697997 | 0.596080 | 0.421675 | 0.267516 |
| *(paper NPA-3)* | 0.895 | 0.698 | 0.596 | 0.422 | 0.268 |
| *(NPA level used)* | 1+AB | 1+AB | 2 | 1+AB | 2 |

### WBC φ = π/6: H_min(A₀B₀|00,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| I_φ (centre) | 1.729981 | 1.069929 | 0.588463 | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 1.730 | 1.070 | 0.588 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| half-way to W₁ | [1.693991, 1.693992] | 0.606132 | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 1.694 | 0.606 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB,2 | 1+AB,2 | 1+AB | 1+AB | 1+AB |
| vertex W₁ | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| CHSH | 0.439524 | 0.409375 | 0.378186 | 0.296221 | 0.187197 |
| *(paper NPA-3)* | 0.440 | 0.409 | 0.378 | 0.296 | 0.187 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB,2 |
| full statistics | [1.759378, 1.759379] | [1.163943, 1.163944] | 0.846721 | 0.475418 | 0.248684 |
| *(paper NPA-3)* | 1.759 | 1.164 | 0.847 | 0.475 | 0.249 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |

### WBC φ = π/6: H_min(A₀|X=0,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| I_φ (centre) | 0.891970 | 0.610278 | 0.359617 | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.892 | 0.610 | 0.360 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| half-way to W₁ | 0.884285 | 0.372425 | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.884 | 0.372 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB,3 | 1+AB,3 | 1+AB | 1+AB | 1+AB |
| vertex W₁ | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| CHSH | 0.356581 | 0.332529 | 0.307668 | 0.242401 | 0.155460 |
| *(paper NPA-3)* | 0.357 | 0.333 | 0.308 | 0.242 | 0.155 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| full statistics | 0.892805 | 0.669448 | 0.541128 | 0.351942 | 0.200250 |
| *(paper NPA-3)* | 0.893 | 0.669 | 0.541 | 0.352 | 0.200 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 2 | 2 |

### WBC φ = π/12 (ancillary): H_min(A₀B₀|00,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| I_φ (centre) | 1.692537 | 0.644473 | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 1.693 | 0.644 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| half-way to W₁ | 1.573891 | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 1.574 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB,2 | 1+AB | 1+AB | 1+AB | 1+AB |
| vertex W₁ | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| CHSH | 0.260635 | 0.241782 | 0.221783 | 0.167054 | 0.090246 |
| *(paper NPA-3)* | 0.261 | 0.242 | 0.222 | 0.167 | 0.090 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB,2 | 1+AB,2 |
| full statistics | [1.707286, 1.707288] | 0.902843 | 0.623863 | 0.303884 | 0.113527 |
| *(paper NPA-3)* | 1.707 | 0.903 | 0.624 | 0.304 | 0.114 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |

### WBC φ = π/12 (ancillary): H_min(A₀|X=0,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| I_φ (centre) | 0.900733 | 0.400217 | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.901 | 0.400 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| half-way to W₁ | 0.848375 | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.848 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB,2 | 1+AB | 1+AB | 1+AB | 1+AB |
| vertex W₁ | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| CHSH | 0.214065 | 0.199042 | 0.183093 | 0.139324 | 0.077203 |
| *(paper NPA-3)* | 0.214 | 0.199 | 0.183 | 0.139 | 0.077 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| full statistics | [0.902934, 0.902935] | 0.598917 | 0.450536 | 0.237619 | 0.097246 |
| *(paper NPA-3)* | 0.903 | 0.599 | 0.451 | 0.238 | 0.097 |
| *(NPA level used)* | 1+AB | 2 | 3 | 3 | 3 |

### Tsirelson point P_T: H_min(A₀|X=0,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| CHSH (centre) | 0.911547 | 0.737568 | 0.642251 | 0.472341 | 0.307179 |
| *(paper NPA-3)* | 0.912 | 0.738 | 0.642 | 0.472 | 0.307 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| half-way to g₀ | 0.908267 | 0.726020 | 0.624817 | 0.437441 | 0.194229 |
| *(paper NPA-3)* | 0.908 | 0.726 | 0.625 | 0.437 | 0.194 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| vertex g₀ | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| full statistics | [0.911632, 0.911633] | 0.739922 | 0.648357 | 0.492569 | 0.336009 |
| *(paper NPA-3)* | 0.912 | 0.740 | 0.648 | 0.493 | 0.336 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 2 |

### Tsirelson point P_T: H_min(A₀B₀|00,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| CHSH (centre) | 1.122306 | 0.911258 | 0.794408 | 0.584108 | 0.377572 |
| *(paper NPA-3)* | 1.122 | 0.911 | 0.794 | 0.584 | 0.378 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| half-way to g₀ | 1.118612 | 0.899307 | 0.777202 | 0.544032 | 0.236652 |
| *(paper NPA-3)* | 1.119 | 0.899 | 0.777 | 0.544 | 0.237 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| vertex g₀ | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| full statistics | [1.133620, 1.133621] | 0.954577 | 0.857586 | 0.646259 | 0.422413 |
| *(paper NPA-3)* | 1.134 | 0.955 | 0.858 | 0.646 | 0.422 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |

### Fig. 2 dense curves (certified points)

| v | (a) tilted | (a) CHSH | (a) full | (b) I_φ centre | (b) CHSH | (b) full |
|---|---|---|---|---|---|---|
| 0.999 | 0.776634 | 0.226659 | [0.835087, 0.835088] | 1.729981 | 0.439524 | [1.759378, 1.759379] |
| 0.995 | 0.540022 | 0.219598 | 0.640180 | 1.388400 | 0.425862 | [1.448192, 1.448193] |
| 0.99 | 0.383682 | 0.210974 | 0.514289 | 1.069929 | 0.409375 | [1.163943, 1.163944] |
| 0.985 | 0.273148 | [0.202561, 0.202562] | 0.433678 | 0.809204 | 0.393498 | [0.977841, 0.977842] |
| 0.98 | 0.185275 | 0.194352 | 0.376273 | 0.588463 | 0.378186 | 0.846721 |
| 0.97 | 0.047068 | 0.178506 | 0.298916 | 0.228099 | 0.349099 | 0.671741 |
| 0.96 | 0 (exact) | 0.163369 | 0.248581 | 0 (exact) | 0.321849 | 0.557727 |
| 0.95 | 0 (exact) | 0.148883 | 0.212696 | 0 (exact) | 0.296221 | 0.475418 |
| 0.94 | 0 (exact) | 0.134999 | 0.185295 | 0 (exact) | 0.272038 | [0.411815, 0.411816] |
| 0.92 | 0 (exact) | 0.108865 | 0.140783 | 0 (exact) | 0.227455 | 0.317337 |
| 0.9 | 0 (exact) | 0.084669 | 0.103973 | 0 (exact) | 0.187197 | 0.248684 |
| 0.88 | 0 (exact) | 0.062172 | 0.072777 | 0 (exact) | 0.150614 | 0.195376 |
| 0.86 | 0 (exact) | 0.041183 | 0.045999 | 0 (exact) | 0.117230 | [0.149035, 0.149036] |
| 0.84 | 0 (exact) | 0.021540 | 0.022925 | 0 (exact) | 0.086696 | 0.107753 |

### von Neumann H(A₁|X=1,E): BFF certified lower bound (../rates_v11) and explicit-attack upper bound (this work)

| θ | v | tilted [BFF, attack] | sup over hexagon [BFF best found, attack for all g ∈ C_θ] | full [BFF, attack] |
|---|---|---|---|---|
| pi/8 | 0.999 | [0.9643, 0.9707] | [0.9643, 0.9707] | [0.9653, 0.9725] |
| pi/8 | 0.99 | [0.6988, 0.7065] | [0.6985, 0.7065] | [0.7619, 0.7729] |
| pi/8 | 0.98 | [0.4075, 0.4130] | [0.4075, 0.4130] | [0.6155, 0.6295] |
| pi/8 | 0.95 | [0.0000, 0.0000] | [–, 0.0000] | [0.3784, 0.3960] |
| 27.76deg | 0.999 | [0.9777, 0.9855] | [0.9775, 0.9855] | [0.9765, 0.9858] |
| 27.76deg | 0.99 | [0.8385, 0.8550] | [0.8385, 0.8550] | [0.8551, 0.8729] |
| 27.76deg | 0.98 | [0.6923, 0.7100] | [0.6922, 0.7100] | [0.7502, 0.7729] |
| 27.76deg | 0.95 | [0.2664, 0.2750] | [0.2660, 0.2750] | [0.5384, 0.5611] |
