# Pseudorapidity dependence of anisotropic flow and its decorrelations using long-range multiparticle correlations in Pb–Pb and Xe–Xe collisions

[paper](https://arxiv.org/pdf/2307.11116v2)

[hepdata](https://www.hepdata.net/record/ins2679248)

## Kinematic cuts
Centrality: V0A ($\eta \in [2.8, 5.1]$) + V0C ($\eta \in [-3.7, -1.7]$)
Particle of interest: $p_T$ in TPC was exptrapolated to 0 GeV with AMPT model.

## Notes
- $v_n\{2, |\Delta \eta| > 0.8\} mid:
  TPC $v_n(\eta < 0)$ correlated with the TPC $v_n(\eta \in [0.8, 1])$
  TPC $v_n(\eta > 0)$ correlated with the TPC $v_n(\eta \in [-1, -0.8])$

- $v_n\{2, |\Delta \eta| > 2.6\} forward:
  FMD $v_n(\eta < 0)$ correlated with the TPC $v_n(\eta \in [0.8, 1])$
  FMD $v_n(\eta > 0)$ correlated with the TPC $v_n(\eta \in [-1, -0.8])$

- $v_n\{2, |\Delta \eta| > 2\} mid:
  TPC $v_n(\eta < 0)$ correlated with the FMD $v_n(\eta \in [2, 5])$
  TPC $v_n(\eta > 0)$ correlated with the FMD $v_n(\eta \in [-3.5, -2])$

- $v_n\{2, |\Delta \eta| > 2\} forward:
  FMD $v_n(\eta > 0)$ correlated with the TPC $v_n(\eta \in [-1, -0.2])$
  FMD $v_n(\eta < 0)$ correlated with the TPC $v_n(\eta \in [0.2, 1])$

- $v_n\{2, |\Delta \eta| > 3.8\} forward:
  FMD $v_n(\eta < 0)$ correlated with the FMD $v_n(\eta \in [2, 5])$
  FMD $v_n(\eta > 0)$ correlated with the FMD $v_n(\eta \in [-3.5, -2])$

- r_{n, n}:
  numerator:    <cos[n(phi_C' - phi_A)] + cos[n(phi_B' - phi_D)]>
  denorminator: <cos[n(phi_C' - phi_D)] + cos[n(phi_B' - phi_A)]>
  FMD: A: (-3.5 < \eta < -1.8), D: (1.8 < \eta < 5)  eta-integrated vn
  TPC: B: (-1 < \eta < 0), C: (0 < \eta < 1)         eta-differential vn
