# Measurements of long-range two-particle correlation over a wide pseudorapidity range in p–Pb collisions at √sNN = 5.02 TeV

[paper](https://arxiv.org/pdf/2308.16590v3)

[hepdata](https://www.hepdata.net/record/ins2693248)

## Kinematic cuts
Centrality: V0A ($\eta \in [-5.1, -2.8]$) (Pb-going side in the simulations)
Particle of interest: $p_T$ in TPC was exptrapolated to 0 GeV.
Measurements used 3x2PC method:
[[ v_n(\eta) = \sqrt{\frac{V_{n, n}(\eta_a, \eta_b) V_{n, n}(\eta_b, \eta_c)}{V_{n, n}(\eta_b, \eta_c)}} ]]

TPC-FMD1,2: ($\eta \in [-0.8, 0.8]$) vs. ($\eta \in [2.9, 3.1]$)
TPC-FMD3: ($\eta \in [-0.8, 0.8]$) vs. ($\eta \in [-3.1, -2.9]$)
FMD1,2-FMD3: ($\eta \in [4.8, 4.6]$) vs. ($\eta \in [-3.1, -2.9]$)

## Notes
ALICE defined the positive pseudorapidity points in the Pb-going direction,
while simulations have Pb setup as a target going in the negative
pseudorapidity direction.
** We need to apply a minus sign for all these cuts for simulations. **

Email exchange with Yuko Sekiguchi (ALICE definition):
- Reference bins for mid-rapidity $(-0.8 < \eta < 0.8)$ are $-3.1 < \eta < -2.9$ and $2.9 < \eta < 3.1$
- Reference bins for forward rapidity $(\eta > 1.7)$ are $-0.4 < \eta < 0$ and $-3.1 < \eta < -2.9$
- Reference bins for backwardrapidity $(\eta < -1.7)$ are $0 < \eta < 0.4$ and $2.9 < \eta < 3.1$.

