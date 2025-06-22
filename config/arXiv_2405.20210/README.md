# The unexpected uses of a bowling pin: anisotropic flow in fixed-target $^{208}$Pb+$^{20}$Ne collisions as a probe of quark-gluon plasma

## eprint: [link](https://arxiv.org/pdf/2405.20210)

## citation

```
@article{Giacalone:2024ixe,
    author = "Giacalone, Giuliano and others",
    title = "{Anisotropic Flow in Fixed-Target Pb208+Ne20 Collisions as a Probe of Quark-Gluon Plasma}",
    eprint = "2405.20210",
    archivePrefix = "arXiv",
    primaryClass = "nucl-th",
    reportNumber = "CERN-TH-2024-074",
    doi = "10.1103/PhysRevLett.134.082301",
    journal = "Phys. Rev. Lett.",
    volume = "134",
    number = "8",
    pages = "082301",
    year = "2025"
}
```

## Setting the parameters in the iEBE-MUSIC framework

One can use the iEBE-MUSIC framework at the tag `arXiv2405.20210` to generate
event simulations.

One can directly use the provided parameter file
`parameters_PbNe_PGCM_cluster.py` to generate event-by-event simulations for
Pb+Ne collisions at 68.5 GeV.

For Pb+O collisions, one can set `Projectile` to "O".

By changing the value of the parameter `lightNucleusOption`, we can choose
nucleus configurations from different theories,

- lightNucleusOption = 2  ->  PGCM clustered
- lightNucleusOption = 3  ->  PGCM uniform
- lightNucleusOption = 4  ->  NLEFT positive weights
- lightNucleusOption = 5  ->  NLEFT negative weights
