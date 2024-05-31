# Unveiling baryon charge carriers through charge stopping in isobar collisions

## eprint: [link](https://arxiv.org/pdf/2405.19439)

## citation

```
@article{Pihan:2024lxw,
    author = {Pihan, Gregoire and Monnai, Akihiko and Schenke, Bj\"orn and Shen, Chun},
    title = "{Unveiling baryon charge carriers through charge stopping in isobar collisions}",
    eprint = "2405.19439",
    archivePrefix = "arXiv",
    primaryClass = "nucl-th",
    month = "5",
    year = "2024"
}
```

## Setting the parameters in the iEBE-MUSIC framework

One can use the iEBE-MUSIC framework at the tag `arXiv2405.19439` to generate
event simulations.

One can directly use the provided parameter file
`parameters_dictRu_LQ0.py` to generate event-by-event simulations for
Ru+Ru collisions at 200 GeV. 

`parameters_dictZr_LQ0.py` to generate event-by-event simulations for
Zr+Zr collisions at 200 GeV. 

By changing the value of the parameter `lambdaQ`, we can adjust the amount
of electric charge stopping at the initial state.

By setting values for `ProjWS_da` and `ProjWS_dR`, we can change the neutron
skin size for the colliding nuclei.

```math
    da = a_n - a_p, dR = R_n - R_p
```
