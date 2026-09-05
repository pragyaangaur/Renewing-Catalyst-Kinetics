# Renewal kinetics

Exact master equation and kinetic Monte Carlo models of SO2 oxidation on a self renewing copper oxide surface.

The motivating question is what happens to a copper oxide catalyst that sulphates in service and recovers by exposing fresh copper underneath. Copper oxide is normally used as a sulphur sorbent rather than as a catalyst, because sulphation converts it to copper sulphate, which is inactive. A surface that renews itself changes the accounting, because sites are no longer a fixed inventory that can only be lost. This repository asks what such a surface can and cannot do, and it tries to separate the answers that depend on assumed rate constants from the answers that do not.

## What is here

The two results worth reading are stated without any dependence on the numerical rate constants. They hold for every strictly positive assignment of rates, so the fact that the barriers in this repository are illustrative does not weaken them.

**A minimum patch size of three sites.** In the adsorbate mechanism modelled here, oxygen arrives as O2 and occupies two neighbouring vacant sites at once, while SO2 adsorbs on a single vacant site, and the two must meet on neighbouring sites to react. On an isolated pair of active sites the oxygen step fills both sites, so SO2 can never adsorb alongside it, and the stationary product rate is exactly zero. Three connected sites are enough to make the rate positive. Every connected graph on three or more sites contains the required path, so the criterion is the patch size and not the shape. This is proved in `src/structure_theory.py` and checked by exhaustive enumeration of every connected graph up to five sites, with six available as an option.

**An exact bound on product per copper atom consumed.** When sulphated sites are recovered only by discarding the oxide layer and exposing fresh copper, the product obtained per copper atom exposed is

```
Y_Cu = (k_d / k_s) * f_sulphate        which is at most  k_d / k_s
```

where `k_d` is product release and `k_s` is sulphation from the same adsorbed intermediate, and `f_sulphate` is the mean sulphated fraction. No patch shape, reset frequency or reoxidation speed can beat that ceiling. Adding a chemical regeneration channel, in which a sulphated site decomposes and returns to service without spending another copper atom, generalises this to

```
Y_Cu = (1 + k_d / k_s) * rho + (k_d / k_s) * f_sulphate
```

where `rho` is the number of regeneration events per copper atom exposed. The first term has no ceiling, which identifies chemical regeneration as the only route past the renewal limit. Both forms are verified numerically to a relative error of about 6e-6 across 200 randomised cases.

A practical consequence is that a surface can be made to produce faster while using its copper less efficiently. Throughput and material efficiency are not the same objective, and the renewal frequency trades one against the other.

![The two rate independent results](figures/structure_results.png)

## What is not here

No experimental validation. No claim that these results are new to the catalysis literature. Poison islands, site ensemble effects and the failure of mean field kinetics to capture spatial correlation are all established subjects, and the site loss and lifetime turnover literature already treats active sites as consumable. A targeted literature search was run and is summarised in [docs/POSITIONING.md](docs/POSITIONING.md). It did not establish priority for anything here, and it was not exhaustive.

The reaction barriers in `src/params.py` are illustrative values chosen to be physically reasonable for oxide supported SO2 oxidation. They are not fitted to data and they are not taken from any published table for CuO. Any result quoted as an absolute rate should be read as a property of the model rather than as a prediction about copper oxide. This is the main reason the two headline results above were developed in a rate independent form.

An earlier and more confident version of this work overstated several conclusions. The corrections are recorded in [docs/CORRECTIONS.md](docs/CORRECTIONS.md) rather than quietly removed, because the way the claims failed is useful information for anyone building on this.

## Layout

| Path | What it is |
| --- | --- |
| `src/structure_theory.py` | The two rate independent results, with proofs in the docstrings and exhaustive checks. |
| `src/structure_figures.py` | Draws the figure above from the saved result files. |
| `src/ensemble_audit.py` | Exact master equation solutions on small patches, with independent Gillespie validation. |
| `src/research_controls.py` | Mechanism controls, uncertainty sweeps and random lattice geometries. |
| `src/kmc_core.py` | Rejection free lattice KMC engine, n-fold way, constant time event sampling. |
| `src/meanfield.py` | The same process list as ordinary differential equations, for comparison. |
| `src/params.py` | Rate constant construction, with the provenance of every barrier stated. |
| `src/thermo.py` | Sulphation phase boundary from tabulated thermochemistry. |
| `src/kmc_geometry.py` | Turnover at matched sulphate coverage, random against clustered. |
| `docs/` | Results, methods, positioning against the literature, corrections and limitations. |
| `results/` | Saved arrays and tables, so every figure regenerates without rerunning the sweeps. |
| `figures/` | Generated figures. |

## Running it

Python 3.10 or newer with numpy, scipy and matplotlib.

```bash
pip install -r requirements.txt
```

The two rate independent results take about half a minute together.

```bash
python3 src/structure_theory.py
```

The exact patch solves and their stochastic validation take a few minutes.

```bash
python3 src/ensemble_audit.py --checks && python3 src/ensemble_audit.py --renewal
```

The lattice KMC scripts are pure Python and take roughly six to eight minutes each. Installing numba is the first thing to do if you want to extend them.

```bash
python3 src/kmc_geometry.py
```

## Contributing and contact

Corrections are more welcome than agreement. If one of the two rate independent results is already in the literature, a pointer to the paper is the most useful thing you can send, and it will be added to `docs/POSITIONING.md` with credit. If the model itself is wrong in a way that matters, an issue describing the missing elementary step or the broken assumption is the right place to start.

## Licence

MIT for the code, CC BY 4.0 for the documentation and figures. See `LICENSE`.
