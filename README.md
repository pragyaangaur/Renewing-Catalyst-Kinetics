# Renewing catalyst kinetics

Exact master equation and kinetic Monte Carlo models of SO2 oxidation on a self renewing copper oxide surface.

The motivating question is what happens to a copper oxide catalyst that sulphates in service and recovers by exposing fresh copper underneath. Copper oxide is normally used as a sulphur sorbent rather than as a catalyst, because sulphation converts it to copper sulphate, which is inactive. A surface that renews itself changes the accounting, because sites are no longer a fixed inventory that can only be lost. This repository asks what such a surface can and cannot do, and it tries to separate the answers that depend on assumed rate constants from the answers that do not.

## Where this came from

This work grew out of SkyCatcher, a project by the author on recovering sulphur from industrial flue gas, supported by Mercedes-Benz and by Emergent Ventures. SkyCatcher raised the question that this repository studies, which is whether a copper oxide surface that sulphates in service can still function as a catalyst rather than as a sorbent once it renews itself.

The repository is a separate piece of work. It contains no device design, no experimental data and nothing specific to that project. Everything here is a model, and the results stand or fall on the model alone.

## What is here

The two results worth reading are stated without any dependence on the numerical rate constants. They hold for every strictly positive assignment of rates, so the fact that the barriers in this repository are illustrative does not weaken them.

**A minimum patch size of three sites.** In the adsorbate mechanism modelled here, oxygen arrives as O2 and occupies two neighbouring vacant sites at once, while SO2 adsorbs on a single vacant site, and the two must meet on neighbouring sites to react. On an isolated pair of active sites the oxygen step fills both sites, so SO2 can never adsorb alongside it, and the stationary product rate is exactly zero. Three connected sites are enough to make the rate positive. Every connected graph on three or more sites contains the required path, so the criterion is the patch size and not the shape. This is proved in `src/structure_theory.py` and checked by exhaustive enumeration of all 996 connected graphs with up to seven sites, with and without adsorbate hops.

**An exact bound on product per copper atom consumed.** When sulphated sites are recovered only by discarding the oxide layer and exposing fresh copper, the product obtained per copper atom exposed is

$$Y_{\mathrm{Cu}} \;=\; \frac{k_d}{k_s}\, f_{\mathrm{sulphate}} \;\leq\; \frac{k_d}{k_s}$$

where $k_d$ is product release and $k_s$ is sulphation from the same adsorbed intermediate, and $f_{\mathrm{sulphate}}$ is the mean sulphated fraction. No patch shape, reset frequency or reoxidation speed can beat that ceiling. Adding a chemical regeneration channel, in which a sulphated site decomposes and returns to service without spending another copper atom, generalises this to

$$Y_{\mathrm{Cu}} \;=\; \left(1 + \frac{k_d}{k_s}\right)\rho \;+\; \frac{k_d}{k_s}\, f_{\mathrm{sulphate}}$$

where $\rho$ is the number of regeneration events per copper atom exposed. The first term has no ceiling, which identifies chemical regeneration as the only route past the renewal limit. Both forms are verified numerically to a relative error of 1.4e-15 across 200 randomised cases.

On the suggestion of Aditya Bhan, Distinguished McKnight University Professor at the University of Minnesota, this ceiling was checked against the cumulative site-loss selectivity of Foley, Johnson and Bhan, [ACS Catalysis 9 (2019) 7065](https://doi.org/10.1021/acscatal.9b01106). It is the same quantity. Counting a discarded copper atom as a lost site makes $Y_{\mathrm{Cu}}$ the inverse cumulative site-loss selectivity exactly, and counting a sulphation event as a lost site makes the ceiling $k_d/k_s$ the inverse site-loss selectivity of the sulphation channel, with $f_{\mathrm{sulphate}}$ as the ratio between the two. The mean sulphated fraction is then the share of discarded copper that had actually deactivated, so the gap below the ceiling is the cost of throwing away copper that was still in service. The check is in `src/site_loss_check.py` and the algebra is written out in [docs/RESULTS.md](docs/RESULTS.md).

**How much of a poisoned surface is dead.** Combining the first result with the cluster statistics of site percolation gives, for a lattice of coordination number $z$ with a fraction $p$ of sites still active, the share of surviving sites that can never turn over,

$$D(p) \;=\; (1-p)^z \;+\; z\,p\,(1-p)^{t_2}$$

where $t_2$ is the perimeter of a neighbouring pair, 4, 6 and 8 on the honeycomb, square and triangular lattices. At 90 percent random deactivation of a square lattice, 86.9 percent of the surviving sites are dead for every rate assignment. The same expansion, with exact rates for each small component shape, gives the rate of a randomly deactivated surface without any sampling, and it shows that a first pass Monte Carlo run which reported a dead surface was wrong. This is in `src/cluster_expansion.py` and `src/firstpass_postmortem.py`.

A practical consequence of the copper bound is that a surface can be made to produce faster while using its copper less efficiently. Throughput and material efficiency are not the same objective, and the renewal frequency trades one against the other.

![The two rate independent results](figures/structure_results.png)

## What is not here

No experimental validation. No claim that these results are new to the catalysis literature. Poison islands, site ensemble effects and the failure of mean field kinetics to capture spatial correlation are all established subjects, and the site loss and lifetime turnover literature already treats active sites as consumable. A targeted literature search was run and is summarised in [docs/POSITIONING.md](docs/POSITIONING.md). It did not establish priority for anything here, and it was not exhaustive.

The reaction barriers in `src/params.py` are illustrative values chosen to be physically reasonable for oxide supported SO2 oxidation. They are not fitted to data and they are not taken from any published table for CuO. Any result quoted as an absolute rate should be read as a property of the model rather than as a prediction about copper oxide. This is the main reason the two headline results above were developed in a rate independent form.

An earlier and more confident version of this work overstated several conclusions. The corrections are recorded in [docs/CORRECTIONS.md](docs/CORRECTIONS.md) rather than quietly removed, because the way the claims failed is useful information for anyone building on this.

## Layout

| Path | What it is |
| --- | --- |
| `run_all.py` | Regenerates every result and figure and writes `results/MANIFEST.json` with checksums. |
| `src/structure_theory.py` | The two rate independent results, with proofs in the docstrings and exhaustive checks. |
| `src/cluster_expansion.py` | Lattice animals, the dead fraction law on three lattices and the exact rate expansion. |
| `src/lattice_validation.py` | The lattice KMC engine against exact component sums on 24 whole lattices. |
| `src/firstpass_postmortem.py` | The first pass zero rate, rebuilt from its seed and explained exactly. |
| `src/paper_figures.py`, `src/paper_numbers.py` | Figures and quoted numbers for the paper, from saved results only. |
| `tests/` | The test suite. |
| `src/structure_figures.py` | Draws the figure above from the saved result files. |
| `src/site_loss_check.py` | Maps Result B onto the cumulative site-loss selectivity of the deactivation literature. |
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

Python 3.10 or newer with numpy, scipy, matplotlib, networkx and mpmath, plus pytest for the tests. The archived results were produced with Python 3.13, numpy 2.5, scipy 1.18, matplotlib 3.11 and networkx 3.6.

```bash
pip install -r requirements.txt
```

Everything, in dependency order, with a checksum manifest at the end. This takes about an hour on a laptop, most of it in the long KMC runs.

```bash
python3 run_all.py
```

The test suite takes under a minute. It checks the solver against 60 digit arithmetic, the enumerations against OEIS and every number quoted in the paper against the saved results.

```bash
python3 -m pytest tests
```

The two rate independent results take about eight minutes together, most of it in the seven site reachability check.

```bash
python3 src/structure_theory.py
```

The comparison with cumulative site-loss selectivity takes about ten seconds.

```bash
python3 src/site_loss_check.py
```

The exact patch solves and their stochastic validation take a few minutes.

```bash
python3 src/ensemble_audit.py --checks && python3 src/ensemble_audit.py --renewal
```

The lattice KMC scripts are pure Python and take roughly six to eight minutes each. Installing numba is the first thing to do if you want to extend them.

```bash
python3 src/kmc_geometry.py
```

## Citing

The archived release is on Zenodo at https://doi.org/10.5281/zenodo.22910821. The paper that describes the results in full is included in that archive and is also published as a preprint on ChemRxiv. `CITATION.cff` carries the same metadata in machine readable form.

## Contributing and contact

Corrections are more welcome than agreement. If one of the two rate independent results is already in the literature, a pointer to the paper is the most useful thing you can send, and it will be added to `docs/POSITIONING.md` with credit. If the model itself is wrong in a way that matters, an issue describing the missing elementary step or the broken assumption is the right place to start.

## Licence

MIT for the code, CC BY 4.0 for the documentation and figures. See `LICENSE`.
