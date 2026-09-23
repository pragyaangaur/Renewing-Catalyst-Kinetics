# Results

This work originated with SkyCatcher, a project on sulphur recovery from industrial flue gas supported by Mercedes-Benz and Emergent Ventures. What follows is a model study and contains no device design or experimental data from it.

Everything below is a property of the model that is defined in this repository. None of it has been validated against experiment. Read [LIMITATIONS.md](LIMITATIONS.md) before quoting any number.

## The model

Each surface site carries one substrate phase and at most one adsorbate. The phase is metallic copper, copper oxide or copper sulphate. The adsorbate is nothing, atomic oxygen, SO2 or SO3.

The elementary steps are the following. Oxygen arrives as O2 and dissociates onto two neighbouring vacant oxide sites, so it needs a pair. SO2 adsorbs on a single vacant oxide site and can desorb again. An adsorbed SO2 reacts with an adsorbed oxygen atom on a neighbouring site to give an adsorbed SO3, which is the productive step. That adsorbed SO3 then either desorbs as product or reacts into the lattice as sulphate, which is the branch point that decides everything. A sulphated site can decompose back to oxide, and freshly exposed metallic copper reoxidises.

The branch between product desorption and sulphation is the whole model. Both channels draw on the same adsorbed SO3 population, so their rates stand in a fixed ratio set by the two rate constants.

![The two rate independent results](../figures/structure_results.png)

## Result A: three sites are the minimum for a nonzero rate

On a connected patch of active oxide sites, with all other surface sites permanently inactive, the stationary rate of product formation is exactly zero for one site and for two sites, and strictly positive for three or more.

The reason is an accessibility argument rather than a rate argument. A single site has no neighbour, so the bimolecular productive step can never fire. On a pair, the only source of adsorbed oxygen is O2, which lands on two empty neighbours and fills both of them. After it lands there is no vacant site left for SO2, and removing an oxygen atom requires the reverse step, which takes both atoms away together. Starting from an empty patch, no reachable configuration holds an oxygen atom and an SO2 at the same time. On three connected sites there is always a path a to b to c, so oxygen can land on the pair a and b, SO2 can then adsorb on the still vacant c, and the reaction fires across b and c.

Every connected graph on three or more vertices contains such a path, because a connected graph in which no vertex has two or more neighbours is at most a single edge. The criterion is therefore the number of sites in the connected component and not its shape.

The positive rate for three or more sites also needs the productive configuration to recur. It does, because the empty patch can be reached from every configuration. SO2 and SO3 can always leave, two neighbouring O atoms can leave together, and an O atom with no O neighbour can be removed by emptying a neighbour, filling it with SO2 and letting the reaction fire. The reachable set is then one communicating class, so every configuration in it has positive stationary probability.

This is checked in two independent ways in `src/structure_theory.py`. The reachability test uses integer logic only and treats every rate as merely possible, so its answer depends on the patch shape alone. It checks both that a productive configuration is reachable and that the reachable set is a single communicating class, with and without adsorbate hops, on all 996 connected graphs with one to seven sites, taken from the Atlas of Graphs. The largest state space has 16,376 configurations. The numerical test solves the stationary master equation for randomised rate constants spanning nine orders of magnitude, on every connected graph up to five sites. All 16 dead cases gave exactly zero, and all 232 live cases gave a strictly positive rate, from 2.3e-21 to 485 per site per second, with relative residuals below 4.1e-17.

One practical caution follows from this. When rate constants are extreme, the stationary rate on a live patch can be tiny, and in a finite stochastic run it can produce no events at all. A general purpose linear solver can also return a tiny rate with a large relative error. None of those is evidence that the true rate is zero. The stationary solver is now the GTH reduction, which is accurate entry by entry, and it matched 60 digit arithmetic to 1.3e-15 in the test suite. Earlier versions of this work made both mistakes, and they are recorded in [CORRECTIONS.md](CORRECTIONS.md).

## Result B: an exact ceiling on product per copper atom consumed

Consider a patch that is periodically discarded and replaced by fresh metallic copper, which then reoxidises. Each reset spends the copper in that layer, including any part of it that had not finished oxidising.

Write $k_d$ for product desorption from adsorbed SO3 and $k_s$ for sulphation from the same intermediate. Both event rates are proportional to the population of that one intermediate, so

$$\frac{R_{\mathrm{product,direct}}}{R_{\mathrm{sulphation}}} \;=\; \frac{k_d}{k_s}$$

At stationary state the sulphate created must equal the sulphate destroyed. With resets at frequency $\gamma$ on a patch of $n$ sites, and a chemical regeneration channel that returns sulphate to oxide without spending copper,

$$R_{\mathrm{sulphation}} \;=\; R_{\mathrm{regeneration}} \;+\; \gamma\, n\, f_{\mathrm{sulphate}}$$

Dividing the total product by the copper spend $\gamma n$, and writing $\rho$ for the regenerations per copper atom exposed, gives

$$Y_{\mathrm{Cu}} \;=\; \left(1 + \frac{k_d}{k_s}\right)\rho \;+\; \frac{k_d}{k_s}\, f_{\mathrm{sulphate}}$$

With no regeneration channel this reduces to $Y_{\mathrm{Cu}} = (k_d/k_s) f_{\mathrm{sulphate}}$, which cannot exceed $k_d/k_s$ because the sulphated fraction cannot exceed one. That ceiling is independent of patch shape, diffusion, reoxidation speed and reset frequency. No arrangement of the surface and no renewal schedule can overcome a poor branching ratio.

With regeneration switched on, the first term has no ceiling and grows as the reset frequency falls. Chemical regeneration is therefore the only route past the renewal limit in this model.

The identity was checked across 200 randomised cases covering four patch shapes, with sulphation and desorption rates drawn over five orders of magnitude and regeneration switched off in a third of them. The maximum relative error was 1.4e-15, which is machine precision. An earlier version reported 5.6e-6, which was the error of a less accurate solver. In the cases without regeneration the largest observed value of $Y_{\mathrm{Cu}}$ divided by $k_d/k_s$ was 0.665, which respects the ceiling. With regeneration the same quantity reached 5.7e6.

Numerical agreement here is a conservation check and a demonstration of the derivation. It is not independent evidence of new physics, because the identity follows from the event balances that the solver is already enforcing.

### Why this matters for a renewing catalyst

Throughput and material efficiency pull in opposite directions. Resetting the surface more often keeps the sulphated fraction low and raises the rate per site, and it also spends copper faster. The table below uses one illustrative barrier set on a square patch with fast reoxidation.

| Reset frequency (1/s) | Product per site per second | Product per copper atom exposed |
| ---: | ---: | ---: |
| 0.01 | 0.443 | 44.28 |
| 0.1 | 2.023 | 20.23 |
| 1 | 5.668 | 5.668 |
| 100 | 262.3 | 2.623 |
| 3162 | 1910 | 0.604 |

The fast end of that table is a mathematical stress test rather than a proposed operating point. At slow reoxidation the picture inverts, because resetting too often leaves mostly bare metal, which is inactive in this scheme. There is therefore an optimum reset rate set by the ratio of reset frequency to reoxidation rate, and a renewal strategy that maximises shedding is not automatically good.

### The relation to cumulative site-loss selectivity

Aditya Bhan, Distinguished McKnight University Professor at the University of Minnesota, pointed out on reading the description of this work that the product obtained per copper atom should be the inverse of the cumulative site-loss selectivity defined by Foley, Johnson and Bhan, [ACS Catalysis 9 (2019) 7065](https://doi.org/10.1021/acscatal.9b01106), and asked for the algebra to be checked. That paper treats active sites as a consumable reactant and defines the cumulative site-loss selectivity as total moles of sites lost divided by total moles of reactant converted to effluent product, so its inverse is the product obtained per site lost. The check was run in `src/site_loss_check.py` over the same 200 randomised cases used above, and the suggestion holds.

The comparison only becomes definite once "a site lost" is tied to an event in this model, and there are two candidates. A reset discards the whole copper layer, including copper that had never sulphated, so the number of copper atoms discarded and the number of sulphation events are different quantities.

$$S_{\mathrm{Cu}} \;=\; \frac{\text{copper atoms discarded}}{\text{product molecules}} \qquad\qquad S_{\mathrm{sulf}} \;=\; \frac{\text{sulphation events}}{\text{product molecules}}$$

Reading the discarded copper atom as the lost site gives the correspondence exactly, to a maximum relative error of 2.2e-16 across the ensemble.

$$Y_{\mathrm{Cu}} \;=\; \frac{1}{S_{\mathrm{Cu}}}$$

That identity is definitional once the mapping is fixed, so it confirms the accounting rather than proving anything new. The second reading is the more informative one. Taking the sulphation event as the lost site, the identity of Result B rearranges to

$$Y_{\mathrm{Cu}} \;=\; \frac{\rho + f_{\mathrm{sulphate}}}{S_{\mathrm{sulf}}}$$

which holds to a maximum relative error of 1.4e-15, the same precision as the identity itself. With renewal only, $\rho$ is zero and the inverse sulphation site-loss selectivity is the ceiling of Result B exactly, $1/S_{\mathrm{sulf}} = k_d/k_s$, verified to 4.2e-16. The gap between the copper efficiency and its ceiling is then exactly $f_{\mathrm{sulphate}}$, which ran from 6.1e-4 to 0.665 over the ensemble.

This gives the ceiling a cleaner reading than the one first written down here. The branching ratio $k_d/k_s$ is the inverse site-loss selectivity of the sulphation channel, so it is the product obtained per site actually deactivated. The mean sulphated fraction is the share of the discarded copper that had in fact deactivated. The distance from the ceiling is therefore a materials accounting cost. It measures the copper that was thrown away while it was still in service.

The 2019 paper also gives a variant of the denominator that counts all reactant consumed rather than only reactant converted to effluent product. Here the extra term is the sulphur that leaves with the discarded layer, and the check confirms

$$\frac{1}{S_{\mathrm{consumed}}} \;=\; Y_{\mathrm{Cu}} \;+\; f_{\mathrm{sulphate}}$$

to a maximum relative error of 2.9e-16. Fifty of the 200 drawn cases are two-site patches, where no product forms at all by Result A and every one of these ratios is undefined. Those cases are reported separately and excluded from the error figures above.

## Result C: patch shape matters, and its ranking is not robust

At equal site count and equal bond count, different patch shapes give different rates. A four site star and a four site chain both have four sites and three bonds, and under one illustrative barrier set their rates differ by a factor of 2.94.

That ranking does not survive a change of mechanism or of rates. A deliberately simple lattice oxygen control, in which SO2 consumes a lattice oxygen atom directly and O2 refills adjacent vacancy pairs, inverts the ordering, giving a star to chain ratio between 0.712 and 0.994 across refill to consumption ratios from 0.01 to 100. Under the original mechanism, a joint sweep of 100 random parameter draws gave star to chain ratios from 0.783 to 272.2, with the star ahead in only 88 of the 100 draws. Diffusion changes the advantage as well.

The conclusion is that shape ranking on its own does not identify the oxygen mechanism. The inversion is a useful hypothesis for designing an experiment and it is not a unique diagnostic. The stricter distinction from Result A, that a pair is dead and three sites are not, does survive the rate and diffusion changes that were tested, although it could still change if additional elementary paths exist.

## Result D: the population of usable patches

At 90 percent random blocking of a square lattice, the fraction of surviving active sites sitting in isolated single sites is $(1-q)^4$ and the fraction in isolated pairs is $4q(1-q)^6$, with $q$ the surviving fraction of 0.1. Together those account for 86.868 percent of the surviving sites, which by Result A can produce nothing. The remaining 13.132 percent sit in components of at least three sites.

Direct simulation on a 128 by 128 lattice gives 13.213 percent with a standard error of 0.215. Moving to a six neighbour lattice raises it to 21.253 percent with a standard error of 0.262, so coordination number changes the usable population substantially.

A finite productive patch does not need a surface spanning cluster, so describing this as a percolation threshold would be misleading. What matters is the size distribution of small components and not the existence of an infinite one.

### The general law

The same count works on any lattice. With coordination number $z$ and a fraction $p$ of sites active, the share of surviving sites that can never turn over is

$$D(p) \;=\; (1-p)^z \;+\; z\,p\,(1-p)^{t_2}$$

where $t_2 = 2z - 2 - c$ is the perimeter of a nearest neighbour pair and $c$ is the number of neighbours the two sites share. That gives $t_2$ = 4 on the honeycomb lattice, 6 on the square lattice and 8 on the triangular lattice, and at 90 percent deactivation $D$ = 0.9258, 0.8687 and 0.7897. The six neighbour sampling above, 21.253 percent usable, agrees with the exact 21.03 percent. Under the lattice oxygen control only isolated sites are dead and $D(p) = (1-p)^z$. The law holds for every positive rate assignment because it combines Result A with pure combinatorics.

`src/cluster_expansion.py` checks it three ways. The lattice animal counts to eight sites match OEIS A001420, A001168 and A001207. The cluster numbers reproduce $p$ exactly through order $p^8$ on all three lattices. Direct sampling of eight 512 by 512 lattices at each of 49 values of $p$ agrees with the closed form to within 3.0, 2.2 and 2.0 standard errors on the honeycomb, square and triangular lattices, and a larger rerun at the worst honeycomb points gave below 1.6.

## Result E: the exact rate of a randomly deactivated surface

Deactivated sites take part in no process, so each component of the surviving sites evolves independently and the lattice rate is a sum over component shapes,

$$\mathrm{TOF}(p) \;=\; \sum_A w_A\, p^{|A|} (1-p)^{t(A)}\, R(A)$$

with $R(A)$ the exact stationary rate of an isolated component shaped like $A$, $t(A)$ its perimeter and $w_A$ one over the sites per unit cell. This is the classical cluster number expansion of site percolation with exact master equation rates as the weights. On the square lattice the 304 animals of three to six sites fall into 18 graph isomorphism classes, each solved once. Every term is non-negative, so a truncated sum is a rigorous lower bound, and the share of active sites it leaves out is known exactly.

With the illustrative 745 K rates of the first pass geometry study and 90 percent random deactivation, the bound is 0.0283 per site per second. The omitted components hold 0.23 percent of the active sites. The first pass reported zero for a surface of exactly this kind. `src/firstpass_postmortem.py` rebuilds that lattice from its seed. Every component in it has five sites or fewer, so it can be solved exactly, and its stationary rate is 0.0300 per site per second. The first pass measured a window of 0.026 seconds that began 0.0127 seconds after an empty start. The transient master equation gives 3.57 expected product events in that window, so a run that saw none had a probability of 0.028. The same engine run for four seconds of simulated time on the same lattice gives 0.0229 plus or minus 0.0048, which is 1.5 standard errors from the exact value.

At 85 percent sulphation the first pass reported 0.157, while the exact stationary rate of that lattice is 0.0751 and a long run gives 0.083 plus or minus 0.016. The first pass value was 2.1 times too high because its window still lay inside the transient that follows the empty start. The exact expectation for that window is 10.9 events, against 5.8 at steady state, and the first pass observed 12. Both first pass numbers were faithful samples of what the engine was simulating. Neither was a measurement of the stationary rate. The numbers are in `results/firstpass_postmortem.json`.

## Result F: the lattice KMC engine against exact answers on whole lattices

The factorisation also gives an exact target for the lattice engine on a full 64 by 64 lattice. `src/lattice_validation.py` builds 24 frozen random lattices at $p$ = 0.1, 0.2 and 0.3, deactivates components larger than six sites so that every component is exactly solvable, and compares the engine with the exact component sum at synthetic rates of order one. All 24 runs passed the acceptance test fixed in advance, five standard errors plus one percent. The pooled chi-square was 21.7 on 24 degrees of freedom, p = 0.60, the largest |z| was 2.14 and the largest relative deviation 0.63 percent.
