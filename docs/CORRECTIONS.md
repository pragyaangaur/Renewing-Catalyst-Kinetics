# Corrections

An earlier and much more confident version of this work reached conclusions that did not hold up. They are listed here rather than deleted, because the way they failed is useful to anyone building on the model.

## Claims that were withdrawn

**A zero rate observed in a finite run is not a zero stationary rate.** The early work reported that a randomly sulphated surface at 90 percent coverage was completely dead, based on lattice KMC runs that produced no product events. A later check on a small patch showed the same pattern and then produced product once the measurement window was extended from 0.1 seconds to 20 seconds. The exact master equation solution for that patch had a nonzero rate all along. Short runs at extreme rate ratios can produce no events without the rate being zero.

**Ratios are not automatically more reliable than absolute values.** The early work argued that quantities measured at matched coverage would cancel most of the parameter uncertainty. A joint parameter sweep showed that the star to chain rate ratio ranges from 0.783 to 272.2 across random draws and that the ordering itself flips under a different oxygen mechanism. The blanket reliability claim is withdrawn. The results that genuinely do not depend on the rates are the ones in `src/structure_theory.py`, and they are stated as accessibility and balance arguments rather than as numerical comparisons.

**Seeding a lattice at the mean field composition does not demonstrate stationarity.** The early comparison between lattice KMC and mean field initialised the lattice at the mean field steady state and treated small drift as evidence that the composition was stationary. The high sulphation runs showed real drift and would need long time analysis before their numbers could be used quantitatively.

**The reactor turnover requirement was a sizing scenario and not a calibration.** A figure of about 0.89 turnovers per site per second was derived from an assumed inlet concentration, an assumed site density and an assumed specific surface area. Those inputs were not measured for the bench observation the number was compared against. One conversion measurement cannot identify a set of rate barriers, and the headroom argument built on that comparison does not stand.

**The plug flow marching routine is not converged.** It uses clipped explicit updates, and its plateaus at complete conversion need numerical convergence testing and a sulphur balance check before they can support any reactor level conclusion. Nothing in the current results depends on it.

**The thermochemistry was stated far too sharply.** The sulphation phase boundary was quoted at 866 K with an implied accuracy of a few kelvin, using constants and a constant heat capacity approximation that were not adequately sourced. That number should not be used as a firm operating limit. The qualitative point, that bulk copper sulphate is the stable phase over much of the temperature range of interest at realistic SO3 partial pressures, is less sensitive to the constants.

**The kinetic network is not thermodynamically closed.** The Langmuir Hinshelwood oxidation step is irreversible in the model, and the sulphation and desulphation rates are not tied to the thermochemical calculation. A kinetic Monte Carlo treatment does not repair uncertain chemistry, and it will happily produce precise numbers from a network that violates detailed balance.

**Relevant DFT literature on this exact material exists and was not consulted early enough.** The barriers in `src/params.py` remain illustrative assumptions rather than literature values. See [POSITIONING.md](POSITIONING.md).

## What survived

The two results in `src/structure_theory.py` survived, and they were reformulated so that they do not depend on the rate constants at all. That reformulation is the main lesson from the corrections above. When the inputs are uncertain by orders of magnitude, the conclusions worth keeping are the ones that hold for every input.

The exact small patch solutions survived, because they solve the master equation directly instead of sampling it, and because each one is checked for a single closed communicating class, a nonnegative stationary distribution, a small residual and closed sulphur and oxygen balances.

## Corrections that came from outside

**Result B is a known quantity under a different name.** Aditya Bhan, Distinguished McKnight University Professor at the University of Minnesota, replied to a question about this work and said that the product obtained per copper atom should be the inverse of the cumulative site-loss selectivity defined in Foley, Johnson and Bhan, [ACS Catalysis 9 (2019) 7065](https://doi.org/10.1021/acscatal.9b01106), and asked for that to be checked. It was checked in `src/site_loss_check.py` and it holds. The identity is therefore the renewal case of an existing deactivation framework rather than an independent result, and the repository now says so in [POSITIONING.md](POSITIONING.md) and in [RESULTS.md](RESULTS.md). The part that is still doing work is the split of the inverse selectivity into the branching ratio and the mean sulphated fraction, which separates the kinetic limit from the cost of discarding copper that had not yet deactivated.
