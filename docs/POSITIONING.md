# Positioning against the literature

This page exists so that nobody has to guess what is claimed as new. The short answer is that nothing here is claimed as new. A targeted literature search was run, it found prior work covering most of the ground, and it was not exhaustive enough to establish priority for the parts it did not cover.

If you know of a paper that already contains either of the two rate independent results, please open an issue or send a pointer. It will be added here with credit.

## What the search found

| Source | Finding relevant here | Consequence |
| --- | --- | --- |
| [Zhdanov, 2005, Pattern formation in heterogeneous catalytic reactions with promoters and poisons](https://doi.org/10.1039/B504324G) | Monte Carlo models already describe attractive interactions and the growth of poison islands. | Making poison islands emerge from the dynamics would not by itself establish novelty. |
| [Kinetic description of site ensembles on catalytic surfaces, 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7923573/) | Explicit site ensembles capture spatial correlations that mean field closure misses, including large errors in rate. | The discrepancy between lattice KMC and mean field found in the early part of this work belongs to a known class of effects. |
| [Liu et al., 2017, A comparative investigation of SO2 oxidative transfer over CuO with a CeO2 surface](https://doi.org/10.1016/j.apsusc.2017.01.088) | DFT and experiment describe SO2 oxidation on CuO through a Mars van Krevelen cycle, with SO3 desorption rate controlling. | The Langmuir Hinshelwood mechanism used here is a hypothesis. The presence of gas phase oxygen does not rule out a lattice oxygen cycle. |
| [Barua and Padak, 2025, SO2 binding on CuO during chemical looping with oxygen uncoupling](https://doi.org/10.1021/acs.energyfuels.5c03924) | DFT treatment of SO2 binding, oxygen vacancies and flue gas components on CuO. | Relevant atomistic work on this exact material exists. No barriers from it have been imported into this model. |
| [A method for assessing catalyst deactivation, 2019](https://doi.org/10.1021/acscatal.9b01106) | Treats active sites as consumable and defines site loss yields and selectivities, connected to lifetime turnovers. | Product per lost catalyst site is an established idea. The copper layer accounting here is an application of it to a spatial renewal model. |

Searches covered kinetic Monte Carlo treatments of poison islands and site blocking, CuO and SO2 oxidation including lattice oxygen pathways, catalyst renewal and site ensembles, chain against branched patch topology, and catalyst lifetime with total turnovers under competing deactivation. The accessed abstracts and available text support the positioning above. No full paper barrier tables were available or extracted, so the numerical barriers in this repository remain illustrative assumptions.

## Where a contribution might sit

The candidate contribution is the combination rather than any single piece. It is a renewal material budget, exact small patch kinetics solved rather than sampled, a branching efficiency bound, and controls that show which design recommendations depend on the oxygen mechanism.

The efficiency bound itself follows from simple event balances and has clear precedent in the site loss and lifetime turnover literature above. The three site minimum is a statement about state space accessibility under a dual site adsorption mechanism, and accessibility arguments of that general kind are standard in lattice gas kinetics. Neither was found stated in this form, and absence from a targeted search is weak evidence.

A manuscript built on this would need to position both results explicitly against the sources above. Presenting poison islands, ensemble effects or the three site minimum as newly discovered CuO physics would not survive review.

## The honest summary

The work is reproducible, it is internally checked, and it makes falsifiable predictions. Its novelty is unestablished. Those are three different things and this repository tries not to confuse them.
