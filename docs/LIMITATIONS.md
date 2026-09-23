# Limitations

## The barriers are illustrative

There is no published density functional study of an in situ grown copper oxide layer under flue gas conditions that this work draws barriers from. The values in `src/params.py` were chosen to be physically reasonable for oxide supported SO2 oxidation, and they are not fitted to data. Absolute rates from this model are properties of the model. This is the reason the two headline results were developed in a form that holds for every positive rate assignment.

## The mechanism is a hypothesis

The model uses a Langmuir Hinshelwood step between adsorbed SO2 and adsorbed atomic oxygen. Published work on SO2 oxidation over CuO describes a Mars van Krevelen cycle using lattice oxygen. A lattice oxygen control is included in `src/research_controls.py` and it changes the minimum productive patch from three sites to two, and it inverts the shape ranking. Which mechanism operates on a real surface is an open experimental question, and the answer changes some of the design conclusions.

## Single layer accounting

The lattice resolves one surface layer. It does not represent the growth of a thick oxide scale, diffusion through that scale, or the consumption of a bulk copper inventory. The renewal model discards one atom per site per reset, which is a material accounting convention rather than a measured oxide thickness or a spallation law. Multilayer removal would change the copper budget and therefore the efficiency numbers.

## The efficiency bound is not universal

The renewal only ceiling assumes that sulphate is removed only by discarding the layer, that the branching ratio is uniform across sites, and that gas phase SO3 readsorption is negligible. Each of those can fail. Chemical regeneration is handled explicitly by the generalised identity, and the others are not. Nonuniform branching, multilayer removal and alternative poisoning routes all change the accounting and must be included before applying any efficiency bound to a real device.

## Finite size and coordination

The exhaustive graph enumeration runs to five sites by default and six as an option. The lattice geometry sampling covers 32, 64 and 128 sites per side on four and six neighbour lattices. Coordination number changes the usable patch population substantially, from 13.2 percent to 21.3 percent at 90 percent blocking, so results should not be carried between lattice types without rechecking.

## Numerical caveats

Rate constants spanning nine orders of magnitude make the stationary balance matrix badly conditioned. The exact solver in `src/structure_theory.py` therefore uses the Grassmann, Taksar and Heyman reduction on the closed communicating class, which never subtracts and is accurate entry by entry. It matched 60 digit arithmetic to about 1e-15 in the tests. The older solver in `src/ensemble_audit.py` uses a sparse direct solve. It is used only with the illustrative rates, where the test suite shows it agrees with the GTH solver to better than 1e-9. Every solve is checked for one closed communicating class, a nonnegative distribution and a small relative residual, and the ensemble audit also checks closed element balances.

Dense GTH costs the cube of the number of states, so the exact solves stop at six sites without renewal, about four thousand states. The rate expansion of `src/cluster_expansion.py` is therefore truncated at components of six sites. The truncated sum is a lower bound, and the share of active sites it omits is reported with it. It is small at low coverage and grows quickly above about 20 percent active sites.

The whole lattice check of the KMC engine uses synthetic rates of order one. The illustrative rates are too stiff for a pure Python engine to sample small rates on a whole lattice in reasonable time, which is how the first pass came to report a zero rate. The long runs in `src/firstpass_postmortem.py` do use the illustrative rates, on two lattices only.

## No experimental validation

Nothing here has been tested against a real catalyst. The discriminating experiments that would be needed are listed at the end of [RESULTS.md](RESULTS.md) in outline, and they centre on a time resolved sulphur balance, a copper inventory measurement over the same run, structural characterisation of the oxide and sulphate, and oxygen switching or isotope work to separate adsorbed oxygen from lattice oxygen pathways.
