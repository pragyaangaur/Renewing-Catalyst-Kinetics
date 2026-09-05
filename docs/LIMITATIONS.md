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

Rate constants spanning nine orders of magnitude make the stationary balance matrix badly conditioned. The solver tries a dense least squares solve and falls back to the Grassmann, Taksar and Heyman reduction, which cannot produce a negative probability through cancellation. Every solve is checked for one closed communicating class, a nonnegative distribution, a small relative residual and closed element balances. Even so, stationary rates on live patches can fall to solver precision under extreme draws, and the reachability test rather than the numerical rate is what establishes Result A.

## No experimental validation

Nothing here has been tested against a real catalyst. The discriminating experiments that would be needed are listed at the end of [RESULTS.md](RESULTS.md) in outline, and they centre on a time resolved sulphur balance, a copper inventory measurement over the same run, structural characterisation of the oxide and sulphate, and oxygen switching or isotope work to separate adsorbed oxygen from lattice oxygen pathways.
