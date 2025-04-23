Reads molecular transition data from the Cologne Database for Molecular Spectroscopy ([CDMS](https://cdms.astro.uni-koeln.de/)) catalogue
to determine the radiative lifetimes of the various states involved in the transitions.
The data must be obtained in the proper format by using the [search and conversion form](https://cdms.astro.uni-koeln.de/cgi-bin/cdmssearch).
The user should ask for A values with the intensities being given as log values.
Units can be in GHz or inverse centimeters.

For now, this only works for one molecule at a time and only for asymmetric top molecules that correspond to the case of Q=23 (see the [CDMS general](https://cdms.astro.uni-koeln.de/classic/general) documentation for more details).

Two sets of output data are produced.
The first set contains the lifetimes for the states with the quantum number F, the second does a statistical average over F.

Docs available [here](https://banana-bred.github.io/CDMSreader/).

## Todo:
- [ ] Add more state types, generalize the types module
- [ ] Organize the types module into submodules for better clarity
- [ ] Add other state types corresponding to different Qs
- [ ] Add the ability to work with different isotopologues/molecules at once so that we can just use one file

## Determining the liftimes of excited states.
Atoms and molecules can be found in excited states in various media, but these excited states have finite lifetimes because they spontaneously decause to lower-lying states.
The ground state is only stable because there is no lower state to which the system could decay on its own.
When an atom/molecule spontaneously decays, it emits radiation as a way of cooling, so to speak.
The Einstein A coefficients $$A_{u\to l}$$ describe the probability per unit time that the atom/molecule in level u will spontaneously decay to level l, emitting a photon of frequency ν.
The probability per unit time that the system decays into any of the available lower states is just given by the sum of the Einstein A coefficients for transitions starting from that state:
$$A_u = \sum\limits_l A_{u\to l}.$$
We can then take the invese of this to get the average lifetime of a specific state, $$\tau_u = 1/A_u$$.

The CDMS [search and conversion form](https://cdms.astro.uni-koeln.de/cgi-bin/cdmssearch) already does the work of determining the Einstein A coefficients.
This code just reads several transitions, determines which states are involved, and then just adds the Einstein A coefficients for the corresponding state to determine the lifetimes.
You can add headers to the file or stream containing the transitions so long as they're commented with `#`.
The format is very delicate, so make sure not to alter the transition lines.

**IMPORANT NOTE**

States that are involved in transitions as lower states but never as upper states will not have an A coefficient, and will therefore have an infinite lifetime.
This may be because it's the ground state, but it might also just be because there are no transitions from that state in the provided input data.

## Installation
This is available as a Fortran Package Manager ([fpm](https://fpm.fortran-lang.org/)) package, so it can just be built with the usual build command in the cloned repository
```
fpm build
```
The package can be built without `fpm`.
On linux, just use the provided `compile` script.
Other compilers than gfortran will probably work just fine.

## Usage
The program expects to read transitions from standard input, so you can just feed it in like so:
```
./CDMSreader < transitions.dat > lifetimes.txt
```
or even with `fpm` if you built the program that way
```
fpm run < transitions.dat > lifetimes.txt
```
