<div class="container-fluid mb-sm-4 mb-xl-2">

<div class="container">

<a href="./index.html" class="navbar-brand">CDMSreader</a>

<div id="navbar" class="navbar-collapse collapse">

- <a href="./lists/files.html" class="nav-link">Source Files</a>
- <a href="./lists/modules.html" class="nav-link">Modules</a>
- <a href="./lists/procedures.html" class="nav-link">Procedures</a>
- <a href="./lists/types.html" class="nav-link">Derived Types</a>

<div class="d-flex align-items-end flex-grow-1">

</div>

</div>

</div>

</div>

<div class="container">

<div id="jumbotron" class="p-5 mb-4 bg-light border rounded-3">

Reads molecular transition data taken from the Cologne Database for
Molecular Spectroscopy ([CDMS](https://cdms.astro.uni-koeln.de/))
catalogue and determines the lifetimes of various quantum states based
on the provided information.

</div>

<div id="text" class="row">

<div class="col-md-12">

# CDMSreader

Reads molecular transition data from the Cologne Database for Molecular
Spectroscopy ([CDMS](https://cdms.astro.uni-koeln.de/)) catalogue to
determine the radiative lifetimes of the various states involved in the
transitions. The data must be obtained in the proper format by using the
[search and conversion
form](https://cdms.astro.uni-koeln.de/cgi-bin/cdmssearch). The user
should ask for $A$ values with the intensities being given as log
values. Units can be in GHz or inverse centimeters.

For now, this only works for one molecule at a time and only for
asymmetric top molecules that correspond to the case of $Q=23$ (see the
[CDMS general](https://cdms.astro.uni-koeln.de/classic/general)
documentation for more details).

Two sets of output data are produced. The first set contains the
lifetimes for the states with the quantum number $F$, the second does a
statistical average over $F$.

## Todo:

- \[ \] Add more state types, generalize the types module
- \[ \] Organize the types module into submodules for better clarity
- \[ \] Add other state types corresponding to different Qs
- \[ \] Add the ability to work with different isotopologues/molecules
  at once so that we can just use one file

## Determining the liftimes of excited states.

Atoms and molecules can be found in excited states in various media, but
these excited states have finite lifetimes because they spontaneously
decause to lower-lying states. The ground state is only stable because
there is no lower state to which the system could decay on its own. When
an atom/molecule spontaneously decays, it emits radiation as a way of
cooling, so to speak. The Einstein A coefficients $A_{u\to l}$ describe
the probability per unit time that the atom/molecule in level $u$ will
spontaneously decay to level $l$, emitting a photon of frequency $\nu$.
The probability per unit time that the system decays into any of the
available lower states is just given by the sum of the Einstein A
coefficients for transitions starting from that state:
$$ A_u = \sum\limits_l A_{u\to l}. $$ We can then take the invese of
this to get the average lifetime of a specific state, $\tau_u = 1/A_u$.

The CDMS [search and conversion
form](https://cdms.astro.uni-koeln.de/cgi-bin/cdmssearch) already does
the work of determining the Einstein A coefficients. This code just
reads several transitions, determines which states are involved, and
then just adds the Einstein A coefficients for the corresponding state
to determine the lifetimes. You can add headers to the file or stream
containing the transitions so long as they're commented with `#`. The
format is very delicate, so make sure not to alter the transition lines.

## Installation

This is available as a Fortran Package Manager
([fpm](https://fpm.fortran-lang.org/)) package, so it can just be built
with the usual build command in the cloned repository

<div class="codehilite">

    fpm build

</div>

The package can be built without `fpm`. On linux, just use the provided
`compile` script. Other compilers than gfortran will probably work just
fine.

## Usage

The program expects to read transitions from standard input, so you can
just feed it in like so:

<div class="codehilite">

    ./CDMSreader < transitions.dat > lifetimes.txt

</div>

or even with `fpm` if you built the program that way

<div class="codehilite">

    fpm run < transitions.dat > lifetimes.txt

</div>

</div>

</div>

<div class="row">

------------------------------------------------------------------------

<div class="col-xs-6 col-sm-3">

<div>

### Source Files

- [CDMSreader\_\_constants.f](sourcefile/cdmsreader__constants.f.html)
- [CDMSreader\_\_readwrite.f](sourcefile/cdmsreader__readwrite.f.html)
- [CDMSreader\_\_system.f](sourcefile/cdmsreader__system.f.html)
- [CDMSreader\_\_types.f](sourcefile/cdmsreader__types.f.html)

</div>

<div>

- [*All source files…*](./lists/files.html)

</div>

</div>

<div class="col-xs-6 col-sm-3">

<div>

### Modules

- [CDMSreader\_\_constants](module/cdmsreader__constants.html)
- [CDMSreader\_\_readwrite](module/cdmsreader__readwrite.html)
- [CDMSreader\_\_system](module/cdmsreader__system.html)
- [CDMSreader\_\_types](module/cdmsreader__types.html)

</div>

<div>

- [*All modules…*](./lists/modules.html)

</div>

</div>

<div class="col-xs-6 col-sm-3">

<div>

### Procedures

- [add_to](interface/add_to.html)
- [CDMS_readline](proc/cdms_readline.html)
- [die](proc/die.html)
- [find_state_number](proc/find_state_number.html)
- [make_asymtop_state](interface/make_asymtop_state.html)
- [sort_last_state](proc/sort_last_state.html)

</div>

<div>

- [*All procedures…*](./lists/procedures.html)

</div>

</div>

<div class="col-xs-6 col-sm-3">

<div>

### Derived Types

- [asymtop_state](type/asymtop_state.html)
- [asymtop_state_hfs](type/asymtop_state_hfs.html)
- [asymtop_state_nohfs](type/asymtop_state_nohfs.html)
- [asymtop_transition](type/asymtop_transition.html)

</div>

<div>

- [*All derived types…*](./lists/types.html)

</div>

</div>

</div>

------------------------------------------------------------------------

</div>

<div class="container">

<div class="row justify-content-between">

<div class="col">

CDMSreader © 2025 gpl3

</div>

<div class="col">

Documentation generated by
[FORD](https://github.com/Fortran-FOSS-Programmers/ford)

</div>

</div>



</div>
