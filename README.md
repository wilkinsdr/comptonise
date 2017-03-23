# comptonise
XSPEC convolution model to Comptonise an input spectrum through a corona. The Kompaneets equation is solved using the input spectrum as the seed photons (for example using the the disc reflection spectrum as the seed photons to simulate the passage of the reflected radiation through the corona, as described in Wilkins & Gallo 2015, MNRAS 448, 703.

XSPEC> model comptonise*[seed_photon_model]

Care needs to be taken when evaluating the Comptonisation model. In particular note that XSPEC only evaluates the input spectrum over the energy range defined by the instrument response. Comptonisation will cause lower energy seed photons to be scattered up into the instrument bandpass. It is therefore necessary to extend the energy response matrix as low as is permitted by the seed photon model to capture all of these photons (at least a factor of 10 below the minium energy of interest).

XSPEC> energies extend low 0.1 100 log

Based on nthcomp; the thermally comptonized continuum model of Zdziarski, Johnson & Magdziarz 1996, MNRAS, 283, 193, as extended by Zycki, Done & Smith 1999, MNRAS 309, 561.

The XSPEC model function is written in C++ (compconv.cxx) but calls Fortran routines to perform the calculation (docompconf.f). The XSPEC model definition is in lmodel_comptonise.dat.

Parameters
----------
1. Gamma    : Photon index of the produced continuum spectrum if >0, the optical depth is calculated automatically to reproduce this photon index from a black body seed spectrum based on the input value of kT. If frozen at -1, the spectrum will be computed for the specific optical depth that is input.
2. Tau_e    : Optical depth due to Thomson scattering through the corona. If Gamma > 0, the input value will be ignored.
3. kT_e     : Electron temperature in the corona
4. Redshift : Applied to the output spectrum

Installation
------------

The model is compiled using the initpackage command in XSPEC:

XSPEC> cd comptonise/model/directory

XSPEC> initpackage comp lmodel_comptonise.dat .

Then can be loaded from the directory it is compiled in (it will need to be loaded every time XSPEC is started):

XSPEC> lmod comp .

Alternatively, the code and model definition can be integrated into a local model package
