# Landscape of fear model

This repository contains an implementation of the ODE and ABM model in our article:

Daniel Kondor, James S. Bennett, Detlef Gronenborn, Peter Turchin (2023).<br>
Landscape of Fear: Indirect Effects of Conflict Can Account for Large-scale Population Declines in Non-state Societies.
_SocArXiv preprint_. https://osf.io/preprints/socarxiv/m7hxw.

The ODE version is self-contained under [ode_model]. The ABM model under [abm_model] requires the preprocessed data about the European landscape that was used in our [previous model](https://github.com/dkondor/neolithic_simulation). Specifically, this means performing the preprocessing steps 3-6 [here](https://github.com/dkondor/neolithic_simulation/tree/main/preprocessing) (note: no need for the climate data set in this case) after downloading the datasets based on the instructions [here](https://github.com/dkondor/neolithic_simulation/tree/main/data) (again, the climate and AgMERRA data can be omitted).


