# code

This folder contains all the code used in generating and maintaining papla-GEM.

* `reconstructionProtocol.m` - the main script that reconstructs the model from template GEMs and other input data in `data/`.
* `newCommit.m` - exports the model (`.txt`, `.xml`, `.yml`) for a commit to a development branch. Run from this folder.
* `newRelease.m` - bumps the model version, exports all model file formats (including binaries), and updates `version.txt`, `history.md` and the stats table in `README.md`. Run from this folder, on `main` only.
* `getEarlierModel.m` - retrieves a model file from an earlier commit or release.

### analysis

Scripts used to analyse and simulate the model, e.g. growth predictions, FBA of lipid production, and identification of overexpression targets.

### curation

Scripts used to curate and clean up the model during reconstruction, e.g. fitting the growth-associated energy costs (GAEC).

### lipidMetabolism

Scripts used to add lipid-related reactions to the model, including the SLIME reactions used to represent lipid species and their scaling.
