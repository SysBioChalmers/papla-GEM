# papla-GEM: Genome-scale metabolic model of _Papiliotrema laurentii_

[![DOI](https://zenodo.org/badge/xxx.svg)](https://zenodo.org/badge/latestdoi/xxx) [![GitHub version](https://badge.fury.io/gh/sysbiochalmers%2Fpapla-gem.svg)](https://badge.fury.io/gh/sysbiochalmers%2Fpapla-gem) 

- Brief model description

This repository contains the current genome-scale metabolic model of _Papiliotrema laurentii UFV-1_, named **papla-gem**. For the latest updated release see [here](https://github.com/SysBioChalmers/papla-gem/releases).

- Abstract

_Papiliotrema laurentii UFV-1_.

- Model keywords

**GEM category:** Species; **Utilisation:** experimental data reconstruction; **Field:** metabolic-network reconstruction; **Type of model:** reconstruction, curated; **Model source:** [rhto-GEM](https://github.com/SysBioChalmers/rhto-GEM); **Omic source:** genomics; **Taxonomy:** _Papiliotrema laurentii_; **Metabolic system:** general metabolism; **Bioreactor**; **Strain:** UFV-1; **Condition:** minimal medium;

- Reference:  

- Last update: 2019-05-29

- Main model descriptors:

| Taxonomy | Template Model | Reactions | Metabolites | Genes |
| ------------- |:-------------:|:-------------:|:-------------:|:-----:|
| _Papiliotrema laurentii_|	[rhto-GEM](https://github.com/SysBioChalmers/rhto-GEM) | xx | xx | xx |

A [Memote](https://memote.readthedocs.io/en/latest/) snapshot report of the most recent release is available [here](https://SysBioChalmers.github.io/papla-gem).

This repository is maintained through a collaboration between Wendel Batista da Silveira (affiliation) and Eduard Kerkhoven ([@edkerk](https://github.com/edkerk/)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.

## Citation

* If you use papla-gem in your research, whether for simulations, data analysis, model reconstruction of other purposes, we ask you to cite the papla-gem paper on [bioRxiv](https://doi.org/10.1101/xxx).
* In addition, it is good practice to cite the specific version of papla-gem that you used, to improve reproducibility. All papla-gem releases are archived in [Zenodo](https://zenodo.org/badge/latestdoi/xx). Find the corresponding DOI for each release [here](https://zenodo.org/search?page=1&size=20&q=conceptrecid:xx&sort=-publication_date&all_versions=True).

## Installation

### Required software

  * This model is recommended to be used with the [**RAVEN toolbox for MATLAB**](https://github.com/SysBioChalmers/RAVEN) (version 2.0).
  * Alternatively, the model can also be directly used with the [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox), [cobrapy](https://github.com/opencobra/cobrapy), or other software or toolboxes that support genome-scale models in SBML L3V1 FBCv2 format.

### Dependencies
* Please see the [RAVEN toolbox](https://github.com/SysBioChalmers/RAVEN) repository for dependencies regarding RAVEN.
* For contribution to development: a [git wrapper](https://github.com/manur/MATLAB-git) added to the search path.

### Installation instructions
* Just want to use the model? Clone it from [`master`](https://github.com/SysBioChalmers/papla-gem) in the Github repo, or just download [the latest release](https://github.com/SysBioChalmers/papla-gem/releases).
* Wish to also contribute? Fork it to your Github account, and create a new branch from [`devel`](https://github.com/SysBioChalmers/papla-gem/tree/devel).

## Model files

The model is available in `.xml`, `.txt`, `.yml`, `.mat` and `.xlsx` (the last 2 extensions only in `master`). Additionally, versions of toolboxes & SBML used for saving the model are tracked in `dependencies.txt`.

### Complementary scripts
### Complementary data

## Contributors

* [Eduard J. Kerkhoven](https://www.chalmers.se/en/staff/Pages/Eduard-Kerkhoven.aspx) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Gothenburg Sweden
