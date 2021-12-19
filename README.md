# papla-GEM: Genome-scale metabolic model of _Papiliotrema laurentii_

## Description

This repository contains the current consensus genome-scale metabolic model of _Papiliotrema laurentii_ UFV-1, named **papla-GEM**. The model distributed on this GitHub repository is continuously updated, with the latest releases available [here](https://github.com/SysBioChalmers/papla-GEM/releases). To get access to the model associated to the Ventorim _et al_. (2022) publication, use [papla-GEM 1.1.0](https://github.com/SysBioChalmers/papla-GEM/releases/tag/1.1.0).

## Citation

* The manuscript has been submitted.

## Keywords:

**Utilisation:** experimental data reconstruction; _in silico_ strain design  
**Field:** metabolic-network reconstruction  
**Type of model:** reconstruction; curated  
**Model source:** [rhto-GEM](https://github.com/SysBioChalmers/rhto-GEM)  
**Omic source:** [genomics](https://doi.org/10.1016/j.fgb.2020.103456  )
**Taxonomic name:** _Papiliotrema laurentii_  
**Taxonomy ID:** [taxonomy:5418](https://identifiers.org/taxonomy:5418)  
**Genome ID:** [insdc.sra:SRR10766837](https://identifiers.org/insdc.sra:SRR10766837)  
**Metabolic system:** general metabolism  
**Strain:** UFV-1  
**Condition:** minimal medium

## Installation & Usage

### **User:**

To obtain papla-GEM, clone it from [`main`](https://github.com/sysbiochalmers/papla-GEM) in the GitHub repository, or just download the [latest release](https://github.com/sysbiochalmers/papla-GEM/releases).

papla-GEM is distributed in SBML L3V1 FBCv1 format (`model/papla-GEM.xml`), and therefore works well with any appropriate constraint-based modelling package, such as [RAVEN Toolbox](https://github.com/sysbiochalmers/raven/), [cobrapy](https://github.com/opencobra/cobrapy), and [COBRA Toolbox](https://github.com/opencobra/cobratoolbox). Installation instructions for each package are provided on their website, after which you can use their default functions for loading and exporting of the models:

***RAVEN Toolbox***
```matlab
model = importModel('papla-GEM.xml')
exportModel(model, 'papla-GEM.xml')
```

***cobrapy***
```python
import cobra
model = cobra.io.read_sbml_model('papla-GEM.xml')
cobra.io.write_sbml_model(model, 'papla-GEM.xml')
```

***COBRA Toolbox*** \*
```matlab
model = readCbModel('papla-GEM.xml')
writeCbModel(model, 'papla-GEM.xml')
```
\* note that some annotation might be lost when exporting the model from COBRA Toolbox.

### **Contributor:**

Development of the model is done via RAVEN, to ensure that model content is retained as much as possible (I/O through other software might result in undesired loss of annotation).

[Fork](https://github.com/sysbiochalmers/papla-GEM/fork) the papla-GEM repository to your own GitHub account, and create a new branch from `devel`.

Load the model in MATLAB using the default code specified [above](#user). Before making a pull-request to the `devel` branch, export the model with the `newCommit` function provided in the repository:
```matlab
cd ./code
newCommit(model);
```
