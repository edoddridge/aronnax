---
title: 'Aronnax: An idealised isopycnal ocean model'
tags:
  - example
  - tags
  - for the paper
authors:
 - name: Arfon M Smith
   orcid: 0000-0000-0000-1234
   affiliation: "1, 2" # (Multiple affiliations must be quoted)
 - name: Mickey Mouse
   orcid: 0000-0000-0000-1234
   affiliation: 2
affiliations:
 - name: Space Telescope Science Institute
   index: 1
 - name: Disney Inc.
   index: 2
date: 14 February 2016
bibliography: paper.bib
---

# Summary

Aronnax is a highly idealised model that approximates the ocean as a number of discrete layers with constant density. The numerical core is written in Fortran to improve performance, and wrapped in Python to improve usability. 

Aronnax serves three distinct purposes. Firstly, many of the studies that use a model like Aronnax do not provide the source code, thereby increasing the likelihood that coding errors go undetected. Aronnax solves this problem by providing an open source, tested, model for the community to use. Secondly, Aronnax furthers the goals of scientific reproducibility, since a simulation can be entirely specified with a set of input files and a version number. Thirdly, Aronnax provides an easy-to-use model that may be included in future modelling hierarchies with minimal effort, thereby enabling new research questions to be addressed.


There are a number of other publicly available ocean models, of these the most abundant are general circulation models and quasigeostrophic models. General circulation models such as [NEMO](https://www.nemo-ocean.eu/), [GOLD](https://www.gfdl.noaa.gov/gold-ocean-model/), [MOM6](https://github.com/NOAA-GFDL/MOM6), the [MITgcm](http://mitgcm.org/) can be coupled with sea ice and atmospheric models to create fully coupled climate models. These are extremely complex models, and require a substantial amount of time from the user to install, setup, and run. The complexity of the models means that there are numerous parameters that must be specified and dependencies that must be satisfied. It is possible to use most of these models in idealised configurations, but doing so requires a substantial investment of time from the user, and non-trivial computing resources. In comparison Aronnax is easy install and cheap to run. The other class of abundant models is quasigeostrophic models. These range in complexity from [QGcm](http://www.q-gcm.org/), which includes the option of a coupled atmosphere, to doubly periodic quasigeostrophic turbulence models such as [PyQG](http://pyqg.readthedocs.io/en/stable/) and [QGModel](https://github.com/joernc/QGModel). While quasigeostrophic models are extremely useful, there are some problems for which they are ill-suited. For example, the adjustment of the ocean circulation often occurs through ageostrophic motions such as boundary waves [@Johnson2002a] which are not represented in quasigeostrophic models. In addition, quasigeostrophic models are limited in their representation of sloping bathymetry. For these reasons it may be preferable to use an idealised non-linear model such as Aronnax.


Aronnax is MIT licensed and can be retrieved from GitHub https://github.com/edoddridge/aronnax


# References
  