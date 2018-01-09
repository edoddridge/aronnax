---
title: 'Aronnax: An idealised isopycnal ocean model'
tags:
  - oceanography
  - isopycnal model
  - python
  - geophysical fluid dynamics
authors:
 - name: Edward W. Doddridge
   orcid: 0000-0002-6097-5729
   affiliation: 1
 - name: Alexey Radul
   affiliation: 2
affiliations:
 - name: Massachusetts Institute of Technology, Earth, Atmospheric and Planetary Sciences
   index: 1
 - name: Massachusetts Institute of Technology
   index: 2
date: 5 December 2017
bibliography: paper.bib
---

# Summary

Aronnax is a highly idealised model for simulating large-scale (100-1000km) flows in the ocean. Aronnax is intended for theoretical and empirical oceanographers, as a (relatively) fast and easy-to-use simulation model, bridging the gap between pencil-and-paper on one hand, and more faithful (and complex) computational models on the other. The numerical core is written in Fortran to improve performance, and wrapped in Python to improve usability.

Aronnax is an _isopyncal_ model: it approximates the ocean as a number of discrete homogeneous layers with spatially variable thicknesses. These layers are stacked vertically and the density difference between neighbouring layers is specified by the user. Other widely used vertical coordinates require solving the equations of motion at specified vertical levels where the density is allowed to vary [@Griffies2000]. Representing the large-scale ocean circulation in layers contributes to Aronnax's speed: one needs only a few layers to achieve the same fidelity as 50 or more fixed vertical levels [@Stewart2017].

Aronnax serves three distinct purposes. Firstly, many of the studies that use a model like Aronnax do not provide the source code, see e.g. [@Davis2014,@Fevrier2007,@Johnson2002a,@Stern1998]. This increases the likelihood that coding errors go undetected, and requires that each research group spend time developing their own idealised model. Aronnax solves these problems by providing an open source, tested model for the community to use. Secondly, Aronnax furthers the goals of scientific reproducibility since a simulation can be entirely specified with a set of input files and a version number. Thirdly, Aronnax provides an easy-to-use model that may be included in future modelling hierarchies with minimal effort, thereby enabling new research questions to be addressed.

There are a number of other publicly available ocean models. Of these the most abundant are general circulation models and quasigeostrophic models. General circulation models such as [NEMO](https://www.nemo-ocean.eu/), [GOLD](https://www.gfdl.noaa.gov/gold-ocean-model/), [MOM6](https://github.com/NOAA-GFDL/MOM6), and the [MITgcm](http://mitgcm.org/) solve a less idealised version of the Navier-Stokes equations and can be coupled with sea ice and atmospheric models to create fully coupled climate models. Because the underlying equations are derived with fewer approximations these models can more faithfully simulate a wider range of flow regimes. However, this comes at a price; general circulation models are extremely complex, with numerous free parameters that must be specified, often prior to compiling the source code. It is possible to use most of these models in idealised configurations, but doing so requires a substantial investment of time from the user, and non-trivial computing resources. In comparison, Aronnax is easy to install and cheap to run.

The other abundant class of models is quasigeostrophic models. Geostrophy is a balance between the Coriolis force and the horizontal gradient of the pressure field; flows in which the Coriolis force and the horizontal pressure gradient _almost_ balance are known as quasigeostrophic. Quasigeostrophic models of the ocean make use of this near balance and a number of other assumptions to simplify the equations of motion from a system of five coupled partial differential equations to a single partial differential equation [@Vallis2006]. Quasigeostrophic models range in complexity from [QGcm](http://www.q-gcm.org/), which includes the option of a coupled atmosphere, to doubly periodic quasigeostrophic turbulence models such as [PyQG](http://pyqg.readthedocs.io/en/stable/) and [QGModel](https://github.com/joernc/QGModel). While quasigeostrophic models are extremely useful, there are some problems for which they are ill-suited. For example, the adjustment of the ocean circulation often occurs through ageostrophic motions such as boundary waves [@Johnson2002a], which are not represented in quasigeostrophic models. In addition, quasigeostrophic models are limited in their representation of sloping bathymetry (depth of the sea floor). For these reasons it may be preferable to use an idealised non-linear model such as Aronnax.

Aronnax is MIT licensed and can be retrieved from GitHub at [https://github.com/edoddridge/aronnax](https://github.com/edoddridge/aronnax).

# References
