Examples
************************

In the following sections we present some simplified examples and use Aronnax to reproduce a number of published results to show how Aronnax can be easily configured for a range of idealised modelling studies.

Canonical examples
===================

These simulations use simple domains and and inputs. They show some of the features available in Aronnax


Gaussian bump on a :math:`\beta`-plane
----------------------------------------

This example is initialised with a Gaussian bump in the layer thickness field. The momentum forcing is set to zero.

1 + 1/2 layers
+++++++++++++++
to do:
- add animated gif of the output


2 layers
+++++++++++
to do:
- add animated gif of the output


Twin gyre on a :math:`\beta`-plane
-------------------------------------

1 + 1/2 layers
+++++++++++++++
to do:
- add animated gif of the output

2 layers
+++++++++++
to do:
- add animated gif of the output

Reproducing published results
===============================

These examples show how Aronnax can be used to reproduce results from the literature.


Davis et al. (2014) - An idealised Beaufort Gyre
-------------------------------------------------
`Davis et al. (2014) <http://dx.doi.org/10.1175/JCLI-D-14-00090.1>`_ used a reduced gravity model to explore the response of an idealised Beaufort Gyre to changes in the seasonal cycle of wind stress. Here we reproduce their control simulation. The domain is set up as a lollipop, with a circular basin for the Beaufort Gyre and a narrow channel connecting it to a region with sponges.

.. figure:: ../reproductions/Davis_et_al_2014/control_final_five/input/wetmask.png
   :width: 45%
   :align: center
   :alt: wetmask defining the domain for Davis et al. (2014)

   The computational domain for Davis et al. (2014). Note: the domain is symmetric, it is the plotting command that makes it look asymmetric.


Over this lollipop basin a wind stress is used to drive an anticyclonic circulation. The magnitude of the wind stress is given as

.. math::
  \frac{1}{r}\int{r \cos^{2}(r)} dr

which is multiplied by :math:`\sin(\theta)` or :math:`-\cos(\theta)` to give the x and y components of the wind stress. Converting the integral into the wind stress requires evaluating :math:`1/r` times the integral as 

.. math::
  \frac{1}{r} \left(\frac{r \sin(2r)}{4} - \frac{\sin^{2}(r)}{4} + \frac{r^{2}}{4}\right)

and normalising the result such that the average wind stress magnitude inside the circular domain is equal to one. This normalised wind stress is then converted into its x and y components.

The y component of the normalised wind stress field is shown on the left, and the time series of wind stress magnitude is on the right.

.. image:: ../reproductions/Davis_et_al_2014/control_final_five/input/tau_y.png
   :width: 37%
.. image:: ../reproductions/Davis_et_al_2014/control_final_five/input/wind_time_series.png
   :width: 62%


After integrating for 40 model years these inputs produce a steady seasonal cycle in velocity and layer thickness. A snap shot is shown on the left, while a time series of the maximum thickness is shown on the right.

.. image:: ../reproductions/Davis_et_al_2014/control_final_five/figures/state_0000155089.png
   :width: 37%
.. image:: ../reproductions/Davis_et_al_2014/control_final_five/figures/h_max.png
   :width: 62%

The seasonal cycle in layer thickness requires a time varying transport through the channel. This is shown below.

.. figure:: ../reproductions/Davis_et_al_2014/control_final_five/figures/transport.png
   :width: 70%
   :align: center
   :alt: time series of transport through the channel

   Time series of transport through the channel due to the seasonal cycle in wind stress.

The paper includes multiple experiments perturbing the seasonal cycle of wind stress. Reproducing the perturbation experiments would require modifying the input variable `wind_mag_time_series_file`.

.. Note:: The configuration used to create these outputs can be found in the reproductions folder of the repository.

Manucharyan and Spall (2016)
-----------------------------
n-layer configuration looking at eddies in the Arctic. (The original experiment was run using a z-level model, but it could also be done in an isopycnal model)


Johnson and Marshall (2002)
----------------------------
Reduced gravity analysis of the adjustment of the MOC to changes in deep water formation rates.
