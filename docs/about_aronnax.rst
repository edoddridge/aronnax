About Aronnax
********************

The Physics
============

Aronnax can run in two different modes:
 * n + 1/2 layer mode, in which the bottom layer is quiescent and infinitely deep
 * n layer mode, in whihc the bottom topography is specified and all layers have finite thickness


n + 1/2 layer mode
-------------------
This is very cheap to run

n layer mode
--------------------
More expensive


Governing Equations
=====================

The model solves the hydrostatic Boussinesq equations within a finite number of discrete isopycnal layers. At the layer interfaces there is a discrete jump in velocity and density, but not in pressure.

Continuity equation
-------------------
.. math::
    \label{eqn:reduced_gravity_layer_continuity} 
    \frac{\partial h_{n}}{\partial t} = \mathbf{\nabla} \cdot \left(h_{n} \mathbf{v_{n}} \right).


in which :math:`\mathbf{v_{n}}` represents the vertically averaged horizontal velocity in layer $n$.

Momentum equations
-------------------
.. math::
    \label{eqn:reduced_grav_layer_1_momentum} 
    \frac{D \mathbf{v_{n}}}{D t} +  \mathbf{f} \times \mathbf{v_{n}} + g^{'}\mathbf{\nabla}h_{n} = \mathbf{F_{n}},



in which :math:`g^{'}` is the reduced gravity given by :math:`{g(\rho_{2} - \rho_{1})}/{\rho_{1}}`. The reduced gravity is dynamically equivalent to gravity, but is scaled to take into account the density difference between the two layers.

This can be rewritten in terms of the Bernoulli Potential to give,

.. math::
    \label{eqn:momentum_Bernoulli_form}
    \frac{\partial\mathbf{v_{n}}}{\partial t} - (f+\zeta_{n}) \times v_{n} + \nabla \Pi_{n} + = \kappa \nabla^{2}v_{n} + \frac{\mathbf{F_{n}}}{\rho_{0}}

where :math:`\Pi_{n}` is the Bernoulli potential, :math:`\left(\mathbf{v_{n}}\cdot\mathbf{v_{n}}\right)/2 + p/\rho_{0}`, and :math:`p` is the hydrostatic pressure. In this form the non-linearity from the material derivative has been moved into the Bernoulli Potential and the vorticity term. 



The model can be used in either reduced gravity mode, with a quiescent abyss, or in n-layer mode with bathymetry. In the n-layer case the model can either be run with a rigid lid, or with a free surface. In simulations with a free surface the following equation is also solved

.. math::
    \label{eqn:}
    \frac{\partial \eta}{\partial t} + \mathbf{\nabla} \cdot (H \mathbf{V}) = 0,

where :math:`H` is the depth from the free-surface to the bathymetry, and :math:`V` is the vertically averaged flow, the barotropic flow. With a rigid lid, the model solves an analogous equation, but it is just for the pressure field to keep the vertically integrated horizontal flow divergence free - the result is not carried from one timestep to the next.