.. iEBE-MUSIC documentation master file, created by
   sphinx-quickstart on Thu Jan  2 12:30:32 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to iEBE-MUSIC's documentation!
######################################

The iEBE-MUSIC is a general framework to simulate the multi-stage
3D dynamics of relativistic heavy-ion collisions event-by-event.

Introduction and Background
***************************

Collective phenomena emerged from many-body Quantum Chromodynamics (QCD)
are studied by colliding ionized heavy nuclei at high energy, relativsitic
heavy-ion collisions. These collisions can create extreme conditions, under
which quarks and gluons are deformed from individual hadrons and form
a new state of matter called the Quark-Gluon Plasma.
The dynamics during the relativistic heavy-ion collisions is extremly short.
During only a few yatoseconds (~10^{-23}), the collision system goes through
a complex multi-stage evolution.
First, it evolves from a pre-equilibrium phase to a strongly coupled 
fluid-dynamic regime. As the system further expands and cools, the
quarks and gluons will recombined into hadrons. These hadrons will scatter
with each other and eventually decoupled and fly freely to the detectors.
Because of this ultra-fast dynamics, we can not directly probe the QGP during
the collisions. Theoretical simulations are essential to rewind the final
measured particle momentum distributions to early stages of the collisions.
These phenomenological studies help us to study the properties of the
strongly-coupled QGP.


.. image:: figs/QGP.jpg
    :width: 500px
    :align: center
    :alt: An illustration of the dynamical evolution of relativistic 
          heavy-ion collisions.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

Quickstart
----------
First, install the dependencies:


User guide
----------
.. toctree::
   :maxdepth: 2

    installation
    docker
    usage
    examples


Attribution
-----------
If you make use of this software in your research, please cite


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
