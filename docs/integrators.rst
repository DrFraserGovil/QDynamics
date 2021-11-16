.. integrators

#######################
Integrators
#######################

An integrator is a section of code which is capable of approximating the evolution x(t) given the initial conditions x(0) and v(0). Different integrators have different strengths and weaknesses. The purpose of Fraser & Schonrich 2021 was to investigate the differences between the integrators provided below.

We provide three classes along an inheritance chain to aid in these calculations. The ``Integrator`` class provides the basic framework of an integrator, including the main evolution loop and the procedures for saving information to file etc. The ``Magi`` integrator is a generalised subclass of ``Integrator`` for the case of rotating systems without nice symmetry properties, and numerically computes a truncated form of the Magnus series. The ``Symi`` integrator in turn subclasses ``Magi`` and overloads the Magnus series computation to the exact analytical expressions for axisymmetric systems.

.. toctree::
   integrator_chassis
   magi
   symi
   :maxdepth: 2
   :caption: Classes:
	
Integrator Enums
--------------------

Both the ``Magi`` and ``Symi`` classes can undergo different step-types, determining how they move from :math:`\mathsf{q}_{i+1}  \to \mathsf{q}_i`. The relevant classes determine which of these to used based on the template argument, UpdateType.

.. doxygenenum:: QDynamics::UpdateType
