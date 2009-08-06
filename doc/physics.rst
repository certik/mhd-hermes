MHD Equations
-------------

The magnetohydrodynamics equations are continuity equation:

.. math::

    \frac{\partial \rho}{\partial t} + \nabla\cdot(\rho {\bf v}) = 0

Navier-Stokes equations (momentum equation) with the Lorentz force on the
right-hand side:

.. math::

    \rho\left(\frac{\partial {\bf v}}{\partial t} + ({\bf v} \cdot \nabla)
     {\bf v} \right) = -\nabla p + {\bf j} \times {\bf B} + \rho {\bf g}

where the current density :math:`{\bf j}` is given by the equations:

.. math::

    {\bf j} = {1\over\mu}\nabla\times{\bf B}

    {1\over\sigma}{\bf j} = {\bf E} + {\bf v}\times{\bf B}
