MHD Equations
=============

Summary
-------

The magnetohydrodynamics (MHD) equations are:

.. math::

    \frac{\partial \rho}{\partial t} + \nabla\cdot(\rho {\bf v}) = 0

    \rho\left(\frac{\partial {\bf v}}{\partial t} + ({\bf v} \cdot \nabla)
     {\bf v} \right) = -\nabla p +
     {1\over\mu}(\nabla\times{\bf B}) \times {\bf B} + \rho {\bf g}

    {\partial {\bf B}\over\partial t}
            = \nabla\times({\bf v}\times{\bf B}) + \eta\nabla^2{\bf B}

    \nabla\cdot{\bf B} = 0

assuming :math:`\eta` is constant.

Derivation
----------

The above equations can easily be derived. We have the continuity equation:

.. math::

    \frac{\partial \rho}{\partial t} + \nabla\cdot(\rho {\bf v}) = 0

Navier-Stokes equations (momentum equation) with the Lorentz force on the
right-hand side:

.. math::

    \rho\left(\frac{\partial {\bf v}}{\partial t} + ({\bf v} \cdot \nabla)
     {\bf v} \right) = -\nabla p + {\bf j} \times {\bf B} + \rho {\bf g}

where the current density :math:`{\bf j}` is given by the Maxwell equation (we
neglect the displacement current :math:`{\partial{\bf E}\over\partial t}`):

.. math::

    {\bf j} = {1\over\mu}\nabla\times{\bf B}

and the Lorentz force:

.. math::

    {1\over\sigma}{\bf j} = {\bf E} + {\bf v}\times{\bf B}

from which we eliminate :math:`{\bf E}`:

.. math::

    {\bf E} = - {\bf v}\times{\bf B} + {1\over\sigma}{\bf j} = 
              - {\bf v}\times{\bf B} + {1\over\sigma\mu}\nabla\times{\bf B}

and put it into the Maxwell equation:

.. math::

    {\partial {\bf B}\over\partial t} = -\nabla\times{\bf E}

so we get:

.. math::

    {\partial {\bf B}\over\partial t} = \nabla\times({\bf v}\times{\bf B})
                - \nabla\times\left({1\over\sigma\mu}\nabla\times{\bf B}\right)

assuming the magnetic diffusivity :math:`\eta={1\over\sigma\mu}` is constant, we
get:

.. math::

    {\partial {\bf B}\over\partial t} = \nabla\times({\bf v}\times{\bf B})
                - \eta\nabla\times\left(\nabla\times{\bf B}\right)
            = \nabla\times({\bf v}\times{\bf B})
                + \eta\left(\nabla^2{\bf B}-\nabla(\nabla\cdot{\bf B})\right)
            = \nabla\times({\bf v}\times{\bf B}) + \eta\nabla^2{\bf B}

where we used the Maxwell equation:

.. math::

    \nabla\cdot{\bf B} = 0
