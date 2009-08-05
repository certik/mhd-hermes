Installation instructions
-------------------------

Install `hermes2d <http://hpfem.org/>`_, so that you can ``import hermes2d``
from Python::

    In [1]: import hermes2d

    In [2]:

Once this works, then just run::

    cmake .
    make

and that's it (cmake will ask the ``hermes2d`` module where all the ``*.h`` and
``*.pxd`` files are).
