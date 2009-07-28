print "Importing..."
from enthought.mayavi import mlab
from enthought.mayavi.filters.warp_scalar import WarpScalar
from enthought.mayavi.modules.surface import Surface
from enthought.mayavi.filters.extract_vector_norm import ExtractVectorNorm
import numpy
from numpy import array
from glob import glob
import re
print "  done."

engine = mlab.get_engine()
vtk_file_reader = engine.open(u'/home/ondrej/repos/mhd-hermes/frame_vec0030.vtk')
extract_vector_norm = ExtractVectorNorm()
engine.add_filter(extract_vector_norm, obj=None)
surface = Surface()
engine.add_module(surface, obj=None)

# show 3D model of the vector norms:
#warp_scalar = WarpScalar()
#engine.add_filter(warp_scalar, obj=None)
#surface1 = Surface()
#engine.add_filter(surface1, warp_scalar)

x, y, z = numpy.mgrid[0:16, 0:6, 0.1:1]
r = numpy.sqrt(x**2 + y**2 + z**4)
u = y*numpy.sin(r)/(r+0.001)
v = -x*numpy.sin(r)/(r+0.001)
w = numpy.zeros_like(z)
mlab.quiver3d(x, y, z, u, v, w, line_width=3, scale_factor=1)

mlab.show()
