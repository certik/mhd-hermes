print "Importing..."
from enthought.mayavi import mlab
from enthought.mayavi.filters.warp_scalar import WarpScalar
from enthought.mayavi.modules.surface import Surface
from enthought.mayavi.filters.extract_vector_norm import ExtractVectorNorm
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

mlab.show()
