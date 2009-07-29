print "Importing..."
from enthought.mayavi import mlab
from enthought.mayavi.filters.warp_scalar import WarpScalar
from enthought.mayavi.modules.surface import Surface
from enthought.mayavi.modules.vectors import Vectors
from enthought.mayavi.filters.extract_vector_norm import ExtractVectorNorm
import numpy
from numpy import array
from glob import glob
import re
print "  done."

engine = mlab.get_engine()
vtk_file_reader = engine.open(u'/home/ondrej/repos/mhd-hermes/frame_vec0030.vtk')

vectors = Vectors()
engine.add_module(vectors, obj=None)

extract_vector_norm = ExtractVectorNorm()
engine.add_filter(extract_vector_norm, obj=None)
surface = Surface()
engine.add_module(surface, obj=None)

# show 3D model of the vector norms:
#warp_scalar = WarpScalar()
#engine.add_filter(warp_scalar, obj=None)
#surface1 = Surface()
#engine.add_filter(surface1, warp_scalar)

def better():
    g = vtk_file_reader.outputs[0]
    points = array(g.points)
    vectors = array(g.point_data.vectors)
    assert len(points) == len(vectors)

    def dist(p, q):
        return (p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2

    def find_nearest_point(x, y, z):
        min = 10**8
        _id = None
        for id, p in enumerate(points):
            d = dist(p, (x, y, z))
            if d < min:
                min = d
                _id = id
        return _id

    x, y, z = numpy.mgrid[0:15.1:0.5, 0:5.1:0.5, 0.1:1]

    u = numpy.sin(x)
    v = numpy.sin(y)
    w = numpy.zeros_like(z)
    print "finding nearest stuff"
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            for k in range(x.shape[2]):
                id = find_nearest_point(x[i, j, k], y[i, j, k], z[i, j, k])
                u[i, j, k] = vectors[id][0]
                v[i, j, k] = vectors[id][1]
                w[i, j, k] = vectors[id][2]
                print i, j, k, id
    print "ok"
    #import IPython; IPython.embed()
    mlab.quiver3d(x, y, z, u, v, w, line_width=1, scale_factor=1, color=(0, 0, 0))

vectors.glyph.color_mode = 'no_coloring'
vectors.actor.property.color = (0, 0, 0)
vectors.glyph.mask_input_points = True
vectors.glyph.mask_points.random_mode = False
vectors.glyph.mask_points.on_ratio = 200

mlab.show()
