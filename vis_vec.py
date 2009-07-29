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

frames = glob("anim/frame_vec*.vtk")
frames.sort()
max_frame_number = int(re.search("(\d+)", frames[-1]).groups()[0])

print "plotting the main frame..."
width = 1000
height = 360
mlab.figure(size=(width, height+32))
mlab.options.offscreen = True

engine = mlab.get_engine()
vtk_file_reader = engine.open(u'anim/frame_vec0000.vtk')

#vectors = Vectors()
#engine.add_module(vectors, obj=None)

extract_vector_norm = ExtractVectorNorm()
engine.add_filter(extract_vector_norm, obj=None)
surface = Surface()
engine.add_module(surface, obj=None)

# show 3D model of the vector norms:
#warp_scalar = WarpScalar()
#engine.add_filter(warp_scalar, obj=None)
#surface1 = Surface()
#engine.add_filter(surface1, warp_scalar)

print "delaunay2d"
field = mlab.pipeline.delaunay2d(vtk_file_reader)
#print "setting tolerance"
#import IPython; IPython.embed()
# this looks ugly:
#field.filter.tolerance = 0.1065
print "grid"
x_g, y_g, z_g = numpy.mgrid[0:15:60j, 0:5:60j, 0:1:1j]
print "scalar_scatter"
sampling_grid = mlab.pipeline.scalar_scatter(x_g, y_g, z_g)

print "ProbeFilter"
filter = mlab.pipeline.user_defined(sampling_grid, filter='ProbeFilter')
print "outputs"
filter.filter.source = field.outputs[0]
# this doesn't look as good:
#filter.filter.source = vtk_file_reader.outputs[0]
#filter.filter.source = extract_vector_norm.outputs[0]
#filter.filter.source = surface.outputs[0]

print "vectors"
vectors = mlab.pipeline.vectors(filter, mode='2darrow')

print "polishing"
vectors.glyph.color_mode = 'no_coloring'
vectors.actor.property.color = (0, 0, 0)
vectors.glyph.glyph.scale_factor = 0.25
vectors.glyph.glyph_source.glyph_source.scale = 1.033
vectors.glyph.glyph_source.glyph_source.center = array([0.5, 0, 0])
print "done"

mlab.view(0.0, 0.0, 10.5, array([ 7.569, 2.5,  0]))
mlab.roll(0)

#mlab.show()

for i in range(max_frame_number+1):
    print "doing:", i
    vtk_file_reader.timestep = i
    vectors.actor.property.color = (0, 0, 0)
    vectors.glyph.glyph.scale_factor = 0.25
    vectors.glyph.glyph_source.glyph_source.scale = 1.033*0.5
    vectors.glyph.glyph_source.glyph_source.center = array([0.5, 0, 0])
    mlab.savefig("output/frame_vec%04d.png" % i)

print "Files saved to output/*"
print """Create the video using:
ffmpeg -i output/frame_vec%04d.png -r 15 -vcodec copy output/output.avi
ffmpeg2theora output/output.avi -o output.ogv

To produce a FLV video, use:
ffmpeg -b 3600k -i output/frame_vec%04d.png -r 15 video.flv
"""
