print "Importing..."
from enthought.mayavi import mlab
from enthought.mayavi.filters.warp_scalar import WarpScalar
from enthought.mayavi.modules.surface import Surface
from numpy import array
from glob import glob
import re
print "  done."

frames = glob("anim/frame_scal*.vtk")
frames.sort()
max_frame_number = int(re.search("(\d+)", frames[-1]).groups()[0])

print "plotting the main frame..."
width = 640
height = 480
mlab.figure(size=(width, height+32))
mlab.options.offscreen = True
engine = mlab.get_engine()
vtk_file_reader = engine.open(u'anim/frame_scal0000.vtk')

warp_scalar = WarpScalar()
engine.add_filter(warp_scalar, vtk_file_reader)
surface = Surface()
engine.add_filter(surface, warp_scalar)
warp_scalar.filter.normal = array([ 0.,  0.,  1.])
warp_scalar.filter.scale_factor = 10.0

module_manager = engine.scenes[0].children[0].children[0].children[0]
module_manager.scalar_lut_manager.use_default_range = False
module_manager.scalar_lut_manager.data_range = array([-1.5, 1.5])
module_manager.scalar_lut_manager.show_scalar_bar = True

mlab.view(-122, 53, 26, array([7.5, 2.5, -0.11]))
mlab.roll(40)
print "  done."

for i in range(max_frame_number+1):
    print "doing:", i
    vtk_file_reader.timestep = i
    mlab.savefig("output/frame_scal%04d.png" % i)

print "Files saved to output/*"
print """Create the video using:
ffmpeg -i output/frame_scal%04d.png -r 15 -vcodec copy output/output.avi
ffmpeg2theora output/output.avi -o output.ogv

To produce a FLV video, use:
ffmpeg -b 3600k -i output/frame_scal%04d.png -r 15 video.flv
"""
