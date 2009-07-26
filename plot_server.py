print "Initializing..."
import cPickle
from SimpleXMLRPCServer import SimpleXMLRPCServer
from threading import Thread, Lock
from time import sleep

import wx
from numpy import zeros
from enthought.mayavi import mlab
from enthought.mayavi.tools import tools
from enthought.mayavi.tools.helper_functions import Mesh, document_pipeline
from enthought.mayavi.tools.sources import MTriangularMeshSource
from enthought.traits.api import Callable
from enthought.pyface.api import GUI

from cPickle import dump

data_source = None

def triangular_mesh_source(x, y, z, triangles, **kwargs):
    #x, y, z, triangles = convert_to_arrays((x, y, z, triangles))

    scalars = kwargs.pop('scalars', None)
    if scalars is None:
        scalars = z

    global data_source
    data_source = MTriangularMeshSource()
    data_source.reset(x=x, y=y, z=z, triangles=triangles, scalars=scalars)

    name = kwargs.pop('name', 'TriangularMeshSource')
    ds = tools.add_dataset(data_source.dataset, name, **kwargs)
    data_source.m_data = ds
    return ds

class TriangularMesh(Mesh):
    _source_function = Callable(triangular_mesh_source)

triangular_mesh = document_pipeline(TriangularMesh())



gui_lock = Lock()

iter = 0
s = None

#f = open("log", "w")

def plot(vert, triangles):
    """
    Makes a mayavi plot.

    Example:

    >>> plot(vert, triangles)

    """
    print "plotting using mayavi..."
    print "unpickling...."
    vert = cPickle.loads(vert)
    triangles = cPickle.loads(triangles)
    print "  done."
    print "converting data..."
    x = vert[:, 0]
    y = vert[:, 1]
    z = zeros(len(y))
    t = vert[:, 2]
    print "  done."
    dump(x, f)
    dump(y, f)
    dump(z, f)
    dump(triangles, f)
    dump(t, f)
    def doit():
        global s
        global iter
        if iter == 0:
            print "plotting the triangular mesh..."
            s = triangular_mesh(x, y, z, triangles, scalars=t)
            #s.mlab_source.reset(x=x, y=y, z=z, triangles=triangles, scalars=t)
            print "  done."
            print "adjusting view..."
            mlab.view(0, 0)
            print "  done."
        else:
            print "changing the source..."
            # produces some messy result, I don't know why:
            #import IPython
            #IPython.ipapi.set_trace()

            # XXX: this segfaults:
            s.mlab_source.reset(x=x, y=y, z=z, triangles=triangles, scalars=t)
            # NOTES:
            # data_source above is equal to s.mlab_source
            # so let's call triangular_mesh again:
            #import IPython
            #IPython.ipapi.set_trace()
            # delete the old plot:
            #scene = #mlab.get_engine().scenes[0]
            #scene = mlab.gcf().scene
            #scene.disable_render = True
            #s = mlab.triangular_mesh(x, y, z, triangles, scalars=t,
            #        representation="surface")
            #mlab.get_engine().scenes[0].children[:1] = []
            #scene.disable_render = False
            print "  done."
        iter += 1
    gui_lock.acquire()
    doit()
    gui_lock.release()

class Server(Thread):

    def run(self):
        server = SimpleXMLRPCServer(("localhost", 8000), allow_none=True)
        server.register_introspection_functions()
        server.register_function(plot)
        print "Listening on port 8000..."
        server.serve_forever()

server = Server()
server.start()
app = wx.GetApp()
assert app is not None
assert wx.Thread_IsMain()
evtloop = wx.EventLoop()
ea = wx.EventLoopActivator(evtloop)
while 1:
    gui_lock.acquire()
    while evtloop.Pending():
        evtloop.Dispatch()
        app.ProcessIdle()
        gui_lock.release()
        gui_lock.acquire()
    gui_lock.release()
    sleep(0.001)
