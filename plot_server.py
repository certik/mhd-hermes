print "Initializing..."
import cPickle
from SimpleXMLRPCServer import SimpleXMLRPCServer
from threading import Thread, Lock
from time import sleep

import wx
from numpy import zeros
from enthought.mayavi import mlab
from enthought.pyface.api import GUI

gui_lock = Lock()

iter = 0
s = None

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
    def doit():
        global s
        global iter
        if iter == 0:
            print "plotting the triangular mesh..."
            s = mlab.triangular_mesh(x, y, z, triangles, scalars=t)
            print "  done."
            print "adjusting view..."
            mlab.view(0, 0)
            print "  done."
        else:
            print "changing the source..."
            # This doesn't work due to a bug in mayavi/VTK:
            # http://github.com/certik/mhd-hermes/issues#issue/1
            #s.mlab_source.reset(x=x, y=y, z=z, triangles=triangles, scalars=t)
            # so until this is fixed, let's call triangular_mesh and delete the
            # old mesh (this is slow but works):
            scene = mlab.gcf().scene
            scene.disable_render = True
            s = mlab.triangular_mesh(x, y, z, triangles, scalars=t)
            mlab.get_engine().scenes[0].children[:1] = []
            scene.disable_render = False
            print "  done."
        iter += 1
    gui_lock.acquire()
    doit()
    gui_lock.release()

def alive():
    return True

class Server(Thread):

    def run(self):
        server = SimpleXMLRPCServer(("localhost", 8000), allow_none=True)
        server.register_introspection_functions()
        server.register_function(plot)
        server.register_function(alive)
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
