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
        print "plotting the triangular mesh..."
        s = mlab.triangular_mesh(x, y, z, triangles, scalars=t)
        print "  done."
        print "adjusting view..."
        mlab.view(0, 0)
        print "  done."
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
    sleep(0.001)
