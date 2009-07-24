import cPickle
from SimpleXMLRPCServer import SimpleXMLRPCServer
from numpy import zeros
from enthought.mayavi import mlab

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
    print "plotting the triangular mesh..."
    s = mlab.triangular_mesh(x, y, z, triangles, scalars=t)
    print "  done."
    print "adjusting view..."
    mlab.view(0, 0)
    print "  done."

server = SimpleXMLRPCServer(("localhost", 8000), allow_none=True)
server.register_introspection_functions()
server.register_function(plot)
print "Listening on port 8000..."
server.serve_forever()
