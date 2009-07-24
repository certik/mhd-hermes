import xmlrpclib
import cPickle

s = xmlrpclib.ServerProxy("http://localhost:8000/", allow_none=True)

def plot(vert, triangles):
    print "plotting using mayavi..."
    v = cPickle.dumps(vert)
    t = cPickle.dumps(triangles)
    s.plot(v, t)
    print "   done."
