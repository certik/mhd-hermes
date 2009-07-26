from cPickle import load
from enthought.mayavi import mlab

f = open("log")
x = load(f)
y = load(f)
z = load(f)
triangles = load(f)
t = load(f)
print len(x), len(y), len(z), len(t)
#print x
#print y
#print z
#print triangles
#print t
s = mlab.triangular_mesh(x, y, z, triangles, scalars=t)
x = load(f)
y = load(f)
z = load(f)
triangles2 = load(f)
t = load(f)
#print x
#print y
#print z
print triangles
print triangles2
print len(triangles)
print len(triangles2)
print len(x), len(y), len(z), len(t)
# XXX: this works, but 2643 segfaults
triangles = triangles2[:2642]
#bad = triangles2[2642]
#bad[0] = 3070
from numpy import array
bad = array([ 500,  504,  508])
print bad
triangles = triangles + bad
s.mlab_source.reset(x=x, y=y, z=z, triangles=triangles, scalars=t)
mlab.show()
