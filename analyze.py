from cPickle import load
from enthought.mayavi import mlab

f = open("log")
x = load(f)
y = load(f)
z = load(f)
triangles = load(f)
t = load(f)
print x
print y
print z
print triangles
print t
s = mlab.triangular_mesh(x, y, z, triangles, scalars=t)
x = load(f)
y = load(f)
z = load(f)
triangles = load(f)
t = load(f)
print x
print y
print z
print triangles
print t
s.mlab_source.reset(x=x, y=y, z=z, triangles=triangles, scalars=t)
