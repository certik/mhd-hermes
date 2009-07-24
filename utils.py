def plot(vert, triangles):
    print "plotting using mayavi..."
    from numpy import zeros
    from enthought.mayavi import mlab
    x = vert[:, 0]
    y = vert[:, 1]
    z = zeros(len(y))
    t = vert[:, 2]
    s = mlab.triangular_mesh(x, y, z, triangles, scalars=t)
    mlab.view(0, 0)
    print "   done."
