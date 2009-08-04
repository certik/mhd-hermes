import sys
import xmlrpclib
import cPickle
import subprocess
from time import sleep
from optparse import OptionParser

from numpy import zeros
import visit_writer


parser = OptionParser()
parser.add_option("--life", "-l", action="store_true", dest="life",
        default=False, help="Show life visualization")
parser.add_option("--vtk", action="store_true", dest="vtk",
        default=False, help="Save to vtk")
options, args = parser.parse_args()
if not options.life and not options.vtk:
    parser.print_help()
    sys.exit()



p = None
s = None

def start_plot_server():
    global p
    if p is None:
        p = subprocess.Popen(["python", "plot_server.py"])

def stop_plot_server():
    if p is not None:
        p.terminate()
        sleep(0.01)
        p.kill()

def plot_server_alive():
    global s
    try:
        s.alive()
    except Exception, e:
        if str(e).endswith("Connection refused"):
            return False
        else:
            raise
    return True


def establish_connection():
    global s
    if s is not None:
        return
    s = xmlrpclib.ServerProxy("http://localhost:8000/", allow_none=True)
    if not plot_server_alive():
        start_plot_server()
        print "waiting for the plot server to start up..."
        while not plot_server_alive():
            sleep(0.05)
        print "  done."

def plot_on_server(vert, triangles):
    establish_connection()
    print "plotting using mayavi..."
    v = cPickle.dumps(vert)
    t = cPickle.dumps(triangles)
    s.plot(v, t)
    print "   done."

def plot(vert, triangles):
    if options.life:
        plot_on_server(vert, triangles)
    if options.vtk:
        save_vtk(vert, triangles)

def plot_vec(vert, triangles):
    if options.life:
        pass
    else:
        save_vtk_vec(vert, triangles)

iter = 0
def save_vtk(vert, triangles):
    global iter
    node_scalars = vert[:, 2]
    node_scalars = tuple(node_scalars)
    pts = vert.copy()
    pts[:, 2] = zeros(vert.shape[0])
    pts = pts.reshape((vert.shape[0]*vert.shape[1], ))
    pts = tuple(pts)
    connectivity = []
    for t in triangles:
        connectivity.append((visit_writer.triangle, int(t[0]), int(t[1]),
            int(t[2])))
    vars = (("pressure", 1, 1, node_scalars), )
    visit_writer.WriteUnstructuredMesh("frame_scal%04d.vtk" % iter, 1, pts,
            connectivity, vars)

def save_vtk_vec(vert, triangles):
    global iter
    vector = vert.copy()
    vector[:, 0] = vector[:, 2]
    vector[:, 1] = vector[:, 3]
    vector[:, 2] = zeros(vert.shape[0])
    vector = vector[:, :3]
    vector = vector.reshape((vert.shape[0]*3, ))
    vector = tuple(vector)
    pts = vert.copy()
    pts[:, 2] = zeros(vert.shape[0])
    pts = pts[:, :3]
    pts = pts.reshape((vert.shape[0]*3, ))
    pts = tuple(pts)
    connectivity = []
    for t in triangles:
        connectivity.append((visit_writer.triangle, int(t[0]), int(t[1]),
            int(t[2])))
    vars = (("velocity", 3, 1, vector), )
    visit_writer.WriteUnstructuredMesh("frame_vec%04d.vtk" % iter, 1, pts,
            connectivity, vars)
    iter += 1
