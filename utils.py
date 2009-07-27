import xmlrpclib
import cPickle
import subprocess
from time import sleep

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

def plot(vert, triangles):
    establish_connection()
    print "plotting using mayavi..."
    v = cPickle.dumps(vert)
    t = cPickle.dumps(triangles)
    s.plot(v, t)
    print "   done."

