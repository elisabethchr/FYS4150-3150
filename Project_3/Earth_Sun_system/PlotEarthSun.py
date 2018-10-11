import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

filename = sys.argv[1]

def Loaddata(filename):
    pos = np.loadtxt(filename)#, skiprows=2)
    print np.shape(pos)
    plt.figure("2D")
    plt.plot(pos[:,0], pos[:,1])
    plt.xlabel("x", size=13)
    plt.ylabel("y", size=13)

    fig = plt.figure("3D")
    ax = fig.gca(projection='3d')
    ax.plot(pos[:,0], pos[:,1], pos[:,2])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.show()


sdklf = Loaddata(filename)
