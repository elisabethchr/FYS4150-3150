import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

filename = sys.argv[1]

def Loaddata(filename):
    pos = np.loadtxt(filename)#, skiprows=2)

    plt.figure("2D")
    plt.plot(0,0,'yo')
    plt.plot(pos[:,0], pos[:,1])
    plt.title('Orbit of Earth in 2D')
    plt.xlabel("x", size=13)
    plt.ylabel("y", size=13)

    fig = plt.figure("3D")
    ax = fig.gca(projection='3d')
    ax.plot(pos[:,0], pos[:,1], pos[:,2])
    ax.scatter(0,0,0, c='y')
    plt.title('Orbit of Earth in 3D')
    ax.set_xlabel("x", size=13)
    ax.set_ylabel("y", size=13)
    ax.set_zlabel("z", size=13)
    plt.show()


sdklf = Loaddata(filename)
