import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

filename = sys.argv[1]
#filename = "C:\Users\Elisabeth\Documents\FYS4150-3150\Project_3\Earth_Sun_system\Euler100yr.txt"
#print(filename)
def Loaddata(filename):
    print filename
    pos = np.loadtxt(filename)#, skiprows=2)
    nyear = len(pos[:,0])/(365.25)+1
    name = filename.split('.')[0]
    print name
    plt.figure("2D")
    plt.plot(0,0,'yo', label='Sun')
    plt.plot(pos[:,0], pos[:,1], '-b', label='Earth')
    plt.title('Orbit of Earth in 2D')
    plt.xlabel("x", size=13)
    plt.ylabel("y", size=13)
    plt.legend(loc=1, fontsize=13)
    plt.axis('equal')
    plt.savefig('%s 2D_%d yr.png'%(name,nyear))
    plt.show()

    fig = plt.figure("3D")
    ax = fig.gca(projection='3d')
    ax.plot(pos[:,0], pos[:,1], pos[:,2], '-b', label='Earth')
    ax.scatter(0,0,0, c='y', label='Sun')
    plt.title('Orbit of Earth in 3D')
    ax.set_xlabel("x", size=13)
    ax.set_ylabel("y", size=13)
    ax.set_zlabel("z", size=13)
    ax.legend(loc=1, fontsize=13)
    ax.set_zlim(-1,1)
    plt.savefig('%s 3D_%d yr.png'%(name,nyear))



Loaddata(filename)
