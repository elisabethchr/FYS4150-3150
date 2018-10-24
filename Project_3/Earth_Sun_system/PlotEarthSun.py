import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys, glob


def Loaddata(name):

    pos = np.loadtxt(name)
    n = len(pos)
    print len(pos)
    plt.figure("2D_%s"%name)
    plt.plot(0,0, 'yo', label='Sun')

    plt.plot(pos[:,0], pos[:,1],'b-', label='Earth')
    plt.plot(pos[0,0], pos[0,1], 'r^', label='Inital position')
    plt.plot(pos[365,0], pos[365,1], 'kx', label='1 year')

    plt.title('Orbit of Earth in 2D')
    plt.xlabel("x [AU]", size=13)
    plt.ylabel("y [AU]", size=13)
    plt.legend(loc=1, fontsize=13)
    plt.axis('equal')
    plt.savefig('Plots/%s2D_.png'%(name))

    fig = plt.figure("3D_%s"%name)
    ax = fig.gca(projection='3d')
    ax.plot(pos[:,0], pos[:,1], pos[:,2],'b-', label='Earth')
    ax.scatter(0,0,0, c='y', label='Sun')
    plt.title('Orbit of Earth in 3D')
    ax.set_xlabel("x [AU]", size=13)
    ax.set_ylabel("y [AU]", size=13)
    ax.set_zlabel("z [AU]", size=13)
    ax.legend(loc=1, fontsize=13)
    ax.set_zlim(-1,1)
    plt.savefig('Plots/%s3D_.png'%(name))


def Plot_dt(Earthfiles):
    i = 0
    pos01 = np.loadtxt(Earthfiles[0])
    pos1 = np.loadtxt(Earthfiles[1])
    pos60 = np.loadtxt(Earthfiles[3])
    N = int(len(pos1[:,0])/2)
    dt = [0.1/356.23, 1/365.25, 50/365.25, 60/365.25]

    plt.figure('time step testing')

    plt.plot(pos60[:,0], pos60[:,1],'-g', label=r'dt=%1.3g $\sim$60d'%dt[3])
    plt.plot(pos1[:,0], pos1[:,1],'-r', label=r'dt=%1.3e$\sim$1d'%dt[1])
    plt.plot(pos01[:,0], pos01[:,1],'--b', label=r'dt=%1.3e$\sim$0.1d'%dt[0])
    plt.xlabel('x [AU]', size=14)
    plt.ylabel('y [AU]', size=14)
    plt.xlim(-3,8)
    plt.ylim(-10,5)
    plt.legend(loc=3, fontsize=12)
    plt.savefig('Plots/Different_dt.png')


def PlotV_escape(txtfiles):

    f5 = np.loadtxt("Vescape5.txt")
    f6 = np.loadtxt("Vescape6.txt")
    f7 = np.loadtxt("Vescape7.txt")
    f8 = np.loadtxt("Vescape8.txt")
    f9 = np.loadtxt("Vescape9.txt")
    f10 = np.loadtxt("Vescape10.txt")
    f11 = np.loadtxt("Vescape11.txt")
    f12 = np.loadtxt("Vescape12.txt")
    print len(f5[:,0])/365
    plt.figure('Different Escape Velocity')
    plt.plot(f5[:,0], f5[:,1],'-b', label='5AU/yr')
    plt.plot(f6[:,0], f6[:,1],'-r', label='6AU/yr')
    plt.plot(f7[:,0], f7[:,1],'-g', label='7AU/yr')
    plt.plot(f8[:,0], f8[:,1],'-y', label='8AU/yr')
    plt.plot(f9[:,0], f9[:,1],'-m', label='9AU/yr')
    plt.plot(f10[:,0], f10[:,1],'-c', label='10AU/yr')
    plt.plot(f11[:,0], f11[:,1],'-k', label='11AU/yr')
    plt.plot(f12[:,0], f12[:,1],'--b', label='12AU/yr')

    plt.xlabel('X [AU]', size=14)
    plt.ylabel('Y [AU]', size=14)
    plt.legend(loc=1, fontsize=12)
    plt.savefig('Plots/EscapeVelocity.png')

def Beta(txtfiles):

    i = 0
    betafiles = []

    for file in txtfiles:
        if file == 'Beta%s.txt'%i:
            betafiles.append(file)
        i+=1
    print betafiles
    pos = np.zeros((len(betafiles)))
    beta = [2,2.25,2.5,2.75,3]
    plt.figure('variating gravity')
    for i, name in enumerate(betafiles):
        pos = np.loadtxt(name)
        plt.plot(pos[:,0], pos[:,1], label=r'$\beta$=%g'%beta[i])

    plt.xlabel('x [AU]', size=14)
    plt.ylabel('y [AU]', size=14)
    plt.xlim(-4,1.2)
    plt.ylim(-2,1.5)
    plt.legend(loc=3, fontsize=12)
    plt.savefig('Plots/Beta.png')

#Eulerfiles = glob.glob('*Euler.txt')
#Verletfile = glob.glob('*Verlet.txt')
#dtfiles_Earth = glob.glob('*E_yr10.txt')        # For different dt of Earth
#dtfiles_Earth = glob.glob('*sun_10yr.txt')     # For different dt for sun
#txtfiles = glob.glob('*.txt')                  # All .txt files.

#Loaddata(Verletfile[2])                        # For verlet files
#Loaddata(Eulerfiles[2])                        # For Euler files
#PlotV_escape(txtfiles)
#Beta(txtfiles)
Plot_dt(dtfiles_Earth)
plt.show()
