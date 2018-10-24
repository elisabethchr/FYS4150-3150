import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys, glob

def PlotOrbits(filenames):
    Sun = np.loadtxt(filenames[0])
    Mercury = np.loadtxt(filenames[1])
    Venus = np.loadtxt(filenames[2])
    Earth = np.loadtxt(filenames[3])
    Mars = np.loadtxt(filenames[4])
    Jupiter = np.loadtxt(filenames[5])
    Saturn = np.loadtxt(filenames[6])
    Uranus = np.loadtxt(filenames[7])
    Neptune = np.loadtxt(filenames[8])
    Pluto = np.loadtxt(filenames[9])
    n = len(Sun[:,0])/40

    plt.figure('inner SS')
    plt.plot(Sun[2:n,0], Sun[2:n,1], '-y', label='Sun')
    plt.plot(Mercury[2:n,0], Mercury[2:n,1], '-k', label='Mercury')
    plt.plot(Venus[2:n,0], Venus[2:n,1], '-g', label='Venus')
    plt.plot(Earth[2:n,0], Earth[2:n,1], '-b', label='Earth')
    plt.plot(Mars[2:n,0], Mars[2:n,1], '-r', label='Mars')
    plt.plot(Jupiter[2:n,0], Jupiter[2:n,1], '-g', label='Jupiter')
    plt.xlabel('x [AU]', size=14)
    plt.ylabel('y [AU]', size=14)
    plt.legend(loc=1, fontsize=12)
    plt.axis('equal')
    plt.grid('on')
    plt.savefig('Plots/inner_SolarSystem.png')
    #plt.savefig('Plots/Jupiter2D_1000.png')


    plt.figure('Orbit 2d outer SS')
    plt.plot(Sun[:n,0], Sun[:n,1], '-y', label='Sun')
    plt.plot(Jupiter[2:n,0], Jupiter[2:n,1], '-g', label='Jupiter')
    plt.plot(Saturn[2:3*n,0], Saturn[2:3*n,1], '-c', label='Saturn')
    plt.plot(Uranus[2:10*n,0], Uranus[2:10*n,1], '-g', label='Uranus')
    plt.plot(Neptune[2:15*n,0], Neptune[2:15*n,1], '-b', label='Neptune')
    plt.plot(Pluto[2:-2,0], Pluto[2:-2,1], '-m', label='Pluto')

    plt.xlabel('x [AU]', size=14)
    plt.ylabel('y [AU]', size=14)
    plt.legend(loc=1, fontsize=12)
    plt.axis('equal')
    plt.grid('on')
    plt.savefig('Plots/outerSolarSystem_verlet.png')


    fig = plt.figure('Orbit 3D')
    ax = fig.gca(projection='3d')
    ax.plot(Sun[2:n,0], Sun[2:n,1], Sun[2:n,2], '-y', label='sun')
    plt.plot(Mercury[:n,0], Mercury[:n,1], Mercury[:n,2], '-k', label='Mercury')
    plt.plot(Venus[:n,0], Venus[:n,1], Venus[:n,2], '-g', label='Venus')
    plt.plot(Earth[2:n,0], Earth[2:n,1], Earth[2:n,2], '-b', label='Earth')
    plt.plot(Mars[:n,0], Mars[:n,1], Mars[:n,2], '-r', label='Mars')
    plt.plot(Jupiter[2:n,0], Jupiter[2:n,1], Jupiter[2:n,2], '-g', label='Jupiter')
    plt.plot(Saturn[2:3*n,0], Saturn[2:3*n,1], Saturn[2:3*n,2], '-c', label='Saturn')
    plt.plot(Uranus[2:10*n,0], Uranus[2:10*n,1], Uranus[2:10*n,2], '-g', label='Uranus')
    plt.plot(Neptune[2:15*n,0], Neptune[2:15*n,1], Neptune[2:15*n,2], '-b', label='Neptune')
    plt.plot(Pluto[:,0], Pluto[:,1], Pluto[:,2], '-m', label='Pluto')

    ax.set_xlabel('x [AU]', size=12)
    ax.set_ylabel('y [AU]', size=12)
    ax.set_zlabel('z [AU]', size=12)
    ax.legend(loc=2, fontsize=10)
    #ax.set_zlim(-1, 1)
    #plt.savefig('Plots/EarthJupiter3D_verlet.png')
    plt.savefig('Plots/SolarSystem3D_verlet.png')
    #plt.savefig('Plots/Jupiter3D_1000.png')


def PrecessionMercury():
    Sun = np.loadtxt('Verlet_mercury0_dt0.000003_yr20.txt')#("Verlet_mercury0_N7340000_yr20.txt")
    Mercury = np.loadtxt('Verlet_mercury1_dt0.000003_yr20.txt')#("Verlet_mercury1_N7340000_yr20.txt")
    N = len(Sun[:,0])
    r = np.zeros(N)

    print (np.linalg.norm(Mercury[0,:] - Sun[0,:]))
    print np.linalg.norm(Mercury[-1,:] - Sun[-1,:])
    r1 = Mercury[0,:] - Sun[0,:]
    print 'r1 = ', r1
    j = 0; k = 0; rend = []
    tan_tp = []; t = []
    for i in range(N):
        r[i] = np.linalg.norm(Mercury[i,:] - Sun[i,:])
        if i >= N-N/100:
            if Mercury[i,0] > 0.25:
                rend.append(np.linalg.norm(Mercury[i,:] - Sun[i,:]))

    rend = np.asarray(rend)
    tan_tp = np.asarray(tan_tp)
    t = np.asarray(t)

    index = np.argwhere(r==np.min(rend))[0][0]
    print Mercury[index,:] - Sun[index,:], index# 14569199

    #rN = Mercury[index,:] - Sun[index,:]
    xp = (Mercury[index,0] - Sun[index,0])# - (Mercury[0,0] - Sun[0,0])
    yp = (Mercury[index,1] - Sun[index,1])# - (Mercury[0,1] - Sun[0,1])
    theta_p = np.arctan(yp/xp)
    print 'theta_p = ', theta_p, 'rad'
    print 'theta_p = ', np.degrees(theta_p)*3600, 'arc sec'
    print 'theta_p = ', np.degrees(theta_p), 'degrees'

    plt.figure('Precession Mercury')
    plt.plot(Sun[:,0], Sun[:,1], 'y.', label='Sun')
    plt.plot(Mercury[:,0], Mercury[:,1], '-k', label='Mercery')
    plt.plot(Mercury[0, 0], Mercury[0, 1], '^r', label='Perihelion 1')
    plt.plot(Mercury[index,0], Mercury[index,1], '^b', label='Perihelion N')
    plt.xlabel('x [AU]', size=14)
    plt.ylabel('y [AU]', size=14)
    plt.legend(loc=1, fontsize=12)
    plt.axis('equal')
    plt.grid('on')
    plt.savefig('Plots/PrecessionMercury.png')



def PlotJupiter(files):
    posSun = np.loadtxt(files[0])
    posEarth = np.loadtxt(files[1])
    posJupiter = np.loadtxt(files[2])
    posS = np.loadtxt(files[3])
    posE = np.loadtxt(files[4])
    posJ = np.loadtxt(files[5])
    N = len(files)
    name = ['Sun', 'Earth', 'Jupiter']
    color = ['y', 'g', 'r']
    n = len(name)
    for i in range(N):
        pos = np.loadtxt(files[i])
        if i < N/2:
            nn = len(pos[:,0])/12

            fig1 = plt.figure('Verlet')
            plt.plot(pos[1000:nn,0], pos[1000:nn,1], c=color[i], label='%s'%name[i])

            plt.axis('equal')
            plt.grid('on')


        if i >= N/2:
            fig2 = plt.figure('Euler')
            plt.plot(pos[:,0], pos[:,1], c=color[i-n], label='%s'%name[i-n])
            plt.axis('equal')
            plt.grid('on')

    plt.xlabel('x [AU]', size=14)
    plt.ylabel('y [AU]', size=14)
    plt.legend(loc=1, fontsize=12)

    fig1.savefig('Plots/JupiterEarth_verlet.png')
    fig2.savefig('Plots/JupiterEarth_euler.png')


##############

filejupiter10 = ["Verlet_jupiter100_dt0.002738_yr25.txt","Verlet_jupiter101_dt0.002738_yr25.txt",\
                "Verlet_jupiter102_dt0.002738_yr25.txt"]

filejupiter1000 = ["Verlet_jupiter10000_dt0.002738_yr10.txt","Verlet_jupiter10001_dt0.002738_yr10.txt",\
                "Verlet_jupiter10002_dt0.002738_yr10.txt"]

FilesSEJ = ['Verlet_pos0.txt', 'Verlet_pos1.txt', 'Verlet_pos2.txt', 'euler_pos0.txt',\
            'euler_pos1.txt', 'euler_pos2.txt']

#"""
filenames = ["verlet_pos0_dt0.002738_yr251.txt","verlet_pos1_dt0.002738_yr251.txt",\
            "verlet_pos2_dt0.002738_yr251.txt","verlet_pos3_dt0.002738_yr251.txt",\
            "verlet_pos4_dt0.002738_yr251.txt","verlet_pos5_dt0.002738_yr251.txt",\
            "verlet_pos6_dt0.002738_yr251.txt","verlet_pos7_dt0.002738_yr251.txt",\
            "verlet_pos8_dt0.002738_yr251.txt","verlet_pos9_dt0.002738_yr251.txt"]
#"""
#PlotJupiter(FilesSEJ)
#PlotOrbits(filenames)
#jup = PlotOrbits(filejupiter1000)
PrecessionMercury()

plt.show()
