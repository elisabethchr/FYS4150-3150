import numpy as np
import matplotlib.pyplot as plt
import sys, glob

txtfiles = glob.glob('*.txt')

def PlotEnergy(files):
    i = 0
    for file in files:
        if file == 'Verlet_energies0.txt':
            Sun = np.loadtxt(file)
        if file == 'Verlet_energies1.txt':
            Earth = np.loadtxt(file)
        if file == 'Verlet_energies2.txt':
            Jupiter = np.loadtxt(file)

    t = Sun[:,0]

    Energies = np.zeros((len(t), 3))

    plt.figure('Energy')
    plt.plot(abs(Sun[:,0]), abs(Sun[:,1]), '-y', label='Sun')
    plt.plot(abs(Earth[:,0]), abs(Earth[:,1]), '-b', label='Earth')
    plt.plot(abs(Jupiter[:,0]), abs(Jupiter[:,1]), '-g', label='Jupiter')

    plt.xlabel('time [yr]', size=14)
    plt.ylabel(r'Energy [J yr$^2$/M$_{\odot} AU^2$]', size=14)
    plt.legend(loc=2, fontsize=12)
    plt.grid('on')
    plt.tight_layout()
    plt.savefig('Energies50yr.png')

def PlotOrbit(txtfiles):
    for file in txtfiles:
        if file == 'Verlet_pos0_dt0.002738_yr50.txt':
            PosSun = np.loadtxt(file)
        if file == 'Verlet_pos1_dt0.002738_yr50.txt':
            PosEarth = np.loadtxt(file)
        if file == 'Verlet_pos2_dt0.002738_yr50.txt':
            PosJupiter = np.loadtxt(file)

    plt.figure('Orbit')
    plt.plot(PosSun[:,0], PosSun[:,1], '-y', label='Sun')
    plt.plot(PosEarth[:,0], PosEarth[:,1], '-b', label='Earth')
    plt.plot(PosJupiter[:,0], PosJupiter[:,1], '-g', label='Jupiter')

    plt.xlabel('x [AU]', size=14)
    plt.ylabel('y [AU]', size=14)
    plt.legend(loc=1, fontsize=12)
    plt.axis('equal')
    plt.savefig('Verlet_JupiterEarth50yr.png')


#PlotOrbit(txtfiles)
PlotEnergy(txtfiles)
plt.show()
