import numpy as np
import matplotlib.pyplot as plt
import sys, glob
#from analytic import Printing
"""
# Load and plot data
filenames = []
for args in range(1,len(sys.argv)):
    filenames.append(sys.argv[args])

filenames = glob.glob('*.txt')
#print
#print filenames
filename1 = 'Ising_nSpin_2_nTemp_1_MC_1000000_.txt'
tot = int(filename1.split('_')[6])
def LoadAndPlot(file):
    data = np.loadtxt(file)
    print data[:, 0]
    #print a

    T = data[:,0]
    E = data[:,1]
    Cv = data[:,2]
    M = data[:,3]
    chi = data[:,4]
    Mabs = data[:,5]


    plt.figure('Suesptibilty')
    plt.plot(E)


#LoadAndPlot(filenames[0])
#LoadAndPlot(filenames[1])
#LoadAndPlot(filenames[2])
#LoadAndPlot(filenames[3])
#LoadAndPlot(filenames[4])
#LoadAndPlot(filenames[5])
LoadAndPlot(filename1)
plt.show()
"""



def PlotMCEnergy(MCs, energy, i, numMc):
    plt.plot(MCs, energy, 'r-')
    plt.xlabel('Monte Carlo cycles')
    plt.ylabel(r'Energy$\cdot$ J')
    plt.legend('T=%3f' %i, loc='best', fontsize=12)
    plt.title(r'$\langle E \rangle$, T = %3f, MC = %d' %i %numMc)
    plt.grid('on')
    plt.savefig()
    plt.show()

def PlotCvTemp(Temp, Cv):
    plt.plot(Temp, Cv, 'b-')
    plt.xlabel('Temperatures')
    plt.ylabel(r'Heat Capacity $C_V$')
    plt.title(r'Heat capacity $C_V$')
    #plt.axis('equal')
    plt.grid('on')
    plt.show()

def PlotMagTemp(Temp, Mag):
    plt.plot(Temp, Mag, 'b-')
    plt.xlabel('Temperatures')
    plt.ylabel(r'Magnetic Moment $\langle M \rangle$')
    plt.title(r'Magnetic Moment $\langle M \rangle$')
    #plt.axis('equal')
    plt.grid('on')
    plt.show()

#filename = 'Ising_nSpin_2_nTemp_1_MC_1000000_.txt'
filename = sys.argv[1]
data = np.genfromtxt(filename)
NumberOfValues = data.shape[0]
NumberOfMc = int(filename.split('_')[6])
NumberOfTemp = int(filename.split('_')[4])


Temps = data[:, 0]
Temperatures = np.linspace(Temps[0], Temps[-1], NumberOfTemp)
print len(Temperatures)
Mcs = data[:, -1]
Energy = data[:, 1]*(-1)
Specific_heat = data[:, 2]
MagneticMoment = data[:, 3]
Susceptibility = data[:, 4]
abs_MagneticMoment = data[:, 5]

for l in xrange(1, NumberOfTemp+1):
    PlotMCEnergy(Mcs[NumberOfMc*(l-1):NumberOfMc*(l)], Energy[NumberOfMc*(l-1):NumberOfMc*(l)], Temps[(NumberOfMc+1)*(l-1)], NumberOfMc)

#for every total number of Monte Carlo cycles for each temperature, extract the final Cv value and insert in Cv1
Cv = np.zeros(NumberOfTemp)
M = np.zeros(NumberOfTemp)

for i in xrange(NumberOfTemp):
    Cv[i] = Specific_heat[i*NumberOfMc]
    M[i] = MagneticMoment[i*NumberOfMc]
#    print Cv[i*NumberOfMc]

print Cv
#PlotCvTemp(Temperatures, Cv)
PlotMagTemp(Temperatures, M)

#
