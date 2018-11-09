import numpy as np
import matplotlib.pyplot as plt
import sys, glob
from analytic import Printing

# Load and plot data
filenames = []
for args in range(1,len(sys.argv)):
    filenames.append(sys.argv[args])

filenames = glob.glob('*.txt')
print
print filenames, len(filenames)

def LoadB(file):
    a = file.split('_')[-1]
    N = int(a.split('.')[0])

    data = np.loadtxt(file)
    MC = data[:,0]
    T = data[:,1]
    E = data[:,2]
    Cv = data[:,3]
    M = data[:,4]
    chi = data[:,5]
    Mabs = data[:,6]
    return data

def LoadC(file):
    a = file.split('_')[-1]
    N = int(a.split('.')[0])

    data = np.loadtxt(file)

    return data


def PlotB1Temp(file, Temp):
    data = LoadB(file)
    outname = file.split('.')[0]
    outname = outname.split('/')[1]
    names = ['temperature', 'Energy', 'Heat capacity', 'Magnetic Moment', \
                'Suseptibility', 'Abs. Magnetization','Monte Carlo cycles']
    units = ['J', 'J/K', ' ', ' ', ' ']

    for ind, name in enumerate(names):
        if ind >= 1 and ind < len(names):
            plt.figure('%s at T=%g'%(name, Temp))
            plt.plot(data[:,-1], data[:,ind], '-b')
            plt.xlabel('Monte Carlo cycles')
            plt.ylabel(r'%s [%s]'%(name, units[ind-2]))
            plt.tight_layout()
            plt.savefig('Plots/Plotsb/%s_%s.png'%(name,outname))

def PlotC1Temp(file, Temp):
    data = LoadC(file)
    outname = file.split('.')[0]
    #outname = outname.split('/')[1]
    names = ['Monte Carlo cycles', 'temperature', 'Energy', 'Heat capacity', 'Magnetic Moment', \
                'Suseptibility', 'Abs. Magnetization', 'Numer of accepted runs']
    units = [' ', 'K', 'J', 'J/K', ' ', ' ', ' ', ' ']

    for ind, name in enumerate(names):
        if ind > 1:

            plt.figure('%s at T=%g'%(name, Temp))
            plt.plot(data[:,0], data[:,ind], '-b')
            plt.xlabel('Monte Carlo cycles')
            plt.ylabel(r'%s [%s]'%(name, units[ind]))
            plt.tight_layout()
            #plt.savefig('Plots/Plotsc/%s_%s.png'%(name,outname))

    plt.figure('Number of accepted runs against temperature')
    plt.plot(data[:,0], data[:,-1], 'g.')
    plt.xlabel(r'Monte Carlo cycles')
    plt.ylabel(r'%s'%(names[-1]))
    plt.tight_layout()
    #plt.savefig('Plots/Plotsc/NumAcc_ofT_%s.png'%(outname))

def PlotAccTemp(file):
    data = LoadC(file)

    plt.figure('Number of accepted runs against temperature')
    plt.plot(data[:,1], data[:,-1], '-b')
    plt.xlabel(r'Temperature, $T$ [K]')
    plt.ylabel(r'Number of accepted runs')
    plt.tight_layout()
    #plt.savefig('Plots/Plotsc/NumAcc_ofT_%s.png'%(outname))


def Probability(filename):
    data = np.loadtxt(filename)
    MC = data[:,0]
    T = data[:,1]
    E = data[:,2]
    Cv = data[:,3]
    M = data[:,4]
    chi = data[:,5]
    Mabs = data[:,6]
    Accepted_runs = data[:,7]

    Tname = np.mean(T)*10
    varE = Cv*T**2

    plt.figure('Probability histogram')
    plt.hist(E, bins=50, facecolor='r')
    plt.plot()
    plt.title('Probability distribution of Energy at T=%g'%(T[-1]))
    plt.xlabel('Energies')
    plt.ylabel('Probability')
    plt.grid('on')
    #plt.savefig('Plots/Plotsd/Probability_%s.png'%Tname)


# Ex P4b:
#PlotB1Temp('Txt_files/4b2x2_10000000.txt', 1.0)

# Ex P4c:
#PlotC1Temp('Txt_files/test4c_T1_1000000.txt', 1.0)
#PlotC1Temp('Txt_files/test4c_randomT1_1000000.txt', 1.0)
#PlotC1Temp('Txt_files/p4c20x20_T24_1000000.txt', 2.4)
#PlotC1Temp('Txt_files/p4c20x20_randomT24_1000000.txt', 2.4)
#PlotC1Temp('p4c20x20_randomT24_1000000.txt', 1)

PlotAccTemp('p4c20x20_diffT_100000.txt')
# Ex P4d:

#Probability('Txt_files/test4c_T24_1000000.txt')
plt.show()
