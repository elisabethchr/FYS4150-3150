import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys, glob


def Load(file):
    data = np.loadtxt(file)
    name = '%s'%data[-1, 0]
    return name, data

def PlotEnergy(file):
    name, data = Load(file)
    name = file.split('.')[0]
    print name
    names = ['T ','E(T)', 'M(T)', 'Cv(T)', '\chi(T)']
    labels = ['Temperature','Energy', 'Magnetization', 'Heat_capacity','Susceptibility']
    units = ['K', 'J', ' ', 'J/K', ' ']
    dT = (data[1,1] - data[0,1])
    #print dT
    for i in range(1,len(data[1,:])):
        plt.figure('%s'%names[i])
        plt.plot(data[:,0], data[:,i], '-b', label=r'$%s$'%names[i])
        plt.title(r'$%s$ with $\Delta T=%.3f$'%(names[i], dT))
        plt.xlabel(r'Temperature, $T$ [K]', size=14)
        plt.ylabel(r'%s, $%s$  [%s]'%(labels[i], names[i], units[i]), size=14)
        plt.grid('on')
        plt.tight_layout()
        #plt.savefig('Plots/Plote/%s_%s.png'%(labels[i], name))


def Plot_TC(files):
    names = ['T ','<E(T)>', '<|M(T)|>', 'Cv(T)', '\chi(T)']
    labels = ['Temperature','Energy', 'Magnetization', 'Heat_capacity','Susceptibility']
    units = ['K', 'J', ' ', 'J/K', ' ']
    print files
    b = files[0].split('_')[0]
    #b = b.split('l')[-1]
    for i in range(1,5):
        plt.figure('%s'%labels[i])
        for file in files:
            data = np.loadtxt(file)
            a = file.split('.')[0]
            name = a.split('_')[-1]
            plt.plot(data[:,0], data[:,i], label='L=%s'%name)

        plt.xlabel(r'Temperature, $k_B T$ [J]', size=14)
        plt.ylabel(r'%s, $%s$ [%s/spin]'%(labels[i], names[i], units[i]), size=14)
        plt.legend(loc=1, fontsize=12)
        plt.grid('on')
        plt.tight_layout()
        plt.savefig('Plots/Plote/comparing_%s_%s.png'%(labels[i], b))


def findTC(files):
    Tc = np.zeros((len(files), 3))
    L = np.zeros((len(files)))
    T_crit = np.zeros(len(files))
    # j=0-> M, j=1-> Cv, j=2-> X
    plt.figure('Critical temperature')
    for ind, file in enumerate(files):
        data = np.loadtxt(file)
        l = file.split('.')[0]

        L[ind] = float(l.split('_')[-1])

        # Tc for M:
        i_m = []
        for j in range(len(data[:,0])):
            if data[j,2] < 0.1:
                i_m.append(j)
        im = i_m[0]
        Tc[ind, 0] = data[im, 0]

        # Tc for Cv:
        i_cv = np.where(data[:,3]==np.max(data[:,3]))
        Tc[ind, 1] = data[i_cv, 0]
        # Tc for chi:
        i_x = np.where(data[:,4]==np.max(data[:,4]))
        Tc[ind, 2] = data[i_x, 0]

        T_crit[ind] = np.sum(Tc[ind,:])/(float(len(Tc[ind,:])))
        plt.plot(L[ind], T_crit[ind], 'ok')

    print Tc, np.mean(Tc)
    print L
    print T_crit
    a = np.zeros(len(L)-1)
    for i in range(1, len(L)):
        for j in range((len(L)-1)):
            if i > j:
                #print i, j
                a[i-1] = (T_crit[i] - T_crit[j])/(L[i] - L[j])
    print sum(a)/len(a), np.mean(a)
    Tc_infty = T_crit - np.mean(a)/L
    print Tc_infty
    print np.mean(Tc_infty)

    plt.axhline(y=2.0/(np.log(1+np.sqrt(2))), linestyle='-', color='r', label=r'Analytic $T_c$')
    plt.axhline(y=np.mean(Tc_infty), linestyle='-', color='b', label=r'Numeric $T_c$')
    plt.xlabel(r'Lattice size, $L$ [spin]', size=14)
    plt.ylabel(r'Temperature, $T$ [k$_B$T/J]', size=14)
    plt.legend(loc=2, fontsize=12)
    plt.grid('on')
    plt.savefig('Plots/Plote/CriticalTemperature.png')


    return 0

files = ['run_100000_40.txt', 'run_100000_60.txt', 'run_100000_80.txt', 'run_100000_100.txt']
#EnergyTC(files)
#Plot_TC(files)
findTC(files)


#PlotEnergy('Txt_files/run1e4_20.txt')
#PlotEnergy('Txt_files/run1e4_40.txt')
#PlotEnergy('Txt_files/run1e4_60.txt')
#PlotEnergy('Txt_files/run1e4_80.txt')
#PlotEnergy('Txt_files/run1e4_100.txt')
#PlotEnergy('Txt_files/run1e5_20.txt')
#PlotEnergy('t10000_20.txt')
plt.show()
