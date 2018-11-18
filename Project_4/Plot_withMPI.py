import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys, glob


def Plot_TC(files):
    names = ['T ','<E(T)>', '<|M(T)|>', 'Cv(T)', '\chi(T)']
    labels = ['Temperature','Energy', 'Magnetization', 'Heat_capacity','Susceptibility']
    units = ['K', 'J', ' ', 'J/K', ' ']
    b = files[0].split('_')[0]

    for i in range(1,5):
        plt.figure('%s'%labels[i])
        for file in files:
            data = np.loadtxt(file)
            a = file.split('.')[0]
            name = a.split('_')[-1]
            plt.plot(data[:,0], data[:,i], label='L=%s'%name)

        plt.xlabel(r'Temperature, $k_B T$ [J]', size=14)
        plt.ylabel(r'%s, $%s$ [%s/spin]'%(labels[i], names[i], units[i]), size=14)
        plt.legend(loc=2, fontsize=12)
        plt.grid('on')
        plt.tight_layout()
        plt.savefig('Plots/Plote/comparing_%s_%s.png'%(labels[i], b))


def findTC(files):
    Tc = np.zeros((len(files), 3))
    L = np.zeros((len(files)))
    T_crit = np.zeros(len(files))


    plt.figure('Critical temperature')
    for ind, file in enumerate(files): # j=0-> M, j=1-> Cv, j=2-> X
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


    a = 0.0; counter =0.0
    for i in range(len(L)):
        for j in range((len(L)-1)):
            if i > j:
                counter += 1.0
                a += (T_crit[i] - T_crit[j])/(L[i] - L[j])
                print i, j, 'a=',a

    Tc_infty = T_crit - a/(counter*L)
    print Tc_infty
    print np.mean(Tc_infty)
    analytic_Tc = 2.0/(np.log(1+np.sqrt(2)))
    plt.axhline(y=analytic_Tc, linestyle='-', color='r', label=r'Analytic $T_c$=%g'%(analytic_Tc))
    plt.axhline(y=np.mean(Tc_infty), linestyle='-', color='b', label=r'Numeric $T_c$=%g'%(np.mean(Tc_infty)))
    plt.xlabel(r'Lattice size, $L$ [spin]', size=14)
    plt.ylabel(r'Temperature, $T$ [k$_B$T/J]', size=14)
    plt.legend(loc=2, fontsize=12)
    plt.grid('on')
    plt.savefig('Plots/Plote/CriticalTemperature.png')


files = ['Txt_files/run_100000_40.txt', 'Txt_files/run_100000_60.txt', \
            'Txt_files/run_100000_80.txt', 'Txt_files/run_100000_100.txt']

Plot_TC(files)
findTC(files)

plt.show()
