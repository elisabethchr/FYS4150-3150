import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys, glob
from analytic import Printing, Z

# Load and plot data
filenames = []
for args in range(1,len(sys.argv)):
    filenames.append(sys.argv[args])

filenames = glob.glob('*.txt')


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
            print 'Mean value for %s = %g'%(names[ind], np.mean(data[:,ind]))

            plt.figure('%s at T=%g'%(name, Temp))
            plt.plot(data[:,-1], data[:,ind], '-b')
            plt.xlabel('Monte Carlo cycles')
            plt.ylabel(r'%s [%s]'%(name, units[ind-2]))
            plt.tight_layout()
            plt.savefig('Plots/Plotsb/%s_%s.png'%(name,outname))


def PlotC1Temp(file, Temp):
    data = LoadC(file)
    outname = file.split('.')[0]
    outname = outname.split('/')[1]
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
            plt.savefig('Plots/Plotsc/%s_%s.png'%(name,outname))



def PlotAcc(files):
    T1order = []; T1random = []; T24order = []; T24random = []

    plt.figure('NumAcc of MC')
    for file in files:
        if file.split('_')[3] == 'mcT1order':
            T1order.append(file)
            mc = file.split('_')[-1]
            mc = mc.split('.')[0]
            data = np.loadtxt(file)
            if len(T1order) == 1:
                plt.semilogx(mc, data[-1,-1]/float(mc), 'xb', label='T=1, spin up')
                #plt.loglog(mc, data[-1,-1]/float(mc), 'xb', label='T=1, spin up')
                #plt.semilogx(mc, data[-1,-1], 'xb', label='T=1, spin up')
            else:
                plt.semilogx(mc, data[-1,-1]/float(mc), 'xb', label='_nolegend_')
                #plt.loglog(mc, data[-1,-1]/float(mc), 'xb', label='_nolegend_')
                #plt.semilogx(mc, data[-1,-1], 'xb', label='_nolegend_')

        if file.split('_')[3] == 'mcT24order':
            T24order.append(file)
            mc = file.split('_')[-1]
            mc = mc.split('.')[0]
            data = np.loadtxt(file)
            if len(T24order) == 1:
                #plt.semilogx(mc, data[-1,-1], '.r', label='T=2.4, spin up')
                plt.semilogx(mc, data[-1,-1]/float(mc), '.r', label='T=2.4, spin up')
                #plt.loglog(mc, data[-1,-1]/float(mc), '.r', label='T=2.4, spin up')
            else:
                #plt.semilogx(mc, data[-1,-1], '.r', label='_nolegend_')
                plt.semilogx(mc, data[-1,-1]/float(mc), '.r', label='_nolegend_')
                #plt.loglog(mc, data[-1,-1]/float(mc), '.r', label='_nolegend_')

        if file.split('_')[3] == 'mcT1random':
            T1random.append(file)
            mc = file.split('_')[-1]
            mc = mc.split('.')[0]
            data = np.loadtxt(file)
            if len(T1random) == 1:
                #plt.semilogx(mc, data[-1,-1], '.g', label='T=1, spin random')
                plt.semilogx(mc, data[-1,-1]/float(mc), 'xg', label='T=1, spin random')
                #plt.loglog(mc, data[-1,-1]/float(mc), 'xg', label='T=1, spin random')
            else:
                #plt.semilogx(mc, data[-1,-1], '.g', label='_nolegend_')
                plt.semilogx(mc, data[-1,-1]/float(mc), 'xg', label='_nolegend_')
                #plt.loglog(mc, data[-1,-1]/float(mc), 'xg', label='_nolegend_')

        if file.split('_')[3] == 'mcT24random':
            T24random.append(file)
            mc = file.split('_')[-1]
            mc = mc.split('.')[0]
            data = np.loadtxt(file)
            if len(T24random) == 1:
                #plt.semilogx(mc, data[-1,-1], '.k', label='T=1, spin random')
                plt.semilogx(mc, data[-1,-1]/float(mc), '.k', label='T=2.4, spin random')
                #plt.loglog(mc, data[-1,-1]/float(mc), '.k', label='T=2.4, spin random')
            else:
                #plt.semilogx(mc, data[-1,-1], '.k')
                plt.semilogx(mc, data[-1,-1]/float(mc), '.k')
                #plt.loglog(mc, data[-1,-1]/float(mc), '.k')


    plt.xlabel('Monte Carlo cycles', size=14)
    plt.ylabel('# accepted runs/ MC cycles', size=14)
    plt.legend(loc='best', ncol=1, fontsize=12)
    plt.savefig('Plots/Plotsc/NumAccRuns_ofMC.png')

def PlotAccTemp(file):
    data = np.loadtxt(file)
    T = data[:,1]

    Index = []; N_accepted = []; Temp = []
    for i in range(1, len(T)):
        if T[i] != T[i-1] or i == len(T)-1:
            Index.append(i-1)
            N_accepted.append(data[i-1,-1])
            Temp.append(T[i-1])
        else:
            pass

    Index = np.asarray(Index); N_accepted = np.asarray(N_accepted); Temp = np.asarray(Temp)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Temp, N_accepted)

    plt.figure('Number of accepted runs against temperature')
    plt.plot(Temp, N_accepted, 'xk', label='Data points')
    plt.plot(Temp, slope*Temp + intercept, '-b', label=r'Lin. Regresion, $\sigma=%.3f$'%std_err)
    plt.xlabel(r'Temperature, $T$ [K]', size=14)
    plt.ylabel(r'Number of accepted runs', size=14)
    plt.legend(loc=4, fontsize=12)
    plt.tight_layout()
    plt.savefig('Plots/Plotsc/NumAcc_ofT.png')


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

    name = filename.split('_')[2]
    Tname = name.split('T')[1]

    varE = Cv*T**2
    print np.mean(E)
    print 'Mean variance of energy: sigma2_E =', np.mean(varE)
    print 'Mean standard deviation of energy: sigma_E =', np.sqrt(np.mean(varE))

    weights = np.ones_like(E)/float(len(E))
    plt.figure('Probability histogram')
    plt.hist(E, bins=50, facecolor='r', weights=weights)
    plt.title('Probability distribution of Energy at T=%g'%(T[-1]))
    plt.xlabel('Energies, [J]', size=14)
    plt.ylabel('Probability', size=14)
    plt.grid('on')
    plt.savefig('Plots/Plotsd/Probability_%s.png'%Tname)

    plt.figure('variance')
    plt.plot(MC, varE, '-b')
    plt.xlabel('Monte Carlo cycles', size=14)
    plt.ylabel(r'Variance, $\sigma_{E}^{2}$ [J$^2$]', size=14)
    plt.title(r'Variance of Energy, $\sigma_E^2$, at T=%g'%(T[-1]))
    plt.grid('on')
    plt.tight_layout()
    plt.savefig('Plots/Plotsd/VarianceE_T%s.png'%Tname)



filesAcc = glob.glob('Txt_files/cAcc/*.txt')

# Ex P4b:
#PlotB1Temp('Txt_files/4b2x2_10000000.txt', 1.0)

# Ex P4c:
#PlotC1Temp('Txt_files/p4c20x20_orderT1_1000000.txt', 1)
#PlotC1Temp('Txt_files/p4c20x20_randomT1_1000000.txt', 1)
#PlotC1Temp('Txt_files/p4c20x20_orderT24_1000000.txt', 2.4)
#PlotC1Temp('Txt_files/p4c20x20_randomT24_1000000.txt', 2.4)

#PlotAcc(filesAcc)
#PlotAccTemp('Txt_files/cAcc/p4c20x20_diffT_100000.txt')

# Ex P4d:
Probability('Txt_files/p4c20x20_orderT1_1000000.txt')
#Probability('Txt_files/p4c20x20_orderT24_1000000.txt')



plt.show()
