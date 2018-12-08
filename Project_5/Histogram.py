import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import rc
from scipy import stats, special
import sys, glob

rc('font',**{'family':'serif'})
#rc('text', usetex=True)
m0 = 100.0
beta = 1.0/m0

def Histogram(file):
    data = np.loadtxt(file)
    Nbins = int(np.mean(data))

    plt.figure('p5a')
    plt.hist(data, bins=Nbins, facecolor='b')
    plt.title('Number of agents with a certain capital after transactions', size=14)
    plt.xlabel(r'Money, $m$', size=14)
    plt.ylabel('Number of agents', size=14)
    plt.savefig('Plots/Histogram.png')

def Gibbs(file):
    data = np.loadtxt(file)

    wm = beta*np.exp(-beta*data)

    plt.figure('p5b')
    plt.semilogy(data, wm, '-b')
    plt.xlabel(r'Money, $m$', size=14)
    plt.ylabel(r'$w_m$', size=14)
    plt.title('Gibbs distribution', size=14)
    plt.savefig('Plots/Gibbs_distribution.png')

    plt.figure('p5b_loglog')
    plt.loglog(data, wm, '-b')
    plt.xlabel(r'Money, $m$', size=14)
    plt.ylabel(r'$w_m$', size=14)
    plt.title('Gibbs distribution', size=14)
    plt.savefig('Plots/Gibbs_distribution_loglog.png')

def WithSaving(files):
    lmbd = [0.0, 0.25, 0.5, 0.9]
    color = ['b', 'r', 'g', 'k']
    for ind, file in enumerate(files):
        data = np.loadtxt(file)
        Nbins = int(np.mean(data))
        name = file.split('.')[0]
        name = name.split('L')[-1]

        x = data/m0
        n = 1.0 + 3.0*lmbd[ind]/(1.0-lmbd[ind])
        an = n**n / special.gamma(n)
        P = an * x**(n-1.0) * np.exp(-n*x)

        values, bin = np.histogram(data/m0, bins=Nbins, normed=True)#, histtype='step', normed=True, color=color[ind],label=r'$\lambda=%g$'%(lmbd[ind]))
        print 'Sum under the graph:', sum(np.diff(bin)*values)
        bin_centers = (bin[:-1] + bin[1:])/2.0

        plt.figure('savings')
        plt.plot(data/m0, P, c=color[ind], label=r'$\lambda=%g$'%(lmbd[ind]))
        plt.errorbar(bin_centers[0::3], values[0::3], yerr=0, fmt='.', color=color[ind])
        plt.xlabel(r'Money, $m/m_0$', size=14)
        plt.ylabel('Number of agents, normalized', size=14)
        plt.title('With savings after transactions', size=14)
        plt.legend(loc=1, fontsize=12)
        plt.xlim(0,4)
        plt.savefig('Plots/Savings.png')

        plt.figure('savings log')
        plt.loglog(data/m0, P, c=color[ind], label=r'$\lambda=%g$'%(lmbd[ind]))
        plt.errorbar(bin_centers[0::3], values[0::3], yerr=0, fmt='.', color=color[ind])
        plt.xlabel(r'Money, $m/m_0$', size=14)
        plt.ylabel('Number of agents, normalized', size=14)
        plt.title('With savings after transactions', size=14)
        plt.legend(loc=1, fontsize=12)
        plt.xlim(0.1,4)
        plt.ylim(0.1,3)
        plt.savefig('Plots/Savings_log.png')

def plot5d(files, label):
    alpha = [0.5, 1.0, 1.5, 2.0]
    color = ['b', 'r', 'g', 'm']

    f1, ax1 = plt.subplots()

    for ind in range(len(files)):

        data = np.loadtxt(files[ind])
        N = len(data)
        Nbins = int(np.mean(data))
        #wm = beta*np.exp(-beta*data)
        w_m, bins = np.histogram(data/m0, bins=N/10, normed=True)
        powlaw = data**(-1-alpha[ind])

        ax1.loglog(data[0::10], w_m/float(N), c=color[ind], linestyle='--', label=r'$\alpha=%g$'%alpha[ind])
        #plt.loglog(data, powlaw/float(N), c=color[ind], linestyle=':', label=r'$m^{-1-\alpha}, \alpha=%g$'%alpha[ind])
        ax1.set_xlabel(r'$m/m_0$', size=14)
        ax1.set_ylabel(r'$w_m/N_{agents}$', size=14)
        ax1.legend(loc=1, fontsize=12)
        ax1.set_title('Nearest neighbourgh interaction no saving, N=%d'%N, size=14)
        #ax1.set_xlim(1e-3,2e1)

    f1.savefig('Plots/Part_d_%s.png'%label)



def Plot5e(files, a):
    """
    Input: list of data files and the value of alpha
    """

    gamma = [0.0, 1.0, 2.0, 3.0, 4.0]
    color = ['b', 'r', 'g', 'm', 'c']
    lmbd = [0, 0.5]
    #alpha = [1.0, 2.0]
    NoSaving = []; j=0
    f1, ax1 = plt.subplots()
    f2, ax2 = plt.subplots()
    for ind,file in enumerate(files):
        name = file.split('_')
        for i in range(len(name)):
            if name[i] == 'L00':
                NoSaving.append(file)
                data = np.loadtxt(file)
                N = len(data)
                Nbins = int(np.mean(data))
                w_m, bins = np.histogram(data/m0, bins=N/10, normed=True)

                ax1.loglog(data[0::10], w_m, c=color[ind], linestyle='--', label=r'$\gamma=%d$'%gamma[ind])
                ax1.set_xlabel(r'$m/m_0$', size=14)
                ax1.set_ylabel(r'$w_m/N_{agents}$', size=14)
                ax1.legend(loc=1, fontsize=12)
                ax1.set_title(r'$\alpha = %g$ and no saving, $\lambda = %g$'%(a, lmbd[0]))
                ax1.grid('on')
                ax1.set_xlim(3e-3, 100)


            if name[i] == 'L05':

                data = np.loadtxt(file)
                N = len(data)
                Nbins = int(np.mean(data))
                w_m, bins = np.histogram(data/m0, bins=N/10, normed=True)

                ax2.loglog(data[0::10], w_m, c=color[j], linestyle='--', label=r'$\gamma=%d$'%gamma[j])
                ax2.set_xlabel(r'$m/m_0$', size=14)
                ax2.set_ylabel(r'$w_m/N_{agents}$', size=14)
                ax2.legend(loc=1, fontsize=12)
                ax2.set_title(r'$\alpha = %g$ and saving, $\lambda = %g$'%(a, lmbd[1]))
                ax2.grid('on')
                ax2.set_xlim(0, 200)

                j += 1

    f1.savefig('Plots/Weldth_e_a1_nosaveing.png')
    f2.savefig('Plots/Weldth_e_a1_saveing.png')



saving_files = glob.glob('Txt_files/part_c/*.txt')

#print files
#file = 'test_sort.txt'
NoSaving1000 = glob.glob('Txt_files/part_d/NoSaving/N1000/*.txt')
NoSaving500 = glob.glob('Txt_files/part_d/NoSaving/N500/*.txt')
Saving1000 = glob.glob('Txt_files/part_d/Saving/N1000/*.txt')
Saving500 = glob.glob('Txt_files/part_d/Saving/N500/*.txt')
gamma_files = glob.glob('Txt_files/part_e/*.txt')

#Histogram(file)
#Gibbs(files[0])
#WithSaving(files)
#plot5d(NoSaving1000, 'no_saving_1000')
#plot5d(NoSaving500, 'no_saving_500')
#plot5d(Saving1000, 'saving_1000')
#plot5d(Saving500, 'saving_500')

Plot5e(gamma_files, 1.0)
plt.show()
