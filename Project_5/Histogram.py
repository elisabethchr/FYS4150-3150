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

        values, bin = np.histogram(data/100.0, bins=Nbins, normed=True)#, histtype='step', normed=True, color=color[ind],label=r'$\lambda=%g$'%(lmbd[ind]))
        print 'Sum under the graph:', sum(np.diff(bin)*values)
        bin_centers = (bin[:-1] + bin[1:])/2.0

        plt.figure('savings')
        plt.plot(data/100.0, P, c=color[ind], label=r'$\lambda=%g$'%(lmbd[ind]))
        plt.errorbar(bin_centers[0::3], values[0::3], yerr=0, fmt='.', color=color[ind])
        plt.xlabel(r'Money, $m$', size=14)
        plt.ylabel('Number of agents, normalized', size=14)
        plt.title('With savings after transactions', size=14)
        plt.legend(loc=1, fontsize=12)
        plt.xlim(0,4)
        plt.savefig('Plots/Savings.png')

        plt.figure('savings log')
        plt.loglog(data/100.0, P, c=color[ind], label=r'$\lambda=%g$'%(lmbd[ind]))
        plt.errorbar(bin_centers[0::3], values[0::3], yerr=0, fmt='.', color=color[ind])
        plt.xlabel(r'Money, $m$', size=14)
        plt.ylabel('Number of agents, normalized', size=14)
        plt.title('With savings after transactions', size=14)
        plt.legend(loc=1, fontsize=12)
        plt.xlim(0.1,4)
        plt.ylim(0.1,3)
        plt.savefig('Plots/Savings_log.png')


files = glob.glob('Txt_files/*.txt')
print files
#file = 'test_sort.txt'

#Histogram(file)
#Gibbs(file)
WithSaving(files)
plt.show()
