import numpy as np
import matplotlib.pyplot as plt
import os, glob, sys

"""
Pick .txt files from directory, then put the approximation and exact values of u''(x) = f(x)
and plot the results for n=10, 100, 1000, 10000.
Plot also the relative error as function of the step length.
"""

filename = glob.glob('*.txt')
print filename
sys.exit()

#filename = ['output1.txt', 'output2.txt','output3.txt','output4.txt','output5.txt','output6.txt']
rel_error = np.zeros(len(filename))
h = np.zeros(len(filename)); # exact = np.zeros(len(filename)); approx = np.zeros(len(filename))
for ind, name in enumerate(filename):
    data = np.loadtxt(name, skiprows=1)
    no_ex = name.split('.')[0]
    rel_error[ind] = np.max(data[:,-2])
    h[ind] = data[0,-1]

print rel_error
#sys.exit()
def plot1():
    """
    Plot the between analytic solutions and numeric solutions. Both guassian elimination
    and LU decomposition.
    """
    for ind in range(len(filename)):
        plt.figure(filename[ind])
        plt.plot(data[:,0], data[:,1], '-r', label='Approximation, v(x)')
        plt.plot(data[:,0], data[:,2], '--b', label='Exact, u(x)')
        plt.title('Comparison of exact and approximation for n=%1.2e'%(ind+1) )
        plt.xlabel('x', size=14); plt.ylabel('Solution', size=14)
        plt.legend(loc='best')
        plt.savefig('Plots/%s_plot%d.png'%(no_ex,ind+1))


def plot2():
    """
    Plot of the dependence of step length in the relative error.
    """
    plt.figure('Error plot')
    plt.semilogx(h, rel_error, '-o')
    plt.xlabel('Step length, log(h)', size=14)
    plt.ylabel(r'Relative error, log($10^{\epsilon}$)', size=14)
    plt.savefig('Plots/Error_plot.png')

#plot1 = plot1()
#plot2 = plot2()
plt.show()
