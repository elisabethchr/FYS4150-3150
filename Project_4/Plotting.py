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
print filenames
def LoadAndPlot(file):
    data = np.loadtxt(file)
    T = data[:,0]
    E = data[:,1]
    Cv = data[:,2]
    M = data[:,3]
    chi = data[:,4]
    Mabs = data[:,5]


    plt.figure('Suesptibilty')
    plt.plot(chi)


#LoadAndPlot(filenames[0])
#LoadAndPlot(filenames[1])
#LoadAndPlot(filenames[2])
#LoadAndPlot(filenames[3])
#LoadAndPlot(filenames[4])
#LoadAndPlot(filenames[5])
LoadAndPlot('test_1000.txt')
plt.show()
