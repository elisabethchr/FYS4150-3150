import numpy as np
import matplotlib.pyplot as plt

J = 1
k_B = 1

def analytic(temps, nSpin):
    nTemp = len(temps)
    energy_analytic = np.zeros(nTemp)
    Cv_analytic = np.zeros(nTemp)
    M_analytic = np.zeros(nTemp)
    for i in xrange(nTemp):
        beta = 1/float(k_B*temps[i])
        beta2 = 1/float(k_B*temps[i]**2)
        Z = 3+np.cosh(8*beta*J)
        energy_analytic[i] = -8*J*np.sinh(8*beta*J)/float(Z)/float(nSpin**2)     #mean energies
        Cv_analytic[i] = beta2*((64*J**2*np.cosh(8*beta*J)/float(Z)) - (64*np.sinh(8*beta*J)**2)/float(Z**2))/float(nSpin*nSpin) #mean Cv
        M_analytic[i] = (2*np.exp(8*beta*J) + 4)/float(Z)/float(nSpin*nSpin)
    return energy_analytic, Cv_analytic, M_analytic

def PlotTempEnergy(temps, energy, nSpin, filename):
    plt.plot(temps, energy, 'r.')
    plt.plot(temps, analytic(temps, nSpin)[0], 'b-')
    plt.xlabel(r'$k_BT$')
    plt.ylabel(r'E/J')
    plt.legend([r'$%d\times%d$ Numeric'%(nSpin, nSpin), r'$%d\times%d$ Analytic'%(nSpin, nSpin)], loc='best', fontsize=12)
    plt.title(r'Energy $\langle E \rangle$')
    plt.grid('on')
#    plt.savefig(filename+'T_'+str(T)+'MC_E.png')
    plt.show()

def PlotTempCv(temps, Cv, nSpin, filename):
    plt.plot(temps, Cv, 'r.')
    plt.plot(temps, analytic(temps, nSpin)[1], 'b-')
    plt.xlabel(r'$k_BT$')
    plt.ylabel(r'$C_V/Jk_B$')
    plt.legend([r'$%d\times%d$ Numeric'%(nSpin, nSpin), r'$%d\times%d$ Analytic'%(nSpin, nSpin)], loc='best', fontsize=12)
    plt.title(r'Heat capacity $C_V$')
    #plt.axis('equal')
    plt.grid('on')
#    plt.savefig(filename+'T_'+str(T)+'MC_Cv.png')
    plt.show()

def PlotTempMag(temps, Mag, nSpin, filename):
    plt.plot(temps, Mag, 'r.')
    plt.plot(temps, analytic(temps, nSpin)[2], 'b-')
    plt.xlabel(r'$k_BT$')
    plt.ylabel(r'$\langle |M| \rangle$')
    plt.legend([r'$%d\times%d$ Numeric'%(nSpin, nSpin), r'$%d\times%d$ Analytic'%(nSpin, nSpin)], loc='best', fontsize=12)
    plt.title(r'Magnetic Moment $\langle |M| \rangle$')
#    plt.savefig(filename+'T_'+str(T)+'MC_Mag.png')
    #plt.axis('equal')
    plt.grid('on')
    plt.show()

def PlotMCAcc_configs(MC, acc, filename):
    plt.plot(MC, acc, 'g.')
    plt.xlabel('Monte Carlo cycles')
    plt.ylabel(r'Accepted states')
    plt.title(r'Accepted states')
    #plt.savefig(filename+'MC_acceptedstates.png')
    #plt.axis('equal')
    plt.grid('on')
    plt.show()





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
