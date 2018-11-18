import numpy as np
import matplotlib.pyplot as plt

J = 1
k_B = 1

def analytic(temps, nSpin):
    nTemp = len(temps)
    energy_analytic = np.zeros(nTemp)
    Cv_analytic = np.zeros(nTemp)
    Mabs_analytic = np.zeros(nTemp)
    chi_analytic = np.zeros(nTemp)
    for i in xrange(nTemp):
        beta = 1/float(k_B*temps[i])
        beta2 = 1/float(k_B*temps[i]**2)
        Z = 3+np.cosh(8*beta*J)
        energy_analytic[i] = -8*J*np.sinh(8*beta*J)/float(Z)/float(nSpin**2)     #mean energies
        Cv_analytic[i] = beta2*((64*J**2*np.cosh(8*beta*J)/float(Z)) - (64*np.sinh(8*beta*J)**2)/float(Z**2))/float(nSpin*nSpin) #mean Cv
        Mabs_analytic[i] = (2*np.exp(8*beta*J) + 4)/float(Z)/float(nSpin*nSpin)
        chi_analytic[i] = beta*((8 * np.exp(8*beta*J) + 8)/float(Z) - (4 * Mabs_analytic[i])**2) / 4
    return energy_analytic, Cv_analytic, Mabs_analytic, chi_analytic

def PlotTempEnergy(temps, energy, nSpin, filename):
    plt.plot(temps, energy, '.', label=r'$%d\times%d$'%(nSpin, nSpin))
#    plt.plot(temps, analytic(temps, nSpin)[0], 'b-')
    plt.xlabel(r'$k_BT$')
    plt.ylabel(r'E/J')
    #plt.legend([r'$%d\times%d$ Numeric'%(nSpin, nSpin), r'$%d\times%d$ Analytic'%(nSpin, nSpin)], loc='best', fontsize=12)
    plt.title(r'Energy $\langle E \rangle$')
    plt.grid('on')
    #plt.savefig(filename+'TempE.png')
    #plt.show()

listCv = []
def PlotTempCv(temps, Cv, nSpin, filename):
    plt.plot(temps, Cv, '.', label=r'$%d\times%d$'%(nSpin, nSpin))
#    plt.plot(temps, analytic(temps, nSpin)[1], 'b-')
    plt.xlabel(r'$k_BT$')
    plt.ylabel(r'$C_V/Jk_B$')
#    plt.legend([r'$%d\times%d$ Numeric'%(nSpin, nSpin)], loc='best', fontsize=12)
    plt.title(r'Heat capacity $C_V$')
    #plt.axis('equal')
    plt.grid('on')
    #plt.savefig(filename+'TempCv.png')
    #plt.show()

def PlotTempMag(temps, Mag, nSpin, filename):
    ax = plt.plot(temps, Mag, '.', label=r'$%d\times%d$'%(nSpin, nSpin))
#    plt.plot(temps, analytic(temps, nSpin)[2], 'b-')
    plt.xlabel(r'$k_BT$')
    plt.ylabel(r'$\langle |M| \rangle$')
#    plt.legend([r'$%d\times%d$ Numeric'%(nSpin, nSpin), r'$%d\times%d$ Analytic'%(nSpin, nSpin)], loc='best', fontsize=12)
    plt.title(r'$\langle |M| \rangle$')
    #plt.axis('equal')
    plt.grid('on')
    #plt.savefig('Plots/Plots_e/' + filename + 'TempMag.png')
    #plt.show()

def PlotTempChi(temps, chi, nSpin, filename):
    plt.plot(temps, chi, '.', label=r'$%d\times%d$'%(nSpin, nSpin))
#    plt.plot(temps, analytic(temps, nSpin)[3], 'b-')
    plt.xlabel(r'$k_BT$')
    plt.ylabel(r'$\chi$')
    #plt.legend([r'$%d\times%d$ Numeric'%(nSpin, nSpin), r'$%d\times%d$ Analytic'%(nSpin, nSpin)], loc='best', fontsize=12)
    plt.title(r'Susceptibility $\chi$')
    #plt.axis('equal')
    plt.grid('on')
    #plt.savefig(filename+'TempChi.png')
    #plt.show()

def PlotMCAcc_configs(MC, acc, filename):
    plt.plot(MC, acc, 'g.')
    plt.xlabel('Monte Carlo cycles')
    plt.ylabel(r'Accepted states')
    plt.title(r'Accepted states')
    #plt.savefig(filename+'MC_acceptedstates.png')
    #plt.axis('equal')
    plt.grid('on')
    plt.show()

def PlotTempAcc(T, acc, nSpin, filename):
    plt.plot(T, acc, 'r-')
    plt.xlabel('Temperatures')
    plt.ylabel('Accepted states')
#    plt.legend([r'T = %.3f' %T], loc='best', fontsize=12)
    plt.title(r'Energy $\langle E \rangle$')
    plt.grid('on')
#    plt.savefig(filename+'T_'+str(T)+'MC_E.png')
    plt.show()


def PlotMCEnergy(MC, energy, T, nSpin, filename):
    plt.plot(MC, energy, 'r-')
    plt.xlabel('Monte Carlo cycles')
    plt.ylabel(r'E/J')
    plt.legend([r'T = %.3f' %T], loc='best', fontsize=12)
    plt.title(r'Energy $\langle E \rangle$')
    plt.grid('on')
    plt.savefig('Plots/Plots_d/'+filename+'T_'+str(T)+'MC_E.png')
    plt.show()

def PlotMCCv(MC, Cv, T, nSpin, filename):
    plt.plot(MC, Cv, 'r-')
    plt.xlabel('Monte Carlo cycles')
    plt.ylabel(r'$C_V$')
    plt.legend([r'T = %.3f' %T], loc='best', fontsize=12)
    plt.title(r'Specific heat C_V')
    plt.grid('on')
    #plt.savefig(filename+'T_'+str(T)+'MC_Cv.png')
    plt.show()

def PlotMCMag(MC, mag, T, nSpin, filename):
    plt.plot(MC, mag, 'r-')
    plt.xlabel(r'$k_BT$')
    plt.ylabel(r'$\langle |M| \rangle$')
    plt.legend([r'T = %.3f' %T], loc='best', fontsize=12)
#    plt.legend([r'$%d\times%d$ Numeric'%(nSpin, nSpin), r'$%d\times%d$ Analytic'%(nSpin, nSpin)], loc='best', fontsize=12)
    plt.title(r'Magnetic Moment $\langle |M| \rangle$')
    plt.grid('on')
    #plt.savefig(filename+'T_'+str(T)+'MC_Mag.png')
    plt.show()

def PlotProbE(energy, temps, filename):
    #c = filename.split('_')[-1]
    #config=c.split('.')[0]
    plt.figure('Probability histogram')
    plt.hist(energy, bins=100, facecolor='r')
    plt.title('Probability distribution of Energy at T = %g' %(temps))
    plt.xlabel('Energies, [J]')#, size=14)
    plt.ylabel('Probability')#, size=14)
    plt.grid('on')
    plt.savefig('Plots/Plots_d/' + filename + 'T_%s'%temps +'probE.png')
    plt.show()

def PlotMCVarE(MCs, varE, T, filename):
    plt.plot(MCs, varE, 'r-', label='T=%g'%T)
    plt.xlabel('Monte Carlo')
    plt.ylabel(r'$\sigma_E^2$')
    plt.title(r'T = %g'%T)
    plt.grid('on')
    plt.savefig('Plots/Plots_d/' + filename + 'T_'+str(T)+'MC_varE.png')
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
