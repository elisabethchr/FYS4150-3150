import numpy as np
import matplotlib.pyplot as plt
import sys, glob
import PlotFunctions as Plot
from scipy import stats
#filename = 'Ising_nSpin_2_nTemp_1_MC_1000000_.txt'
plt.figure()

Tc = []
L_inv = []
# Gather files containing different data results in list
filenames = []
for args in range(1,len(sys.argv)):
    filenames.append(sys.argv[args])

# Loop over each filename and obtain values
for file in xrange(len(filenames)):
    filename = filenames[file]
    data = np.genfromtxt(filename)
    NumberOfValues = data.shape[0]
    NumberOfSpins = int(filename.split('_')[2])
    NumberOfTemp = int(filename.split('_')[4])
    NumberOfMC = int(int(filename.split('_')[6]))#/10)
    spins = []
    spins.append(NumberOfSpins)

    # Decide number of values per temperature depending on which file size we're looking at
    if 100000 >= NumberOfMC > 10000:
        steps = int(0.0001*NumberOfMC)
        print steps
        if NumberOfSpins > 20:
            NumberOfSteps = 900     #If excluding first 10% of data values
        else:
            NumberOfSteps = 10000
    elif NumberOfMC > 100000:
        if NumberOfSpins >= 20:
            NumberOfSteps = 9000    #If excluding first 10% of data values
        else:
            NumberOfSteps = 10000   #If including first 10% of data values
    else:
        steps = 1
        NumberOfSteps = 900
        #NumberOfSteps = 900
    #NumberOfSteps=10000
    NumberOfSteps = 1000

    print "Number of steps = ", NumberOfSteps
    # Create list containing each temperature in a given file
    Temps = data[:, 1]
    if NumberOfValues > NumberOfSteps:  #If several temperatures in a file
#        print Temps[NumberOfSteps*(NumberOfTemp)-1], Temps[NumberOfSteps*(NumberOfTemp-1)-1]
        temp_step = Temps[NumberOfSteps*(NumberOfTemp)-1] - Temps[NumberOfSteps*(NumberOfTemp-1)-1]
        print "temp_step = %.3f" %temp_step
    else:
        temp_step = 1   #If only one temperature in a file
    if NumberOfTemp > 2:
        Temperatures = np.arange(Temps[0], Temps[-1], temp_step)
    elif NumberOfTemp == 2:
        Temperatures = [1.0, 2.4]

    #Obtain values from each coulumn in file
    NumberOfTemp = len(Temperatures)
    MCs = data[:, 0]
    Energy = data[:, 2]
    Specific_heat = data[:, 3]
    MagneticMoment = data[:, 4]
    Susceptibility = data[:, 5]
    abs_MagneticMoment = data[:, 6]
    varE = np.zeros(NumberOfSteps)
    #acc_configs = data[:, -1]


    #Plot for different variables as a function Monte Carlo cycles
    for l in xrange(0, NumberOfTemp):
        #print "MCs[%.3f] = %.3f" %(NumberOfSteps*(l)+l, MCs[NumberOfSteps*(l)+l])
        #print "T = %3f" %Temps[NumberOfSteps*(l)+l]
        #print "MCs[%.3f] = %.3f" %(NumberOfSteps*(l+1)+l, MCs[NumberOfSteps*(l+1)+l])
        #print "T = %3f \n" %Temps[NumberOfSteps*(l+1)+l]
        #Plot.PlotMCVarE(MCs[NumberOfSteps*(l)+l:NumberOfSteps*(l+1)+l], Specific_heat[NumberOfSteps*(l)+l:NumberOfSteps*(l+1)+l]*Temperatures[l]**2, Temps[(NumberOfSteps+1)*(l)], filename)
        #Plot.PlotProbE(Energy[NumberOfSteps*(l)+l:NumberOfSteps*(l+1)+l], Temperatures[l], filename)
        #Plot.PlotProbE(Energy[NumberOfSteps*(l):NumberOfSteps*(l+1)-l], Temps[(NumberOfSteps+1)*(l)], filename)
        #Plot.PlotMCEnergy(MCs[NumberOfSteps*(l)+l:NumberOfSteps*(l+1)+l], Energy[NumberOfSteps*(l)+l:NumberOfSteps*(l+1)+l], Temperatures[l], NumberOfSpins, filename)
        #Plot.PlotMCCv(MCs[NumberOfSteps*(l)+l:NumberOfSteps*(l+1)+l], Specific_heat[NumberOfSteps*(l)+l:NumberOfSteps*(l+1)+l], Temps[(NumberOfSteps+1)*(l)], NumberOfSpins, filename)
        Plot.PlotMCMag(MCs[NumberOfSteps*(l)+l:NumberOfSteps*(l+1)+l], abs_MagneticMoment[NumberOfSteps*(l)+l:NumberOfSteps*(l+1)+l], Temps[(NumberOfSteps+1)*(l)], NumberOfSpins, filename)
    #for l in xrange(1, NumberOfTemp+1):
        #Plot.PlotMCAcc_configs(MCs[NumberOfMc*(l-1):NumberOfMc*(l)], acc_configs[NumberOfMc*(l-1):NumberOfMc*(l)], filename)
        #print NumberOfMc*(l-1), NumberOfMc*(l)


    # Plot for different variables as a function of temperature
    accepted_temp = np.zeros(NumberOfTemp)
    Cv_mean = np.zeros(NumberOfTemp)
    Mabs_mean = np.zeros(NumberOfTemp)
    E_mean = np.zeros(NumberOfTemp)
    Chi_mean = np.zeros(NumberOfTemp)
    varE_mean = np.zeros(NumberOfTemp)
    M_mean = np.zeros(NumberOfTemp)

    Tc_possible = []
    for i in xrange(NumberOfTemp):
        Cv_mean[i] = sum(Specific_heat[i*NumberOfSteps:(i+1)*NumberOfSteps])/NumberOfSteps
        Mabs_mean[i] = sum(abs_MagneticMoment[i*NumberOfSteps:(i+1)*NumberOfSteps])/NumberOfSteps
        E_mean[i] = sum(Energy[i*NumberOfSteps:(i+1)*NumberOfSteps])/NumberOfSteps
        Chi_mean[i] = sum(Susceptibility[i*NumberOfSteps:(i+1)*NumberOfSteps])/NumberOfSteps
        M_mean[i] = sum(MagneticMoment[i*NumberOfSteps:(i+1)*NumberOfSteps])/NumberOfSteps
        varE_mean[i] = Cv_mean[i]*Temperatures[i]**2
    #Plot.PlotTempAcc(Temperatures, accepted_ordered, NumberOfSpins, filename)
    #Plot.PlotTempAcc(Temperatures, accepted_random, NumberOfSpins, filename)
    #Plot.PlotTempEnergy(Temperatures, E_mean, NumberOfSpins, filename)
    Plot.PlotTempCv(Temperatures, Cv_mean, NumberOfSpins, filename)
    #Plot.PlotTempMag(Temperatures, Mabs_mean, NumberOfSpins, filename)
    #Plot.PlotTempMag(Temperatures, M_mean, NumberOfSpins, filename)
    #Plot.PlotTempVarE(Temperatures, varE_mean)
    #Plot.PlotTempChi(Temperatures, Chi_mean, NumberOfSpins, filename)

    # Get possible critical temperatures for each lattice if |M| < 0


    for k in xrange(NumberOfValues):
        if MagneticMoment[k] < 0:
            Tc_possible.append(Temps[k])

    # Append the lowest critical temperature for each lattice to list Tc
    if len(Tc_possible) > 0:
        Tc.append(Tc_possible[0])
        #L.append(NumberOfSpins)
        L_inv.append(1/float(NumberOfSpins**2))

    print "# Tc_possible = ", len(Tc_possible)
#print "Tc_possible = ", Tc_possible
#print "T_c =", Tc

plt.legend(loc='best', fontsize=12)
#plt.savefig('Plots/Plots_e/Temp_chi_multiple_lattices')
plt.show()
# Plot critical temperature for each lattice
#slope, intercept, r_value, p_value, std_err = stats.linregress(L_inv, Tc)
#print slope
#print "T_c = %.3f" %intercept
#a = slope*L + intercept
#print a
"""
plt.plot(L_inv, slope*np.array(L_inv) + intercept, 'r-')
plt.plot(L_inv, Tc, 'b.')
plt.title(r'Critical temperature $T_c$')
plt.legend([r'$T_C = %g$' %intercept])
plt.xlabel(r'$1/L^2$')
plt.ylabel(r'$T_C [J/k_B]$')
plt.savefig('Plots/Plots_f/Critical_temperature')
plt.show()
"""
#If NumberOfSteps=10000:
"""
for l in xrange(0, NumberOfTemp):
    print "MCs[%.3f] = %.3f" %(NumberOfSteps*(l)-l, MCs[NumberOfSteps*(l)-l])
    print "T = %3f" %Temps[NumberOfSteps*(l)-l]
    print "MCs[%.3f] = %.3f" %(NumberOfSteps*(l+1)-l, MCs[NumberOfSteps*(l+1)-l])
    print "T = %3f \n" %Temps[NumberOfSteps*(l+1)-l]

    #Plot.PlotMCVarE(MCs[NumberOfSteps*(l)+l:NumberOfSteps*(l+1)+l], Specific_heat[NumberOfSteps*(l)+l:NumberOfSteps*(l+1)+l], Temps[(NumberOfSteps+1)*(l)])
    #Plot.PlotProbE(Energy[NumberOfSteps*(l)+l:NumberOfSteps*(l+1)+l], Temperatures[l], filename)
    #Plot.PlotMCEnergy(MCs[NumberOfSteps*(l)+l:NumberOfSteps*(l+1)+l], Energy[NumberOfSteps*(l)+l:NumberOfSteps*(l+1)+l], Temperatures[l], NumberOfSpins, filename)
    #Plot.PlotMCCv(MCs[NumberOfSteps*(l):NumberOfSteps*(l+1)-1], Specific_heat[NumberOfSteps*(l):NumberOfSteps*(l+1)-1], Temps[(NumberOfSteps+1)*(l)], NumberOfSpins, filename)
    #Plot.PlotMCMag(MCs[NumberOfSteps*(l):NumberOfSteps*(l+1)-1], abs_MagneticMoment[NumberOfSteps*(l):NumberOfSteps*(l+1)-1], Temps[(NumberOfSteps+1)*(l)], NumberOfSpins, filename)
"""
