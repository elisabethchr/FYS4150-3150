import numpy as np
import matplotlib.pyplot as plt
import sys, glob
import PlotFunctions as Plot

#filename = 'Ising_nSpin_2_nTemp_1_MC_1000000_.txt'
filename = sys.argv[1]
#filename = glob.glob('*.txt')
data = np.genfromtxt(filename)
NumberOfValues = data.shape[0]
NumberOfSpins = int(filename.split('_')[2])
NumberOfTemp = int(filename.split('_')[4])
NumberOfMC = int(int(filename.split('_')[6]))#/10)

if NumberOfMC >= 10000:
    steps = int(0.0001*NumberOfMC)
    NumberOfSteps = 9000
else:
    steps = 1
    NumberOfSteps = NumberOfMC

NumberOfMc = int(NumberOfMC/steps)
#NumberOfMc = int(NumberOfValues/steps)
#print NumberOfMc
#NumberOfTemp = int(filename.split('_')[4])


MCs = data[:, 0]
Temps = data[:, 1]
#print Temps[NumberOfMc*1-NumberOfMc]
#Temperatures = np.arange(Temps[0], Temps[-1], temp_step)
Temperatures = np.arange(Temps[0], Temps[-1], 0.2)
Energy = data[:, 2]
Specific_heat = data[:, 3]
MagneticMoment = data[:, 4]
Susceptibility = data[:, 5]
abs_MagneticMoment = data[:, 6]
if data.shape[1] >= 7:
    acc_configs = data[:, -1]
    #for l in xrange(1, NumberOfTemp+1):
        #Plot.PlotMCAcc_configs(MCs[NumberOfMc*(l-1):NumberOfMc*(l)], acc_configs[NumberOfMc*(l-1):NumberOfMc*(l)], filename)
        #print NumberOfMc*(l-1), NumberOfMc*(l)
"""
for l in xrange(1, NumberOfTemp+1):
    Plot.PlotMCEnergy(MCs[NumberOfMc*(l-1):NumberOfMc*(l)], Energy[NumberOfMc*(l-1):NumberOfMc*(l)], Temps[(NumberOfMc+1)*(l-1)], NumberOfMc, filename)
    Plot.PlotMCCv(MCs[NumberOfMc*(l-1):NumberOfMc*(l)], Specific_heat[NumberOfMc*(l-1):NumberOfMc*(l)], Temps[(NumberOfMc+1)*(l-1)], filename)
    Plot.PlotMCMag(MCs[NumberOfMc*(l-1):NumberOfMc*(l)], abs_MagneticMoment[NumberOfMc*(l-1):NumberOfMc*(l)], Temps[(NumberOfMc+1)*(l-1)], filename)
"""
"""
for l in xrange(10, NumberOfTemp+1):

    print "MCs[%.3f] = %.3f" %(NumberOfSteps*(l-1), MCs[NumberOfSteps*(l-1)])
    print "T = %3f" %Temps[NumberOfSteps*(l-1)]
    print "MCs[%.3f] = %.3f" %(NumberOfSteps*(l), MCs[NumberOfSteps*(l)-1])
    print "T = %3f" %Temps[NumberOfSteps*(l)-1]0

    Plot.PlotMCEnergy(MCs[NumberOfSteps*(l-1):NumberOfSteps*(l)-1], Energy[NumberOfSteps*(l-1):NumberOfSteps*(l)-1], Temps[(NumberOfSteps+1)*(l-1)], NumberOfMc, filename)
    Plot.PlotMCCv(MCs[NumberOfSteps*(l-1):NumberOfSteps*(l)], Specific_heat[NumberOfSteps*(l-1):NumberOfSteps*(l)], Temps[(NumberOfMc+1)*(l-1)], filename)
    Plot.PlotMCMag(MCs[NumberOfSteps*(l-1):NumberOfSteps*(l)], abs_MagneticMoment[NumberOfSteps*(l-1):NumberOfSteps*(l)], Temps[(NumberOfMc+1)*(l-1)], filename)
"""
"""
for l in xrange(1, NumberOfTemp+1):
    Plot.PlotMCEnergy(MCs[NumberOfMc*(l-1):NumberOfMc*(l)], Energy[NumberOfMc*(l-1):NumberOfMc*(l)], Temps[(NumberOfMc+1)*(l-1)], NumberOfMc, filename)
    Plot.PlotMCCv(MCs[NumberOfMc*(l-1):NumberOfMc*(l)], Specific_heat[NumberOfMc*(l-1):NumberOfMc*(l)], Temps[(NumberOfMc+1)*(l-1)], filename)
    Plot.PlotMCMag(MCs[NumberOfMc*(l-1):NumberOfMc*(l)], abs_MagneticMoment[NumberOfMc*(l-1):NumberOfMc*(l)], Temps[(NumberOfMc+1)*(l-1)], filename)
"""

#for every total number of Monte Carlo cycles for each temperature, extract the final Cv value and insert in Cv1
Cv_mean = np.zeros(NumberOfTemp)
M_mean = np.zeros(NumberOfTemp)
accepted = np.zeros(NumberOfTemp)
E_mean = np.zeros(NumberOfTemp)

for i in xrange(NumberOfTemp):
    Cv_mean[i] = sum(Specific_heat[i*NumberOfSteps:(i+1)*NumberOfSteps])/NumberOfSteps
    M_mean[i] = sum(abs_MagneticMoment[i*NumberOfSteps:(i+1)*NumberOfSteps])/NumberOfSteps
    E_mean[i] = sum(Energy[i*NumberOfSteps:(i+1)*NumberOfSteps])/NumberOfSteps
#    accepted[]
#    print Cv[i*NumberOfMc]
#PlotCvTemp(Temperatures, Cv)
Plot.PlotTempEnergy(Temperatures, E_mean, NumberOfSpins, filename)
Plot.PlotTempCv(Temperatures, Cv_mean, NumberOfSpins, filename)
Plot.PlotTempMag(Temperatures, M_mean, NumberOfSpins, filename)
#
