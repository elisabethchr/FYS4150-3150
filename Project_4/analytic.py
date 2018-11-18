import numpy as np
import matplotlib.pyplot as plt

J = 1.0
E = 8.0*J
kB = 1.0
T = 1.0
beta = 1.0/(kB*T)
spin = 16

def Z():
    return 12.0 + 4.0*(np.exp(E*beta) + np.exp(-E*beta))

def ExpEnergy():
    Energy = E*(np.exp(-E*beta) - 2*np.exp(E*beta))/Z()
    return Energy/2.0

def ExpEnergy2():
    E2 = E*E*(np.exp(-E*beta) + 2*np.exp(E*beta))/Z()
    return E2/2.0

def MagneticMoment():
    return (8*np.exp(E*beta) + 4)/(Z()*2.0)#12.0

def MagneticMoment2():
    return 8*(4*np.exp(E*beta) + 1)/(Z()*2.0)

def Cv():
    return beta*(ExpEnergy2() - ExpEnergy()*ExpEnergy())/(T)

def chi():
    return beta*(MagneticMoment2() - MagneticMoment()*MagneticMoment())

def Printing():
    E = ExpEnergy()
    EE = E*E
    E2 = ExpEnergy2()
    print
    print '=========================================================='
    print '             Analytic values for 2x2 lattice              '
    print '=========================================================='
    print 'Partition function:           Z = %g'%(Z())
    print '----------------------------------------------------------'
    print 'Energies:       <E> = %5.3g,  <E>^2 = %5.3g, <E^2> = %5.3g'%(ExpEnergy(),EE,E2)
    print 'Magnetization: <|M|> = %3.3g, <|M|>^2 = %3.3g, <M^2> = %5.3g'%(MagneticMoment(), MagneticMoment()**2, MagneticMoment2())
    print '----------------------------------------------------------'
    print 'Heat capacity and Suseptibility:'
    print 'Cv = %5.5g,           chi = %5.5g'%(Cv(), chi())
    print '=========================================================='

Printing()
