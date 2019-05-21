import numpy as np
import matplotlib.pyplot as plt
import csv

path = 'C:\Users\Elisabeth\Documents\FYS3150\Project_1\build-Project_1-Desktop_Qt_5_11_1_MinGW_32bit-Debug\debug'

files = ['output1.txt', 'output2.txt', 'output3.txt', 'output4.txt', 'output5.txt']

min_RelError = np.zeros(len(files))
h = np.zeros(len(files))

for i in xrange(len(files)):
    a = np.loadtxt(files[i], skiprows=1, usecols=[3, -1])
    min_RelError[i] = min(a[:, 0])
    h[i] = a[0][1]

b = np.loadtxt(files[1], skiprows=1, usecols=[0, 1, 2])
print b

x = b[:, 0]
approx = b[:, 1]
exact = b[:, 2]

plt.plot(x, approx)
plt.plot(x, exact)
plt.show()

print "RelError: ", min_RelError
print "h: ", h

plt.semilogx(h, min_RelError, 'r-')
plt.show()
