import numpy as np
from astroquery.jplhorizons import conf
conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'
from astroquery.jplhorizons import Horizons


def get_object_info(obj):
    obj = Horizons(id=obj, location='500@10',epochs={'start':'2018-12-31','stop':'2019-01-02','step':'1d'}, id_type='majorbody')
    table = obj.vectors()

    initial_position = np.array([table['x'][0], table['y'][0], table['z'][0]])    # AU
    initial_velocity = np.array([table['vx'][0], table['vy'][0], table['vz'][0]])*365.25 # AU/yr

    return table['targetname'][0], initial_position, initial_velocity



def writePostionFile(file, planets):

    planets.append(get_object_info('999'))

    outfile = open(file, 'w')
    for obj in range(len(planets)):
        name = get_object_info(obj)[0].split(' ')[0]
        outfile.write('%s'%(name))
        for i in range(3):
            outfile.write('  %g'%(get_object_info(obj)[1][i]))
        outfile.write('\n')



def WriteVelocityFile(file, planets):

    planets.append(get_object_info('999'))

    outfile = open(file, 'w')
    for obj in range(len(planets)):
        name = get_object_info(obj)[0].split(' ')[0]
        outfile.write('%s'%(name))
        for i in range(3):
            outfile.write('  %g'%(get_object_info(obj)[2][i]))
        outfile.write('\n')

def WriteMassFile(file, masses):
    outfile = open(file, 'w')
    m = np.zeros((len(masses), 3))
    m[:,0] = masses
    print (m)
    #outfile.write(m)
    for obj in range(len(masses)):
        outfile.write('Planet%d  '%obj)
        for j in range(3):
            outfile.write('%11g'%(m[obj, j]))
        outfile.write('\n')

# calls:
planets = ['199', '299', '399', '499', '599', '699', '799', '899', '999']
Msun = 1.99e30
planetmasses = np.array([1e-24*Msun, 0.33, 4.9, 6.0, 0.64,1900,570,87, 103, 0.0005])*1e24/Msun   # M/Msun
writePostionFile('Initialposition.txt', planets)
WriteVelocityFile('Initialvelocity.txt', planets)
WriteMassFile('masses.txt', planetmasses)
