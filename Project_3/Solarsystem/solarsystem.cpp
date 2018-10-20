#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "include_headers.h"
/*
#include "solarsystem.h"
#include "gravitationalforce.h"
//#include "Readfile.h"
#include "readfile_test.h"
#include "celestialobject.h"
*/

SolarSystem::SolarSystem(std::string positionfile, std::string velocityfile, std::string massfile) :
    m_kineticEnergy(0),
    m_potentialEnergy(0)
{
    Readfile_test readfile;
    arma::mat pos0 = readfile.Readfile_(positionfile);    // AU
    arma::mat vel0 = readfile.Readfile_(velocityfile);
    arma::mat masses = readfile.Readfile_(massfile);   // AU/yr
    vec3 pos0_CoM; vec3 vel0_CoM;
    vec3 pos0_Mercury; vec3 vel0_Mercury;
    vec3 pos0_Venus; vec3 vel0_Venus;
    vec3 pos0_Earth; vec3 vel0_Earth;
    vec3 pos0_Mars; vec3 vel0_Mars;
    vec3 pos0_Jupiter; vec3 vel0_Jupiter;
    vec3 pos0_Saturn; vec3 vel0_Saturn;
    vec3 pos0_Uranus; vec3 vel0_Uranus;
    vec3 pos0_Neptune; vec3 vel0_Neptune;
    vec3 pos0_Pluto; vec3 vel0_Pluto;

    double mass_CoM = masses(0, 0), mass_Mercury = masses(1, 0), mass_Venus = masses(2, 0), mass_Earth = masses(3, 0), mass_Mars = masses(4, 0);
    double mass_Jupiter = masses(5, 0), mass_Saturn= masses(6, 0), mass_Uranus= masses(7, 0), mass_Neptune= masses(8, 0), mass_Pluto= masses(9, 0);

    for (int j=0; j<3;j++){
        pos0_CoM(j) = pos0(0,j); vel0_CoM(j) = vel0(0,j);
        pos0_Mercury(j) = pos0(1,j); vel0_Mercury(j)=vel0(1,j);
        pos0_Venus(j) = pos0(2,j); vel0_Venus(j) = vel0(2,j);
        pos0_Earth(j) = pos0(3,j); vel0_Earth(j) = vel0(3,j);
        pos0_Mars(j) = pos0(4,j); vel0_Mars(j) = vel0(4,j);
        pos0_Jupiter(j) = pos0(5,j); vel0_Jupiter(j) = vel0(5,j);
        pos0_Saturn(j) = pos0(6,j); vel0_Saturn(j) = vel0(6,j);
        pos0_Uranus(j) = pos0(7,j); vel0_Uranus(j) = vel0(7,j);
        pos0_Neptune(j) = pos0(8,j); vel0_Neptune(j) = vel0(8,j);
        pos0_Pluto(j) = pos0(9,j); vel0_Pluto(j) = vel0(9,j);
    }

    CelestialObject &sun = createCelestialObject(pos0_CoM, vel0_CoM, mass_CoM);
    CelestialObject &earth = createCelestialObject(pos0_Earth, vel0_Earth, mass_Earth);
    CelestialObject &jupiter = createCelestialObject(pos0_Jupiter, vel0_Jupiter, mass_Jupiter);
    /*
    CelestialObject &mercury = createCelestialObject(pos0_Mercury, vel0_Mercury, mass_Mercury);
    CelestialObject &venus = createCelestialObject(pos0_Venus, vel0_Venus, mass_Venus);
    CelestialObject &earth = createCelestialObject(pos0_Earth, vel0_Earth, mass_Earth);
    CelestialObject &mars = createCelestialObject(pos0_Mars, vel0_Mars, mass_Mars);
    CelestialObject &jupiter = createCelestialObject(pos0_Jupiter, vel0_Jupiter, mass_Jupiter);
    CelestialObject &saturn = createCelestialObject(pos0_Saturn, vel0_Saturn, mass_Saturn);
    CelestialObject &uranus = createCelestialObject(pos0_Uranus, vel0_Uranus, mass_Uranus);
    CelestialObject &neptune = createCelestialObject(pos0_Neptune, vel0_Neptune, mass_Neptune);
    CelestialObject &pluto = createCelestialObject(pos0_Pluto, vel0_Pluto, mass_Pluto);
*/
}

CelestialObject &SolarSystem::createCelestialObject(vec3 position, vec3 velocity, double mass){
    m_objects.push_back( CelestialObject(position, velocity, mass) );  //.push_bak = add element at end of vector
    return m_objects.back(); // Return reference to the newest added celstial body
}

void SolarSystem::calculateForcesAndEnergy()
{
    double pi = acos(-1);
    double G = 4*pi*pi;

    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
    m_angularMomentum.zeros();


    /*
 * Calculating the gravitational forces as well as potential and kinetic energies:
 */
    for(int i=0; i<numberOfObjects(); i++){
        CelestialObject &object1 = m_objects[i];
        for(int j=i+1; j<numberOfObjects(); j++) {
            CelestialObject &object2 = m_objects[j];
            vec3 dRvec = object2.position - object1.position;
            double dr = dRvec.length();
            //std::cout << "object1_pos = " << object1.position << std::endl;
            //std::cout << "object2_pos = " << object2.position << std::endl;
            //std::cout << dr << " " << dRvec << std::endl;
            vec3 F1 = G*object2.mass*object1.mass*dRvec/(dr*dr*dr);
            vec3 F2 = (-1)*F1;
            object1.addForce(F1);
            object2.addForce(F2);
            //            std::cout << "F1 = " << F1 << std::endl;
            m_kineticEnergy += 0.5*object1.mass*object1.velocity.lengthSquared();
            m_potentialEnergy += (-G*object1.mass*object2.mass)/dr;
        }

        //        m_kineticEnergy += 0.5*object1.mass*object1.velocity.lengthSquared();
        //        m_potentialEnergy += (-G*object1.mass*object2.mass)
    }
}

int SolarSystem::numberOfObjects() const
{
    return m_objects.size();
}

double SolarSystem::totalEnergy() const
{
    return m_kineticEnergy + m_potentialEnergy;
}

double SolarSystem::potentialEnergy() const
{
    return m_potentialEnergy;
}

double SolarSystem::kineticEnergy() const
{
    return m_kineticEnergy;
}

std::vector<CelestialObject> &SolarSystem::objects()
{
    return m_objects;
}


/*
void addNewPlanet(CelestialObject ){

}
*/
