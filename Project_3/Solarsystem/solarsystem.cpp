#include "solarsystem.h"
#include "gravitationalforce.h"
#include <iostream>
//using namespace std;

SolarSystem::SolarSystem() :
    m_kineticEnergy(0),
    m_potentialEnergy(0)
{
}

CelestialObject& SolarSystem::createCelestialObject(vec3 position, vec3 velocity, double mass) {
    m_objects.push_back( CelestialObject(position, velocity, mass) );  //.push_bak = add element at end of vector
    return m_objects.back(); // Return reference to the newest added celstial body
}

void SolarSystem::calculateForcesAndEnergy()
{
    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
    m_angularMomentum.zeros();
/*
    for(CelestialObject &object : m_objects) {
        // Reset forces on all bodies
        object.force.zeros();
    }
*/

    for(int i=0; i<numberOfObjects(); i++) {
        CelestialObject &object1 = m_objects[i];
        for(int j=i+1; j<numberOfObjects(); j++) {
            //for (int k=0; k<)
            CelestialObject &object2 = m_objects[j];
            vec3 deltaRVector = object1.position - object2.position;
            double dr = deltaRVector.length();
            // Calculate the force and potential energy here

            GravitationalForce::Gravity(&object1, &object2);
        }

        m_kineticEnergy += 0.5*object1.mass*object1.velocity.lengthSquared();
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
/*
void SolarSystem::writeToFile(string filename)
{
    if(!m_file.good()) {
        m_file.open(filename.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }

    m_file << numberOfObjects() << endl;
//    m_file << "Comment line that needs to be here" << endl;
    for(CelestialObject &object : m_objects) {
        m_file << " " << object.position.x() << " " << object.position.y() << " " << object.position.z() << "\n";
    }
}
*/
/*
vec3 SolarSystem::angularMomentum() const
{
    return m_angularMomentum;
}
*/
std::vector<CelestialObject> &SolarSystem::objects()
{
    return m_objects;
}
