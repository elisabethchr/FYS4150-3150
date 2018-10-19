#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "celestialobject.h"
#include "gravitationalforce.h"
#include <vector>
#include <string>
#include <fstream>

class SolarSystem
{
public:
    SolarSystem(std::string filename1, std::string filename2, std::string filename3);
    CelestialObject &createCelestialObject(vec3 position, vec3 velocity, double mass);
    void calculateForcesAndEnergy();
    int numberOfObjects() const;

    double totalEnergy() const;
    double potentialEnergy() const;
    double kineticEnergy() const;
    void writeToFile(std::string filename);
    vec3 angularMomentum() const;
    std::vector<CelestialObject> &objects();
    void addNewPlanet();

   // GravitationalForce *gravForce = nullptr;
    std::vector<CelestialObject> m_objects;

private:
//    std::vector<CelestialObject> m_objects;
    vec3 m_angularMomentum;
    std::ofstream m_file;
    double m_kineticEnergy;
    double m_potentialEnergy;
};

#endif // SOLARSYSTEM_H
