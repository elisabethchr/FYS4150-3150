#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "celestialobject.h"
#include <vector>
#include <string>
#include <fstream>

class SolarSystem
{
public:
    SolarSystem();
    CelestialObject &createCelestialObject(mat position, mat velocity, double mass);
    void calculateForcesAndEnergy();
    int numberOfObjects() const;

    double totalEnergy() const;
    double potentialEnergy() const;
    double kineticEnergy() const;
    void writeToFile(std::string filename);
    vec3 angularMomentum() const;
    std::vector<CelestialObject> &objects();

private:
    std::vector<CelestialObject> m_objects;
    vec3 m_angularMomentum;
    std::ofstream m_file;
    double m_kineticEnergy;
    double m_potentialEnergy;
};

#endif // SOLARSYSTEM_H
