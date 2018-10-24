# Project 3
FYS4150/FYS3150: Project 3 by Eirik Bratli and Elisabeth Christensen

#### Project description:
This project aims at the use of the forward Euler method and the velocity Verlet algorithm , through object-oriented code, to calculate the orbits of astronomical objects within our Solar system. By testing their accuracy through both conservation of energy and angular momentum we find that the velocity Verlet algorithm is preferrable over the forward Euler method due to it not being asymmetric in time such as the forward Euler is.

#### Files included:
##### Earth_Sun_system (non-object oriented)
  
  - celestialobject.h/cpp
  - Earth_Sun.cpp (main file)
  - gravitationalforce.h/cpp
  - initialize.h/cpp
  - solarsystem.h/cpp
  - velocityerlet.h/cpp
  - Method_Earth_Sun.hpp

##### Solarsystem (object oriented)
  - celestialobject.h/cpp
  - forwardeuler.h/cpp
  - MainSS.cpp
  - readfile_test.h/cpp
  - solarsystem.h/cpp
  - velocityverlet.h/cpp
  - writetofile.h/cpp
  - GetPlanetValues.py
