TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += \
        C:\Armadillo\include \
        C:\Users\Elisabeth\Documents\FYS3150\Project_3

SOURCES += \
        Earth_Sun.cpp \
    forwardeuler.cpp \
    initialize.cpp \
    gravitationalforce.cpp \
    celestialobject.cpp \
    velocityverlet.cpp
    solarsystem.cpp

HEADERS += \
        Method_Earth_sun.hpp \
        EnergyTest.hpp\
    forwardeuler.h\
    initialize.h \
    gravitationalforce.h \
    celestialobject.h \
    velocityverlet.h
    solarsystem.h
    WriteToFile.hpp
