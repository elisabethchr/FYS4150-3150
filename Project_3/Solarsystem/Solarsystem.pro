TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += \
        C:\Armadillo\include \
        C:\Users\Elisabeth\Documents\FYS3150\Project_3

SOURCES += \
    celestialobject.cpp \
    solarsystem.cpp \
    MainSS.cpp \
    vec3.cpp \
    forwardeuler.cpp \
    readfile_test.cpp
    velocityverlet.cpp
    gravitationalforce.cpp

HEADERS += \
    initialize.h \
    celestialobject.h \
    solarsystem.h \
    forwardeuler.h \
    vec3.h \
    readfile_test.h
    Readfile.h

    velocityverlet.h
    WriteToFile.hpp
    gravitationalforce.h
    EnergyTest.hpp
