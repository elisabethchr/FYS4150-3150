TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += \
        C:\Armadillo\include \
        C:\Users\Elisabeth\Documents\FYS3150\Project_3

SOURCES += \
    forwardeuler.cpp \
    gravitationalforce.cpp \
    celestialobject.cpp \
    velocityverlet.cpp \
    solarsystem.cpp \
    MainSS.cpp

HEADERS += \
    forwardeuler.h\
    initialize.h \
    gravitationalforce.h \
    celestialobject.h \
    velocityverlet.h \
    solarsystem.h \
    Readfile.h
    WriteToFile.hpp
    EnergyTest.hpp
