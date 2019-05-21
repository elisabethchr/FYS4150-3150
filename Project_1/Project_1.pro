TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += C:\Armadillo\include\
x64: LIBS += -LC:\Armadillo\examples\lib_win64
    -lblas_win64_MT \
    -llapack_win64_MT

SOURCES += \
    Gaussian_elimination.cpp
