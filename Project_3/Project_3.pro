TEMPLATE = subdirs
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += C:\Armadillo\include

SOURCES += \


x64: LIBS += \
    -L\C:\Users\Elisabeth\Documents\cpp1\armadillo\examples\lib_win64\
    -lblas_win64_MT \
    -llapack_win64_MT

SUBDIRS += \
    Earth_Sun_system
