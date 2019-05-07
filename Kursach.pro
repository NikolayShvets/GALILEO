TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    custom.cpp \
    integrator.cpp \
    model.cpp \
    generator.cpp \
    tmatrix.cpp \
    tvector.cpp \
    consumer.cpp \
    mnk.cpp

HEADERS += \
    custom.h \
    integrator.h \
    model.h \
    generator.h \
    tmatrix.h \
    tvector.h \
    consumer.h \
    mnk.h
