TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    custom.cpp \
    integrator.cpp \
    model.cpp \
    generator.cpp \
    consumer.cpp

HEADERS += \
    custom.h \
    integrator.h \
    model.h \
    generator.h \
    consumer.h
LIBS += -llapack -lblas -larmadillo
