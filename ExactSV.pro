QT += core
QT -= gui

CONFIG += c++11

TARGET = ExactSV
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    SureMap.cpp \
    LocalAligner.cpp

HEADERS += \
    Header.h \
    ExtraTools.h \
    DataStructures.h \
    SureMap.h \
    LocalAligner.h

LIBS += -pthread
