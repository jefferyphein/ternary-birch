bin_PROGRAMS = birch
lib_LTLIBRARIES = libbirch.la

AM_CXXFLAGS = -g -O3 -Wall -Werror -std=c++11
AM_CXXFLAGS += -fvar-tracking-assignments-toggle
AM_CXXFLAGS += -fomit-frame-pointer
AM_CXXFLAGS += -funroll-all-loops
AM_LDFLAGS = -lm -lgmp -lgmpxx

SOURCES  = birch.h
SOURCES += birch_util.cpp
SOURCES += birch_util.h
SOURCES += Eigenvector.h
SOURCES += Fp.cpp
SOURCES += Fp.h
SOURCES += Genus.h
SOURCES += HashMap.h
SOURCES += Isometry.cpp
SOURCES += Isometry.h
SOURCES += IsometrySequence.h
SOURCES += Math.cpp
SOURCES += Math.h
SOURCES += NeighborManager.h
SOURCES += QuadForm.cpp
SOURCES += QuadForm.h
SOURCES += SetCover.cpp
SOURCES += SetCover.h
SOURCES += Spinor.h

birch_SOURCES = birch.cpp $(SOURCES)
birch_CXXFLAGS = $(AM_CXXFLAGS)
libbirch_la_SOURCES = $(SOURCES)
libbirch_la_CXXFLAGS = -shared -fPIC $(AM_CXXFLAGS)
