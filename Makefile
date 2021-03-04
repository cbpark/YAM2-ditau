PKGNAME  := YAM2-ditau
SRCDIR   := src
BINDIR   := bin
CXXFLAGS := -g -O2 -Wall -Wextra -std=c++17 -pedantic -I$(SRCDIR) $(CXXFLAGS)
LIBS     := -lm

# Targets
EXE    	 := $(BINDIR)/m2ditau
EXESRC 	 := $(patsubst $(BINDIR)/%,$(SRCDIR)/%.cc,$(EXE))
EXEOBJ 	 := $(EXESRC:.cc=.o)
LIBSRC 	 := $(filter-out $(EXESRC),$(wildcard $(SRCDIR)/*.cc))
LIBOBJ 	 := $(LIBSRC:.cc=.o)

# HepMC3 (https://gitlab.cern.ch/hepmc/HepMC3)
CXXFLAGS += $(shell HepMC3-config --cflags)
LIBS     += $(shell HepMC3-config --libs)

# YAM2 (https://github.com/cbpark/YAM2)
YAM2     ?= /usr/local
CXXFLAGS += -I$(YAM2)/include
LIBS     += -L$(YAM2)/lib -lYAM2 -Wl,-rpath $(YAM2)/lib

# NLopt (https://nlopt.readthedocs.io/
NLOPT    ?= /usr
LIBS     += -L$(NLOPT)/lib -lnlopt -Wl,-rpath $(NLOPT)/lib

.PHONY: all build clean

all: $(EXE)

$(EXE): $(EXEOBJ) $(LIBOBJ) build
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBOBJ) $(LIBS)

build:
	mkdir -p $(BINDIR)

clean::
	rm -f $(EXE) $(EXEOBJ) $(LIBOBJ)
	rmdir $(BINDIR)
