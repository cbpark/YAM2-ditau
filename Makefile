PKGNAME  := YAM2-ditau
SRCDIR   := src
BINDIR   := bin
CXXFLAGS := -g -O0 -Wall -Wextra -std=c++17 -pedantic -I$(SRCDIR) $(CXXFLAGS)
LIBS     :=

# Targets
EXE    	 := $(BINDIR)/m2ditau
EXESRC 	 := $(patsubst $(BINDIR)/%,$(SRCDIR)/%.cc,$(EXE))
EXEOBJ 	 := $(EXESRC:.cc=.o)
LIBSRC 	 := $(filter-out $(EXESRC),$(wildcard $(SRCDIR)/*.cc))
LIBOBJ 	 := $(LIBSRC:.cc=.o)

# HepMC3 (https://gitlab.cern.ch/hepmc/HepMC3)
CXXFLAGS += $(shell HepMC3-config --cflags)
LIBS     += $(shell HepMC3-config --libs)

.PHONY: all build clean

all: $(EXE)

$(EXE): $(EXEOBJ) $(LIBOBJ) build
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBOBJ) $(LIBS)

build:
	mkdir -p $(BINDIR)

clean::
	rm -f $(EXE) $(EXEOBJ) $(LIBOBJ)
	rmdir $(BINDIR)
