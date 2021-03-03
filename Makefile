SRCDIR   := src
BINDIR   := bin
CXXFLAGS := -g -O0 -Wall -Wextra -std=c++17 -pedantic -I$(SRCDIR) $(CXXFLAGS)
MKDIR    := mkdir -p
RM       := rm -f

EXE    	 := $(BINDIR)/m2ditau
EXESRC 	 := $(patsubst $(BINDIR)/%,$(SRCDIR)/%.cc,$(EXE))
EXEOBJ 	 := $(EXESRC:.cc=.o)

.PHONY: all build clean

all: $(EXE)

$(EXE): $(EXEOBJ) build
	$(CXX) $(LDFLAGS) -o $@ $<

build:
	$(MKDIR) $(BINDIR)

clean::
	$(RM) $(EXE) $(EXEOBJ)
	$(RM) -r $(BINDIR)
