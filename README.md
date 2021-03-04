YAM2-ditau
==========
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

An example analysis calculating the M<sub>2</sub> collider observable by using the [YAM2](https://github.com/cbpark/YAM2) library. The process under consideration is e+ e- --> tau+ tau- (ditau) with known center-of-mass energy and longitudinal momentum. It can handle both one- and three-prong decays of tau to leptons and charged pions.

## Dependencies

* [HepMC3](https://gitlab.cern.ch/hepmc/HepMC3) for [`HepMC3/LHEF.h`](http://home.thep.lu.se/~leif/LHEF/): See [`README.md`](https://gitlab.cern.ch/hepmc/HepMC3/-/blob/master/README.md). For Ubuntu, install [`libhepmc3-dev`](https://launchpad.net/ubuntu/+source/hepmc3). After HepMC3 has been installed, make sure that `HepMC3-config` is in the `PATH`.

``` no-hightlight
$ HepMC3-config --version
3.02.03

$ HepMC3-config --prefix
/usr
```

* [NLopt](https://nlopt.readthedocs.io): See [NLopt installation](https://nlopt.readthedocs.io/en/latest/NLopt_Installation/). For some Linux distributions, it can be installed by using the system package manager. For example,

```
# For Arch Linux/Manjaro:
sudo pacman -S nlopt

# For CentOS/Fedora:
sudo dnf install NLopt-devel

# For Debian/Ubuntu:
sudo apt-get install libnlopt-cxx-dev

# For openSUSE:
sudo zypper install nlopt
```

For Ubuntu prior to 20.04 LTS (Focal Fossa), install `libnlopt-dev`. (Check whether there exists `/usr/include/nlopt.hpp`.) In CentOS, the [EPEL](https://fedoraproject.org/wiki/EPEL) repository must be installed and enabled. In macOS, you can install NLopt using [Homebrew](https://brew.sh/).

``` no-hightlight
brew install nlopt
```

* [YAM2](https://github.com/cbpark/YAM2): See [How to build](https://github.com/cbpark/YAM2/blob/master/README.md). We assume that the installation path for YAM2 is `/usr/local`:

``` no-hightlight
cd YAM2
make
sudo DESTDIR=/usr/local make install
```

## Build and run

Before building, see [`Makefile`](./Makefile) to check whether the paths are all correct. The main routine is in [`src/m2ditau.cc`](./src/m2ditau.cc). Then, run

``` no-hightlight
make
```

If the build is successful, the executable `m2ditau` will be created in the `bin` directory.

``` no-hightlight
$ ./bin/m2ditau
usage: ./bin/m2ditau <event.lhe> <output.dat> [mInvisible]
  <event.lhe>: input file in LHEF (required).
  <output.dat>: output file to store the result (required).
  [mInvisible]: the input mass for invisible particles (optional, default = 0)

$ ./bin/m2ditau event.lhe output.dat 0.0
m2ditau: input LHE file: event.lhe
m2ditau: the output will be stored in output.dat
m2ditau: the invisible mass is 0
m2ditau: the number of events processed: 3

$ head output.dat
# M2, k1x, k1y, k1z, k2x, k2y, k2z
   1.70758    -2.15427     0.577034      1.88708     1.88566    -0.912675    0.0649369
   1.72208    -2.96835     -1.12938      3.05195     1.97229      1.12215     -1.12123
   1.09666     1.40419      1.33507      2.26728   -0.230774    -0.320558     -0.62984
```
