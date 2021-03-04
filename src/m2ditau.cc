/*
 *  Copyright 2021 Chan Beom Park <cbpark@gmail.com>
 *
 *  This calculate the M2 for the e+ e- --> tau+ tau- process. using
 *  YAM2 (https://github.com/cbpark/YAM2). The decay topology is
 *
 *  X --> Y1 + Y2 --> vis1(p1) inv1(k1) + vis2(p2) inv2(k2)
 *
 *  where X corresponds to sqrt(s) and Y is decaying particle (tau).
 *  vis and inv are visible and invisible particles, respectively.
 */

#include <cstdlib>        // std::atof
#include <exception>      // std::exception
#include <fstream>        // std:ifstream, std::ofstream
#include <iomanip>        // std::setw
#include <iostream>       // std::cout ,std::cerr
#include <set>            // std::set
#include <utility>        // std::pair
#include <vector>         // std::vector
#include "HepMC3/LHEF.h"  // LHEF::Reader
#include "YAM2/yam2.h"
#include "lhef_helper.h"

// defined in `lhef_helper.h`
using analysis::Particles;
// defined in `lhef_helper.h`
using analysis::FourMomentum;
using std::cerr;
using std::cout;
using std::setw;

/// the IDs of invisible particles (neutrinos, ...)
const std::set<int> INVISIBLES = {12, -12, 14, -14, 16, -16, 40, 3000};

/// the IDs of visible particles (electron, muon, charged pions).
const std::set<int> VISIBLES = {11, -11, 13, -13, 211, -211};

/// to check the event contains ditau (defined below main).
bool ditauEvent(const Particles &ps);

/// to get visible and invisible particle momenta (defined below main).
std::pair<FourMomentum, FourMomentum> getVisInvis(const Particles &ps);

int main(int argc, char *argv[]) {
    const auto appname{"m2ditau"};
    if (!(argc == 3 || argc == 4)) {
        cout << "usage: ./bin/" << appname
             << " <event.lhe> <output.dat> [mInvisible]\n"
             << "  <event.lhe>: input file in LHEF (required).\n"
             << "  <output.dat>: output file to store the result (required).\n"
             << "  [mInvisible]: the input mass ifor invisible particles "
                "(optional, default = 0)\n";
        return 1;
    }

    cout << appname << ": input LHE file: " << argv[1] << '\n';
    std::ifstream fin{argv[1]};
    if (!fin) {
        cerr << appname << ": failed to open the input file!\n";
        return 1;
    }

    try {
        // check whether the LHE file is valid.
        LHEF::Reader lhef{fin};

        cout << appname << ": the output will be stored in " << argv[2] << '\n';
        std::ofstream fout{argv[2]};
        fout << "# M2, k1x, k1y, k1z, k2x, k2y, k2z\n";

        // invisible particle mass.
        yam2::Mass m_inv;
        if (argc == 3) {
            m_inv = yam2::Mass{0.0};
        } else {
            m_inv = yam2::Mass{std::atof(argv[3])};
        }
        cout << appname << ": the invisible mass is " << m_inv.value << '\n';

        // zero momentum for convenience to calculate M2.
        const auto zero = yam2::FourMomentum();

        // --------------------------------------------------------------------------
        // event loop begins
        int nev = 0;
        while (lhef.readEvent() && ++nev) {
            const auto event = lhef.hepeup;

            // sqrt(s) of the particle collision: sqrt( (7 + 4)^2 - (7 - 4)^2 )
            const double sqrt_s = 10.583;
            // const double sqrt_s = event.SCALUP;
#ifdef DEBUG
            event.print(cout);
            cout << "sqrt(s) = " << sqrt_s << '\n';
#endif

            // get all the particles entries from the event.
            const auto ps = analysis::getParticlesAll(event);
            // parse the initial states.
            const auto taus = analysis::initialStates(ps);
            if (!ditauEvent(taus)) {
                cerr << appname << ": this isn't ditau event! (" << nev
                     << ")\n";
                continue;
            }

#ifdef DEBUG
            cout << "m_inv = " << m_inv.value << '\n';
#endif
            if (m_inv.value > taus[0].mass() || m_inv.value > taus[1].mass()) {
                cerr << appname << ": the input mass is unphysical! (" << nev
                     << ")\n";
                continue;
            }

            // the longitudinal momentum of total system: 7 - 4 GeV.
            // const oduble pz_ditau = 3.0;
            const double pz_ditau =
                taus[0].four_momentum()[2] + taus[1].four_momentum()[2];
#ifdef DEBUG
            cout << "pz_ditau = " << pz_ditau << '\n';
#endif

            // the decay products of taus.
            std::vector<Particles> from_tau;
            for (const auto &tau : taus) {
                from_tau.push_back(analysis::finalStatesOf(tau.linenum(), ps));
            }
            // to make sure that we have two decay chains.
            if (from_tau.size() != 2) { continue; }

            const auto &[vis1, inv1] = getVisInvis(from_tau[0]);
            const auto &[vis2, inv2] = getVisInvis(from_tau[1]);
#ifdef DEBUG
            cout << "vis1: " << vis1 << "\ninv1: " << inv1 << '\n';
            cout << "vis2: " << vis2 << "\ninv2: " << inv2 << '\n';
#endif

            // missing transverse momentum.
            const yam2::TransverseMomentum ptmiss{inv1.px() + inv2.px(),
                                                  inv1.py() + inv2.py()};
            // input kinematic configuration for M2.
            const auto input = yam2::mkInput(
                {{vis1.e(), vis1.px(), vis1.py(), vis1.pz()},
                 {vis2.e(), vis2.px(), vis2.py(), vis2.pz()}},
                {zero, zero}, ptmiss, m_inv, {}, sqrt_s, {pz_ditau});
            if (!input) {
                cerr << appname << ": invalid input! (" << nev << ")\n";
                continue;
            }

            // calculate M2. the latter arguments are tolerance and maximal
            // evaluation number.
            const auto m2sol = yam2::m2Cons(input.value(), 1.0e-6, 1000);
            if (!m2sol) {
                cerr << appname << ": failed to find minimum! (" << nev
                     << ")\n";
                continue;
            } else {
                // print the M2 value.
                fout << setw(10) << m2sol.value().m2();

                // the invisible momentum solution of M2.
                const auto k1sol = m2sol.value().k1();
                const auto k2sol = m2sol.value().k2();
                // print the momentum solution.
                fout << std::right << setw(12) << k1sol.px() << ' ' << setw(12)
                     << k1sol.py() << ' ' << setw(12) << k1sol.pz();
                fout << std::right << setw(12) << k2sol.px() << ' ' << setw(12)
                     << k2sol.py() << ' ' << setw(12) << k2sol.pz() << '\n';
            }

            if ((nev < 10000 && nev % 1000 == 0) || nev % 10000 == 0) {
                cout << appname << ": ... processed " << nev << " events.\n";
            }
        }
        // event loop ends
        // --------------------------------------------------------------------------

        // close the input LHE file.
        fin.close();

        // close the output file.
        fout.close();

        // finished.
        cout << appname << ": the number of events processed: " << nev << '\n';
    } catch (std::exception &err) {
        // invalid LHE file.
        cerr << appname << ": " << err.what() << '\n';
        fin.close();
        return 1;
    }
}

bool ditauEvent(const Particles &ps) {
    if (ps.size() != 2) { return false; }

    if (ps[0].id() == 15 && ps[1].id() == -15) {
        return true;
    } else if (ps[0].id() == -15 && ps[1].id() == 15) {
        return true;
    } else {
        return false;
    }
}

std::pair<FourMomentum, FourMomentum> getVisInvis(const Particles &ps) {
    const auto visibles = particlesOf(VISIBLES, ps);
    const auto p_vis = analysis::sum(visibles);
    const auto invisibles = particlesOf(INVISIBLES, ps);
    const auto p_invis = analysis::sum(invisibles);
    return {p_vis, p_invis};
}
