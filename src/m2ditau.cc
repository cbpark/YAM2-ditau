/*
 *  Copyright 2021 Chan Beom Park <cbpark@gmail.com>
 */

#include <array>
#include <fstream>
#include <iostream>
#include <set>
#include <utility>
#include <vector>
#include "HepMC3/LHEF.h"  // LHEF::Reader
#include "lhef_helper.h"

// defined in `lhef_helper.h`
using analysis::Particles;
// defined in `lhef_helper.h`
using analysis::FourMomentum;
using std::cout;

/// the IDs of invisible particles.
const std::set<int> INVISIBLES = {12, -12, 14, -14, 16, -16, 40, 3000};

/// the IDs of visible particles.
const std::set<int> VISIBLES = {11, -11, 13, -13, 211, -211};

/// to check the event contains ditau.
bool ditau_event(const Particles &ps);

/// to get visible and invisible particle momenta.
std::pair<FourMomentum, FourMomentum> get_vis_invis(const Particles &ps);

int main(int, char *argv[]) {
    const std::string appname{"m2ditau"};

    cout << appname << ": input LHE file: " << argv[1] << '\n';
    std::ifstream fin{argv[1]};
    LHEF::Reader lhef{fin};

    cout << appname << ": the output will be stored in " << argv[2] << '\n';
    std::ofstream fout{argv[2]};

    int nev = 0;
    // -- event loop begins -- //
    while (lhef.readEvent() && ++nev) {
        const auto event = lhef.hepeup;

        // sqrt(s) of the particle collision.
        const double sqrts = event.SCALUP;
#ifdef DEBUG
        event.print(cout);
        cout << "sqrt(s) = " << sqrts << '\n';
#endif

        // get all the particles entries from the event.
        const auto ps = analysis::get_particles_all(event);
        // parse the initial states.
        const auto taus = analysis::initial_states(ps);
        if (!ditau_event(taus)) {
            cout << appname << ": this isn't ditau event! (" << nev << ")\n";
            continue;
        }

        // the longitudinal momentum of total system.
        const double pz_ditau =
            taus[0].four_momentum()[2] + taus[1].four_momentum()[2];
#ifdef DEBUG
        cout << "pz_ditau = " << pz_ditau << '\n';
#endif

        // the decay products of taus.
        std::vector<Particles> from_tau;
        for (const auto &tau : taus) {
            from_tau.push_back(analysis::final_states_of(tau.linenum(), ps));
        }
        // to make sure that we have two decay chains.
        if (from_tau.size() != 2) { continue; }

        auto [vis1, inv1] = get_vis_invis(from_tau[0]);
        auto [vis2, inv2] = get_vis_invis(from_tau[1]);

#ifdef DEBUG
        cout << "vis1: " << vis1 << ", inv1: " << inv1 << '\n';
        cout << "vis2: " << vis2 << ", inv2: " << inv2 << '\n';
#endif
    }
    // -- event loop ends -- //
    fin.close();

    cout << appname << ": the number of events processed: " << nev << '\n';
}

bool ditau_event(const Particles &ps) {
    if (ps.size() != 2) { return false; }

    if (ps[0].id() == 15 && ps[1].id() == -15) {
        return true;
    } else if (ps[0].id() == -15 && ps[1].id() == 15) {
        return true;
    }

    return false;
}

std::pair<FourMomentum, FourMomentum> get_vis_invis(const Particles &ps) {
    const auto visibles = particles_of(VISIBLES, ps);
    const auto p_vis = analysis::sum(visibles);
    const auto invisibles = particles_of(INVISIBLES, ps);
    const auto p_invis = analysis::sum(invisibles);
    return {p_vis, p_invis};
}
