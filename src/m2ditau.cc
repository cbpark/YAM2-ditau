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
using std::cout;

/// the IDs of invisible particles.
const std::set<long> INVISIBLES = {12, -12, 14, -14, 16, -16, 40, 3000};

/// the IDs of visible particles.
const std::set<long> VISIBLES = {11, -11, 13, -13, 211, -211};

/// to check the event contains ditau.
bool ditau_event(const Particles &ps);

int main(int, char *argv[]) {
    const std::string appname{"m2ditau"};

    cout << appname << ": input LHE file: " << argv[1] << '\n';
    std::ifstream fin{argv[1]};
    LHEF::Reader lhef{fin};

    cout << appname << ": the output will be stored in " << argv[2] << '\n';
    std::ofstream fout{argv[2]};

    long nev = 0;
    while (lhef.readEvent() && ++nev) {
        // get the particles entries from the event.
        const auto ps = analysis::get_particles(lhef.hepeup);
        // parse the initial states.
        const auto taus = analysis::initial_states(ps);
        if (!ditau_event(taus)) {
            cout << appname << ": this isn't an ditau event! (" << nev << ")\n";
            continue;
        }

        // the decay products of taus.
        std::vector<Particles> from_tau;
        for (const auto &tau : taus) {
            from_tau.push_back(analysis::final_states_of(tau.linenum(), ps));
        }
        // to make sure that we have two decay chains.
        if (from_tau.size() != 2) { continue; }

    }  // event loop
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
