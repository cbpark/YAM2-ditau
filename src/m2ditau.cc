#include <fstream>
#include <iostream>
#include "HepMC3/LHEF.h"
#include "lhef_helper.h"

using std::cout;

int main(int, char *argv[]) {
    const std::string appname{"m2ditau"};

    cout << appname << ": input LHE file: " << argv[1] << '\n';
    std::ifstream fin{argv[1]};
    LHEF::Reader lhef{fin};

    cout << appname << ": the output will be stored in " << argv[2] << '\n';
    std::ofstream fout{argv[2]};

    long nev = 0;
    while (lhef.readEvent()) {
        const auto event = lhef.hepeup;
        // event.print(cout);

        const auto ps = analysis::get_particles(event);
        const auto from_tau1 = analysis::final_states_of(3, ps);
        for (const auto &p : from_tau1) { cout << p.linenum() << '\n'; }
        ++nev;
    }
    fin.close();

    cout << appname << ": the number of events processed: " << nev << '\n';
}
