#include "lhef_helper.h"
#include <algorithm>
#include <vector>
#include "HepMC3/LHEF.h"

namespace analysis {
bool Particle::produced_from(int parent,
                             const std::vector<Particle> &ps) const {
    const int parent_line = parents_.first;
    if (parent_line == 1 || parent_line == 2) {  // omit the initial states.
        return false;
    } else if (parent_line == parent) {  // we've found the parent.
        return true;
    }

    // the line of the particle in the right previous decay step.
    const auto direct_parent = this->parent(ps);
    if (!direct_parent) {  // orphan?
        return false;
    } else {
        return direct_parent->produced_from(parent, ps);
    }
}

Particles get_particles(const LHEF::HEPEUP &event) {
    Particles ps;
    ps.reserve(event.NUP);
    for (auto iptl = 0; iptl < event.NUP; ++iptl) {
        ps.emplace_back(iptl + 1, event.IDUP[iptl], event.ISTUP[iptl],
                        event.MOTHUP[iptl], event.PUP[iptl]);
    }
    return ps;
}

Particles final_states_of(int parent, const Particles &ps) {
    Particles final_states;
    for (const auto &p : ps) {
        if (p.status() == 1) { final_states.push_back(p); }
    }

    auto pos = std::remove_if(final_states.begin(), final_states.end(),
                              [parent, &ps](const Particle &p) {
                                  return !p.produced_from(parent, ps);
                              });
    final_states.erase(pos, final_states.end());
    return final_states;
}
}  // namespace analysis
