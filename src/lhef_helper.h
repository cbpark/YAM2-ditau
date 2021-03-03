#ifndef SRC_LHEF_HELPER_H_
#define SRC_LHEF_HELPER_H_

#include <array>
#include <optional>
#include <utility>
#include <vector>
#include "HepMC3/LHEF.h"

namespace analysis {
class Particle {
private:
    int linenum_;
    long id_;
    int status_;
    std::pair<int, int> parents_;
    /// Lab frame momentum (Px, Py, Pz, E) of particle in GeV.
    std::array<double, 4> p_;
    double mass_;

public:
    Particle() = delete;
    Particle(int linenum, long id, int status,
             const std::pair<int, int> &parents, const std::vector<double> &p)
        : linenum_(linenum),
          id_(id),
          status_(status),
          parents_(parents),
          p_({{p[0], p[1], p[2], p[3]}}),
          mass_(p[4]) {}

    int linenum() const { return linenum_; }
    double status() const { return status_; }
    double mass() const { return mass_; }

    std::optional<Particle> parent(const std::vector<Particle> &ps) const {
        for (const auto &p : ps) {
            if (p.linenum() == parents_.first) { return {p}; }
        }
        return {};
    }

    bool produced_from(int parent, const std::vector<Particle> &ps) const;
};

using Particles = std::vector<Particle>;

Particles get_particles(const LHEF::HEPEUP &event);

Particles final_states_of(int parent, const Particles &ps);
}  // namespace analysis

#endif  // SRC_LHEF_HELPER_H_
