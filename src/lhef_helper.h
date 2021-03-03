#ifndef SRC_LHEF_HELPER_H_
#define SRC_LHEF_HELPER_H_

#include <array>
#include <optional>
#include <utility>
#include <vector>
#include "HepMC3/LHEF.h"

namespace analysis {
class Particle;
using Particles = std::vector<Particle>;

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
    long id() const { return id_; }
    double status() const { return status_; }
    std::array<double, 4> four_momentum() const { return p_; }
    double mass() const { return mass_; }

    std::optional<Particle> parent(const std::vector<Particle> &ps) const {
        for (const auto &p : ps) {
            if (p.linenum() == parents_.first) { return {p}; }
        }
        return {};
    }

    bool produced_from(int parent, const Particles &ps) const;

    friend Particles initial_states(const Particles &event);
};

inline Particles initial_states(const Particles &ps) {
    Particles initial_states;
    for (const auto &p : ps) {
        if (p.parents_.first == 1) { initial_states.push_back(p); }
    }
    return initial_states;
}

Particles get_particles(const LHEF::HEPEUP &event);

Particles final_states_of(int parent, const Particles &ps);

Particles particles_of(std::set<long> pid, const Particles &ps);

class FourMomentum {
private:
    std::array<double, 4> p_;

public:
    FourMomentum() = delete;
    FourMomentum(double x, double y, double z, double t) : p_({{x, y, z, t}}) {}
    FourMomentum(const Particle &p) : p_(p.four_momentum()) {}

    double px() const { return p_[0]; }
    double py() const { return p_[1]; }
    double pz() const { return p_[2]; }
    double e() const { return p_[3]; }

    FourMomentum &operator+=(const FourMomentum &p) {
        this->p_[0] += p.px();
        this->p_[1] += p.py();
        this->p_[2] += p.pz();
        this->p_[3] += p.e();
        return *this;
    }

    friend FourMomentum operator+(FourMomentum p1, const FourMomentum &p2) {
        p1 += p2;
        return p1;
    }
};

/// returns the four-momentum (Px, Py, Pz, E) summed over the particles (`ps`).
inline FourMomentum sum(const Particles &ps) {
    FourMomentum psum{0.0, 0.0, 0.0, 0.0};
    for (const auto &p : ps) { psum += FourMomentum{p}; }
    return psum;
}
}  // namespace analysis

#endif  // SRC_LHEF_HELPER_H_
