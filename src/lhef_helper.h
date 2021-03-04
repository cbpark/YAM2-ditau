/*
 *  Copyright 2021 Chan Beom Park <cbpark@gmail.com>
 *
 *  Types and functions for analyzing events in LHEF.
 */

#ifndef SRC_LHEF_HELPER_H_
#define SRC_LHEF_HELPER_H_

#include <array>          // std::array
#include <cmath>          // std::sqrt
#include <optional>       // std::optional
#include <ostream>        // std::ostream
#include <set>            // std::set
#include <utility>        // std::pair
#include <vector>         // std::vector
#include "HepMC3/LHEF.h"  // LHEF::HEPEUP

namespace analysis {
class Particle;
using Particles = std::vector<Particle>;

class Particle {
private:
    /// the line number in the LHEF event block (starting from 1, not 0).
    int linenum_;
    int id_;
    int status_;
    std::pair<int, int> parents_;
    /// Lab frame momentum (Px, Py, Pz, E) of particle in GeV.
    std::array<double, 4> p_;
    double mass_;

public:
    Particle() = delete;
    Particle(int linenum, int id, int status,
             const std::pair<int, int> &parents, const std::vector<double> &p)
        : linenum_(linenum),
          id_(id),
          status_(status),
          parents_(parents),
          p_({{p[0], p[1], p[2], p[3]}}),
          mass_(p[4]) {}

    int linenum() const { return linenum_; }
    int id() const { return id_; }
    double status() const { return status_; }
    std::array<double, 4> four_momentum() const { return p_; }
    double mass() const { return mass_; }

    std::optional<Particle> parent(const std::vector<Particle> &ps) const {
        for (const auto &p : ps) {
            if (p.linenum() == parents_.first) { return {p}; }
        }
        return {};
    }

    /// is the particle produced from `parent`?
    bool produced_from(int parent, const Particles &ps) const;

    friend Particles initialStates(const Particles &event);
};

/// initial state particles.
inline Particles initialStates(const Particles &ps) {
    Particles initial_states;
    for (const auto &p : ps) {
        if (p.parents_.first == 1) { initial_states.push_back(p); }
    }
    return initial_states;
}

/// collect all the particle entries in LHEF event block.
Particles getParticlesAll(const LHEF::HEPEUP &event);

/// the final state particles in the decay process of the particle in `parent`.
Particles finalStatesOf(int parent, const Particles &ps);

/// collect particles having ID in `pid`.
Particles particlesOf(std::set<int> pid, const Particles &ps);

class FourMomentum {
private:
    std::array<double, 4> p_;

public:
    FourMomentum() = delete;
    FourMomentum(double x, double y, double z, double t) : p_({{x, y, z, t}}) {}
    explicit FourMomentum(const Particle &p) : p_(p.four_momentum()) {}

    double px() const { return p_[0]; }
    double py() const { return p_[1]; }
    double pz() const { return p_[2]; }
    double e() const { return p_[3]; }
    std::optional<double> mass() const {
        double m_sq =
            p_[3] * p_[3] - p_[0] * p_[0] - p_[1] * p_[1] - p_[2] * p_[2];
        if (m_sq < 0.0) {
            return {};
        } else {
            return std::sqrt(m_sq);
        }
    }

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

    friend std::ostream &operator<<(std::ostream &os, const FourMomentum &p);
};

/// returns the four-momentum (Px, Py, Pz, E) summed over the particles (`ps`).
inline FourMomentum sum(const Particles &ps) {
    FourMomentum psum{0.0, 0.0, 0.0, 0.0};
    for (const auto &p : ps) { psum += FourMomentum{p}; }
    return psum;
}
}  // namespace analysis

#endif  // SRC_LHEF_HELPER_H_
