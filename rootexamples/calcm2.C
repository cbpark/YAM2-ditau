R__LOAD_LIBRARY(/usr/local/lib/libYAM2.so)
R__LOAD_LIBRARY(/usr/lib/libnlopt.so)

#include "TLorentzVector.h"
#include "TVector2.h"
#include "YAM2/yam2.h"

void calcm2() {
    // For the include path of YAM2, where YAM2/yam2.h exists.
    gSystem->AddIncludePath("/usr/local/include");

    // note that in the interactive ROOT environment, you should also run
    // gSystem->Load("/usr/local/lib/libYAM2.so");
    // gSystem->Load("/usr/lib/libnlopt.so");

    /*
     *  sqrt(s) --> v1(p1) chi1(k1) + v2(p2) chi2(k2).
     */
    // (p1x, p1y, p1z, e1)
    TLorentzVector v1{-3.3470, -0.2686, 1.8677, 3.8437};
    // (p2x, p2y, p2z, e2)
    TLorentzVector v2{1.0035, -0.2445, 0.340724, 1.0927};
    // (ptmiss_x, ptmiss_y): (k1x + k2x, k1y + k2y)
    TVector2 ptmiss{2.3435, 0.5131};
    // invisible particle mass
    double m_invis = 0.0;
    // longitudinal momentum of the total system: 7 - 4 = 3
    double ptot_z = 3.0;
    // sqrt(s) of the total system: sqrt( (7 + 4)^2 - (7 - 4)^2 )
    double sqrt_s = 10.583;

    auto zero = yam2::FourMomentum();

    auto input = yam2::mkInput({{v1.E(), v1.Px(), v1.Py(), v1.Pz()},
                                {v2.E(), v2.Px(), v2.Py(), v2.Pz()}},
                               {zero, zero}, {ptmiss.Px(), ptmiss.Py()},
                               yam2::Mass{m_invis}, {}, sqrt_s, ptot_z);

    if (!input) {
        std::cerr << "Invalid input!\n";
        return 1;
    }

    const auto m2sol = yam2::m2Cons(input, 1.0e-6, 1000);
    if (!m2sol) {
        std::cerr << "Failed to find minimum!\n";
        return 1;
    }

    std::cout << "M2 = " << m2sol.value().m2() << '\n'
              << "where \n"
              << "  k1x: " << m2sol.value().k1().px()
              << ", k1y: " << m2sol.value().k1().py()
              << ", k1z: " << m2sol.value().k1().pz() << '\n'
              << "  k2x: " << m2sol.value().k2().px()
              << ", k2y: " << m2sol.value().k2().py()
              << ", k2z: " << m2sol.value().k2().pz() << '\n';
}
