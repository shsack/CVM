#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <complex>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
using namespace itensor;
#include "cvm.h"

template<typename sites_t>
MPO heisenberg(sites_t const& sites, int N) {

    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j) {

        ampo += "Sz",j,"Sz",j+1;
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
    }
    return MPO(ampo);
}

template<typename sites_t>
MPO Sx(int i, sites_t const& sites) {
    auto ampo = AutoMPO(sites);
    ampo += "Sx", i;
    return MPO(ampo);
}


double ground_state(MPS & psi, MPO H) {

    auto sweeps = Sweeps(5);
    // sweeps.maxm() = 5, 10; //gradually increase states kept
    // sweeps.cutoff() = cut; //desired truncation error
    return dmrg(psi, H, sweeps, "Quiet");
}


int main(int argc, char* argv[]) {

    const double omega = std::atof(argv[1]);
    const int N = std::atoi(argv[2]);
    const double eta = std::atof(argv[3]);
    const unsigned int max_iter = std::atoi(argv[4]);
    const double tol = std::atof(argv[5]);
    const unsigned int i = std::atoi(argv[6]);
    const unsigned int j = std::atoi(argv[7]);
    const int maxm = std::atoi(argv[8]);
    const double cut = std::atof(argv[9]);

    auto sites = SpinHalf(N);
    MPS psi = MPS(sites);
    MPO H = heisenberg(sites, N);
    double energy = ground_state(psi, H);

    double A = spectral_function(psi, H, Sx(i, sites), Sx(j, sites), omega, eta, energy, tol, max_iter, maxm, cut, sites);

    std::ofstream myfile;
    myfile.open("data/heisenberg_" + std::to_string(omega) + ".txt");

    myfile << A << std::endl;
    myfile.close();

    return 0;
}


