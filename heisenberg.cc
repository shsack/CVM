#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <complex>
#include <stdlib.h>
#include <math.h>
#include <iomanip>

using namespace itensor;
#include "cvm.h"

const int N = 7;
// const int maxm = 10;
const double cut = 1E-6;
// auto sites = Hubbard(N, {"ConserveNf", false}); // Define Hubbard model

auto sites = SpinHalf(N);
// auto args = Args({"Maxm", maxm, "Cutoff", cut});




MPO Hamiltonian() {
/*    double U = 1.0;
    double t = 1.0;

    auto ampo = AutoMPO (sites);

    for(int i = 1; i <= N; ++i)
    {
        ampo += U,"Nupdn",i;
    }
    for(int b = 1; b < N; ++b)
    {
        ampo += -t,"Cdagup",b,"Cup",b+1;
        ampo += -t,"Cdagup",b+1,"Cup",b;
        ampo += -t,"Cdagdn",b,"Cdn",b+1;
        ampo += -t,"Cdagdn",b+1,"Cdn",b;
    }
    return MPO(ampo);*/

    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j) {

        ampo += "Sz",j,"Sz",j+1;
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
    }
    return MPO(ampo);
}

MPO Sx(int i) {
    auto ampo = AutoMPO(sites);
    ampo += "Sx", i;
    return MPO(ampo);
}


double ground_state(MPS & psi, MPO H) {

    auto sweeps = Sweeps(5);
    sweeps.maxm() = 5, 10; //gradually increase states kept
    sweeps.cutoff() = cut; //desired truncation error

    return dmrg(psi, H, sweeps, "Quiet");
}



int main(int argc, char* argv[]) {

    auto psi = MPS(sites);

/*
    InitState state(sites,"Up");
    //Now set every other spin to be Dn
    for(int j = 2; j <= N; j += 2)
    {
        state.set(j,"Dn");
    }

    MPS psi(state);
*/


    auto H = Hamiltonian();
    auto energy = ground_state(psi, H);

    const double omega = std::atof(argv[1]);
    const double eta = std::atof(argv[2]);
    const unsigned int max_iter = std::atoi(argv[3]);
    const double tol = std::atof(argv[4]);
    const unsigned int i = std::atoi(argv[5]);
    const unsigned int j = std::atoi(argv[6]);
    const int maxm = std::atoi(argv[7]);
    const double cut = std::atof(argv[8]);

    double A = spectral_function(psi, H, Sx(i), Sx(j), omega, eta, energy, tol, max_iter, maxm, cut, sites);

    std::ofstream myfile;
    myfile.open("data/heisenberg_" + std::to_string(omega) + ".txt");

    myfile << A << std::endl;
    myfile.close();

    return 0;
}


