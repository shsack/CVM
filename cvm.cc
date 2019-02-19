#include "itensor/all.h"
#include "cvmclass.h"
#include <iostream>
#include <fstream>
#include <complex>
#include <stdlib.h>
#include <math.h>
#include <iomanip>


using namespace itensor;

const int N = 3;
const int maxm = 1000;
const double cut = 1E-13;
// auto sites = Hubbard(N, {"ConserveNf", false}); // Define Hubbard model

auto sites = SpinHalf(N);


auto args = Args({"Maxm", maxm, "Cutoff", cut});


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


MPO I() {

    auto ampo = AutoMPO(sites);
    for(int i = 1; i <= N; ++i){

        ampo += "Id", i;

    }

    return (1./N) * MPO(ampo);

}


double ground_state(MPS & psi, MPO H) {

    auto sweeps = Sweeps(5);
    sweeps.maxm() = 10, 20, 50, 70; //gradually increase states kept
    sweeps.cutoff() = cut; //desired truncation error

    return dmrg(psi, H, sweeps, "Quiet");
}

MPS conjugate_(MPS psi){

    auto tmp = psi;

    for(int i = 1; i <= N; ++i){

        tmp.Aref(i).conj();

    }

    return tmp;
}


MPS conjugate_gradient_squared(MPO A, MPS b, int iter) {

    auto x = b; // initialize |x> with |b> to iteratively solve A|x> = |b>
    MPS r_old = sum(b, -1 * applyMPO(A, x, args));
    MPS r_ = r_old;
    MPS p = r_old;
    MPS u = r_old;
    MPS q;
    MPS r_new;
    MPS Ap;
    std::complex<double> alpha;
    std::complex<double> beta;

    //
    //alpha = 0.0;
    //auto alpha = std::complex<long double>(std::real(0), std::imag(0));

    MPS u_q;

    // Print(overlapC(r_old, r_));

    for(int i = 0; i < iter; ++i){

        Ap = applyMPO(A, p, args);
        alpha = overlapC(conjugate_(r_old), r_) / overlapC(conjugate_(Ap), r_); // print

        // Print(alpha);

        //Print(overlapC(conjugate_(Ap), r_));

        std::cout << alpha << std::endl;

        q = sum(u, -alpha * Ap, args);
        u_q = sum(u, q, args);
        x = sum(x, alpha * u_q, args);
        r_new = sum(r_old, -alpha * applyMPO(A, u_q, args), args);
        beta = overlapC(conjugate_(r_new), r_) / overlapC(conjugate_(r_old), r_); // prin


        // Print(beta);

        // Print(overlapC(conjugate_(r_old), r_));

        std::cout << beta << std::endl;


        u = sum(r_new, beta * q, args);
        p = sum(u, beta * sum(q, beta * p, args), args);
        r_old = r_new;

    }

    std::cout << "\nResidue = " << norm(r_old) << std::endl;

    return x;
}

double spectral_function(MPS psi, MPO H, MPO (*Sx)(int), double omega, double eta, double energy, int iter, int i, int j) {

    const std::complex<double> z(omega + energy, eta);
    auto A = sum(z * I(), -1. * H, args);
    auto b =  applyMPO(Sx(j), psi, args);
    auto x = conjugate_gradient_squared(A, b, iter);

    std::complex<double> G = overlapC(psi, Sx(i), x);

    return -G.imag() / M_PI;
};


int main(int argc, char* argv[]) {

    // auto psi = MPS(sites);

    InitState state(sites,"Up");
//Now set every other spin to be Dn
    for(int j = 2; j <= N; j += 2)
    {
        state.set(j,"Dn");
    }

    MPS psi(state);


    auto H = Hamiltonian();
    auto energy = ground_state(psi, H);

    const double omega = std::atof(argv[1]);
    const double eta = std::atof(argv[2]);
    const unsigned int iter = std::atoi(argv[3]);
    const unsigned int i = std::atoi(argv[4]);
    const unsigned int j = std::atoi(argv[5]);

    double A = spectral_function(psi, H, Sx, omega, eta, energy, iter, i, j);
/*
    std::cout << "\nomega = " << omega << std::endl;
    std::cout << "\neta = " << eta << std::endl;
    std::cout << "\nA = " << A << std::endl;

*/

/*
    std::string str = std::to_string (omega);
    str.erase ( str.find('.') + 2, std::string::npos );
*/


    // Print(omega);
    std::ofstream myfile;
    myfile.open("data/cvm_" + std::to_string(omega) + ".txt");

    myfile << A << std::endl;
    myfile.close();

    return 0;
}


