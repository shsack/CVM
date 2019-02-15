#include "itensor/all.h"
#include "cvmclass.h"
#include <iostream>
#include <complex>
#include <stdlib.h>
#include <math.h>

using namespace itensor;
using namespace std;

const int N = 7;
const int maxm = 100;
const double cut = 1E-6;
// auto sites = Hubbard(N, {"ConserveNf", false}); // Define Hubbard model

auto sites = SpinHalf(N);


auto args = Args({"Maxm", maxm, "Cutoff", cut});


MPO Hamiltonian()
{
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
    for(int j = 1; j < N; ++j)
    {
        ampo += "Sz",j,"Sz",j+1;
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
    }
    return MPO(ampo);
}

MPO Sz(int i)
{
    auto ampo = AutoMPO(sites);
    ampo += "Sz", i;
    return MPO(ampo);
}



double ground_state(MPS & psi, MPO H)
{
    auto sweeps = Sweeps(10);
    sweeps.maxm() = 10, 20, 50, 70; //gradually increase states kept
    sweeps.cutoff() = cut; //desired truncation error

    return dmrg(psi, H, sweeps, "Quiet");
}


MPS conjugate_gradient_squared(MPO H, complex<double> z, MPS b)
{

    auto x = b; // initialize |x> with |b> to iteratively solve A|x> = |b>
    MPS r_old = sum(b, -1 * sum(applyMPO(H, x, args), -z * x, args));
    MPS r_ = r_old;
    MPS p = r_old;
    MPS u = r_old;
    MPS q;
    MPS r_new;
    MPS Ap;
    complex<double> alpha;
    complex<double> beta;
    MPS u_q;

    // Print(overlapC(r_old, r_));

    for(int i = 0; i < 50; ++i){

        Ap = sum(applyMPO(H, p, args), -z * p, args);
        alpha = overlapC(r_old, r_) / overlapC(Ap, r_);
        q = sum(u, -alpha * Ap, args);
        u_q = sum(u, q, args);
        x = sum(x, alpha * u_q, args);
        r_new = sum(r_old, -alpha * sum(applyMPO(H, u_q, args), -z * u_q, args), args);
        beta = overlapC(r_new, r_) / overlapC(r_old, r_);
        u = sum(r_new, beta * q, args);
        p = sum(u, beta * sum(q, beta * p, args), args);
        r_old = r_new;

    }

    cout << "\nResidue = " << norm(r_old) << endl;

    return x;
}

double spectral_function(MPS psi, MPO H, MPO (*Sz)(int), double omega, double eta, complex<double> z, int i, int j)
{
    auto x = conjugate_gradient_squared(H, z, applyMPO(Sz(j), psi, args));
    complex<double> G = overlapC(psi, applyMPO(Sz(i), x, args));

    return - G.imag() / M_PI;
};




int main()
{
    auto psi = MPS(sites);
    auto H = Hamiltonian();
    auto energy = ground_state(psi, H);

    const double omega = 0.;
    const double eta = 0.1;
    const complex<double> z(omega + energy, -eta);

    int i = 2;
    int j = 3;

    double A = spectral_function(psi, H, Sz, omega, eta, z, i, j);
    cout << "\nA = " << A << endl;

    return 0;
}


