#include "itensor/all.h"
#include "cvmclass.h"
#include <iostream>
#include <complex>
#include <stdlib.h>

using namespace itensor;
using namespace std;

const int N = 5;
const int maxm = 50;
const double cut = 1E-6;
auto sites = Hubbard(N, {"ConserveNf", false}); // Define Hubbard model


auto args = Args({"Maxm", maxm, "Cutoff", cut});


MPO Hamiltonian()
{
    double U = 1.0;
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
    return MPO(ampo);
}

double ground_state(MPS & psi, MPO H)
{
    auto sweeps = Sweeps(5);
    sweeps.maxm() = 10, 20, 50; //gradually increase states kept
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

    Print(overlapC(r_old, r_));

    for(int i = 0; i < 10; ++i){

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

        Print(sqrt(overlapC(r_old, r_old)));

    }//while(abs(norm(r_old) - 1) > 1E-3);

    return x;
}


int main()
{
    auto psi = MPS(sites);
    auto H = Hamiltonian();
    auto energy = ground_state(psi, H);

    Print(energy);

    const double omega = 0.5;
    const double eta = 0.1;
    const complex<double> z(energy + omega, eta);

//  Solve (H - z)|x> = c|Gs> by iteratively solving A|x> = b

    int i = 1;
    int j = 1;

    // c_dagger |Gs>
    auto b = psi;
    b.position(i);
    b.setA(i, (b.A(i) * sites.op("Cdagup", i)).noprime());

    auto x = conjugate_gradient_squared(H, z, b);

    // c |x>
    auto c_x = x;
    c_x.position(j);
    c_x.setA(j, (c_x.A(j) * sites.op("Cdn", j)).noprime());
    double G = overlap(psi, c_x);

    Print(G);

    return 0;
}


