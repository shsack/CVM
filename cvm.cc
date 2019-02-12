#include "itensor/all.h"
#include "cvmclass.h"
#include <iostream>
#include <complex>
#include <stdlib.h>

using namespace itensor;
using namespace std;

const int N = 10;
const int maxm = 100;
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


MPS conj_gradient(MPO H, complex<double> z, MPS b)
{
    auto x = b; // initialize |x> with |b> to iteratively solve A|x> = |b>
    MPS r_old = sum(b, -1 * sum(applyMPO(H, x, args), -z * x, args));
    MPS p = r_old;
    MPS r_new;
    MPS Ap = sum(applyMPO(H, p, args), -z * p, args);
    complex<double> alpha;
    complex<double> beta;
    MPS tmp1;
    MPS tmp2;

    for(int i = 0; i < 10; ++i){

        tmp1 = sum(applyMPO(H, r_old, args),  -z * r_old, args);
        alpha = overlapC(r_old, tmp1) / overlapC(Ap, Ap);
        x = sum(x, alpha * p, args);
        r_new = sum(r_old, -alpha * Ap, args);
        tmp2 = sum(applyMPO(H, r_new, args), -z * r_new, args);
        beta = overlapC(r_new, tmp2) / overlapC(r_old, tmp1);
        p = sum(r_new, beta * p, args);
        Ap = sum(tmp2, beta * Ap, args);
        r_old = r_new;

        Print(norm(r_old));

    }//while(abs(norm(r_old) - 1) > 1E-3);


    return x;
}


int main()
{
    auto psi = MPS(sites);
    auto H = Hamiltonian();
    auto energy = ground_state(psi, H);

    const double omega = 0.0;
    const double eta = 0.01;
    const complex<double> z(energy + omega, eta);

//  Solve (H - z)|x> = c|Gs> by iteratively solving A|x> = b

    int i = 5;
    int j = 5;


    auto b = psi;
    b.position(i);
    b.setA(i, (b.A(i) * sites.op("Cdagup", i)).noprime());
    auto x = conj_gradient(H, z, b);
    auto c_x = x;
    c_x.position(j);
    c_x.setA(j, (c_x.A(j) * sites.op("Cdn", j)).noprime());
    double G = overlap(psi, c_x);

    Print(G);

    return 0;
}


