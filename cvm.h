template<typename sites_t>
MPO Iden(sites_t sites) {
    auto ampo = AutoMPO(sites);
    // AutoMPO ampo(sites);
    ampo += "Id", 1;
    return MPO(ampo);
};

MPS conjMPS(MPS psi){
    auto tmp = psi*1.;
    for(int i = 0; i < tmp.N(); ++i){
        tmp.Aref(i+1).conj();
    }
    return tmp;
}


// CVM solver
MPS bicgstab(MPO A, MPS b, double tol, int max_it, Args const& args){

    MPS x = b;
    MPS r_old = sum(b, -1 * applyMPO(A, x, args));
    MPS r_new;
    MPS r_ = r_old;
    MPS p = r_old;
    MPS s;
    MPS Ap;
    MPS As;
    std::complex<double> alpha;
    std::complex<double> beta;
    std::complex<double> w;
    double res;
    int k = 0;

    while(k < max_it){

        Ap = applyMPO(A, p, args);
        alpha = overlapC(conjMPS(r_old), r_) / overlapC(conjMPS(Ap), r_);
        s = sum(r_old, -alpha * Ap, args);
        As = applyMPO(A, s, args);
        w = overlapC(conjMPS(As), s) / overlapC(conjMPS(As), As);
        x = sum(x, sum(alpha * p, w * s, args), args);
        r_new = sum(s, -w * As, args);
        res = sqrt(abs(overlapC(conjMPS(r_new), r_new).real()));

        if(res <= tol){
            std::cout << "Residue = " << res << std::endl;
            break;
        }

        beta = (alpha / w) * overlapC(conjMPS(r_new), r_) / overlapC(conjMPS(r_old), r_);
        p = sum(r_new, beta * sum(p, -w * Ap, args), args);
        r_old = r_new;
        k++;

    };

    return x;
}

// main CVM function
template <typename sites_t>
double spectral_function(MPS psi, MPO H, MPO S1, MPO S2, double omega,
                                 double eta, double energy, double tol, int max_it,
                                 int maxm, double cut, sites_t sites) {

    auto args = Args({"Maxm", maxm, "Cutoff", cut});
    const std::complex<double> z(omega + energy, eta);
    MPO A = sum(z * Iden(sites), -1. * H, args);
    MPS b =  applyMPO(S2, psi, args);
    MPS x = bicgstab(A, b, tol, max_it, args);
    std::complex<double> G = overlapC(psi, S1, x);

    return -G.imag() / M_PI;
};
