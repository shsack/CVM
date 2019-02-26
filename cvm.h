MPO I() {
    auto ampo = AutoMPO(sites);
    ampo += "Id", 0;
    return MPO(ampo);
}

MPS conjMPS(MPS psi){
    auto tmp = psi;
    for(int i = 1; i <= N; ++i){
        tmp.Aref(i).conj();
    }
    return tmp;
}

MPS bicstab(MPO A, MPS b, double tol, auto args){

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
    int n = 100;
    int k = 0;

    while(k < n){

        Ap = applyMPO(A, p, args);
        alpha = overlapC(conjMPS(r_old), r_) / overlapC(conjMPS(Ap), r_);
        s = sum(r_old, -alpha * Ap, args);
        As = applyMPO(A, s, args);
        w = overlapC(conjMPS(As), s) / overlapC(conjMPS(As), As);
        x = sum(x, sum(alpha * p, w * s, args), args);
        r_new = sum(s, -w * As, args);
        res = overlapC(conjMPS(r_new), r_new).real();

        if(abs(res) <= tol * tol){
            std::cout << "Residue = " << res << std::endl;
            break;
        }

        beta = (alpha / w) * overlapC(conjMPS(r_new), r_) / overlapC(conjMPS(r_old), r_);
        p = sum(r_new, beta * sum(p, -w * Ap, args), args);
        r_old = r_new;
        k++

    };

    return x;
}

double spectral_function(MPS psi, MPO H, MPO S1, MPO S2, double omega, double eta, double energy, double tol, int i, int j) {

    int maxm = get_int_value("Maxm");
    double cut = get_float_value("Cutoff");
    auto args = Args({"Maxm", maxm, "Cutoff", cut});
    const std::complex<double> z(omega + energy, eta);
    auto A = sum(z * I(), -1. * H, args);
    auto b =  applyMPO(S2j, psi, args);
    auto x = bicstab(A, b, tol, args);
    std::complex<double> G = overlapC(psi, S1, x);

    return -G.imag() / M_PI;
};
