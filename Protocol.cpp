#include <KZG.hpp>
#include <NTT.hpp>
#include <iostream>

//univariate zero test PIOP
void UZPIOP();
//univariate sum check PIOP
void USCPIOP();

mcl::Fr evalPoly(const std::vector<mcl::Fr>&, const mcl::Fr&);
mcl::Fr evalZx(size_t, const mcl::Fr&);

int main()
{
    initPairing(mcl::BN_SNARK1);
    NTT_initialize();
    UZPIOP();
    return 0;
}


void UZPIOP()
{
    //initialize H
    size_t l = 2;
    std::vector<mcl::Fr> H(l);
    H[0] = 1; H[1] = -1;

    //initialize f in ZeroTest
    std::vector<mcl::Fr> f(4);
    f[0] = 0; f[1] = -1; f[2] = 0; f[3] = 1;
    std::cout << "Univariate ZeroTest PIOP\n\n";

    //setup KZG
    PK pk = Setup(128, 3);

    //P commits to f
    mcl::G1 PCS_f = Commit(pk, f);
   std:: cout << "Prover's commitment:\nf: " << PCS_f.getStr() << "\n";

    //P commits to q (f = Z * q)
    std::vector<mcl::Fr> q(2);
    q[0] = 0; q[1] = 1;
    mcl::G1 PCS_q = Commit(pk, q);
    std::cout << "q: " << PCS_q.getStr() << "\n\n";

    //V sends challenge x
    mcl::Fr x; x.setRand();
    std::cout << "Verifier's chanllege: " << x.getStr() << "\n\n";

    //P sends evaluations and witness
    mcl::Fr fx, qx, Zx;
    fx = evalPoly(f, x);
    qx = evalPoly(q, x);
    Witness fw = CreateWitness(pk, f, x);
    Witness qw = CreateWitness(pk, q, x);
    std::cout << "Prover:\nf(x) = " << fx << "\nWitness: " << fw.witness.getStr() << "\n"
        << "q(x) = " << qx << "\nWitness: " << qw.witness.getStr() << "\n\n";

    //V verifies evaluations
    bool check = true;
    check = check && VerifyEval(pk, PCS_f, fw);
    check = check && VerifyEval(pk, PCS_q, qw);
    Zx = evalZx(H.size(), x);
    check = check && (fx == qx * Zx);
    if (check)
        std::cout << "Verifier accepts\n";
    else std::cout << "Verifier rejects\n";

    std::cout << "Protocol ends\n";
}

void USCPIOP()
{
    //initialize H
    size_t l = 2;
    std::vector<mcl::Fr> H(l);
    H[0] = 1; H[1] = -1;

    //initialize f in SumCheck
    std::vector<mcl::Fr> f(4);
    f[0] = 0; f[1] = -2; f[2] = 0; f[3] = 1;
    std::cout << "Univariate SumCheck PIOP\n\n";

    //setup KZG
    PK pk = Setup(128, 3);

    //P commits to f
    mcl::G1 PCS_f = Commit(pk, f);
    std::cout << "Prover's commitment:\nf: " << PCS_f.getStr() << "\n";

    //P commits to q (f = Z * q + x * p)
    std::vector<mcl::Fr> q(2);
    q[0] = 0; q[1] = 1;
    mcl::G1 PCS_q = Commit(pk, q);
    std::cout << "q: " << PCS_q.getStr() << "\n";
    
    //P commits to p (f = Z * q + x * p)
    std::vector<mcl::Fr> p(1);
    p[0] = -1;
    mcl::G1 PCS_p = Commit(pk, p);
    std::cout << "p: " << PCS_p.getStr() << "\n\n";

    //V sends challenge x
    mcl::Fr x; x.setRand();
    std::cout << "Verifier's chanllege: " << x.getStr() << "\n\n";

    //P sends evaluations and witness
    mcl::Fr fx, qx, px, Zx;
    fx = evalPoly(f, x);
    qx = evalPoly(q, x);
    px = evalPoly(p, x);
    Witness fw = CreateWitness(pk, f, x);
    Witness qw = CreateWitness(pk, q, x);
    Witness pw = CreateWitness(pk, p, x);
    std::cout << "Prover:\nf(x) = " << fx << "\nWitness: " << fw.witness.getStr() << "\n"
        << "q(x) = " << qx << "\nWitness: " << qw.witness.getStr() << "\n"
        << "p(x) = " << px << "\nWitness: " << pw.witness.getStr() << "\n\n";

    //V verifies evaluations
    bool check = true;
    check = check && VerifyEval(pk, PCS_f, fw);
    check = check && VerifyEval(pk, PCS_q, qw);
    check = check && VerifyEval(pk, PCS_p, pw);
    Zx = evalZx(H.size(), x);
    check = check && (fx == qx * Zx + x * px);
    if (check)
        std::cout << "Verifier accepts\n";
    else std::cout << "Verifier rejects\n";

    std::cout << "Protocol ends\n";
}

mcl::Fr evalPoly(const std::vector<mcl::Fr>& poly, const mcl::Fr& i)
{
    mcl::Fr x = 1;
    mcl::Fr sum = 0;
    size_t n = poly.size();
    for (size_t j = 0; j < n; j++)
    {
        sum += x * poly[j];
        x *= i;
    }
    return sum;
}

mcl::Fr evalZx(size_t order, const mcl::Fr& x)
{
    mcl::Fr p;
    mcl::Fr::pow(p, x, order);
    return p - 1;
}