#include <KZG.hpp>
#include <NTT.hpp>

static std::vector<mcl::Fr> calc_B_poly(const std::vector<mcl::Fr>&, size_t, size_t);
static void poly_div(const std::vector<mcl::Fr>& f, const std::vector<mcl::Fr>& Z, std::vector<mcl::Fr>& Q, std::vector<mcl::Fr>& R);

PK Setup(int k, int t)
{
    if (k > 128)
    {
        std::cout << "Unable to provide required security level with BN_SNARK1!";
        exit(-1);
    }
    mcl::G1 g; mcl::Fr alpha;
    mcl::Fp r; r.setByCSPRNG();
    mapToG1(g, r);    
    mcl::G2 h; mcl::hashAndMapToG2(h, "yeah");
    alpha.setByCSPRNG(); if (alpha == 0) alpha = alpha + 1;
    std::vector<mcl::G1> v1(t + 1);
    std::vector<mcl::G2> v2(t + 1);
    mcl::Fr p = 1;
    for (int i = 0; i <= t; i++)
    {
        v1[i] = g * p;
        v2[i] = h * p;
        p = p * alpha;
    }

    return PK(v1, v2);
}

mcl::G1 Commit(const PK& pk, const std::vector<mcl::Fr>& poly)
{
    if (poly.size() > pk.srs.size())
    {
        std::cout << "Invalid input polynomial!";
        exit(-1);
    }
    mcl::G1 commit = (pk.srs[0] * poly[0]);
    for (int i = 1; i < poly.size(); i++)
        commit = commit + (pk.srs[i] * poly[i]);
    return commit;
}

Witness CreateWitness(const PK& pk, const std::vector<mcl::Fr>& poly, const mcl::Fr& i)
{
    size_t t = poly.size() - 1;
    mcl::Fr x = 1;
    mcl::Fr phi_i = 0;
    for (int j = 0; j <= t; j++)
    {
        phi_i = phi_i + x * poly[j];
        x = x * i;
    }
    if (t == 0) return Witness(i, phi_i, pk.srs[0] * 0);
    std::vector<mcl::Fr> ps1(t);
    ps1[t - 1] = poly[t];
    for (int j = t - 2; j >= 0; j--)
        ps1[j] = poly[j + 1] + i * ps1[j + 1];
    mcl::G1 witness = pk.srs[0] * ps1[0];
    for (int i = 1; i < t; i++)
        witness = witness + (pk.srs[i] * ps1[i]);
    return Witness(i, phi_i, witness);
}

bool VerifyEval(const PK& pk, const mcl::G1& c, const Witness& w)
{
    mcl::Fp12 e1;
    mcl::pairing(e1, w.witness, pk.hs[1] - (pk.hs[0] * w.i));
    mcl::Fp12 e2;
    mcl::pairing(e2, c - (pk.srs[0] * w.phi_i), pk.hs[0]);
    return e1 == e2; 
}

Witness_B CreateWitnessBatch(const PK& pk, const std::vector<mcl::Fr>& poly, const std::vector<mcl::Fr>& B)
{
    std::vector<mcl::Fr> B_poly = calc_B_poly(B, 0, B.size());
    std::vector<mcl::Fr> Q, r;
    if (B_poly.size() > poly.size())
    {
        std::vector<mcl::Fr> r;
        r.assign(poly.begin(), poly.end());
        return Witness_B(B, r, pk.srs[0] * 0);
    }
    poly_div(poly, B_poly, Q, r);
    mcl::G1 witness = pk.srs[0] * Q[0];
    for (int i = 1; i < Q.size(); i++)
        witness = witness + (pk.srs[i] * Q[i]);
    return Witness_B(B, r, witness);
}

bool VerifyEvalBatch(const PK& pk, mcl::G1& c, const Witness_B& w)
{
    mcl::Fp12 e1, e2, e3;
    mcl::pairing(e1, c, pk.hs[0]);
    mcl::G2 r_alpha = pk.hs[0] * w.r[0];
    for (int i = 1; i < w.r.size(); i++)
        r_alpha = r_alpha + (pk.hs[i] * w.r[i]);
    mcl::pairing(e2, pk.srs[0], r_alpha);
    std::vector<mcl::Fr> Z = calc_B_poly(w.B, 0, w.B.size());
    mcl::G2 h_pi = pk.hs[0] * Z[0];
    for (int i = 1; i < Z.size(); i++)
        h_pi = h_pi + (pk.hs[i] * Z[i]);
    mcl::pairing(e3, w.witness, h_pi);
    return e1 == e2 * e3;
}

static std::vector<mcl::Fr> calc_B_poly(const std::vector<mcl::Fr>& B, size_t begin, size_t end)
{
    if (end - begin == 1)
    {
        std::vector<mcl::Fr> v(2);
        v[0] = -1 * B[begin]; v[1] = 1;
        return v;
    }
    else
    {
        size_t mid = begin + (end - begin) / 2;
        std::vector<mcl::Fr> v1, v2;
        v1 = calc_B_poly(B, begin, mid);
        v2 = calc_B_poly(B, mid, end);
        v1 = poly_mult(v1, v2);
        return v1;
    }
}

static void poly_div(const std::vector<mcl::Fr>& f, const std::vector<mcl::Fr>& Z, std::vector<mcl::Fr>& Q, std::vector<mcl::Fr>& r)
{
    size_t n = f.size(), m = Z.size();
    r.assign(f.begin(), f.end());
    Q = std::vector<mcl::Fr>(n - m + 1);
    for (int i = n - m; i >= 0; --i) {
        mcl::Fr t = r[i + m - 1] / Z[m - 1];
        Q[i] = t;
        for (int j = i; j < m; ++j)
            r[j] = r[j] - t * Z[j - i];
    }
    r.resize(m - 1);
}