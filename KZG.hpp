#pragma once

#include <mcl/bls12_381.hpp>
#include <vector>

struct Witness
{
    mcl::Fr i;
    mcl::Fr phi_i;
    mcl::G1 witness;

    Witness(mcl::Fr i, mcl::Fr phi_i, mcl::G1 witness)
    {
        this->i = i;
        this->phi_i = phi_i;
        this->witness = witness;
    }
};

struct PK
{
    std::vector<mcl::G1> srs;
    std::vector<mcl::G2> hs;
    
    PK(std::vector<mcl::G1> srs, std::vector<mcl::G2> hs)
    {
        this->srs = srs;
        this->hs = hs;
    }
};

struct Witness_B
{
    std::vector<mcl::Fr> B;
    std::vector<mcl::Fr> r;
    mcl::G1 witness;

    Witness_B(std::vector<mcl::Fr> B, std::vector<mcl::Fr> r, mcl::G1 witness)
    {
        this->B = B;
        this->r = r;
        this->witness = witness;
    }
};

PK Setup(int k, int t);
mcl::G1 Commit(const PK& pk, const std::vector<mcl::Fr>& poly);
Witness CreateWitness(const PK& pk, const std::vector<mcl::Fr>& poly, const mcl::Fr& i);
bool VerifyEval(const PK& pk, const mcl::G1& c, const Witness& w);
Witness_B CreateWitnessBatch(const PK& pk, const std::vector<mcl::Fr>& poly, const std::vector<mcl::Fr>& B);
bool VerifyEvalBatch(const PK& pk, mcl::G1& c, const Witness_B& w);