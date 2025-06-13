//Order: 2^16

#pragma once

#include <mcl/bls12_381.hpp>
#include <vector>

static const size_t order = 1 << 16;
static mcl::Fr g, inv_g;

std::vector<mcl::Fr> iter_NTT(const std::vector<mcl::Fr>&);
std::vector<mcl::Fr> inv_NTT(const std::vector<mcl::Fr>&);
std::vector<mcl::Fr> poly_mult(const std::vector<mcl::Fr>&, const std::vector<mcl::Fr>&);
void NTT_initialize();