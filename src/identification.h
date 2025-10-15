#pragma once

#include "algebra.h"
#include "random.h"

#include <optional>
#include <variant>
#include <vector>

namespace fasttreeid {

struct IdentificationResult {
	std::vector<std::optional<std::variant<Fraction, Cycle, Path>>> identification;
	uint64_t seed;
	int64_t prime;

	IdentificationResult (size_t n, uint64_t seed_, int64_t prime_) : identification(n), seed(seed_), prime(prime_) {}
};

IdentificationResult identify (const std::vector<std::pair<size_t, size_t>> &bidirected, const std::vector<size_t> &directed, std::optional<uint64_t> seed, std::optional<int64_t> prime, bool (*)(int64_t));

}
