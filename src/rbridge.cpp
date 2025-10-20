#include <Rcpp.h>
#include <stdexcept>

#include "identification.h"

namespace {

std::vector<std::pair<size_t, size_t>> convert_bidirected (SEXP bidirectedSEXP, size_t n) {
	Rcpp::IntegerMatrix bidirectedMat(bidirectedSEXP);
	std::vector<std::pair<size_t, size_t>> bidirected;
	bidirected.reserve(bidirectedMat.nrow());
	for (int i = 0; i < bidirectedMat.nrow(); i++) {
		int u = bidirectedMat(i, 0), v = bidirectedMat(i, 1);
		u--, v--;
		if (!(0 <= u && u < (int)n && 0 <= v && v < (int)n && u != v)) {
			throw std::invalid_argument("bidirected should be a matrix with rows (u, v) where 1<=u,v<=n and u != v");
		}
		bidirected.emplace_back(static_cast<size_t>(u), static_cast<size_t>(v));
	}
	return bidirected;
}

std::vector<size_t> convert_directed (SEXP directedSEXP) {
	/* Everything is 1-indexed in R and 0-indexed in C++ */
	Rcpp::IntegerVector directedVec(directedSEXP);
	std::vector<size_t> directed;
	directed.reserve(directedVec.size());
	for (R_xlen_t i = 0; i < directedVec.size(); i++) {
		int x = directedVec[i];
		x--;
		if (!(0 <= x && x <= i)) {
			throw std::invalid_argument("directed should be a vector containing n-1 integers p_2, …, p_n with 1<=p_i<i, the directed parents of 2, …, n");
		}
		directed.push_back(static_cast<size_t>(x));
	}
	return directed;
}

std::optional<uint64_t> convert_seedOrPrime (SEXP seedSEXP) {
	if (Rf_isNull(seedSEXP)) {
		return std::nullopt;
	}
	if (TYPEOF(seedSEXP) != STRSXP || Rf_length(seedSEXP) != 1) {
		throw std::invalid_argument("seed and prime should be (decimal) strings or NULL");
	}
	const char *seedStr = CHAR(STRING_ELT(seedSEXP, 0));
	if (!seedStr) {
		throw std::invalid_argument("seed and prime should be (decimal) strings or NULL");
	}

	errno = 0;
	char *end = nullptr;
	unsigned long long seed = std::strtoull(seedStr, &end, 10); /* 10 is the base */
	if (errno == ERANGE) {
		throw std::invalid_argument("Out of range: seed is limited to 64 bits and prime should be between 2^58 and 2^59");
	}
	while (end && *end == ' ') {
		end++; /* Allow trailing spaces */
	}
	if (end && *end != '\0') {
		throw std::invalid_argument("seed and prime should be pure decimal strings or NULL");
	}
	return static_cast<uint64_t>(seed);
}

SEXP convert_output (fasttreeid::IdentificationResult identificationResult) {
	return Rcpp::wrap(static_cast<int>(identificationResult.identification.size()));
	// TODO
}

bool isProbablyPrime ([[maybe_unused]] int64_t candidate) {
	Rcpp::Function requireNamespace = Rcpp::Function("requireNamespace");

	if (Rcpp::as<bool>(requireNamespace("gmp", Rcpp::Named("quietly") = true))) {
		Rcpp::Environment gmp = Rcpp::Environment::namespace_env("gmp");
		Rcpp::Function as_bigz = gmp["as.bigz"];
		Rcpp::Function isprime = gmp["isprime"];
		Rcpp::CharacterVector candidateStr(1);
		candidateStr[0] = std::to_string(candidate);
		SEXP candidateBigz = as_bigz(candidateStr);

		return Rcpp::as<Rcpp::IntegerVector>(isprime(candidateBigz))[0] >= 1;
	} else {
		throw std::logic_error("The package gmp should be installed if no prime is given");
	}
}

}

// [[Rcpp::export(name = ".fasttreeid_identify_bridge")]]
SEXP fasttreeid_identify_bridge (SEXP bidirectedSEXP, SEXP directedSEXP, SEXP seedSEXP, SEXP primeSEXP) {
	auto directed = convert_directed(directedSEXP);
	auto bidirected = convert_bidirected(bidirectedSEXP, directed.size() + 1);
	auto seed = convert_seedOrPrime(seedSEXP);
	auto prime = static_cast<std::optional<int64_t>>(convert_seedOrPrime(primeSEXP));

	auto identificationResult = fasttreeid::identify(bidirected, directed, seed, prime, isProbablyPrime);

	return convert_output(identificationResult);
}
