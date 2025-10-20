#include "../src/algebra.h"
#include "../src/identification.h"

#include <cstring>
#include <iostream>
#include <stdexcept>

#if defined(HAVE_GMP)
#include <gmp.h>
#elif defined(HAVE_OPENSSL)
#include <openssl/bn.h>
#endif

namespace {

enum Verbosity {
	minimal,
	normal
};

#if defined(HAVE_GMP)

bool isProbablyPrimeGMP (int64_t candidate) {
	mpz_t z;
	mpz_init(z);
	mpz_import(z, 1, -1, sizeof(candidate), 0, 0, &candidate);
	int result = mpz_probab_prime_p(z, 40);
	mpz_clear(z);
	return result > 0;
}

#elif defined(HAVE_OPENSSL)

bool isProbablyPrimeOPENSSL (int64_t candidate) {
	unsigned char bytes[8];
	for (size_t i = 0; i < 8; i++) {
		bytes[i] = static_cast<unsigned char>((candidate >> (8 * i)) & 0xFF);
	}
	BIGNUM *bn = BN_lebin2bn(bytes, 8, nullptr);
	if (!bn) {
		throw std::runtime_error("OpenSSL: BN_lebin2bn failed");
	}
	BN_CTX *ctx = BN_CTX_new();
	if (!ctx) {
		BN_free(bn);
		throw std::runtime_error("OpenSSL: BN_CTX_new failed");
	}

#if OPENSSL_VERSION_NUMBER >= 0x30000000L
	int result = BN_check_prime(bn, ctx, nullptr);
#else
	int result = BN_is_prime_ex(bn, 40, ctx, nullptr);
#endif

	BN_CTX_free(ctx);
	BN_free(bn);
	if (result < 0) {
		throw std::runtime_error("OpenSSL prime check failed with error code < 0");
	}
	return result > 0;
}

#else

bool isProbablyPrimeError ([[maybe_unused]] int64_t candidate) {
	throw std::logic_error("Neither HAVE_GMP nor HAVE_OPENSSL defined but no prime given");
}

#endif

}

int main (int argc, char **argv) {
	std::optional<uint64_t> seed;
	std::optional<int64_t> prime;
	Verbosity verbosity = normal; 

	for (int i = 1; i < argc; i++) {
		if (std::strcmp(argv[i], "--seed") == 0) {
			if (! (i + 1 < argc)) {
				throw std::invalid_argument("When using --seed SEED, SEED cannot be omitted");
			}
			seed = std::stoull(argv[i + 1]);
			i++;
		} else if (std::strcmp(argv[i], "--prime") == 0) {
			if (! (i + 1 < argc)) {
				throw std::invalid_argument("When using --prime PRIME, PRIME cannot be omitted");
			}
			prime = std::stoll(argv[i + 1]);
			if (! ((1ll << 58) <= prime && prime < (1ll << 59))) {
				throw std::invalid_argument("When using --prime PRIME, PRIME must be between 2^58 and 2^59");
			}
			i++;
		} else if (std::strcmp(argv[i], "--minimal") == 0) {
			verbosity = minimal;
		} else if (std::strcmp(argv[i], "--help") == 0) {
			std::cout << "Usage: " << argv[0] << " [OPTION]..." << std::endl
					  << "Reads the input from stdin and writes to stdout" << std::endl << std::endl
					  << "Input format:" << std::endl
					  << "n\np_1 … p_{n-1}\nm\nu_0 v_0\n…\nu_{m-1} v_{m-1}\n" << std::endl
					  << "where n is the number of nodes, m is the number of bidirected edges," << std::endl
					  << "the nodes are numbered 0, …, n-1," << std::endl
					  << "p_1, …, p_{n-1} are the directed parents of 1, …, n-1" << std::endl
					  << "and {u_0, v_0}, …, {u_{m-1},v_{m-1}} are the bidirected edges" << std::endl << std::endl
					  << "Options:" << std::endl
					  << "   --seed SEED        Seed the random with SEED. Generates a random seed by default." << std::endl
					  << "   --prime PRIME      Use PRIME as the prime. PRIME should be between 2^58 and 2^59." << std::endl
					  << "                      Selects a random prime in this range by default." << std::endl
					  << "   --minimal          Only output n-1 integers i_1 … i_{n-1}, i_j∈{0,1,2}" << std::endl
					  << "                      0 means unidentifiable, 1 means 1-identifiable, 2 means 2-identifiable" << std::endl
					  << "   --help             Display this help and exit" << std::endl;
			return 0;
		} else {
			throw std::invalid_argument(std::string("Unknown command line option: ") + argv[i]);
		}
	}

	size_t n, m;
	if (! ((std::cin >> n) && n > 0)) {
		throw std::invalid_argument("The input should first contain a positive integer n, the number of nodes");
	}
	std::vector<size_t> directed(n - 1);
	for (size_t i = 0; i < n - 1; i++) {
		if (! ((std::cin >> directed[i]) && directed[i] <= i)) {
			throw std::invalid_argument("The second line of the input should contain n-1 integers p_1, …, p_{n-1} with 0<=p_i<i, the directed parents of 1, …, n-1");
		}
	}
	if (! ((std::cin >> m))) {
		throw std::invalid_argument("The third line of the input should contain a non-negative integer m, the number of bidirected edges");
	}
	std::vector<std::pair<size_t, size_t>> bidirected(m);
	for (size_t i = 0; i < m; i++) {
		size_t u, v;
		if (! ((std::cin >> u >> v) && u != v && u < n && v < n)) {
			throw std::invalid_argument("Lines 4 to 4+m of the input should each contain two integers u and v with 0<=u,v<n and u!=v, indicating that there is a bidirected edge between u and v");
		}
		if (u > v) std::swap(u, v);
		bidirected[i] = {u, v};
	}

#if defined(HAVE_OPENSSL)
	auto identificationResult = fasttreeid::identify(bidirected, directed, seed, prime, isProbablyPrimeOPENSSL);
#elif defined(HAVE_GMP)
	auto identificationResult = fasttreeid::identify(bidirected, directed, seed, prime, isProbablyPrimeGMP);
#else
	auto identificationResult = fasttreeid::identify(bidirected, directed, seed, prime, isProbablyPrimeError);
#endif

	if (verbosity != minimal) {
		std::cout << "Identification completed with seed = " << identificationResult.seed << " and prime = " << identificationResult.prime << std::endl << "Result:" << std::endl;
	}
	for (size_t i = 1; i < n; i++) {
		if (identificationResult.identification[i]) {
			if (std::holds_alternative<fasttreeid::Fraction>(*identificationResult.identification[i])) {
				if (verbosity != minimal) {
					std::cout << i << " is 1-identifiable as "
						<< std::get<fasttreeid::Fraction>(*identificationResult.identification[i]) << std::endl;
				} else {
					std::cout << "1 ";
				}
			} else if (std::holds_alternative<fasttreeid::Cycle>(*identificationResult.identification[i])) {
				if (std::get<fasttreeid::Cycle>(*identificationResult.identification[i]).isTwoIdentifiable()) {
					if (verbosity != minimal) {
						std::cout << i << " is 2-identifiable";
					} else {
						std::cout << "2 ";
					}
				} else {
					if (verbosity != minimal) {
						std::cout << i << " is 1-identifiable";
					} else {
						std::cout << "1 ";
					}
				}
				if (verbosity != minimal) {
					std::cout << " by the cycle " << std::get<fasttreeid::Cycle>(*identificationResult.identification[i]) << std::endl;
				}
			} else {
				size_t identifiedBy = std::get<fasttreeid::Path>(*identificationResult.identification[i]).getBack();
				if (std::holds_alternative<fasttreeid::Fraction>(*identificationResult.identification[identifiedBy])
						|| !std::get<fasttreeid::Cycle>(*identificationResult.identification[identifiedBy]).isTwoIdentifiable()) {
					if (verbosity != minimal) {
						std::cout << i << " is 1-identifiable";
					} else {
						std::cout << "1 ";
					}
				} else {
					if (verbosity != minimal) {
						std::cout << i << " is 2-identifiable";
					} else {
						std::cout << "2 ";
					}
				}
				if (verbosity != minimal) {
					std::cout << " by using the rank 2 edges on the path " << i
						<< std::get<fasttreeid::Path>(*identificationResult.identification[i]) << std::endl;
				}
			}
		} else {
			if (verbosity != minimal) {
				std::cout << i << " is unidentifiable" << std::endl;
			} else {
				std::cout << "0 ";
			}
		}
	}
	if (verbosity == minimal) {
		std::cout << std::endl;
	}
}
