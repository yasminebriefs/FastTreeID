#include "random.h"

#include <stdexcept>

#if defined(HAVE_GMP)
#include <gmp.h>
#elif defined(HAVE_OPENSSL)
#include <openssl/bn.h>
#endif

namespace fasttreeid {

uint64_t deriveGeneralSeed (uint64_t seed) {
	std::mt19937_64 gen(seed);
	return gen();
}

uint64_t derivePrimeSeed (uint64_t seed) {
	std::mt19937_64 gen(seed);
	gen();
	return gen();
}

bool isProbablyPrime ([[maybe_unused]] int64_t candidate, [[maybe_unused]] int numchecks = 25) { /* TODO: number of checks? */

#if defined(HAVE_GMP)

	mpz_t z;
	mpz_init(z);
	mpz_import(z, 1, -1, sizeof(candidate), 0, 0, &candidate);
	int result = mpz_probab_prime_p(z, numchecks);
	mpz_clear(z);
	return result > 0;

#elif defined(HAVE_OPENSSL)

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
	int result = BN_is_prime_ex(bn, numchecks, ctx, nullptr);
#endif

	BN_CTX_free(ctx);
	BN_free(bn);
	if (result < 0) {
		throw std::runtime_error("OpenSSL prime check failed with error code < 0");
	}
	return result > 0;

#else

	throw std::logic_error("Neither HAVE_GMP nor HAVE_OPENSSL defined but no prime given");

#endif
}

int64_t generateCandidate (std::mt19937_64 &gen) {
	return std::uniform_int_distribution<int64_t>(1ll << 58, (1ll << 59) - 1)(gen) | 1ll;
}

int64_t generatePrime (uint64_t prime_seed) {
	std::mt19937_64 gen(prime_seed);

	int64_t candidate = generateCandidate(gen);
	while (!isProbablyPrime(candidate)) {
		candidate = generateCandidate(gen);
	}
	return candidate;
}

RandomTool::RandomTool (uint64_t seed, std::optional<int64_t> prime_) : rng(deriveGeneralSeed(seed)), prime(prime_.has_value() ? prime_.value() : generatePrime(derivePrimeSeed(seed))) {}

int64_t RandomTool::getRandom (int64_t l, int64_t r) {
	if (l > r) {
		throw std::logic_error("getRandom: Call with invalid bounds. Please report this error.");
	}
	return std::uniform_int_distribution<int64_t>(l, r)(this->rng);
}

int64_t RandomTool::getRandomSmallerPrime () {
	return this->getRandom(1, prime - 1);
}

int64_t RandomTool::getPrime () const noexcept {
	return this->prime;
}

}
