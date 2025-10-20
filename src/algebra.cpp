#include "algebra.h"

namespace fasttreeid {

std::ostream& operator<< (std::ostream &os, const Variable &var) {
	return var.print(os);
}

int64_t multiply (int64_t a, int64_t b, int64_t prime) noexcept {
#if defined(__SIZEOF_INT128__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"

	__int128 A = a, B = b;
	return (int64_t)(A * B % prime);

#pragma GCC diagnostic pop
#else

	auto multiplyReduce = [](int64_t left, int64_t right, int64_t prime) {
		/* Assumes 0 <= left, right < 2^30 and 2^58 < prime < 2^59. Returns left * right % prime */
		int64_t ans = left * right;
		if (ans >= prime) ans -= prime;
		if (ans >= prime) ans -= prime;
		if (ans >= prime) ans -= prime;
		return ans;
	};

	auto double30Times = [](int64_t num, int64_t prime) {
		/* Assumes 0 <= num < prime. Returns num * 2^30 % prime */
		for (size_t i = 0; i < 30; i++) {
			num *= 2;
			if (num >= prime) num -= prime;
		}
		return num;
	};

	int64_t a0 = a & ((1 << 30) - 1), a1 = a >> 30;
	int64_t b0 = b & ((1 << 30) - 1), b1 = b >> 30;

	int64_t tmp = multiplyReduce(a0, b1, prime) + multiplyReduce(a1, b0, prime) + double30Times(multiplyReduce(a1, b1, prime), prime);
	if (tmp >= prime) tmp -= prime;
	if (tmp >= prime) tmp -= prime;

	int64_t res = multiplyReduce(a0, b0, prime) + double30Times(tmp, prime);
	if (res >= prime) res -= prime;

	return res;

#endif
}

Lambda::Lambda (const std::shared_ptr<RandomTool> &randomTool, size_t x_, size_t y_) : x(x_), y(y_), value(randomTool->getRandomSmallerPrime()) {}

Omega::Omega (const std::shared_ptr<RandomTool> &randomTool, size_t x_, size_t y_) : x(x_), y(y_), value(randomTool->getRandomSmallerPrime()) {}

Sigma::Sigma (size_t x_, size_t y_) : x(x_), y(y_), inTermsOfLambdaAndOmega(0) {}

std::ostream& Lambda::print (std::ostream &os) const {
	return os << "λ_" << this->x << "," << this->y;
}

std::ostream& Omega::print (std::ostream &os) const {
	return os << "ω_" << this->x << "," << this->y;
}

std::ostream& Sigma::print (std::ostream &os) const {
	return os << "σ_" << this->x << "," << this->y;
}

int64_t Lambda::getValue () const noexcept {
	return this->value;
}

int64_t Omega::getValue () const noexcept {
	return this->value;
}

void Sigma::setInTermsOfLambdaAndOmega (int64_t value) noexcept {
	this->inTermsOfLambdaAndOmega = value;
}

int64_t Sigma::getInTermsOfLambdaAndOmega () const noexcept {
	return this->inTermsOfLambdaAndOmega;
}

size_t Sigma::getX () const noexcept {
	return this->x;
}

size_t Sigma::getY () const noexcept {
	return this->y;
}

FASTP::FASTP () noexcept {}

FASTP::FASTP (int64_t p_, int64_t q_, int64_t r_, int64_t s_, int64_t t_) noexcept
	: p(p_), q(q_), r(r_), s(s_), t(t_) {}

int64_t FASTP::getP () const noexcept {
	return this->p;
}

int64_t FASTP::getQ () const noexcept {
	return this->q;
}

int64_t FASTP::getR () const noexcept {
	return this->r;
}

int64_t FASTP::getS () const noexcept {
	return this->s;
}

int64_t FASTP::getT () const noexcept {
	return this->t;
}

void FASTP::setP (int64_t p_) noexcept {
	this->p = p_;
}

void FASTP::setQ (int64_t q_) noexcept {
	this->q = q_;
}

void FASTP::setR (int64_t r_) noexcept {
	this->r = r_;
}

void FASTP::setS (int64_t s_) noexcept {
	this->s = s_;
}

void FASTP::setT (int64_t t_) noexcept {
	this->t = t_;
}

std::ostream& operator<< (std::ostream &os, const FASTP &fastp) {
	os << "(" << fastp.p << " + " << fastp.q << " sqrt(" << fastp.s << "))";
	os << "/(" << fastp.r << " + " << fastp.t << " sqrt(" << fastp.s << "))";
	return os;
}

int64_t getA (const FASTP &i, const FASTP &j, int64_t prime,
		int64_t sigma_pq, int64_t sigma_iq, int64_t sigma_pj, int64_t sigma_ij) noexcept {

	int64_t p = i.getP(), q = i.getQ(), r = i.getR(), s = i.getS(), t = i.getT();
	int64_t P = j.getP(), Q = j.getQ(), R = j.getR(), T = j.getT();

	int64_t result = multiply(multiply(p, P, prime), sigma_pq, prime);
	result -= multiply(multiply(p, R, prime), sigma_pj, prime);
	if (result < 0) result += prime;
	result -= multiply(multiply(r, P, prime), sigma_iq, prime);
	if (result < 0) result += prime;
	result += multiply(multiply(r, R, prime), sigma_ij, prime);
	if (result >= prime) result -= prime;
	result += multiply(multiply(multiply(q, Q, prime), s, prime), sigma_pq, prime);
	if (result >= prime) result -= prime;
	result -= multiply(multiply(multiply(q, T, prime), s, prime), sigma_pj, prime);
	if (result < 0) result += prime;
	result -= multiply(multiply(multiply(t, Q, prime), s, prime), sigma_iq, prime);
	if (result < 0) result += prime;
	result +=  multiply(multiply(multiply(t, T, prime), s, prime), sigma_ij, prime);
	if (result >= prime) result -= prime;

	return result;
}

int64_t getB (const FASTP &i, const FASTP &j, int64_t prime,
		int64_t sigma_pq, int64_t sigma_iq, int64_t sigma_pj, int64_t sigma_ij) noexcept {
	
	int64_t p = i.getP(), q = i.getQ(), r = i.getR(), t = i.getT();
	int64_t P = j.getP(), Q = j.getQ(), R = j.getR(), T = j.getT();

	int64_t result = multiply(multiply(p, Q, prime), sigma_pq, prime);
	result += multiply(multiply(P, q, prime), sigma_pq, prime);
	if (result >= prime) result -= prime;
	result -= multiply(multiply(p, T, prime), sigma_pj, prime);
	if (result < 0) result += prime;
	result -= multiply(multiply(q, R, prime), sigma_pj, prime);
	if (result < 0) result += prime;
	result -= multiply(multiply(t, P, prime), sigma_iq, prime);
	if (result < 0) result += prime;
	result -= multiply(multiply(r, Q, prime), sigma_iq, prime);
	if (result < 0) result += prime;
	result += multiply(multiply(r, T, prime), sigma_ij, prime);
	if (result >= prime) result -= prime;
	result += multiply(multiply(R, t, prime), sigma_ij, prime);
	if (result >= prime) result -= prime;

	return result;
}

Fraction::Fraction (std::shared_ptr<Sigma> p_, std::shared_ptr<Sigma> q_) noexcept : p(std::move(p_)), q(std::move(q_)) {}

std::shared_ptr<Sigma> Fraction::getP () const noexcept {
	return this->p;
}

std::shared_ptr<Sigma> Fraction::getQ () const noexcept {
	return this->q;
}

std::ostream& operator<< (std::ostream &os, const Fraction &fraction) {
	return os << "(" << *fraction.p << ")/(" << *fraction.q << ")";
}

Cycle::Cycle (const std::vector<size_t> &nodes_) : nodes(nodes_), identifiability(Cycle::twoIdentifiable) {}

void Cycle::setTwoIdentifiable () noexcept {
	this->identifiability = Cycle::twoIdentifiable;
}

void Cycle::setOneIdentifiableAZero (size_t reasonI_, size_t reasonJ_) noexcept {
	this->identifiability = Cycle::oneIdentifiableAZero;
	this->reasonI = reasonI_;
	this->reasonJ = reasonJ_;
}

void Cycle::setOneIdentifiableDiscriminantZero () noexcept {
	this->identifiability = Cycle::oneIdentifiableDiscriminantZero;
}

void Cycle::setOneIdentifiableOneOption (size_t reasonI_, size_t reasonJ_) noexcept {
	this->identifiability = Cycle::oneIdentifiableOneOption;
	this->reasonI = reasonI_;
	this->reasonJ = reasonJ_;
}

bool Cycle::isTwoIdentifiable () const noexcept {
	return this->identifiability == Cycle::twoIdentifiable;
}

std::vector<size_t> Cycle::getNodes () const {
	return this->nodes;
}

Cycle::Identifiability Cycle::getIdentifiability () const noexcept {
	return this->identifiability;
}

size_t Cycle::getReasonI () const {
	if (this->identifiability != Cycle::oneIdentifiableAZero && this->identifiability != Cycle::oneIdentifiableOneOption) {
		throw std::logic_error("getReasonI: value not set");
	}
	return this->reasonI;
}

size_t Cycle::getReasonJ () const {
	if (this->identifiability != Cycle::oneIdentifiableAZero && this->identifiability != Cycle::oneIdentifiableOneOption) {
		throw std::logic_error("getReasonJ: value not set");
	}
	return this->reasonJ;
}

std::ostream& operator<< (std::ostream &os, const Cycle &cycle) {
	for (size_t v : cycle.nodes) {
		os << v << "-";
	}
	os << cycle.nodes[0];
	if (!cycle.isTwoIdentifiable()) {
		os << ". Reason: ";
		switch (cycle.identifiability) {
			case Cycle::oneIdentifiableAZero:
				os << "a = 0 for λ_" << cycle.reasonI << "," << cycle.reasonJ;
				break;
			case Cycle::oneIdentifiableDiscriminantZero:
				os << "Δ = 0";
				break;
			case Cycle::oneIdentifiableOneOption:
				os << "Equation for missing edge {" << cycle.reasonI << ", " << cycle.reasonJ
					<< "} only satisfied by one of the options";
				break;
			default:
				break;
		}
	}
	return os;
}

Path::Path (const std::vector<size_t> &nodes_) : nodes(nodes_) {}

size_t Path::getBack () const noexcept {
	return this->nodes.back();
}

std::vector<size_t> Path::getNodes () const {
	return this->nodes;
}

std::ostream& operator<< (std::ostream &os, const Path &path) {
	for (size_t v : path.nodes) {
		os << "->" << v;
	}
	return os;
}

}
