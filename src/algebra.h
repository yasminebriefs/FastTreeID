#pragma once

#include "random.h"

#include <memory>
#include <ostream>
#include <vector>

namespace fasttreeid {

int64_t multiply (int64_t, int64_t, int64_t) noexcept;

class Variable {
public:
	friend std::ostream& operator<< (std::ostream &, const Variable &);

private:
	virtual std::ostream& print (std::ostream &) const = 0;
};

class Lambda : public Variable {
public:
	Lambda (const std::shared_ptr<RandomTool> &, size_t, size_t);

	int64_t getValue () const noexcept;

private:
	size_t x, y;
	const int64_t value;

	std::ostream& print (std::ostream &) const override;
};

class Omega : public Variable {
public:
	Omega (const std::shared_ptr<RandomTool> &, size_t, size_t);

	int64_t getValue () const noexcept;

private:
	size_t x, y;
	const int64_t value;

	std::ostream& print (std::ostream &) const override;
};

class Sigma : public Variable {
public:
	Sigma (size_t, size_t);

	void setInTermsOfLambdaAndOmega (int64_t) noexcept;

	int64_t getInTermsOfLambdaAndOmega () const noexcept;

	size_t getX() const noexcept;
	size_t getY() const noexcept;

private:
	size_t x, y;
	int64_t inTermsOfLambdaAndOmega;

	std::ostream& print (std::ostream &) const override;
};

class FASTP {
public:
	FASTP () noexcept;
	FASTP (int64_t, int64_t, int64_t, int64_t, int64_t) noexcept;

	int64_t getP () const noexcept;
	int64_t getQ () const noexcept;
	int64_t getR () const noexcept;
	int64_t getS () const noexcept;
	int64_t getT () const noexcept;

	void setP (int64_t) noexcept;
	void setQ (int64_t) noexcept;
	void setR (int64_t) noexcept;
	void setS (int64_t) noexcept;
	void setT (int64_t) noexcept;

	friend std::ostream& operator<< (std::ostream &, const FASTP &);

	friend int64_t getA (const FASTP &, const FASTP&, int64_t, int64_t, int64_t, int64_t, int64_t) noexcept;
	friend int64_t getB (const FASTP &, const FASTP&, int64_t, int64_t, int64_t, int64_t, int64_t) noexcept;

private:
	int64_t p, q, r, s, t;
};

class Fraction {
public:
	Fraction (std::shared_ptr<Sigma>, std::shared_ptr<Sigma>) noexcept;

	std::shared_ptr<Sigma> getP () const noexcept;
	std::shared_ptr<Sigma> getQ () const noexcept;

	friend std::ostream& operator<< (std::ostream &, const Fraction &);

private:
	std::shared_ptr<Sigma> p, q;
};

class Cycle {
public:
	enum Identifiability {
		twoIdentifiable,
		oneIdentifiableAZero,
		oneIdentifiableDiscriminantZero,
		oneIdentifiableOneOption
	};

	Cycle (const std::vector<size_t> &);

	void setTwoIdentifiable () noexcept;
	void setOneIdentifiableAZero (size_t, size_t) noexcept;
	void setOneIdentifiableDiscriminantZero () noexcept;
	void setOneIdentifiableOneOption (size_t, size_t) noexcept;

	bool isTwoIdentifiable () const noexcept;

	std::vector<size_t> getNodes () const;
	Identifiability getIdentifiability () const noexcept;
	size_t getReasonI () const;
	size_t getReasonJ () const;

	friend std::ostream& operator<< (std::ostream &, const Cycle &);

private:

	std::vector<size_t> nodes;
	Identifiability identifiability;
	size_t reasonI, reasonJ;
};

class Path {
public:
	Path (const std::vector<size_t> &);

	size_t getBack () const noexcept;
	std::vector<size_t> getNodes () const;

	friend std::ostream& operator<< (std::ostream &, const Path &);

private:
	std::vector<size_t> nodes;
};

}
