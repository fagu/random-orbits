// C++ wrapper of part of the Arb library

#ifndef RANDOM_ORBIT_ARBXX_H
#define RANDOM_ORBIT_ARBXX_H

#include <cassert>
#include <random>
#include "fmpzxx.h"
#if __FLINT_VERSION == 2
#include <arb.h>
#else
#include <flint/arb.h>
#endif

using namespace std;

// Macros for operations in ball arithmetic that use the variable "prec" as the precision argument.
#define ADD(a,b) add(a, b, prec)
#define SUB(a,b) sub(a, b, prec)
#define MUL(a,b) mul(a, b, prec)
#define DIV(a,b) div(a, b, prec)
#define POW_UI(a,b) pow_ui(a, b, prec)
#define POW_SI(a,b) pow_si(a, b, prec)
#define ROOT(a,b) root(a, b, prec)
#define INV(a) inv(a, prec)
#define SQRT(a) sqrt(a, prec)
#define LOG(a) log(a, prec)
#define FLOOR(a) floor(a, prec)
#define CEIL(a) ceil(a, prec)

// Exception indicating that the precision was insufficient to decide a predicate, find the unique integer in an interval, etc.
struct insufficient_precision {};

class arbxx {
public:
	// bool destroyed = false;
	arb_t inner;
	arbxx() {
		arb_init(inner);
	}
	arbxx(slong x) {
		arb_init(inner);
		arb_set_si(inner, x);
	}
	arbxx(fmpzxx x) {
		arb_init(inner);
		arb_set_fmpz(inner, x.inner);
	}
	arbxx(const arbxx& a) {
		arb_init(inner);
		arb_set(inner, a.inner);
	}
	arbxx& operator=(const arbxx& a) {
		arb_set(inner, a.inner);
		return *this;
	}
	~arbxx() {
		arb_clear(inner);
	}
	friend arbxx operator-(const arbxx& a) {
		arbxx res;
		arb_neg(res.inner, a.inner);
		return res;
	}
	// Whether a < b. May throw insufficient_precision.
	friend bool operator<(const arbxx& a, const arbxx& b) {
		if (arb_lt(a.inner, b.inner))
			return true;
		if (arb_ge(a.inner, b.inner))
			return false;
		throw insufficient_precision();
	}
	// The integer in the given interval. Throws insufficient_precision if there is not *exactly* one such integer.
	friend fmpzxx to_int(const arbxx& a) {
		fmpzxx res;
		if (!arb_get_unique_fmpz(res.inner, a.inner))
			throw insufficient_precision();
		return res;
	}
};

inline arbxx add(const arbxx& a, const arbxx& b, slong prec) {
	arbxx res;
	arb_add(res.inner, a.inner, b.inner, prec);
	return res;
}

inline arbxx sub(const arbxx& a, const arbxx& b, slong prec) {
	arbxx res;
	arb_sub(res.inner, a.inner, b.inner, prec);
	return res;
}

inline arbxx mul(const arbxx& a, const arbxx& b, slong prec) {
	arbxx res;
	arb_mul(res.inner, a.inner, b.inner, prec);
	return res;
}

inline arbxx div(const arbxx& a, const arbxx& b, slong prec) {
	arbxx res;
	arb_div(res.inner, a.inner, b.inner, prec);
	return res;
}

inline arbxx inv(const arbxx& a, slong prec) {
	arbxx res;
	arb_inv(res.inner, a.inner, prec);
	return res;
}

inline arbxx sqrt(const arbxx& a, slong prec) {
	arbxx res;
	arb_sqrt(res.inner, a.inner, prec);
	return res;
}

inline arbxx root(const arbxx& a, ulong k, slong prec) {
	arbxx res;
	arb_root_ui(res.inner, a.inner, k, prec);
	return res;
}

inline arbxx pow_ui(const arbxx& a, ulong e, slong prec) {
	arbxx res;
	arb_pow_ui(res.inner, a.inner, e, prec);
	return res;
}

inline arbxx pow_si(const arbxx& a, slong e, slong prec) {
	if (e == 0)
		return 1;
	else if (e < 0)
		return pow_ui(inv(a, prec), -e, prec);
	else
		return pow_ui(a, e, prec);
}

inline arbxx floor(const arbxx& a, slong prec) {
	arbxx res;
	arb_floor(res.inner, a.inner, prec);
	return res;
}

inline arbxx ceil(const arbxx& a, slong prec) {
	arbxx res;
	arb_ceil(res.inner, a.inner, prec);
	return res;
}

// A random real number in the interval [0,1) with given precision.
class rand_real {
private:
	fmpzxx numerator;
	slong prec;
public:
	rand_real() : numerator(0), prec(0) {}
	// Using the given random number generator, generates more bits of the random number until exactly new_prec bits (after the decimal) are fixed.
	template<class Generator>
	void refine(slong new_prec, Generator& gen) {
		assert(new_prec >= prec);
		numerator *= fmpzxx::exp2(new_prec - prec);
		for (slong i = new_prec - prec - 1; i >= 0; i--)
			if (uniform_int_distribution<slong>(0,1)(gen))
				numerator.set_bit(i);
		prec = new_prec;
	}
	// Returns the interval containing the random number with the current precision.
	arbxx get() const {
		fmpzxx num_f(2*numerator+1);
		fmpzxx prec_f(-(prec+1));
		arbxx res;
		arb_set_fmpz_2exp(res.inner, num_f.inner, prec_f.inner);
		arb_add_error_2exp_si(res.inner, -(prec+1));
		return res;
	}
};

#endif
