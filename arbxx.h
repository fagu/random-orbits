#include <cassert>
#include <random>
#include "fmpzxx.h"
#if __FLINT_VERSION == 2
#include <arb.h>
#else
#include <flint/arb.h>
#endif
using namespace std;

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
		arb_set_fmpz(inner, x._fmpz());
	}
	arbxx(const arbxx& a) {
		arb_init(inner);
		arb_set(inner, a.inner);
	}
	// arbxx(arbxx&& a) {
	// 	inner[0] = a.inner[0];
	// 	a.destroyed = true;
	// }
	// 2^e
	static arbxx pow2(fmpzxx e) {
		arbxx res;
		fmpzxx one(1);
		arb_set_fmpz_2exp(res.inner, one._fmpz(), e._fmpz());
		return res;
	}
	arbxx& operator=(const arbxx& a) {
		arb_set(inner, a.inner);
		return *this;
	}
	// arbxx& operator=(arbxx&& a) {
	// 	if (!destroyed)
	// 		arb_clear(inner);
	// 	inner[0] = a.inner[0];
	// 	a.destroyed = true;
	// 	return *this;
	// }
	~arbxx() {
		// if (!destroyed)
			arb_clear(inner);
	}
};

inline arbxx neg(const arbxx& a) {
	arbxx res;
	arb_neg(res.inner, a.inner);
	return res;
}

inline arbxx add(const arbxx& a, const arbxx& b, slong prec) {
	arbxx res;
	arb_add(res.inner, a.inner, b.inner, prec);
	return res;
}

// inline arbxx add(arbxx&& a, const arbxx& b, slong prec) {
// 	// assert(false);
// 	arbxx res(std::move(a));
// 	arb_add(res.inner, res.inner, b.inner, prec);
// 	return res;
// }

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

inline arbxx log(const arbxx& a, slong prec) {
	arbxx res;
	arb_log(res.inner, a.inner, prec);
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

inline bool certainly_negative(const arbxx& a) {
	return arb_is_negative(a.inner);
}

inline bool certainly_positive(const arbxx& a) {
	return arb_is_positive(a.inner);
}

inline bool lt(const arbxx& a, const arbxx& b) {
	if (arb_lt(a.inner, b.inner))
		return true;
	if (arb_ge(a.inner, b.inner))
		return false;
	throw insufficient_precision();
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

inline fmpzxx to_int(const arbxx& a) {
	fmpzxx res;
	if (!arb_get_unique_fmpz(res._fmpz(), a.inner))
		throw insufficient_precision();
	return res;
}

inline string to_str(const arbxx& a, slong digits) {
	char* s = arb_get_str(a.inner, digits, 0);
	string res(s);
	free(s);
	return res;
}

inline arbxx unionn(const arbxx& a, const arbxx& b, slong prec) {
	arbxx res;
	arb_union(res.inner, a.inner, b.inner, prec);
	return res;
}

// A random real number in the interval [0,1) with given precision.
class rand_real {
private:
	fmpzxx numerator;
	slong prec;
public:
	rand_real() : numerator(0), prec(0) {}
	// Using the given random number generator, generates bits of the random number until it exactly new_prec bits (after the decimal) are fixed.
	template<class Generator>
	void refine(slong new_prec, Generator& gen) {
		assert(new_prec >= prec);
		numerator *= fmpzxx::exp2(new_prec - prec);
		for (slong i = new_prec - prec - 1; i >= 0; i--)
			if (uniform_int_distribution<slong>(0,1)(gen))
				numerator.set_bit(i);
		prec = new_prec;
	}
	// Returns the interval with current precision contain the random number.
	arbxx get() const {
		fmpzxx num_f(2*numerator+1);
		fmpzxx prec_f(-(prec+1));
		arbxx res;
		arb_set_fmpz_2exp(res.inner, num_f._fmpz(), prec_f._fmpz());
		arb_add_error_2exp_si(res.inner, -(prec+1));
		return res;
	}
};
