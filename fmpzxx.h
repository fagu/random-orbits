#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>

class fmpzxx {
public:
	fmpz_t inner;
	fmpzxx() {
		fmpz_init(inner);
	}
	fmpzxx(slong val) {
		fmpz_init(inner);
		fmpz_set_si(inner, val);
	}
	fmpzxx(fmpzxx& x) {
		fmpz_init(inner);
		fmpz_set(inner, x.inner);
	}
	fmpzxx& operator=(const fmpzxx& x) {
		fmpz_set(inner, x.inner);
		return *this;
	}
	~fmpzxx() {
		fmpz_clear(inner);
	}
	fmpz_t& _fmpz() {
		return inner;
	}
	// 2^e
	static fmpzxx exp2(ulong e) {
		fmpzxx res = 1;
		fmpz_mul_2exp(res.inner, res.inner, e);
		return res;
	}
	friend fmpzxx& operator*=(fmpzxx& a, const fmpzxx& b) {
		fmpz_mul(a.inner, a.inner, b.inner);
		return a;
	}
	void set_bit(ulong i) {
		fmpz_setbit(inner, i);
	}
	friend fmpzxx operator+(const fmpzxx& a, const fmpzxx& b) {
		fmpzxx res;
		fmpz_add(res.inner, a.inner, b.inner);
		return res;
	}
	friend fmpzxx operator-(const fmpzxx& a, const fmpzxx& b) {
		fmpzxx res;
		fmpz_sub(res.inner, a.inner, b.inner);
		return res;
	}
	friend fmpzxx operator*(const fmpzxx& a, const fmpzxx& b) {
		fmpzxx res;
		fmpz_mul(res.inner, a.inner, b.inner);
		return res;
	}
	friend fmpzxx operator-(const fmpzxx& a) {
		fmpzxx res;
		fmpz_neg(res.inner, a.inner);
		return res;
	}
	friend fmpzxx cdiv_q(const fmpzxx& a, const fmpzxx& b) {
		fmpzxx res;
		fmpz_cdiv_q(res.inner, a.inner, b.inner);
		return res;
	}
	fmpzxx abs() const {
		fmpzxx res;
		fmpz_abs(res.inner, inner);
		return res;
	}
	friend bool operator==(const fmpzxx& a, const fmpzxx& b) {
		return fmpz_equal(a.inner, b.inner);
	}
	friend bool operator!=(const fmpzxx& a, const fmpzxx& b) {
		return !fmpz_equal(a.inner, b.inner);
	}
	friend bool operator>(const fmpzxx& a, const fmpzxx& b) {
		return fmpz_cmp(a.inner, b.inner) > 0;
	}
	friend bool operator<=(const fmpzxx& a, const fmpzxx& b) {
		return fmpz_cmp(a.inner, b.inner) <= 0;
	}
	friend bool operator<(const fmpzxx& a, const fmpzxx& b) {
		return fmpz_cmp(a.inner, b.inner) < 0;
	}
	int sgn() const {
		return fmpz_sgn(inner);
	}
	bool is_square() const {
		return fmpz_is_square(inner);
	}
};

class fmpz_polyxx {
public:
	fmpz_poly_t inner;
	fmpz_polyxx() {
		fmpz_poly_init(inner);
	}
	fmpz_polyxx(const fmpz_polyxx& x) {
		fmpz_poly_init(inner);
		fmpz_poly_set(inner, x.inner);
	}
	fmpz_polyxx& operator=(const fmpz_polyxx& x) {
		fmpz_poly_set(inner, x.inner);
		return *this;
	}
	~fmpz_polyxx() {
		fmpz_poly_clear(inner);
	}
	void set_coeff(slong n, const fmpzxx& c) {
		fmpz_poly_set_coeff_fmpz(inner, n, c.inner);
	}
	fmpzxx coeff(slong n) {
		fmpzxx res;
		fmpz_poly_get_coeff_fmpz(res.inner, inner, n);
		return res;
	}
	fmpzxx disc() const {
		fmpzxx res;
		fmpz_poly_discriminant(res.inner, inner);
		return res;
	}
	fmpzxx content() const {
		fmpzxx res;
		fmpz_poly_content(res.inner, inner);
		return res;
	}
	// Returns whether a given nonconstant polynomial is irreducible over Q.
	bool irreducible_over_Q() const {
		fmpz_poly_factor_t pfac;
		fmpz_poly_factor_init(pfac);
		fmpz_poly_factor(pfac, inner);
		// FLINT removes the content of the polynomial before factorization, so this is
		// effectively a factorization over Q.
		// Hence, f is irreducible over Q if and only if the number of factors is 1.
		bool irred = pfac->num == 1;
		fmpz_poly_factor_clear(pfac);
		return irred;
	}
};
